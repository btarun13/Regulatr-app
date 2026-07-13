#!/usr/bin/env python3
"""
regulatr — cell-type-aware non-coding variant prioritization (webapp backend).

Layer-1 deterministic tools (per SPEC.md §6) + a thin HTTP API that composes them.
Zero third-party dependencies: stdlib only, so it runs anywhere Python 3.8+ exists.

Pipeline for a single ENCODE accession:
    resolve_accession -> fetch_peaks -> fetch_variants -> overlap
                      -> map_to_gene -> opentargets_score -> rank

Scope guard (SPEC.md §0): accessibility (ATAC peaks) tells us WHERE regulatory DNA
is open; variants come from ClinVar. We never call variants from ATAC reads.
Internal coordinate standard is hg38 / GRCh38 (SPEC.md §1).
"""

import gzip
import io
import json
import os
import re
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
STATIC = os.path.join(HERE, "static")
CACHE = os.path.join(HERE, ".cache")
os.makedirs(CACHE, exist_ok=True)

ENCODE = "https://www.encodeproject.org"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL = "https://rest.ensembl.org"
OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"
GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

# Locked worked example (SPEC.md §14): SORT1/CELSR2/PSRC1 locus.
DEFAULT_REGION = "chr1:109200000-109350000"
MAX_VARIANTS = 1000  # cap on ClinVar records pulled per region (keeps the request responsive)
ENRICH_MAX_GENES = 500  # cap on the ordered gene list sent to g:Profiler
ENRICH_SOURCES = ["GO:BP", "REAC", "KEGG"]  # gene-set collections queried for enrichment

# Region discovery: window built around each accessible hotspot, and how many to offer.
REGION_HALF = 75000        # ± bp around a peak summit -> ~150 kb candidate window
REGION_MIN_SEP = 250000    # min bp between two chosen hotspots (so windows don't pile up)
MAX_REGIONS = 3000         # safety cap on how many hotspots to surface (default: all)
GENE_TILE = 4_000_000      # bp window for one Ensembl gene fetch (< Ensembl's 5 Mb limit)
GENE_LABEL_BUDGET = 60     # max gene-tile fetches per discovery (strongest hotspots first)
STD_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}

# Gene tiles are immutable — cache them across requests so repeat discoveries are fast.
_GENE_TILE_CACHE = {}

# Preference order for which ENCODE peak file to use when an experiment has several.
PEAK_PREF = [
    "conservative IDR thresholded peaks",
    "IDR thresholded peaks",
    "pseudoreplicated peaks",
    "replicated peaks",
    "stable peaks",
    "peaks",
]

UA = {"User-Agent": "regulatr/0.1 (hackathon demo)"}


# --------------------------------------------------------------------------- #
# small HTTP helpers
# --------------------------------------------------------------------------- #
def _get_json(url, timeout=30):
    req = urllib.request.Request(url, headers={**UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def _post_json(url, payload, timeout=30):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        url, data=data, headers={**UA, "Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def parse_region(region):
    """'chr1:109200000-109350000' -> (chrom_no_prefix, chrom_with_prefix, start, end)."""
    m = re.match(r"(?:chr)?([\w]+):([\d,]+)-([\d,]+)", (region or "").strip())
    if not m:
        raise ValueError(f"Bad region '{region}'. Expected e.g. chr1:109200000-109350000")
    chrom = m.group(1)
    start = int(m.group(2).replace(",", ""))
    end = int(m.group(3).replace(",", ""))
    return chrom, "chr" + chrom, start, end


# --------------------------------------------------------------------------- #
# Layer-1 tool: resolve_accession  (SPEC.md §6.1/6.2)
# --------------------------------------------------------------------------- #
def resolve_accession(code):
    """ENCODE experiment (ENCSR) or file (ENCFF) accession -> chosen hg38 peak file + provenance."""
    code = code.strip().upper()
    if code.startswith("ENCSR"):
        exp = _get_json(f"{ENCODE}/experiments/{code}/?format=json")
        fields = "&".join(
            f"field={f}" for f in ("accession", "output_type", "file_size", "href", "assembly")
        )
        search = _get_json(
            f"{ENCODE}/search/?type=File&dataset=/experiments/{code}/"
            f"&file_format=bed&assembly=GRCh38&status=released&format=json&limit=all&{fields}"
        )
        files = search.get("@graph", [])
        if not files:
            raise ValueError(f"{code}: no released GRCh38 bed peak files found.")
        chosen = _pick_peak_file(files)
        return _provenance(chosen, exp, code)

    if code.startswith("ENCFF"):
        f = _get_json(f"{ENCODE}/files/{code}/?format=json")
        if f.get("file_format") != "bed":
            raise ValueError(f"{code} is not a bed file (got {f.get('file_format')}).")
        if f.get("assembly") != "GRCh38":
            raise ValueError(f"{code} is {f.get('assembly')}, not GRCh38 — build guard (SPEC §1).")
        exp_id = f.get("dataset", "")
        exp = _get_json(f"{ENCODE}{exp_id}?format=json") if exp_id else {}
        return _provenance(f, exp, exp.get("accession", exp_id.strip("/").split("/")[-1]))

    raise ValueError(f"'{code}' is not an ENCODE accession (expected ENCSR… or ENCFF…).")


def _pick_peak_file(files):
    def rank(f):
        ot = f.get("output_type", "")
        return (PEAK_PREF.index(ot) if ot in PEAK_PREF else len(PEAK_PREF), -(f.get("file_size") or 0))
    return sorted(files, key=rank)[0]


def _provenance(f, exp, exp_acc):
    onto = (exp or {}).get("biosample_ontology") or {}
    return {
        "file_accession": f.get("accession"),
        "href": f.get("href"),
        "output_type": f.get("output_type"),
        "file_size": f.get("file_size"),
        "genome": f.get("assembly", "GRCh38"),
        "experiment": exp_acc,
        "assay": (exp or {}).get("assay_title") or (exp or {}).get("assay_term_name"),
        "biosample": (exp or {}).get("biosample_summary")
        or onto.get("term_name")
        or "unknown biosample",
        "lab": ((exp or {}).get("lab") or {}).get("title"),
    }


# --------------------------------------------------------------------------- #
# Layer-1 tool: fetch_peaks  (SPEC.md §6.2)
# --------------------------------------------------------------------------- #
def _peak_path(href, file_accession):
    """Download the peak bed.gz once and cache it; return the local path."""
    path = os.path.join(CACHE, f"{file_accession}.bed.gz")
    if not os.path.exists(path):
        req = urllib.request.Request(ENCODE + href, headers=UA)
        with urllib.request.urlopen(req, timeout=120) as r, open(path, "wb") as out:
            out.write(r.read())
    return path


# Different ENCODE peak beds use different schemas. narrowPeak (headerless, 10 col)
# puts signalValue in col 7; some density beds carry a #header and put a genomic
# COORDINATE there — so we resolve the signal column from the header when present,
# preferring named signal fields and never a coordinate column.
_SIGNAL_NAMES = ["signalvalue", "smoothed_peak_height", "peak_height",
                 "max_density", "summit_density", "density", "signal", "fold_enrichment"]


def _signal_col_from_header(header_cols):
    names = [h.strip().lstrip("#").strip().lower() for h in header_cols]
    for want in _SIGNAL_NAMES:
        for i, n in enumerate(names):
            if n == want:
                return i
    return None


def _iter_bed(path):
    """Yield (chrom, start, end, signal, score) from a bed.gz, schema-aware."""
    signal_col = None
    header_seen = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == "#":
                if not header_seen:
                    header_seen = True
                    signal_col = _signal_col_from_header(line.rstrip("\n").split("\t"))
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            try:
                start, end = int(c[1]), int(c[2])
            except ValueError:
                continue
            score = _num(c[4]) if len(c) > 4 else 0.0
            if signal_col is not None and len(c) > signal_col:
                signal = _num(c[signal_col])          # named signal column from header
            elif len(c) > 6:
                signal = _num(c[6])                   # narrowPeak signalValue
            else:
                signal = score
            yield c[0], start, end, signal, score


def fetch_peaks(href, file_accession, chrom_p, start, end):
    """Return cached peaks overlapping the region (downloads the bed.gz on first use)."""
    path = _peak_path(href, file_accession)
    peaks = []
    for chrom, p_start, p_end, signal, score in _iter_bed(path):
        if chrom != chrom_p or p_end <= start or p_start >= end:
            continue
        peaks.append({"chrom": chrom, "start": p_start, "end": p_end,
                      "score": score, "signal": signal})
    peaks.sort(key=lambda p: p["start"])
    return peaks


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


# --------------------------------------------------------------------------- #
# Region discovery: derive the dropdown from the accession itself.
# The "right" regions for a cell type are where THAT experiment has the strongest
# open chromatin — not a generic gene list. Scan the peak file genome-wide, pick the
# top hotspots (spaced apart), and label each by its nearest gene.
# --------------------------------------------------------------------------- #
def _load_all_peaks(href, file_accession):
    """Parse the whole cached bed.gz (standard chromosomes only)."""
    path = _peak_path(href, file_accession)
    peaks = []
    for chrom, start, end, signal, score in _iter_bed(path):
        if chrom not in STD_CHROMS:
            continue
        peaks.append({"chrom": chrom, "start": start, "end": end,
                      "signal": signal, "score": score})
    return peaks


def _nearest_symbol(genes, center):
    """Nearest (protein-coding preferred) gene symbol to a position, from a gene list."""
    if not genes:
        return None
    coding = [g for g in genes if g.get("biotype") == "protein_coding"] or genes
    best, best_dist = None, None
    for g in coding:
        if g["start"] <= center <= g["end"]:
            return g["symbol"]
        d = g["start"] - center if center < g["start"] else center - g["end"]
        if best_dist is None or d < best_dist:
            best, best_dist = g, d
    return best["symbol"] if best else None


def _fetch_tile(key):
    """(chrom_no, tile_idx) -> genes in that GENE_TILE window (for concurrent prefetch)."""
    chrom_no, tile_idx = key
    ts = tile_idx * GENE_TILE
    return fetch_genes(chrom_no, max(1, ts), ts + GENE_TILE, [])


def _prefetch_gene_tiles(tile_keys):
    """Fetch the needed gene tiles concurrently (Ensembl is the slow part), cache them."""
    todo = [k for k in tile_keys if k not in _GENE_TILE_CACHE]
    if not todo:
        return
    with ThreadPoolExecutor(max_workers=8) as ex:
        for key, genes in zip(todo, ex.map(_fetch_tile, todo)):
            _GENE_TILE_CACHE[key] = genes


def discover_regions(accession, limit=None):
    """accession -> its accessible regions (all by default), labelled by nearest gene.

    Hotspots are the strongest peaks, spaced >= REGION_MIN_SEP apart so windows don't pile
    up. Gene labels come from Ensembl gene tiles (fetched concurrently, cached across
    requests); the strongest hotspots' tiles are fetched first, so if a huge file exceeds
    the fetch budget only the weak tail is coordinate-only.
    """
    prov = resolve_accession(accession)
    peaks = _load_all_peaks(prov["href"], prov["file_accession"])
    peaks.sort(key=lambda p: p["signal"], reverse=True)

    hard_cap = min(limit or MAX_REGIONS, MAX_REGIONS)
    occupied, chosen = {}, []
    for p in peaks:
        c = (p["start"] + p["end"]) // 2
        b = c // REGION_MIN_SEP
        occ = occupied.setdefault(p["chrom"], set())
        if b in occ or (b - 1) in occ or (b + 1) in occ:   # keep hotspots spaced apart
            continue
        occ.add(b)
        chosen.append((p, c))
        if len(chosen) >= hard_cap:
            break

    # distinct gene tiles, strongest-first, up to the label budget -> prefetch concurrently
    want_tiles, seen = [], set()
    for p, c in chosen:
        chrom_no = p["chrom"][3:] if p["chrom"].startswith("chr") else p["chrom"]
        key = (chrom_no, c // GENE_TILE)
        if key not in seen:
            seen.add(key)
            want_tiles.append(key)
        if len(want_tiles) >= GENE_LABEL_BUDGET:
            break
    _prefetch_gene_tiles(want_tiles)

    labelled = 0
    regions = []
    for p, c in chosen:                                     # strongest-first
        chrom = p["chrom"]
        chrom_no = chrom[3:] if chrom.startswith("chr") else chrom
        rs, re = max(0, c - REGION_HALF), c + REGION_HALF
        gene = _nearest_symbol(_GENE_TILE_CACHE.get((chrom_no, c // GENE_TILE)), c)
        if gene:
            labelled += 1
        label = f"{gene or chrom}  ·  {chrom}:{rs:,}-{re:,}  ·  signal {p['signal']:.0f}"
        regions.append({
            "label": label, "region": f"{chrom}:{rs}-{re}",
            "gene": gene, "signal": round(p["signal"], 1),
            "chrom": chrom, "start": rs, "end": re,
        })

    return {
        "biosample": prov.get("biosample"),
        "file_accession": prov.get("file_accession"),
        "assay": prov.get("assay"),
        "total_peaks": len(peaks),
        "region_count": len(regions),
        "gene_labelled": labelled,
        "capped": len(chosen) >= MAX_REGIONS,
        "regions": regions,
    }


# --------------------------------------------------------------------------- #
# Layer-1 tool: fetch_variants  (SPEC.md §6.3) — ClinVar pathogenic/likely in region
# --------------------------------------------------------------------------- #
def fetch_variants(chrom, start, end, warnings):
    # All ClinVar records in the region (not just pathogenic) — significance becomes a
    # ranking term, not a hard filter, so the table is as full as the locus allows.
    term = f"{chrom}[chr] AND {start}:{end}[chrpos38]"
    try:
        es = _get_json(
            f"{EUTILS}/esearch.fcgi?db=clinvar&retmax={MAX_VARIANTS}&retmode=json"
            f"&term={urllib.parse.quote(term)}"
        )
        result = es.get("esearchresult", {})
        ids = result.get("idlist", [])
        total = int(result.get("count", len(ids)))
        if total > len(ids):
            warnings.append(f"{total} ClinVar records in region; showing first {len(ids)} (cap {MAX_VARIANTS}).")
        if not ids:
            return []
        variants = []
        for chunk in _chunks(ids, 100):
            summ = _get_json(
                f"{EUTILS}/esummary.fcgi?db=clinvar&retmode=json&id={','.join(chunk)}"
            )
            res = summ.get("result", {})
            for uid in res.get("uids", []):
                v = _parse_clinvar(res.get(uid, {}), chrom, start, end)
                if v:
                    variants.append(v)
        return variants
    except Exception as e:  # graceful degradation (SPEC §13)
        warnings.append(f"ClinVar query failed ({e}); variant list may be incomplete.")
        return []


def _parse_clinvar(rec, chrom, start, end):
    vset = rec.get("variation_set") or []
    if not vset:
        return None
    vs = vset[0]
    loc = next((l for l in vs.get("variation_loc", []) if l.get("assembly_name") == "GRCh38"), None)
    if not loc or not loc.get("start"):
        return None
    pos = int(loc["start"])
    if loc.get("chr") != str(chrom) or not (start <= pos <= end):
        return None
    rsid = None
    for x in vs.get("variation_xrefs", []):
        if x.get("db_source") == "dbSNP":
            rsid = "rs" + str(x.get("db_id"))
            break
    af, af_source = _best_gnomad_af(vs.get("allele_freq_set") or [])
    clnsig = (
        (rec.get("germline_classification") or {}).get("description")
        or (rec.get("clinical_significance") or {}).get("description")
        or rec.get("clinical_significance_description")
        or "unknown"
    )
    genes = rec.get("genes") or []
    gene = genes[0].get("symbol") if genes else None
    # canonical_spdi ("NC_000001.11:109231759:A:T") is the only place esummary carries
    # ref/alt — needed to join this variant to gnomAD (which keys on chrom-pos-ref-alt).
    ref, alt = _spdi_ref_alt(vs.get("canonical_spdi"))
    return {
        "rsid": rsid,
        "clinvar_id": rec.get("accession"),
        "chrom": "chr" + str(chrom),
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "title": vs.get("variation_name") or rec.get("title"),
        "clnsig": clnsig,
        "gnomad_af": af,
        "gnomad_af_source": af_source,
        "clinvar_gene": gene,
    }


def _spdi_ref_alt(spdi):
    """'NC_000001.11:109231759:A:T' -> ('A', 'T'). Reliable for SNVs; ('', '') otherwise."""
    if not spdi:
        return "", ""
    parts = spdi.split(":")
    if len(parts) != 4:
        return "", ""
    ref, alt = parts[2].strip().upper(), parts[3].strip().upper()
    return ref, alt


def _best_gnomad_af(freqs):
    """Best gnomAD allele frequency from a ClinVar allele_freq_set, preferring genomes.

    Bug this fixes: the raw list interleaves several population sources — TOPMed, ExAC,
    1000G and *two* gnomAD entries ('...(gnomAD)' = genomes, '...(gnomAD), exomes'). Taking
    the first 'gnomAD' match could grab the exomes value, which for the non-coding variants
    this tool cares about is usually absent or unrepresentative. gnomAD genomes is the right
    subset for regulatory DNA, so we prefer it and fall back to exomes only if genomes is
    missing. (gnomAD's own region API can't serve our 150 kb windows in one call — it 502s
    on ~70k variants — so ClinVar's embedded gnomAD frequency is the reliable source here.)
    """
    genomes = exomes = None
    for a in freqs:
        src = (a.get("source") or "")
        if "gnomad" not in src.lower():
            continue
        val = a.get("value")
        if val in (None, ""):
            continue
        try:
            f = float(val)
        except (TypeError, ValueError):
            continue
        if "exome" in src.lower():
            if exomes is None:
                exomes = f
        elif genomes is None:
            genomes = f
    if genomes is not None:
        return genomes, "gnomAD genomes"
    if exomes is not None:
        return exomes, "gnomAD exomes"
    return None, None


def _chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


# --------------------------------------------------------------------------- #
# Layer-1 tool: overlap  (SPEC.md §6.4)
# --------------------------------------------------------------------------- #
def overlap(variants, peaks):
    for v in variants:
        pos = v["pos"]
        hit, nearest = None, None
        best_dist = None
        for p in peaks:
            if p["start"] <= pos < p["end"]:
                hit = p
                break
            d = p["start"] - pos if pos < p["start"] else pos - p["end"]
            if best_dist is None or d < best_dist:
                best_dist, nearest = d, p
        v["in_peak"] = hit is not None
        v["peak_signal"] = (hit or {}).get("signal", 0.0) if hit else 0.0
        v["peak_score"] = (hit or {}).get("score", 0.0) if hit else 0.0
        v["dist_to_peak"] = 0 if hit else (best_dist if best_dist is not None else None)
    return variants


# --------------------------------------------------------------------------- #
# Layer-1 tool: map_to_gene  (SPEC.md §6.5) — nearest gene via Ensembl
# --------------------------------------------------------------------------- #
def fetch_genes(chrom, start, end, warnings):
    try:
        genes = _get_json(
            f"{ENSEMBL}/overlap/region/human/{chrom}:{start}-{end}"
            f"?feature=gene;content-type=application/json"
        )
        return [
            {
                "ensembl_id": g.get("gene_id") or g.get("id"),
                "symbol": g.get("external_name"),
                "start": g.get("start"),
                "end": g.get("end"),
                "biotype": g.get("biotype"),
            }
            for g in genes
            if (g.get("gene_id") or g.get("id"))
        ]
    except Exception as e:
        warnings.append(f"Ensembl gene lookup failed ({e}).")
        return []


def assign_nearest_gene(variant, genes):
    pos = variant["pos"]
    best, best_dist = None, None
    for g in genes:
        if g["start"] <= pos <= g["end"]:
            best, best_dist = g, 0
            break
        d = g["start"] - pos if pos < g["start"] else pos - g["end"]
        if best_dist is None or d < best_dist:
            best, best_dist = g, d
    if best:
        variant["gene_symbol"] = best["symbol"]
        variant["ensembl_id"] = best["ensembl_id"]
        variant["gene_distance"] = best_dist
        variant["link_source"] = "nearest-gene (Ensembl)"
    return variant


# --------------------------------------------------------------------------- #
# Layer-1 tool: opentargets_score  (SPEC.md §6.6)
# --------------------------------------------------------------------------- #
# No disease id from the user: ask Open Targets for the gene's single strongest
# disease association and surface that (disease is inferred from the locus, not typed).
_OT_QUERY = """
query top($ensg: String!) {
  target(ensemblId: $ensg) {
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 1}, orderByScore: "score") {
      rows { score datatypeScores { id score } disease { id name } }
    }
  }
}"""


def opentargets_top_disease(ensembl_id, cache, warnings):
    if not ensembl_id:
        return None
    if ensembl_id in cache:
        return cache[ensembl_id]
    try:
        r = _post_json(OPENTARGETS, {"query": _OT_QUERY, "variables": {"ensg": ensembl_id}})
        rows = (((r.get("data") or {}).get("target") or {}).get("associatedDiseases") or {}).get("rows") or []
        if not rows:
            cache[ensembl_id] = None
            return None
        row = rows[0]
        out = {
            "overall_score": row.get("score"),
            "disease": (row.get("disease") or {}).get("name"),
            "disease_id": (row.get("disease") or {}).get("id"),
            "datatype_scores": {d["id"]: round(d["score"], 3) for d in row.get("datatypeScores", [])},
        }
        cache[ensembl_id] = out
        return out
    except Exception as e:
        warnings.append(f"Open Targets query failed for {ensembl_id} ({e}).")
        cache[ensembl_id] = None
        return None


# --------------------------------------------------------------------------- #
# Layer-1 tool: enrich  (SPEC.md §6.8) — gene-set enrichment over the accessible genes
# The input is the nearest gene of each open-chromatin hotspot, ORDERED strongest-signal
# first, so g:Profiler runs a ranked query: it scans down the list for the cutoff that
# maximises enrichment. An ordered query is a ranked-list statistic (no foreground/
# background split), so the standard annotated-genome domain is the correct choice here —
# the §11 "background = accessible genes" rule applies to over-representation tests, not
# to this ranked scan (a few hundred nearest genes is also too small a custom background
# to be stable — it returns nothing). The ranking itself is the cell-type signal.
# --------------------------------------------------------------------------- #
def enrich(ordered_genes, warnings):
    genes = ordered_genes[:ENRICH_MAX_GENES]
    if not genes:
        return {"terms": [], "genes": [], "sources": ENRICH_SOURCES}
    try:
        r = _post_json(
            GPROFILER,
            {"organism": "hsapiens", "query": genes, "ordered": True,
             "sources": ENRICH_SOURCES, "user_threshold": 0.05,
             "significance_threshold_method": "g_SCS", "no_evidences": True},
            timeout=90,
        )
        res = sorted(r.get("result") or [], key=lambda t: t.get("p_value", 1.0))
        terms = [{
            "source": t.get("source"),
            "id": t.get("native"),
            "name": t.get("name"),
            "p_value": t.get("p_value"),
            "intersection_size": t.get("intersection_size"),
            "term_size": t.get("term_size"),
        } for t in res]
        return {"terms": terms, "genes": genes, "sources": ENRICH_SOURCES}
    except Exception as e:  # graceful degradation (SPEC §13)
        warnings.append(f"g:Profiler enrichment failed ({e}).")
        return {"terms": [], "genes": genes, "sources": ENRICH_SOURCES, "error": str(e)}


def _ordered_accessible_genes(discovery):
    """Ordered, de-duplicated nearest-gene list from region discovery (strongest signal first)."""
    seen, ordered = set(), []
    for r in discovery.get("regions", []):     # discover_regions already sorts by signal desc
        g = r.get("gene")
        if g and g not in seen:
            seen.add(g)
            ordered.append(g)
    return ordered


# --------------------------------------------------------------------------- #
# Layer-1 tool: rank  (SPEC.md §6.9 / §7 composite score — available terms only)
# --------------------------------------------------------------------------- #
_SIG_WEIGHT = [
    ("likely pathogenic", 0.75), ("pathogenic", 1.0),   # order: check 'likely' before 'pathogenic'
    ("likely benign", 0.2), ("benign", 0.1),
    ("risk factor", 0.6), ("association", 0.6), ("drug response", 0.55), ("protective", 0.5),
    ("conflicting", 0.5), ("uncertain", 0.4),
]


def _sig_weight(clnsig):
    s = (clnsig or "").strip().lower()
    for key, w in _SIG_WEIGHT:
        if key in s:
            return w
    return 0.3


def rank(variants, max_signal):
    for v in variants:
        acc = (v["peak_signal"] / max_signal) if (v["in_peak"] and max_signal) else 0.0
        path = _sig_weight(v["clnsig"])
        ot = ((v.get("opentargets") or {}).get("overall_score")) or 0.0
        # equal-weighted over the three terms we can compute standalone (SPEC §7 start-equal).
        v["score_terms"] = {"accessibility": round(acc, 3),
                            "clinvar_pathogenicity": round(path, 3),
                            "opentargets_gene_score": round(ot, 3)}
        v["score"] = round((acc + path + ot) / 3.0, 4)
    variants.sort(key=lambda x: (x["in_peak"], x["score"]), reverse=True)
    return variants


# --------------------------------------------------------------------------- #
# Orchestration (Layer-2 lite): compose the tools + write a templated narration.
# --------------------------------------------------------------------------- #
def analyze(accession, region):
    warnings = []
    chrom, chrom_p, start, end = parse_region(region or DEFAULT_REGION)

    prov = resolve_accession(accession)
    prov["region"] = f"{chrom_p}:{start}-{end}"

    peaks = fetch_peaks(prov["href"], prov["file_accession"], chrom_p, start, end)
    max_signal = max((p["signal"] for p in peaks), default=0.0)

    variants = fetch_variants(chrom, start, end, warnings)
    overlap(variants, peaks)

    genes = fetch_genes(chrom, start, end, warnings)
    ot_cache = {}
    for v in variants:
        assign_nearest_gene(v, genes)
        # disease is inferred from the gene (its strongest Open Targets association), not typed.
        v["opentargets"] = opentargets_top_disease(v.get("ensembl_id"), ot_cache, warnings)

    rank(variants, max_signal)

    stats = {
        "peaks_in_region": len(peaks),
        "variants_in_region": len(variants),
        "variants_in_open_chromatin": sum(1 for v in variants if v["in_peak"]),
        "variants_with_gnomad_af": sum(1 for v in variants if v.get("gnomad_af") is not None),
        "region_bp": end - start,
        "max_peak_signal": round(max_signal, 2),
    }
    return {
        "provenance": prov,
        "stats": stats,
        "peaks": peaks,
        "variants": variants,
        "warnings": warnings,
        "narration": _narrate(prov, stats, variants),
    }


def _narrate(prov, stats, variants):
    """Deterministic, data-driven summary. In the agent version, Claude writes this (SPEC §4.3)."""
    open_v = [v for v in variants if v["in_peak"]]
    lines = [
        f"Accessibility from {prov['file_accession']} "
        f"({prov.get('output_type')}) in {prov.get('biosample')} — "
        f"{prov.get('assay') or 'ATAC/DNase'}, {prov.get('genome')}.",
        f"{stats['peaks_in_region']} accessible peaks and "
        f"{stats['variants_in_region']} ClinVar variants "
        f"in {prov['region']}; {stats['variants_in_open_chromatin']} of those variants "
        f"fall inside open chromatin in this cell type.",
    ]
    if open_v:
        top = open_v[0]
        ot = top.get("opentargets") or {}
        lines.append(
            f"Top candidate: {top.get('rsid') or top.get('clinvar_id')} "
            f"(pos {top['pos']}, {top['clnsig']}) sits in an open peak "
            f"(signal {top['peak_signal']:.1f}) and maps to "
            f"{top.get('gene_symbol') or '?'}"
            + (f", whose strongest disease association is {ot.get('disease')} "
               f"(Open Targets score {ot['overall_score']:.2f})"
               if ot.get("overall_score") is not None else "")
            + ". This is the cell-type-specific mechanistic hypothesis to probe first."
        )
    else:
        lines.append(
            "No pathogenic ClinVar variants fall inside open chromatin here for this cell type — "
            "consistent with this locus being regulatorily inert in this biosample (the contrast the tool is built to show)."
        )
    return " ".join(lines)


# --------------------------------------------------------------------------- #
# HTTP server
# --------------------------------------------------------------------------- #
class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):  # quieter console
        sys.stderr.write("· " + (a[0] % a[1:]) + "\n")

    def _send(self, code, body, ctype="application/json"):
        if isinstance(body, (dict, list)):
            body = json.dumps(body).encode()
        elif isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _file(self, name, ctype):
        try:
            with open(os.path.join(STATIC, name), "rb") as f:
                self._send(200, f.read(), ctype)
        except FileNotFoundError:
            self._send(404, {"error": "not found"})

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        if u.path in ("/", "/index.html"):
            return self._file("index.html", "text/html; charset=utf-8")
        if u.path == "/app.js":
            return self._file("app.js", "application/javascript")
        if u.path == "/style.css":
            return self._file("style.css", "text/css")
        if u.path == "/api/health":
            return self._send(200, {"ok": True})
        if u.path == "/api/regions":
            q = urllib.parse.parse_qs(u.query)
            accession = (q.get("accession") or [""])[0].strip()
            if not accession:
                return self._send(400, {"error": "accession is required"})
            try:
                return self._send(200, discover_regions(accession))
            except Exception as e:
                return self._send(502, {"error": str(e)})
        if u.path == "/api/enrichment":
            q = urllib.parse.parse_qs(u.query)
            accession = (q.get("accession") or [""])[0].strip()
            if not accession:
                return self._send(400, {"error": "accession is required"})
            try:
                disc = discover_regions(accession)
                ordered = _ordered_accessible_genes(disc)
                warnings = []
                out = enrich(ordered, warnings)
                out.update({
                    "biosample": disc.get("biosample"),
                    "file_accession": disc.get("file_accession"),
                    "gene_count": len(ordered),
                    "warnings": warnings,
                })
                return self._send(200, out)
            except Exception as e:
                return self._send(502, {"error": str(e)})
        if u.path == "/api/analyze":
            q = urllib.parse.parse_qs(u.query)
            accession = (q.get("accession") or [""])[0].strip()
            if not accession:
                return self._send(400, {"error": "accession is required"})
            try:
                result = analyze(accession, (q.get("region") or [DEFAULT_REGION])[0])
                return self._send(200, result)
            except Exception as e:
                return self._send(502, {"error": str(e)})
        self._send(404, {"error": "not found"})


def main():
    # Bind 0.0.0.0 so a host platform (Render/Railway/Fly) can route to us; PORT is
    # supplied by the platform. Locally this is still reachable at 127.0.0.1.
    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", "8000"))
    srv = ThreadingHTTPServer((host, port), Handler)
    print(f"regulatr running →  http://{host}:{port}")
    print("  default region", DEFAULT_REGION, "· disease inferred per gene from Open Targets")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")


if __name__ == "__main__":
    main()
