# regulatr — webapp (accession → analysis)

A single-accession vertical slice of the SPEC: enter an ENCODE accession, get the
prioritized non-coding variants that sit in open chromatin **in that cell type**,
mapped to genes and scored against a disease — populated live on the page, with a
loading overlay while it runs and JSON/CSV download.

## Run

```bash
cd app
./run.sh                 # http://127.0.0.1:8765   (or PORT=9000 ./run.sh)
```

No third-party packages — Python 3.8+ stdlib only. All data is fetched live.

## Input

- **Accession** (required): an ATAC-seq experiment `ENCSR…` (the app picks the best
  GRCh38 peak file — conservative IDR > IDR > pseudoreplicated > …) or a peak file
  `ENCFF…` directly. This is the only thing you must supply.
- **Region**: *derived from the accession.* Hit **Load regions** and the dropdown fills
  with **all** the loci where **this** experiment has open chromatin — the peak file is
  scanned genome-wide, hotspots are picked (spaced ≥250 kb apart, strongest first, capped
  at 3,000) and each is labelled by its nearest gene (e.g.
  `LONRF1 · chr8:12,526,526-12,676,526 · signal 43`). Gene labels for the strongest
  hotspots come from Ensembl gene tiles fetched concurrently and cached; the weak tail may
  be coordinate-only. A different cell type surfaces different regions — that contrast is
  the point. A **Custom region…** option remains for any `chrN:start-end`.
- **Disease**: *not entered.* For each target gene the app pulls its single strongest
  Open Targets association and shows that disease + score. Disease is inferred from the
  locus, not typed.

Try `ENCSR042AWH` (HepG2 / liver ATAC) vs `ENCSR000EMT` (GM12878 / B-cell) — Load regions
for each and compare which loci top the list. Deep-linkable:
`/?accession=ENCSR042AWH&region=chr8:12526526-12676526` loads regions then auto-runs.

BED schemas differ across ENCODE files (headerless 10-col narrowPeak vs. density beds
with a `#header`); the peak signal column is resolved per file so ranking never reads a
coordinate column as "signal".

## Pipeline (Layer-1 tools, `server.py`)

`resolve_accession` → `fetch_peaks` → `fetch_variants` (all ClinVar in region, cap 500) →
`overlap` → `map_to_gene` (Ensembl nearest) → `opentargets_top_disease` → `rank` (§7
composite, equal weights over the terms computable standalone: accessibility strength,
ClinVar significance weight, Open Targets gene–disease score). Variants inside open
chromatin sort to the top.

Data sources: ENCODE REST, NCBI ClinVar E-utilities, Ensembl REST, Open Targets GraphQL.
hg38 is the only internal coordinate standard. The first call for a given accession
downloads and caches its peak `bed.gz` under `app/.cache/`; later calls reuse it.

## Pathway enrichment (g:Profiler)

After **Load regions**, a **Run enrichment** panel appears. It takes the nearest gene of
each accessible hotspot, ordered **strongest-signal first**, and runs a g:Profiler
*ordered* query (GO:BP / Reactome / KEGG) — g:Profiler scans down the ranked list for the
cutoff that maximises enrichment, so the accessibility ranking itself is the signal. Results
are shown as an ordered list (most significant first) with source badge, intersection/term
size, p-value, and a −log10(p) bar. A different cell type surfaces different pathways.

Because an ordered query is a ranked-list statistic (not a foreground/background
over-representation test), it uses g:Profiler's standard annotated-genome domain; the
§11 "background = accessible genes" rule targets the over-representation case (and a few
hundred nearest genes is too small a custom background to be stable).

## gnomAD allele frequency

AF comes from ClinVar's embedded population frequencies, preferring **gnomAD genomes** over
**gnomAD exomes** — for the non-coding variants this tool targets, exomes is usually absent
or unrepresentative, so genome frequencies are the right subset. The column shows the value
(hover for the source) and renders ClinVar's rounded-to-zero ultra-rare variants as `<1e-5`
rather than a misleading `0`. gnomAD's own region API can't serve 150 kb windows in one call
(it 502s on ~70k variants), so ClinVar's embedded gnomAD frequency is the reliable source.

## Ranked-variants table

Every ClinVar variant in the region is pulled (cap 1,000) and ranked. The table is
**sortable** (click any column header; click again to reverse) and **paginated** (25 / 50 /
100 / All rows per page, with first/prev/next/last controls). Default order is the rank:
in-open-chromatin variants first, then composite score. Downloads (JSON/CSV) always contain
the full set, not just the visible page.

## Visualization

An embedded **igv.js** genome browser (hg38) with two feature tracks — accessibility
peaks and ClinVar variants (green = in an open peak) — over RefSeq gene models. If
igv.js can't load, the panel falls back to a self-contained static SVG track.

## Not yet built (SPEC features beyond this slice)

Cell-type contrast view (§8), NL cell-type resolution (§4.1),
Borzoi effect prediction (§6.7), user-uploaded BED. The narration is currently a
deterministic template; in the agent version Claude authors it (§4.3). g:Profiler
enrichment (§6.8) is now built (see above).
