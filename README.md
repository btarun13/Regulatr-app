# Regulatr

**Cell-type-aware non-coding variant prioritization.**

Give it a single ENCODE chromatin-accessibility experiment and it tells you which
disease-associated non-coding variants land in open regulatory DNA *in that specific cell
type* — then maps them to their likely target genes and diseases and ranks them as
mechanistic hypotheses.

Python standard library only (no third-party packages). All data is fetched live from
public genomics resources at query time.

![Demonstration of the project's videos](https://github.com/btarun13/Regulatr-app/blob/main/regulatr_teaser.mp4)

---

## Why

Most disease-associated genetic variation is **non-coding**, where its effect — if any — is
*regulatory*: it changes whether and where a gene is switched on. But regulatory DNA is
**cell-type-specific**: an enhancer variant only matters in a cell type where that enhancer
is active. So the useful question isn't "is this variant functional?" but:

> **"Is this variant in open, active regulatory DNA *in the relevant cell type*, and which
> gene does it plausibly regulate?"**

regulatr uses **chromatin accessibility** (ATAC-seq / DNase peaks from ENCODE) as the readout
of where regulatory DNA is open in a given biosample, and intersects that with known variants
from ClinVar. Swap the accession for a different cell type and different loci rise to the top
— that contrast is the point.

**Scope guard:** accessibility tells us *where* regulatory DNA is open; it never *calls*
variants. Variants come from ClinVar. The single internal coordinate standard is **hg38 /
GRCh38**, and build mismatches are rejected rather than silently lifted.

---

## Features

- **Accession resolution** — enter an ENCODE experiment (`ENCSR…`) or peak file (`ENCFF…`);
  the tool picks the best GRCh38 peak file and records full provenance (assay, biosample, lab).
- **Region discovery** — scans the peak file genome-wide and surfaces the loci where *this*
  experiment has the strongest open chromatin, each labelled by its nearest gene.
- **Variant retrieval** — every ClinVar record in the region, with clinical significance,
  rsID, gnomAD allele frequency, and ref/alt.
- **Accessibility overlap** — marks which variants fall inside an open-chromatin peak in this
  cell type, and the peak's signal.
- **Gene mapping** — assigns each variant to its nearest gene via Ensembl.
- **Disease inference** — each target gene's single strongest Open Targets association
  (disease is inferred from the locus, never typed).
- **Composite ranking** — accessibility strength + ClinVar significance + Open Targets gene
  score, equally weighted; in-peak variants sort to the top.
- **Pathway enrichment** — an ordered g:Profiler query (GO:BP / Reactome / KEGG) over the
  cell type's most-accessible genes.
- **Genome-browser view** — embedded igv.js (hg38) with a self-contained SVG fallback.
- **Export** — full result set as JSON or CSV.

---

## Quick start

```bash
cd app
./run.sh                 # http://127.0.0.1:8765   (or PORT=9000 ./run.sh)
```

No third-party packages — Python 3.8+ stdlib only.

Try `ENCSR637OPZ` (CD8 T cell) or `ENCSR000EMT` (GM12878 / B-cell): hit **Load regions** for
each and compare which loci top the list. Deep-linkable:
`/?accession=ENCSR637OPZ&region=chr1:109200000-109350000` loads regions, then auto-runs.

---

## How it works

For a single accession + region the backend composes seven deterministic tools:

```
resolve_accession → fetch_peaks → fetch_variants → overlap
                 → map_to_gene → opentargets_score → rank
```

Each variant's composite score is the mean of three terms in `[0, 1]`: **accessibility**
(peak signal normalized to the strongest peak in the region), **clinvar_pathogenicity** (a
weight from ClinVar significance), and **opentargets_gene_score** (the target gene's strongest
association). In-peak status is the primary sort key, so cell-type accessibility dominates
ordering.

The first call for a given accession downloads and caches its peak `bed.gz` under
`app/.cache/`; later calls reuse it.

---

## Project structure

```
app/
  server.py          # stdlib HTTP server + the Layer-1 tools and pipeline
  run.sh             # launch script
  static/            # index.html, app.js, style.css (frontend)
  README.md          # detailed backend/pipeline notes
LICENSE              # MIT
```

## Data sources

ENCODE REST · NCBI ClinVar E-utilities · Ensembl REST · Open Targets GraphQL · g:Profiler.
Genome build: hg38 / GRCh38 throughout.

---

## Limitations

- **Accessibility ≠ causality.** A variant in an open peak is a *hypothesis*, not a proven
  regulatory effect.
- **Nearest-gene mapping is a heuristic** — not chromatin-contact evidence; enhancers can
  regulate distal genes.
- **Coverage is bounded by the inputs** — ClinVar is capped per region, and the gene-label
  pass during region discovery has a fetch budget, so the weakest hotspots may be
  coordinate-only.

---

## License

[MIT](LICENSE) © 2026 Tarun Naithani
