"use strict";

const $ = (id) => document.getElementById(id);
let lastResult = null;
let igvBrowser = null;

// The region dropdown is derived from the accession (top accessible loci in that
// cell type), fetched from /api/regions. Until then it holds only placeholders.
const CUSTOM = "__custom__";
let regionsLoaded = false;

initControls();

function setSelectPlaceholder(text) {
  const sel = $("region-select");
  sel.innerHTML =
    `<option value="" disabled selected>${esc(text)}</option>` +
    `<option value="${CUSTOM}">Custom region…</option>`;
  sel.disabled = false;
}

function populateRegions(regions) {
  const sel = $("region-select");
  const opts = regions.map(
    (r) => `<option value="${esc(r.region)}">${esc(r.label)}</option>`
  );
  opts.push(`<option value="${CUSTOM}">Custom region…</option>`);
  sel.innerHTML = opts.join("");
  sel.disabled = false;
  if (regions.length) sel.selectedIndex = 0; // preselect the strongest hotspot
  syncCustomField();
  regionsLoaded = true;
  updateRunEnabled();
  // enrichment runs on the accession's accessible genes — offer it once regions are loaded
  $("enrich-panel").classList.remove("hidden");
  $("enrichbtn").disabled = false;
  $("enrich-list").innerHTML = "";
  setEnrichStatus("Which pathways do this cell type's most-accessible genes converge on? Run it to find out.");
}

function syncCustomField() {
  $("custom-region-field").classList.toggle("hidden", $("region-select").value !== CUSTOM);
}

function updateRunEnabled() {
  $("run").disabled = !currentRegion();
}

async function loadRegions() {
  const accession = $("accession").value.trim();
  if (!accession) {
    setRegionStatus("Enter an ENCODE accession first.", true);
    return false;
  }
  regionsLoaded = false;
  $("run").disabled = true;
  $("loadregions").disabled = true;
  $("region-select").disabled = true;
  $("enrich-panel").classList.add("hidden"); // hide stale enrichment for the previous accession
  setSelectPlaceholder("finding accessible regions…");
  setRegionStatus("⏳ Scanning " + accession + " genome-wide for its strongest open chromatin — first load downloads the peak file (~15s).");
  try {
    const res = await fetch("/api/regions?" + new URLSearchParams({ accession }));
    const data = await res.json();
    if (!res.ok || data.error) throw new Error(data.error || `HTTP ${res.status}`);
    if (!data.regions || !data.regions.length) {
      setSelectPlaceholder("no accessible peaks found — use Custom region…");
      setRegionStatus("No standard-chromosome peaks found for " + accession + ".", true);
      return false;
    }
    populateRegions(data.regions);
    const top = data.regions[0];
    const capped = data.capped ? ` (capped at ${data.regions.length.toLocaleString()})` : "";
    setRegionStatus(
      `${data.regions.length.toLocaleString()} accessible regions${capped} in ` +
      `${data.biosample || accession} (${(data.total_peaks || 0).toLocaleString()} peaks genome-wide), ` +
      `sorted by accessibility. Strongest: ${top.gene || top.chrom} (signal ${top.signal}). ` +
      `Pick one and hit Analyze.`
    );
    return true;
  } catch (err) {
    setSelectPlaceholder("region load failed — use Custom region…");
    setRegionStatus("Could not load regions: " + err.message, true);
    return false;
  } finally {
    $("loadregions").disabled = false;
  }
}

function setRegionStatus(msg, isError) {
  const el = $("region-status");
  el.textContent = msg;
  el.classList.toggle("err", !!isError);
}

// ---- pathway enrichment (g:Profiler ordered query over top-signal nearest genes) ----
async function runEnrichment() {
  const accession = $("accession").value.trim();
  if (!accession) return;
  $("enrichbtn").disabled = true;
  $("enrich-list").innerHTML = "";
  setEnrichStatus("⏳ Ranking this cell type's most-accessible genes and testing them against GO / Reactome / KEGG (ordered query)…");
  try {
    const res = await fetch("/api/enrichment?" + new URLSearchParams({ accession }));
    const data = await res.json();
    if (!res.ok || data.error) throw new Error(data.error || `HTTP ${res.status}`);
    renderEnrichment(data);
  } catch (err) {
    setEnrichStatus("Enrichment failed: " + err.message, true);
  } finally {
    $("enrichbtn").disabled = false;
  }
}

function setEnrichStatus(msg, isError) {
  const el = $("enrich-status");
  el.textContent = msg;
  el.classList.toggle("err", !!isError);
}

function renderEnrichment(d) {
  const terms = d.terms || [];
  const who = d.biosample || "this cell type";
  if (!terms.length) {
    $("enrich-list").innerHTML = "";
    setEnrichStatus(
      `No pathways pass the g:SCS threshold among the ${d.gene_count || 0} most-accessible genes in ` +
      `${who} — the strongest peaks may sit at generic promoters rather than a coherent gene set.`);
    return;
  }
  setEnrichStatus(
    `${terms.length} enriched term${terms.length > 1 ? "s" : ""} (GO / Reactome / KEGG) among the ` +
    `${d.gene_count} most-accessible genes in ${who}, ordered strongest-signal first. Showing the top ` +
    `${Math.min(25, terms.length)} by significance.`);
  const maxNL = Math.max(1, ...terms.map((t) => -Math.log10(t.p_value || 1)));
  $("enrich-list").innerHTML = terms.slice(0, 25).map((t, i) => {
    const nl = -Math.log10(t.p_value || 1);
    const w = Math.max(4, Math.round((nl / maxNL) * 100));
    const cls = String(t.source || "").replace(/[^a-z]/gi, "");
    return `<li>
      <div class="erow">
        <span class="rank">${i + 1}</span>
        <span class="src-badge ${cls}">${esc(t.source)}</span>
        <span class="term-name">${esc(t.name)}</span>
        <span class="term-meta">${t.intersection_size}/${t.term_size} · p=${fmtP(t.p_value)}</span>
      </div>
      <div class="term-bar-wrap"><span class="term-bar" style="width:${w}%"></span></div>
    </li>`;
  }).join("");
}

function fmtP(p) {
  if (p == null) return "—";
  return p < 1e-3 ? p.toExponential(1) : p.toFixed(3);
}

function initControls() {
  setSelectPlaceholder("load an accession to list its accessible regions");

  const sel = $("region-select");
  sel.addEventListener("change", () => { syncCustomField(); updateRunEnabled(); });
  $("region-custom").addEventListener("input", updateRunEnabled);

  $("loadregions").addEventListener("click", loadRegions);
  $("enrichbtn").addEventListener("click", runEnrichment);
  $("form").addEventListener("submit", (e) => { e.preventDefault(); run(); });
  document.querySelectorAll(".chip").forEach((c) =>
    c.addEventListener("click", () => {
      $("accession").value = c.dataset.acc;
      loadRegions();
    })
  );
  $("dljson").addEventListener("click", () => download("json"));
  $("dlcsv").addEventListener("click", () => download("csv"));

  // ranked-variants table: page size + pager
  $("pagesize").addEventListener("change", () => {
    const v = $("pagesize").value;
    tableState.pageSize = v === "all" ? "all" : parseInt(v, 10);
    tableState.page = 0;
    drawTable();
  });
  $("firstpage").addEventListener("click", () => gotoPage(0));
  $("prevpage").addEventListener("click", () => gotoPage(tableState.page - 1));
  $("nextpage").addEventListener("click", () => gotoPage(tableState.page + 1));
  $("lastpage").addEventListener("click", () => gotoPage(Infinity));

  // deep link: /?accession=ENCSR042AWH&region=chr8:...  loads regions then runs
  const qs = new URLSearchParams(location.search);
  if (qs.get("accession")) {
    $("accession").value = qs.get("accession").trim();
    loadRegions().then((ok) => {
      const r = qs.get("region");
      if (r) {
        const opt = [...sel.options].find((o) => o.value === r);
        if (opt) sel.value = r;
        else { sel.value = CUSTOM; $("region-custom").value = r; }
        syncCustomField();
        updateRunEnabled();
      }
      if (ok || (r && currentRegion())) run();
    });
  }
}

function currentRegion() {
  const sel = $("region-select");
  return sel.value === CUSTOM ? $("region-custom").value.trim() : sel.value;
}

async function run() {
  const accession = $("accession").value.trim();
  if (!accession) return;
  const region = currentRegion();
  if (!region) {
    showError("Enter a custom region or pick one from the list.");
    return;
  }

  $("results").classList.remove("hidden");
  $("errbox").classList.add("hidden");
  $("overlay").classList.remove("hidden"); // buffer icon over table + figure
  $("run").disabled = true;
  $("dljson").disabled = $("dlcsv").disabled = true;

  try {
    const q = new URLSearchParams({ accession, region });
    const res = await fetch("/api/analyze?" + q.toString());
    const data = await res.json();
    if (!res.ok || data.error) throw new Error(data.error || `HTTP ${res.status}`);
    lastResult = data;
    await render(data); // await so the spinner stays until igv is ready
    $("dljson").disabled = $("dlcsv").disabled = false;
  } catch (err) {
    showError("Analysis failed: " + err.message);
  } finally {
    $("overlay").classList.add("hidden");
    $("run").disabled = false;
  }
}

function showError(msg) {
  const e = $("errbox");
  e.textContent = msg;
  e.classList.remove("hidden");
}

async function render(d) {
  renderProvenance(d.provenance);
  renderStats(d.stats);
  $("narration").textContent = d.narration || "";
  renderTable(d.variants);
  $("warnings").textContent = (d.warnings || []).join("  ·  ");
  await renderTrack(d);
}

function renderProvenance(p) {
  const rows = [
    ["Cell type", `<b>${esc(p.biosample)}</b>`],
    ["Accession used", `${esc(p.file_accession)} <small>(${esc(p.output_type || "")})</small>`],
    ["Experiment", esc(p.experiment || "—")],
    ["Assay / build", `${esc(p.assay || "—")} · ${esc(p.genome)}`],
    ["Region", esc(p.region)],
    ["Lab", esc(p.lab || "—")],
  ];
  $("prov").innerHTML = rows
    .map(([k, v]) => `<div class="k">${k}</div><div class="v">${v}</div>`)
    .join("");
}

function renderStats(s) {
  const cells = [
    ["Peaks in region", s.peaks_in_region, false],
    ["ClinVar variants", s.variants_in_region, false],
    ["In open chromatin", s.variants_in_open_chromatin, true],
    ["Max peak signal", s.max_peak_signal, false],
    ["Region size (bp)", (s.region_bp || 0).toLocaleString(), false],
  ];
  $("stats").innerHTML = cells
    .map(([lbl, num, hi]) =>
      `<div class="stat ${hi ? "hi" : ""}"><div class="num">${num}</div><div class="lbl">${lbl}</div></div>`)
    .join("");
}

// ---- sortable + paginated ranked-variants table ----
// Each column: how to sort it (num/str) and how to render its cell.
const COLS = [
  { key: "variant", label: "Variant", type: "str",
    sort: (v) => (v.rsid || v.clinvar_id || "").toLowerCase(),
    cell: (v) => esc(v.rsid || v.clinvar_id || "—") +
      (v.ref && v.alt ? ` <small class="allele">${esc(v.ref)}&gt;${esc(v.alt)}</small>` : "") },
  { key: "pos", label: "Position", type: "num",
    sort: (v) => v.pos || 0,
    cell: (v) => `${v.chrom}:${(v.pos || 0).toLocaleString()}` },
  { key: "clnsig", label: "ClinVar", type: "str",
    sort: (v) => (v.clnsig || "").toLowerCase(),
    cell: (v) => esc(v.clnsig || "—") },
  { key: "af", label: "gnomAD AF", type: "num",
    sort: (v) => (v.gnomad_af == null ? -1 : v.gnomad_af),
    cell: (v) => {
      if (v.gnomad_af == null) return "—";
      // ClinVar rounds AF to 5 decimals, so ultra-rare variants (true AF < 5e-6) arrive as 0.
      const txt = v.gnomad_af > 0 ? v.gnomad_af.toExponential(1) : "&lt;1e-5";
      const src = v.gnomad_af_source ? ` title="${esc(v.gnomad_af_source)}"` : "";
      return `<span${src}>${txt}</span>`;
    } },
  { key: "inpeak", label: "In open peak", type: "num",
    sort: (v) => (v.in_peak ? 1 : 0),
    cell: (v) => `<span class="badge ${v.in_peak ? "yes" : "no"}">${v.in_peak ? "yes" : "no"}</span>` },
  { key: "signal", label: "Peak signal", type: "num",
    sort: (v) => (v.in_peak ? (v.peak_signal || 0) : -1),
    cell: (v) => (v.in_peak ? (v.peak_signal || 0).toFixed(1) : "—") },
  { key: "gene", label: "Target gene", type: "str",
    sort: (v) => (v.gene_symbol || "").toLowerCase(),
    cell: (v) => esc(v.gene_symbol || "—") },
  { key: "disease", label: "Top disease (Open Targets)", type: "num",
    sort: (v) => ((v.opentargets || {}).overall_score ?? -1),
    cell: (v) => {
      const ot = v.opentargets || {};
      return ot.overall_score != null
        ? `${esc(ot.disease || "")} <b>${ot.overall_score.toFixed(2)}</b>` : "—";
    } },
  { key: "score", label: "Score", type: "num",
    sort: (v) => v.score || 0,
    cell: (v, maxScore) => {
      const barW = Math.round(((v.score || 0) / maxScore) * 60);
      return `<span class="scorebar" style="width:${barW}px"></span>${(v.score || 0).toFixed(3)}`;
    } },
];

const tableState = { rows: [], sortKey: null, sortDir: -1, page: 0, pageSize: 25 };

function renderTable(vs) {
  tableState.rows = vs || [];
  tableState.page = 0;
  tableState.sortKey = null; // start in backend rank order (in-peak first, then score)
  buildTableHead();
  drawTable();
}

function buildTableHead() {
  $("vhead").innerHTML =
    "<tr>" +
    COLS.map((c) => {
      const active = tableState.sortKey === c.key;
      const arrow = active ? (tableState.sortDir < 0 ? " ▼" : " ▲") : "";
      return `<th data-key="${c.key}" class="${active ? "sorted" : ""}">${esc(c.label)}<span class="arrow">${arrow}</span></th>`;
    }).join("") +
    "</tr>";
  $("vhead").querySelectorAll("th").forEach((th) =>
    th.addEventListener("click", () => sortBy(th.dataset.key))
  );
}

function sortBy(key) {
  if (tableState.sortKey === key) {
    tableState.sortDir = -tableState.sortDir;
  } else {
    tableState.sortKey = key;
    tableState.sortDir = -1; // first click: descending (biggest / Z-A first)
  }
  tableState.page = 0;
  buildTableHead();
  drawTable();
}

function sortedRows() {
  if (!tableState.sortKey) return tableState.rows;
  const col = COLS.find((c) => c.key === tableState.sortKey);
  const dir = tableState.sortDir;
  return [...tableState.rows].sort((a, b) => {
    const x = col.sort(a), y = col.sort(b);
    if (x < y) return -dir;
    if (x > y) return dir;
    return 0;
  });
}

function pageCount(total) {
  if (tableState.pageSize === "all") return 1;
  return Math.max(1, Math.ceil(total / tableState.pageSize));
}

function drawTable() {
  const rows = sortedRows();
  const total = rows.length;
  const pages = pageCount(total);
  if (tableState.page >= pages) tableState.page = pages - 1;

  const size = tableState.pageSize === "all" ? total : tableState.pageSize;
  const start = tableState.pageSize === "all" ? 0 : tableState.page * size;
  const slice = rows.slice(start, start + (size || total));

  const maxScore = Math.max(0.001, ...tableState.rows.map((v) => v.score || 0));
  $("vbody").innerHTML = slice
    .map((v) => `<tr class="${v.in_peak ? "inpeak" : ""}">` +
      COLS.map((c) => `<td>${c.cell(v, maxScore)}</td>`).join("") + "</tr>")
    .join("");

  const shownFrom = total ? start + 1 : 0;
  const shownTo = Math.min(start + slice.length, total);
  $("table-count").textContent = total
    ? `${shownFrom.toLocaleString()}–${shownTo.toLocaleString()} of ${total.toLocaleString()} variants`
    : "no variants";
  $("pageinfo").textContent = `Page ${tableState.page + 1} / ${pages}`;
  const atFirst = tableState.page === 0, atLast = tableState.page >= pages - 1;
  $("firstpage").disabled = $("prevpage").disabled = atFirst;
  $("nextpage").disabled = $("lastpage").disabled = atLast;
}

function gotoPage(p) {
  const pages = pageCount(sortedRows().length);
  tableState.page = Math.max(0, Math.min(p, pages - 1));
  drawTable();
}

// ---- igv.js genome browser, with a self-contained SVG fallback ----
async function renderTrack(d) {
  const div = $("track");
  const peaks = d.peaks.map((p) => ({
    chr: p.chrom, start: p.start, end: p.end,
    name: `peak · signal ${(p.signal || 0).toFixed(1)}`,
  }));
  const vars = d.variants.map((v) => ({
    chr: v.chrom, start: Math.max(0, v.pos - 1), end: v.pos,
    name: `${v.rsid || v.clinvar_id || ""} · ${v.clnsig || ""}${v.in_peak ? " · in open peak" : ""}`,
    color: v.in_peak ? "#37d39b" : "#9aa3b2",
  }));

  try {
    if (typeof igv === "undefined") throw new Error("igv unavailable");
    igv.removeAllBrowsers();
    div.innerHTML = "";
    const create = igv.createBrowser(div, {
      genome: "hg38",
      locus: d.provenance.region,
      tracks: [
        { name: `Accessible peaks — ${d.provenance.biosample || ""}`, type: "annotation",
          displayMode: "COLLAPSED", color: "#3a5f8a", height: 60, features: peaks },
        { name: "ClinVar variants (green = in open peak)", type: "annotation",
          displayMode: "COLLAPSED", height: 70, features: vars },
      ],
    });
    igvBrowser = await withTimeout(create, 20000);
  } catch (e) {
    renderTrackSVG(d, div); // fallback so the panel is never blank
  }
}

function withTimeout(promise, ms) {
  return Promise.race([
    promise,
    new Promise((_, rej) => setTimeout(() => rej(new Error("igv load timed out")), ms)),
  ]);
}

function renderTrackSVG(d, div) {
  const [, span] = d.provenance.region.split(":");
  const [rs, re] = span.split("-").map(Number);
  const W = Math.max(820, d.peaks.length * 3), H = 150, padL = 8, padR = 8, baseY = 108, topY = 34, axisY = 120;
  const x = (pos) => padL + ((pos - rs) / (re - rs)) * (W - padL - padR);
  const maxSig = Math.max(1, ...d.peaks.map((p) => p.signal || 0));
  let svg = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}">`;
  svg += `<line x1="${padL}" y1="${baseY}" x2="${W - padR}" y2="${baseY}" stroke="#2a2f3a"/>`;
  for (let i = 0; i <= 4; i++) {
    const pos = rs + ((re - rs) * i) / 4, xx = x(pos);
    svg += `<text x="${xx}" y="${axisY + 10}" text-anchor="middle">${(pos / 1e6).toFixed(3)} Mb</text>`;
  }
  for (const p of d.peaks) {
    const x1 = x(p.start), x2 = Math.max(x(p.end), x1 + 1.5);
    const h = 6 + ((p.signal || 0) / maxSig) * (baseY - topY - 6);
    svg += `<rect x="${x1}" y="${baseY - h}" width="${x2 - x1}" height="${h}" rx="1.5" fill="#3a5f8a" opacity="0.85"/>`;
  }
  for (const v of d.variants) {
    const xx = x(v.pos), color = v.in_peak ? "#37d39b" : "#6b7280";
    svg += `<line x1="${xx}" y1="${topY - 6}" x2="${xx}" y2="${baseY}" stroke="${color}" stroke-width="1.2" opacity="0.85"/>`;
    svg += `<circle cx="${xx}" cy="${topY - 8}" r="3.5" fill="${color}"/>`;
  }
  svg += `</svg>`;
  div.innerHTML = `<p class="warnings">igv.js unavailable — showing static fallback track.</p>` + svg;
}

// ---- downloads ----
function download(kind) {
  if (!lastResult) return;
  let blob, name;
  if (kind === "json") {
    blob = new Blob([JSON.stringify(lastResult, null, 2)], { type: "application/json" });
    name = `regulatr_${lastResult.provenance.file_accession}.json`;
  } else {
    blob = new Blob([toCSV(lastResult.variants)], { type: "text/csv" });
    name = `regulatr_${lastResult.provenance.file_accession}.csv`;
  }
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = name; a.click();
  URL.revokeObjectURL(url);
}

function toCSV(vs) {
  const cols = ["rsid", "clinvar_id", "chrom", "pos", "ref", "alt", "clnsig",
    "gnomad_af", "gnomad_af_source", "in_peak", "peak_signal", "gene_symbol",
    "ensembl_id", "score"];
  const rows = vs.map((v) => {
    const ot = v.opentargets || {};
    return [...cols.map((c) => csv(v[c])), csv(ot.disease), csv(ot.overall_score)].join(",");
  });
  return [cols.join(",") + ",top_disease,opentargets_score", ...rows].join("\n");
}

function csv(x) {
  if (x == null) return "";
  const s = String(x);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

function esc(s) {
  return String(s == null ? "" : s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}
