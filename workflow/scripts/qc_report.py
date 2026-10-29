#!/usr/bin/env python3
"""Write one self-contained HTML QC report for all pools (no dependencies: charts are inline SVG)."""

import argparse
import csv
import glob
import html
import json
from datetime import datetime
from pathlib import Path


# Below these, the result is still written but flagged for a closer look.
MIN_CELLS_PER_DONOR = 200
MAX_DOUBLET_FRACTION = 0.30
MAX_UNASSIGNED_FRACTION = 0.10
MIN_HLA_READS = 1000

SEQUENTIAL = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
              "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument(
        "--pools", nargs="+", required=True, help="POOL=N_DONORS for every pool"
    )
    parser.add_argument("--matching-genes", default="A,B,C")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def esc(value):
    return html.escape(str(value))


def pct(part, total):
    return 100.0 * part / total if total else 0.0


def read_csv(path, delimiter=","):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def to_int(value):
    try:
        return int(float(str(value).replace(",", "")))
    except ValueError:
        return None


# --------------------------------------------------------------------------
# Data of one pool
# --------------------------------------------------------------------------

def cellranger_metrics(pool_dir):
    paths = glob.glob(str(pool_dir / "cellranger" / "outs" / "per_sample_outs" / "*" / "metrics_summary.csv"))
    if not paths:
        return {}
    metrics = {}
    for row in read_csv(paths[0]):
        if row.get("Category") == "Cells" and not row.get("Group Name"):
            metrics[(row["Library Type"], row["Metric Name"])] = row["Metric Value"]
    return metrics


def load_pool(results_dir, pool, n_expected, matching_genes):
    pool_dir = Path(results_dir) / pool
    manifest = read_csv(pool_dir / "demultiplexed" / "manifest.tsv", delimiter="\t")
    scores = read_csv(pool_dir / "reports" / "donor_match_scores.csv")

    calls = {"doublet": 0, "unassigned": 0}
    with (pool_dir / "demultiplexed" / "cell_to_donor.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["souporcell_status"] in calls:
                calls[row["souporcell_status"]] += 1

    donors = []
    for row in manifest:
        genes_json = pool_dir / "arcasHLA" / "genotype" / f"{row['cluster_donor']}_filtered.genes.json"
        hla_reads = {}
        if genes_json.exists():
            data = json.loads(genes_json.read_text())
            hla_reads = {gene: int(data[gene][0]) for gene in matching_genes if gene in data}
        donors.append(
            {
                "label": row["donor"],
                "cluster": row["cluster_donor"].replace("donor", ""),
                "status": row["match_status"],
                "best": row["best_reference_donor"],
                "n_cells": to_int(row["n_cells"]) or 0,
                "n_tcr_cells": to_int(row["n_tcr_cells"]),
                "hla_reads": hla_reads,
                "alleles": {gene: (row.get(f"{gene}_1", ""), row.get(f"{gene}_2", "")) for gene in matching_genes},
            }
        )
    donors.sort(key=lambda donor: int(donor["cluster"]) if donor["cluster"].isdigit() else 0)

    return {
        "pool": pool,
        "n_expected": n_expected,
        "donors": donors,
        "calls": calls,
        "scores": scores,
        "has_tcr": any(donor["n_tcr_cells"] is not None for donor in donors),
        "metrics": cellranger_metrics(pool_dir),
        "pool_dir": pool_dir,
    }


def pool_checks(data, matching_genes):
    """List of (ok, text) answering: did this pool run as it should?"""
    donors, calls = data["donors"], data["calls"]
    singlets = sum(donor["n_cells"] for donor in donors)
    total = singlets + calls["doublet"] + calls["unassigned"]
    checks = []

    n_found = len(donors)
    checks.append((
        n_found == data["n_expected"],
        f"souporcell found {n_found} patients with singlet cells (expected {data['n_expected']})",
    ))

    small = [donor["label"] for donor in donors if donor["n_cells"] < MIN_CELLS_PER_DONOR]
    checks.append((
        not small,
        f"every patient has at least {MIN_CELLS_PER_DONOR} cells"
        if not small else f"few cells (<{MIN_CELLS_PER_DONOR}) for {', '.join(small)}: may be noise, or n_donors too high",
    ))

    doublet_fraction = calls["doublet"] / total if total else 0
    checks.append((
        doublet_fraction <= MAX_DOUBLET_FRACTION,
        f"{pct(calls['doublet'], total):.1f}% doublets (flagged above {MAX_DOUBLET_FRACTION:.0%})",
    ))
    unassigned_fraction = calls["unassigned"] / total if total else 0
    checks.append((
        unassigned_fraction <= MAX_UNASSIGNED_FRACTION,
        f"{pct(calls['unassigned'], total):.1f}% unassigned cells (flagged above {MAX_UNASSIGNED_FRACTION:.0%})",
    ))

    low = [
        f"{donor['label']} ({', '.join(f'{gene} {reads:,}' for gene, reads in donor['hla_reads'].items() if reads < MIN_HLA_READS)})"
        for donor in donors
        if any(reads < MIN_HLA_READS for reads in donor["hla_reads"].values())
    ]
    checks.append((
        not low,
        f"every patient has at least {MIN_HLA_READS:,} HLA reads per gene ({', '.join(matching_genes)})"
        if not low else f"low HLA read support, check these calls: {'; '.join(low)}",
    ))

    if not data["scores"]:
        checks.append((None, "no HLA table given: patients keep the names donor0, donor1, ..."))
    else:
        matched = [donor for donor in donors if donor["status"] == "matched"]
        not_matched = [f"{donor['label']} ({donor['status']})" for donor in donors if donor["status"] != "matched"]
        checks.append((
            not not_matched,
            f"all {len(donors)} patients matched to the HLA table"
            if not not_matched else f"{len(matched)} of {len(donors)} matched; not matched: {', '.join(not_matched)}",
        ))
        table_donors = sorted({row["reference_donor"] for row in data["scores"]})
        unused = [name for name in table_donors if name not in {donor["label"] for donor in matched}]
        if unused:
            checks.append((False, f"in the HLA table but not found in the pool: {', '.join(unused)}"))

    return checks


# --------------------------------------------------------------------------
# Charts (inline SVG)
# --------------------------------------------------------------------------

def bar_chart(rows, value_format):
    """Horizontal bars. rows = [(label, value, muted, tooltip)]."""
    label_width, bar_width, row_height, right = 150, 440, 30, 90
    width = label_width + bar_width + right
    height = row_height * len(rows) + 8
    top = max((value for _, value, _, _ in rows), default=0) or 1

    parts = [f'<svg class="chart" viewBox="0 0 {width} {height}" style="max-width:{width}px" role="img">']
    for index, (label, value, muted, tooltip) in enumerate(rows):
        y = 4 + index * row_height
        length = max(2.0, bar_width * value / top) if value else 0
        fill = "var(--muted-bar)" if muted else "var(--series-1)"
        parts.append(
            f'<g class="mark"><title>{esc(tooltip)}</title>'
            f'<rect x="0" y="{y}" width="{width}" height="{row_height}" fill="transparent"/>'
            f'<text x="{label_width - 10}" y="{y + 19}" text-anchor="end" class="label">{esc(label)}</text>'
            f'<rect x="{label_width}" y="{y + 7}" width="{length:.1f}" height="16" rx="3" fill="{fill}"/>'
            f'<text x="{label_width + length + 8:.1f}" y="{y + 19}" class="value">{esc(value_format(value))}</text></g>'
        )
    parts.append(
        f'<line x1="{label_width}" x2="{label_width}" y1="0" y2="{height}" class="axis"/></svg>'
    )
    return "".join(parts)


def heatmap(data):
    """Clusters (rows) x patients of the HLA table (columns): shared / compared alleles."""
    references = sorted({row["reference_donor"] for row in data["scores"]})
    by_pair = {(row["cluster_donor"], row["reference_donor"]): row for row in data["scores"]}
    assigned = {f"donor{donor['cluster']}": donor["label"] for donor in data["donors"] if donor["status"] == "matched"}

    label_width, cell_w, cell_h, head = 190, 96, 40, 34
    width = label_width + cell_w * len(references) + 4
    height = head + cell_h * len(data["donors"]) + 4
    parts = [f'<svg class="chart" viewBox="0 0 {width} {height}" style="max-width:{width}px" role="img">']
    for col, reference in enumerate(references):
        x = label_width + col * cell_w + cell_w / 2
        parts.append(f'<text x="{x}" y="{head - 12}" text-anchor="middle" class="label">{esc(reference)}</text>')

    for row_index, donor in enumerate(data["donors"]):
        cluster = f"donor{donor['cluster']}"
        y = head + row_index * cell_h
        parts.append(
            f'<text x="{label_width - 10}" y="{y + cell_h / 2 + 5}" text-anchor="end" class="label">'
            f'cluster {esc(donor["cluster"])} → {esc(donor["label"])}</text>'
        )
        for col, reference in enumerate(references):
            score = by_pair.get((cluster, reference))
            if not score:
                continue
            shared, compared = int(score["exact_allele_matches"]), int(score["compared_alleles"])
            fraction = shared / compared if compared else 0
            fill = SEQUENTIAL[round(fraction * (len(SEQUENTIAL) - 1))] if shared else "var(--empty-cell)"
            # ink follows the fill, not the theme (light-blue cells are light in both modes)
            ink = "#ffffff" if fraction > 0.55 else ("#0b0b0b" if shared else "var(--text-primary)")
            is_assigned = assigned.get(cluster) == reference
            x = label_width + col * cell_w
            tooltip = (
                f"cluster {donor['cluster']} vs {reference}: {shared} of {compared} alleles shared"
                + (f" ({score['matched_alleles']})" if score["matched_alleles"] else "")
                + (" - assigned" if is_assigned else "")
            )
            parts.append(
                f'<g class="mark"><title>{esc(tooltip)}</title>'
                f'<rect x="{x + 1}" y="{y + 1}" width="{cell_w - 2}" height="{cell_h - 2}" rx="3" fill="{fill}"'
                + (' class="assigned"' if is_assigned else "")
                + f'/><text x="{x + cell_w / 2}" y="{y + cell_h / 2 + 5}" text-anchor="middle" fill="{ink}" class="cell">'
                f'{"✓ " if is_assigned else ""}{shared}/{compared}</text></g>'
            )
    parts.append("</svg>")
    return "".join(parts)


# --------------------------------------------------------------------------
# Page
# --------------------------------------------------------------------------

def status_badge(ok):
    if ok is None:
        return '<span class="badge info">• Note</span>'
    if ok:
        return '<span class="badge good">✓ OK</span>'
    return '<span class="badge warn">⚠ Check</span>'


def pool_section(data, matching_genes):
    pool, donors, calls = data["pool"], data["donors"], data["calls"]
    checks = pool_checks(data, matching_genes)
    singlets = sum(donor["n_cells"] for donor in donors)
    total = singlets + calls["doublet"] + calls["unassigned"]
    n_warn = sum(1 for ok, _ in checks if ok is False)

    out = [f'<section id="{esc(pool)}"><h2>{esc(pool)}</h2>']
    out.append(
        '<p class="verdict">' + (
            status_badge(True) + " All checks passed" if not n_warn
            else status_badge(False) + f" {n_warn} thing{'s' if n_warn > 1 else ''} to check"
        ) + f' <span class="muted">· {"GEX + TCR" if data["has_tcr"] else "GEX only"}</span></p>'
    )
    out.append('<ul class="checks">' + "".join(
        f"<li>{status_badge(ok)} {esc(text)}</li>" for ok, text in checks
    ) + "</ul>")

    metrics = data["metrics"]
    if metrics:
        items = [
            ("Cells (Cell Ranger)", metrics.get(("Gene Expression", "Cells"))),
            ("Median genes per cell", metrics.get(("Gene Expression", "Median genes per cell"))),
            ("Mean reads per cell", metrics.get(("Gene Expression", "Mean reads per cell"))),
            ("Reads mapped in cells", metrics.get(("Gene Expression", "Confidently mapped reads in cells"))),
            ("TCR cells", metrics.get(("VDJ T", "VDJ cells"))),
            ("Productive TRA+TRB pair", metrics.get(("VDJ T", "Cells with productive V-J spanning pair"))),
        ]
        out.append('<div class="tiles">' + "".join(
            f'<div class="tile"><div class="tile-value">{esc(value)}</div><div class="tile-label">{esc(label)}</div></div>'
            for label, value in items if value
        ) + "</div>")
        out.append(f'<p class="muted small">Full Cell Ranger report: {esc(data["pool_dir"] / "cellranger" / "outs" / "per_sample_outs")}/*/web_summary.html</p>')

    # 1. cells per patient
    rows = [
        (donor["label"], donor["n_cells"], False,
         f"{donor['label']} (cluster {donor['cluster']}): {donor['n_cells']:,} cells, {pct(donor['n_cells'], total):.1f}% of all cells")
        for donor in donors
    ] + [
        ("doublets", calls["doublet"], True, f"{calls['doublet']:,} cells holding two patients' RNA"),
        ("unassigned", calls["unassigned"], True, f"{calls['unassigned']:,} cells souporcell could not place"),
    ]
    out.append(
        f'<figure><h3>{singlets:,} of {total:,} cells assigned to {len(donors)} patients</h3>'
        f'<p class="sub">Singlet cells per patient, from souporcell. Doublets ({pct(calls["doublet"], total):.1f}%) '
        f'and unassigned cells ({pct(calls["unassigned"], total):.1f}%) are left out of the patient folders.</p>'
        + bar_chart(rows, lambda value: f"{value:,}") + "</figure>"
    )

    # 2. HLA matching
    if data["scores"]:
        matched = sum(1 for donor in donors if donor["status"] == "matched")
        out.append(
            f'<figure><h3>{matched} of {len(donors)} clusters matched a patient of the HLA table</h3>'
            f'<p class="sub">Each cell: HLA alleles shared between a souporcell cluster and a patient, out of those compared '
            f'({", ".join(matching_genes)}, two-field). Darker = more shared. ✓ = the patient the cluster was named after.</p>'
            + heatmap(data) + "</figure>"
        )

    # 3. TCR
    if data["has_tcr"]:
        fractions = [pct(donor["n_tcr_cells"] or 0, donor["n_cells"]) for donor in donors]
        tcr_cells = sum(donor["n_tcr_cells"] or 0 for donor in donors)
        rows = [
            (donor["label"], fraction, False,
             f"{donor['label']}: {donor['n_tcr_cells']:,} of {donor['n_cells']:,} cells have a TCR contig")
            for donor, fraction in zip(donors, fractions)
        ]
        out.append(
            f'<figure><h3>TCR found in {pct(tcr_cells, singlets):.0f}% of assigned cells '
            f'({min(fractions):.0f}–{max(fractions):.0f}% per patient)</h3>'
            '<p class="sub">Share of each patient\'s cells with at least one TCR contig (Cell Ranger filtered contigs). '
            'Non-T cells have none, so this also reflects how many T cells each patient has.</p>'
            + bar_chart(rows, lambda value: f"{value:.0f}%") + "</figure>"
        )

    # 4. HLA calls table
    head = "".join(f"<th>HLA-{esc(gene)}</th>" for gene in matching_genes)
    reads_head = "".join(f'<th class="num">{esc(gene)} reads</th>' for gene in matching_genes)
    body = []
    for donor in donors:
        cells = "".join(f"<td>{esc(' / '.join(a for a in donor['alleles'][gene] if a))}</td>" for gene in matching_genes)
        reads = "".join(
            f'<td class="num">{"⚠ " if donor["hla_reads"].get(gene, 0) < MIN_HLA_READS else ""}{donor["hla_reads"].get(gene, 0):,}</td>'
            for gene in matching_genes
        )
        body.append(
            f"<tr><td><b>{esc(donor['label'])}</b></td><td>{esc(donor['cluster'])}</td>"
            f"<td>{esc(donor['status'])}</td>{cells}{reads}</tr>"
        )
    out.append(
        '<figure><h3>HLA alleles called per patient</h3>'
        '<p class="sub">arcasHLA calls on each cluster\'s reads, with the number of reads behind each gene. '
        f'Fewer than {MIN_HLA_READS:,} reads (⚠) makes a call less reliable.</p>'
        f'<div class="table-wrap"><table><thead><tr><th>Patient</th><th>Cluster</th><th>Match</th>{head}{reads_head}</tr></thead>'
        f'<tbody>{"".join(body)}</tbody></table></div></figure>'
    )
    out.append("</section>")
    return "".join(out), n_warn, singlets, total


CSS = """
:root { color-scheme: light;
  --page: #f9f9f7; --surface-1: #fcfcfb; --text-primary: #0b0b0b; --text-secondary: #52514e;
  --muted: #898781; --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --series-1: #2a78d6; --muted-bar: #c3c2b7; --empty-cell: #f0efec;
  --good: #0ca30c; --warn: #fab219; --info: #898781; }
@media (prefers-color-scheme: dark) { :root:where(:not([data-theme="light"])) { color-scheme: dark;
  --page: #0d0d0d; --surface-1: #1a1a19; --text-primary: #f0efec; --text-secondary: #c3c2b7;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --series-1: #3987e5; --muted-bar: #52514e; --empty-cell: #2c2c2a; } }
:root[data-theme="dark"] { color-scheme: dark;
  --page: #0d0d0d; --surface-1: #1a1a19; --text-primary: #f0efec; --text-secondary: #c3c2b7;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --series-1: #3987e5; --muted-bar: #52514e; --empty-cell: #2c2c2a; }
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--text-primary);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 900px; margin: 0 auto; padding: 32px 16px 64px; }
h1 { font-size: 26px; margin: 0 0 4px; } h2 { font-size: 21px; margin: 0 0 8px; }
h3 { font-size: 16px; margin: 0 0 2px; }
section, .overview { background: var(--surface-1); border: 1px solid var(--border); border-radius: 10px;
  padding: 20px; margin-top: 20px; }
figure { margin: 24px 0 0; } .sub, .muted { color: var(--text-secondary); } .sub { margin: 0 0 10px; font-size: 13.5px; }
.small { font-size: 12.5px; word-break: break-all; }
.chart { width: 100%; height: auto; display: block; }
.chart .label { font-size: 13px; fill: var(--text-primary); }
.chart .value { font-size: 12.5px; fill: var(--text-secondary); font-variant-numeric: tabular-nums; }
.chart .cell { font-size: 13px; font-variant-numeric: tabular-nums; }
.chart .axis { stroke: var(--axis); stroke-width: 1; }
.chart .assigned { stroke: var(--text-primary); stroke-width: 2; }
.chart .mark:hover rect:not([fill="transparent"]) { opacity: 0.85; }
.checks { list-style: none; padding: 0; margin: 8px 0 0; }
.checks li { padding: 4px 0; border-top: 1px solid var(--grid); font-size: 14px; }
.badge { display: inline-block; min-width: 68px; font-weight: 600; font-size: 12.5px; }
.badge.good { color: var(--good); } .badge.warn { color: var(--warn); } .badge.info { color: var(--info); }
.verdict { font-size: 15px; margin: 0 0 6px; }
.tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(130px, 1fr)); gap: 10px; margin-top: 16px; }
.tile { border: 1px solid var(--border); border-radius: 8px; padding: 10px; }
.tile-value { font-size: 19px; font-weight: 600; } .tile-label { font-size: 12.5px; color: var(--text-secondary); }
.table-wrap { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { text-align: left; padding: 6px 8px; border-bottom: 1px solid var(--grid); vertical-align: top; }
.num, th { white-space: nowrap; }
th { color: var(--text-secondary); font-weight: 600; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
a { color: var(--series-1); }
"""


def main():
    args = parse_args()
    matching_genes = [gene.strip().upper() for gene in args.matching_genes.split(",") if gene.strip()]

    sections, overview = [], []
    for item in args.pools:
        pool, n_expected = item.rsplit("=", 1)
        data = load_pool(args.results_dir, pool, int(n_expected), matching_genes)
        section, n_warn, singlets, total = pool_section(data, matching_genes)
        sections.append(section)
        matched = sum(1 for donor in data["donors"] if donor["status"] == "matched")
        overview.append(
            f'<tr><td><a href="#{esc(pool)}">{esc(pool)}</a></td>'
            f'<td>{"GEX + TCR" if data["has_tcr"] else "GEX"}</td>'
            f'<td class="num">{len(data["donors"])}</td>'
            f'<td class="num">{(str(matched) if data["scores"] else "–")}</td>'
            f'<td class="num">{singlets:,}</td><td class="num">{pct(singlets, total):.0f}%</td>'
            f'<td>{status_badge(not n_warn)}{"" if not n_warn else f" {n_warn}"}</td></tr>'
        )

    n_ok = sum(1 for row in overview if "✓ OK" in row)
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>scHLA-demux QC</title><style>{CSS}</style></head>
<body><main>
<h1>scHLA-demux QC report</h1>
<p class="muted">{n_ok} of {len(overview)} pools passed every check · generated {datetime.now():%Y-%m-%d %H:%M}</p>
<div class="overview"><h2>Pools</h2><div class="table-wrap"><table>
<thead><tr><th>Pool</th><th>Data</th><th class="num">Patients</th><th class="num">Matched to HLA table</th>
<th class="num">Assigned cells</th><th class="num">% of cells</th><th>Status</th></tr></thead>
<tbody>{"".join(overview)}</tbody></table></div>
<p class="muted small">⚠ Check means a result is unusual and worth a look, not that the run failed. Hover a bar or cell for its exact numbers.</p></div>
{"".join(sections)}
</main></body></html>
"""
    Path(args.output).write_text(page)
    print(f"Wrote QC report to {args.output}")


if __name__ == "__main__":
    main()
