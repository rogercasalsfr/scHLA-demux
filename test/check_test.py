#!/usr/bin/env python3
"""Check the test run against the published ground truth.

    python3 test/check_test.py        (after the test run; ./test/run_test.sh runs it for you)

The test data (CellBench) is a 10x 3' library. 3' reads cover the end of the HLA
genes where the alleles hardly differ, so clusters often cannot be named from the
HLA table (HLA naming needs 5' data). The test therefore checks:

1. Clustering: souporcell finds the 3 cell lines, each cluster holds one line
   (published demuxlet calls, data/sc_10X.best), and at least 95% of the cells
   called a singlet by both get the same line.
2. Naming is safe: the pipeline reports the 3' chemistry, and no cluster is ever
   named after the wrong cell line (unnamed clusters are fine with 3' data; with
   5' data every cluster must be named).
3. Elimination: the matching step is run again with one cell line's HLA typing
   removed. Still no cluster may get a wrong name (with 5' data, that line must be
   named, as matched_by_elimination).
"""

import csv
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent
REPO = TEST_DIR.parent
sys.path.insert(0, str(REPO / "workflow" / "scripts"))
import chemistry  # noqa: E402

POOL_DIR = TEST_DIR / "results" / "CellBench"
TRUTH = TEST_DIR / "data" / "sc_10X.best"
HLA_TABLE = TEST_DIR / "config" / "hla.csv"
QC_REPORT = TEST_DIR / "results" / "qc_report.html"
MIN_CONCORDANCE = 0.95
MIN_PURITY = 0.95
UNTYPED_FOR_ELIMINATION = "H1975"
NAMED = {"matched", "matched_by_elimination", "manual"}


def read_csv(path, delimiter=","):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def check(ok, text, results):
    results.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {text}")


def check_naming(assignments, truth_names, five_prime, results, untyped=None):
    """Named clusters must carry the right name; with 5' data every cluster must be named."""
    check(set(assignments) == set(truth_names),
          f"the naming covers the souporcell clusters ({', '.join(sorted(assignments))})", results)
    for cluster, row in sorted(assignments.items()):
        expected = truth_names.get(cluster, "?")
        named = row["match_status"] in NAMED
        evidence = f", read fit {row['read_fit']}" if row.get("read_fit") else ""
        text = (f"{cluster}: named {row['assigned_label']} ({row['match_status']}{evidence}), "
                f"ground truth {expected}")
        if five_prime:
            wanted = "matched_by_elimination" if expected == untyped else "matched"
            check(named and row["assigned_label"] == expected and row["match_status"] == wanted,
                  f"{text} (expected {wanted})", results)
        else:
            check(not named or row["assigned_label"] == expected,
                  text + ("" if named else " (not named: fine with 3' data)"), results)


def main():
    for path in (TRUTH, POOL_DIR / "demultiplexed" / "cell_to_donor.csv"):
        if not path.exists():
            raise SystemExit(f"Not found: {path}\nRun ./test/run_test.sh first.")

    # singlet calls for the 3 cell lines (demuxlet was also given lines that are not in this library)
    lines = {row["donor"] for row in read_csv(HLA_TABLE)}
    truth = {
        row["BARCODE"]: row["BEST"][4:]
        for row in read_csv(TRUTH, delimiter="\t")
        if row["BEST"].startswith("SNG-") and row["BEST"][4:] in lines
    }
    cells = read_csv(POOL_DIR / "demultiplexed" / "cell_to_donor.csv")
    assignments = {
        row["cluster_donor"]: row
        for row in read_csv(POOL_DIR / "reports" / "cluster_to_donor_assignments.csv")
    }
    library = chemistry.gene_expression_chemistry(POOL_DIR / "cellranger" / "outs")
    five_prime = chemistry.is_five_prime(library)
    results = []

    print("Clustering (ground truth = published calls of the cluster's cells):")
    by_cluster = {}
    for cell in cells:
        if cell["souporcell_status"] == "singlet":
            line = truth.get(cell["barcode"].split("-")[0])
            if line:
                by_cluster.setdefault(f"donor{cell['souporcell_cluster']}", Counter())[line] += 1
    truth_names = {}
    for cluster, counts in sorted(by_cluster.items()):
        expected, n = counts.most_common(1)[0]
        truth_names[cluster] = expected
        purity = n / sum(counts.values())
        check(purity >= MIN_PURITY,
              f"{cluster}: {n} of {sum(counts.values())} cells are {expected} ({purity:.1%}, need {MIN_PURITY:.0%})",
              results)
    check(len(by_cluster) == 3, f"{len(by_cluster)} clusters with ground-truth cells (expected 3)", results)
    check(len(set(truth_names.values())) == len(truth_names),
          f"each cluster is a different cell line ({', '.join(sorted(truth_names.values()))})", results)

    both = [
        (truth_names.get(f"donor{cell['souporcell_cluster']}"), truth[cell["barcode"].split("-")[0]])
        for cell in cells
        if cell["souporcell_status"] == "singlet" and cell["barcode"].split("-")[0] in truth
    ]
    agree = sum(1 for ours, theirs in both if ours == theirs)
    concordance = agree / len(both) if both else 0
    check(
        concordance >= MIN_CONCORDANCE,
        f"{agree} of {len(both)} cells called a singlet by both are in the cluster of their cell line "
        f"({concordance:.1%}, need {MIN_CONCORDANCE:.0%})",
        results,
    )

    print(f"Naming ({library or 'chemistry unknown'}):")
    if five_prime is False:
        unnamed = [row for row in assignments.values() if row["match_status"] not in NAMED]
        report = QC_REPORT.read_text().replace("&#x27;", "'") if QC_REPORT.exists() else ""
        warned = all("needs 5' data" in row["notes"] for row in unnamed) and "needs 5' data" in report
        check(warned, "the 3' chemistry is reported (QC report, and the notes of the unnamed clusters)",
              results)
    check_naming(assignments, truth_names, five_prime, results)

    print(f"Elimination ({UNTYPED_FOR_ELIMINATION} HLA typing removed, matching run again):")
    with tempfile.TemporaryDirectory() as tmp:
        rows = read_csv(HLA_TABLE)
        for row in rows:
            if row["donor"] == UNTYPED_FOR_ELIMINATION:
                row.update({key: "" for key in row if key[:1] in "ABC" and key[1:2] == "_"})
        table = Path(tmp) / "hla_untyped.csv"
        with table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        # same matching method as the run: read fits when the run made them
        read_fit = POOL_DIR / "reports" / "donor_read_fit.csv"
        by_reads = read_fit.exists()
        summary = "cluster_hla_from_reads.csv" if by_reads else "souporcell_cluster_hla_summary.csv"
        subprocess.run(
            [
                sys.executable, str(REPO / "workflow" / "scripts" / "build_donor_assignments.py"),
                "--genotype-summary", str(POOL_DIR / "reports" / summary),
                "--reference-csv", str(table),
                "--assignment-csv", str(Path(tmp) / "assignments.csv"),
                "--score-csv", str(Path(tmp) / "scores.csv"),
            ] + (["--read-fit-csv", str(read_fit)] if by_reads else []),
            check=True,
            capture_output=True,
        )
        rerun = {row["cluster_donor"]: row for row in read_csv(Path(tmp) / "assignments.csv")}
    check_naming(rerun, truth_names, five_prime, results, untyped=UNTYPED_FOR_ELIMINATION)

    print()
    if all(results):
        print(f"TEST PASSED ({len(results)} checks)")
    else:
        print(f"TEST FAILED ({results.count(False)} of {len(results)} checks failed)")
        sys.exit(1)


if __name__ == "__main__":
    main()
