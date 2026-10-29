#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


DONOR_COLUMNS = [
    "pool",
    "donor",
    "match_status",
    "n_cells",
    "n_tcr_cells",
    "cluster_donor",
    "best_reference_donor",
    "exact_allele_matches",
    "compared_alleles",
    "matched_alleles",
    "notes",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine the per-pool donor manifests and cell tables into one file each."
    )
    parser.add_argument("--manifests", nargs="+", required=True)
    parser.add_argument("--cell-tables", nargs="+", required=True)
    parser.add_argument("--donor-summary", required=True)
    parser.add_argument("--cell-summary", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    rows = []
    hla_columns = []
    for manifest in args.manifests:
        with Path(manifest).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for column in reader.fieldnames:
                if column[-2:] in ("_1", "_2") and column not in hla_columns:
                    hla_columns.append(column)
            rows.extend(reader)

    fieldnames = DONOR_COLUMNS + hla_columns + ["output_dir"]
    with Path(args.donor_summary).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    with Path(args.cell_summary).open("w", newline="") as out:
        header_written = False
        for cell_table in args.cell_tables:
            with Path(cell_table).open() as handle:
                header = handle.readline()
                if not header_written:
                    out.write(header)
                    header_written = True
                out.writelines(handle)

    print(f"Wrote {len(rows)} donors from {len(args.manifests)} pools to {args.donor_summary}")
    print(f"Wrote all cells to {args.cell_summary}")


if __name__ == "__main__":
    main()
