#!/usr/bin/env python3

import argparse
import csv
import json
import re
from pathlib import Path


PREFERRED_GENE_ORDER = ["A", "B", "C", "DRB1", "DRB2"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert arcasHLA genotype JSON files into one donor summary CSV."
    )
    parser.add_argument(
        "--input-dir",
        default="genotype",
        help="Folder containing donor*_filtered.genotype.json files.",
    )
    parser.add_argument(
        "--output-csv",
        default="genotype/donor_hla_summary.csv",
        help="Path for the summary CSV file.",
    )
    return parser.parse_args()


def donor_sort_key(path):
    match = re.search(r"donor(\d+)_filtered\.genotype\.json$", path.name)
    if match:
        return int(match.group(1))
    return path.name


def extract_donor_name(path):
    match = re.search(r"(donor\d+)_filtered\.genotype\.json$", path.name)
    if match:
        return match.group(1)
    return path.stem


def normalize_alleles(value):
    if isinstance(value, list):
        return [str(item) for item in value]
    if value in (None, ""):
        return []
    return [str(value)]


def ordered_genes(all_genes):
    ordered = [gene for gene in PREFERRED_GENE_ORDER if gene in all_genes]
    extras = sorted(gene for gene in all_genes if gene not in PREFERRED_GENE_ORDER)
    return ordered + extras


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_csv = Path(args.output_csv)
    genotype_files = sorted(
        input_dir.glob("donor*_filtered.genotype.json"),
        key=donor_sort_key,
    )

    if not genotype_files:
        raise SystemExit(f"No genotype JSON files found in: {input_dir}")

    genotype_data = []
    all_genes = set()

    for json_file in genotype_files:
        with json_file.open() as handle:
            data = json.load(handle)

        genotype_data.append((json_file, data))
        all_genes.update(data.keys())

    genes = ordered_genes(all_genes)
    fieldnames = ["donor"] + [f"{gene}_alleles" for gene in genes] + ["all_alleles"]

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for json_file, data in genotype_data:
            row = {"donor": extract_donor_name(json_file)}
            combined_alleles = []

            for gene in genes:
                alleles = normalize_alleles(data.get(gene, []))
                row[f"{gene}_alleles"] = " | ".join(alleles)
                combined_alleles.extend(alleles)

            row["all_alleles"] = " ; ".join(combined_alleles)
            writer.writerow(row)

    print(f"Wrote {len(genotype_files)} donor summaries to {output_csv}")


if __name__ == "__main__":
    main()
