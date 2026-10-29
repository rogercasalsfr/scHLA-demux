#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path

from hla_table import ordered_genes, write_hla_table


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert arcasHLA genotype JSON files into one HLA table "
            "(donor,A_1,A_2,B_1,B_2,...)."
        )
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
    donors = {}

    for json_file, data in genotype_data:
        profile = {}
        for gene in genes:
            alleles = normalize_alleles(data.get(gene, []))
            if len(alleles) == 1:
                # one reported allele = homozygous, same convention as the input table
                alleles = alleles * 2
            profile[gene] = alleles[:2]
        donors[extract_donor_name(json_file)] = profile

    write_hla_table(output_csv, genes, donors)

    print(f"Wrote {len(genotype_files)} donor summaries to {output_csv}")


if __name__ == "__main__":
    main()
