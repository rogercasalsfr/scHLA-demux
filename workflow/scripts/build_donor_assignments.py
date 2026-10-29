#!/usr/bin/env python3

import argparse
import csv
from collections import Counter
from pathlib import Path

import hla_table


ASSIGNMENT_COLUMNS = [
    "cluster_donor",
    "assigned_label",
    "assigned_reference_donor",
    "best_reference_donor",
    "match_status",
    "exact_allele_matches",
    "compared_alleles",
    "overlap_fraction",
    "shared_gene_count",
    "matched_genes",
    "matched_alleles",
    "tied_reference_donors",
    "notes",
]

SCORE_COLUMNS = [
    "cluster_donor",
    "reference_donor",
    "exact_allele_matches",
    "compared_alleles",
    "overlap_fraction",
    "shared_gene_count",
    "matched_genes",
    "matched_alleles",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Assign souporcell donor clusters to reference donor HLA annotations."
        )
    )
    parser.add_argument("--genotype-summary", required=True)
    parser.add_argument("--reference-csv")
    parser.add_argument("--assignment-csv", required=True)
    parser.add_argument("--score-csv", required=True)
    parser.add_argument(
        "--genes",
        default="A,B,C",
        help="Comma-separated list of HLA genes to use for matching.",
    )
    parser.add_argument(
        "--two-field-resolution",
        dest="two_field_resolution",
        action="store_true",
        help="Compare alleles at two-field resolution.",
    )
    parser.add_argument(
        "--no-two-field-resolution",
        dest="two_field_resolution",
        action="store_false",
        help="Compare alleles using the full allele string.",
    )
    parser.set_defaults(two_field_resolution=True)
    return parser.parse_args()


def write_csv(path, fieldnames, rows):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_gene_list(value):
    return [gene.strip() for gene in str(value).split(",") if gene.strip()]


def normalize_allele(allele, two_field_resolution):
    allele = allele.strip()
    if not allele or not two_field_resolution or "*" not in allele:
        return allele

    gene, remainder = allele.split("*", 1)
    fields = [field for field in remainder.split(":") if field]
    if not fields:
        return allele

    return f"{gene}*{':'.join(fields[:2])}"


def build_profile(alleles_by_gene, genes, two_field_resolution):
    profile = {}
    for gene in genes:
        alleles = {
            normalize_allele(allele, two_field_resolution)
            for allele in alleles_by_gene.get(gene, [])
        }
        alleles.discard("")
        profile[gene] = alleles
    return profile


def score_pair(cluster_name, cluster_profile, reference_name, reference_profile, genes):
    matched_genes = []
    matched_alleles = []
    exact_allele_matches = 0
    compared_alleles = 0
    shared_gene_count = 0

    for gene in genes:
        cluster_alleles = cluster_profile.get(gene, set())
        reference_alleles = reference_profile.get(gene, set())
        if not cluster_alleles or not reference_alleles:
            continue

        shared_gene_count += 1
        overlap = sorted(cluster_alleles & reference_alleles)
        exact_allele_matches += len(overlap)
        compared_alleles += min(len(cluster_alleles), len(reference_alleles))

        if overlap:
            matched_genes.append(gene)
            matched_alleles.extend(overlap)

    overlap_fraction = (
        exact_allele_matches / compared_alleles if compared_alleles else 0.0
    )

    return {
        "cluster_donor": cluster_name,
        "reference_donor": reference_name,
        "exact_allele_matches": exact_allele_matches,
        "compared_alleles": compared_alleles,
        "overlap_fraction": f"{overlap_fraction:.3f}",
        "shared_gene_count": shared_gene_count,
        "matched_genes": " | ".join(matched_genes),
        "matched_alleles": " | ".join(matched_alleles),
    }


def default_assignment(cluster_name, status, note=""):
    return {
        "cluster_donor": cluster_name,
        "assigned_label": cluster_name,
        "assigned_reference_donor": "",
        "best_reference_donor": "",
        "match_status": status,
        "exact_allele_matches": 0,
        "compared_alleles": 0,
        "overlap_fraction": "0.000",
        "shared_gene_count": 0,
        "matched_genes": "",
        "matched_alleles": "",
        "tied_reference_donors": "",
        "notes": note,
    }


def resolve_assignments(score_rows):
    by_cluster = {}
    for row in score_rows:
        by_cluster.setdefault(row["cluster_donor"], []).append(row)

    assignments = []
    for cluster_name, rows in sorted(by_cluster.items()):
        rows = sorted(
            rows,
            key=lambda row: (
                -int(row["exact_allele_matches"]),
                -float(row["overlap_fraction"]),
                row["reference_donor"],
            ),
        )

        best = rows[0]
        best_exact = int(best["exact_allele_matches"])
        best_fraction = float(best["overlap_fraction"])
        tied = [
            row["reference_donor"]
            for row in rows
            if int(row["exact_allele_matches"]) == best_exact
            and float(row["overlap_fraction"]) == best_fraction
        ]

        if best_exact == 0:
            assignment = default_assignment(
                cluster_name,
                status="unmatched",
                note="No exact HLA allele overlap with the reference donor table.",
            )
            assignment["best_reference_donor"] = best["reference_donor"]
        elif len(tied) > 1:
            assignment = default_assignment(
                cluster_name,
                status="ambiguous",
                note="Multiple reference donors share the same top HLA match score.",
            )
            assignment["best_reference_donor"] = best["reference_donor"]
            assignment["exact_allele_matches"] = best_exact
            assignment["compared_alleles"] = best["compared_alleles"]
            assignment["overlap_fraction"] = best["overlap_fraction"]
            assignment["shared_gene_count"] = best["shared_gene_count"]
            assignment["matched_genes"] = best["matched_genes"]
            assignment["matched_alleles"] = best["matched_alleles"]
            assignment["tied_reference_donors"] = " | ".join(tied)
        else:
            assignment = {
                "cluster_donor": cluster_name,
                "assigned_label": best["reference_donor"],
                "assigned_reference_donor": best["reference_donor"],
                "best_reference_donor": best["reference_donor"],
                "match_status": "matched",
                "exact_allele_matches": best_exact,
                "compared_alleles": best["compared_alleles"],
                "overlap_fraction": best["overlap_fraction"],
                "shared_gene_count": best["shared_gene_count"],
                "matched_genes": best["matched_genes"],
                "matched_alleles": best["matched_alleles"],
                "tied_reference_donors": "",
                "notes": "",
            }

        assignments.append(assignment)

    counts = Counter(
        row["assigned_reference_donor"]
        for row in assignments
        if row["match_status"] == "matched" and row["assigned_reference_donor"]
    )

    for row in assignments:
        assigned_reference = row["assigned_reference_donor"]
        if assigned_reference and counts[assigned_reference] > 1:
            row["assigned_label"] = row["cluster_donor"]
            row["assigned_reference_donor"] = ""
            row["match_status"] = "conflicting_match"
            row["notes"] = (
                f"Reference donor {assigned_reference} was assigned to multiple "
                "souporcell clusters."
            )

    return assignments


def read_table(path, label):
    try:
        return hla_table.read_hla_table(path)
    except hla_table.HLATableError as error:
        raise SystemExit(f"Could not read the {label}: {error}")


def main():
    args = parse_args()

    genotype = read_table(args.genotype_summary, "genotype summary")

    requested_genes = [gene.upper() for gene in parse_gene_list(args.genes)]
    cluster_genes = [gene for gene in genotype.genes if gene in requested_genes]

    if not cluster_genes:
        raise SystemExit(
            "None of the requested matching genes were found in the genotype summary."
        )

    cluster_profiles = {
        cluster_name: build_profile(alleles, cluster_genes, args.two_field_resolution)
        for cluster_name, alleles in genotype.donors.items()
    }

    score_rows = []

    if not args.reference_csv:
        assignments = [
            default_assignment(
                cluster_name,
                status="no_reference_annotation",
                note="No reference donor HLA annotation table was provided.",
            )
            for cluster_name in sorted(cluster_profiles)
        ]
    else:
        reference = read_table(args.reference_csv, "donor HLA table")
        shared_genes = [gene for gene in cluster_genes if gene in reference.genes]

        if not shared_genes:
            assignments = [
                default_assignment(
                    cluster_name,
                    status="no_shared_genes",
                    note=(
                        "The reference donor HLA table does not share any requested "
                        "HLA gene columns with the genotype summary."
                    ),
                )
                for cluster_name in sorted(cluster_profiles)
            ]
        else:
            reference_profiles = {
                donor_name: build_profile(alleles, shared_genes, args.two_field_resolution)
                for donor_name, alleles in reference.donors.items()
            }

            for cluster_name, cluster_profile in sorted(cluster_profiles.items()):
                for reference_name, reference_profile in sorted(reference_profiles.items()):
                    score_rows.append(
                        score_pair(
                            cluster_name,
                            cluster_profile,
                            reference_name,
                            reference_profile,
                            shared_genes,
                        )
                    )

            assignments = resolve_assignments(score_rows)

    write_csv(args.assignment_csv, ASSIGNMENT_COLUMNS, assignments)
    write_csv(args.score_csv, SCORE_COLUMNS, score_rows)

    print(f"Wrote donor assignments to {args.assignment_csv}")
    print(f"Wrote pairwise match scores to {args.score_csv}")


if __name__ == "__main__":
    main()
