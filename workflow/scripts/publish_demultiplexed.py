#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path

import hla_table


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Collect donor-specific barcode, BAM, FASTQ, genotype and TCR outputs "
            "into one demultiplexed folder per donor, plus a cell-to-donor table."
        )
    )
    parser.add_argument("--pool", required=True)
    parser.add_argument("--clusters", required=True, help="souporcell clusters.tsv")
    parser.add_argument(
        "--tcr-contigs",
        help="Cell Ranger vdj_t/filtered_contig_annotations.csv (TCR runs only).",
    )
    parser.add_argument("--assignment-csv", required=True)
    parser.add_argument("--genotype-summary", required=True)
    parser.add_argument("--barcode-dir", required=True)
    parser.add_argument("--bam-dir", required=True)
    parser.add_argument("--fastq-dir", required=True)
    parser.add_argument("--genotype-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--cell-table", required=True)
    return parser.parse_args()


def read_csv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"CSV file has no header: {path}")

        rows = []
        for row in reader:
            cleaned = {}
            for key, value in row.items():
                cleaned[key] = value.strip() if isinstance(value, str) else value
            rows.append(cleaned)

    return reader.fieldnames, rows


def safe_label(value):
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    cleaned = cleaned.strip("._")
    return cleaned or "unknown_donor"


def unique_label(base_label, cluster_donor, used_labels):
    candidate = safe_label(base_label)
    if candidate not in used_labels:
        used_labels.add(candidate)
        return candidate

    fallback = f"{candidate}__{safe_label(cluster_donor)}"
    counter = 2
    unique_candidate = fallback
    while unique_candidate in used_labels:
        unique_candidate = f"{fallback}_{counter}"
        counter += 1

    used_labels.add(unique_candidate)
    return unique_candidate


def link_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    destination.symlink_to(source.resolve())


def expected_artifacts(cluster_donor, barcode_dir, bam_dir, fastq_dir, genotype_dir):
    return {
        "barcodes.txt": Path(barcode_dir) / f"{cluster_donor}_barcodes.txt",
        "sample_alignments.bam": Path(bam_dir) / f"{cluster_donor}_filtered.bam",
        "sample_alignments.bam.bai": Path(bam_dir)
        / f"{cluster_donor}_filtered.bam.bai",
        "reads.extracted.fq.gz": Path(fastq_dir)
        / f"{cluster_donor}_filtered.extracted.fq.gz",
        "hla.alignment.p": Path(genotype_dir) / f"{cluster_donor}_filtered.alignment.p",
        "hla.genes.json": Path(genotype_dir) / f"{cluster_donor}_filtered.genes.json",
        "hla.genotype.json": Path(genotype_dir)
        / f"{cluster_donor}_filtered.genotype.json",
        "hla.genotype.log": Path(genotype_dir) / f"{cluster_donor}_filtered.genotype.log",
    }


def read_barcodes(path):
    with Path(path).open() as handle:
        return {line.strip() for line in handle if line.strip()}


def read_tcr_contigs(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def write_tcr_contigs(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_cell_table(path, pool, clusters_path, label_by_cluster):
    """One row per cell barcode with its souporcell call and published donor label."""
    with Path(clusters_path).open() as handle, Path(path).open("w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(
            ["pool", "barcode", "donor", "souporcell_status", "souporcell_cluster"]
        )
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] == "barcode":
                continue
            barcode, status, cluster = fields[0], fields[1], fields[2]
            if status == "singlet":
                donor = label_by_cluster.get(f"donor{cluster}", f"donor{cluster}")
            else:
                donor = status  # "doublet" or "unassigned"
            writer.writerow([pool, barcode, donor, status, cluster])


def main():
    args = parse_args()

    assignment_fields, assignment_rows = read_csv(args.assignment_csv)
    try:
        summary = hla_table.read_hla_table(args.genotype_summary)
    except hla_table.HLATableError as error:
        raise SystemExit(f"Could not read the genotype summary: {error}")

    if not assignment_rows:
        raise SystemExit(f"No rows found in assignment CSV: {args.assignment_csv}")

    allele_columns = hla_table.hla_columns(summary.genes)

    tcr_fields, tcr_rows = (None, None)
    if args.tcr_contigs:
        tcr_fields, tcr_rows = read_tcr_contigs(args.tcr_contigs)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    used_labels = set()
    manifest_rows = []

    for row in assignment_rows:
        cluster_donor = row["cluster_donor"]
        if cluster_donor not in summary.donors:
            raise SystemExit(
                f"Cluster {cluster_donor} was present in the assignment CSV but not "
                "in the genotype summary."
            )

        preferred_label = row.get("assigned_label") or cluster_donor
        published_label = unique_label(preferred_label, cluster_donor, used_labels)
        donor_dir = output_dir / published_label

        artifacts = expected_artifacts(
            cluster_donor=cluster_donor,
            barcode_dir=args.barcode_dir,
            bam_dir=args.bam_dir,
            fastq_dir=args.fastq_dir,
            genotype_dir=args.genotype_dir,
        )

        for link_name, source_path in artifacts.items():
            if not source_path.exists():
                raise SystemExit(
                    f"Expected donor artifact is missing for {cluster_donor}: {source_path}"
                )
            link_file(source_path, donor_dir / link_name)

        barcodes = read_barcodes(artifacts["barcodes.txt"])
        n_tcr_cells = ""
        tcr_path = ""
        if tcr_rows is not None:
            donor_contigs = [contig for contig in tcr_rows if contig.get("barcode") in barcodes]
            tcr_path = donor_dir / "tcr_filtered_contig_annotations.csv"
            write_tcr_contigs(tcr_path, tcr_fields, donor_contigs)
            n_tcr_cells = len({contig["barcode"] for contig in donor_contigs})

        manifest_row = {
            "pool": args.pool,
            "donor": published_label,
            "n_cells": len(barcodes),
            "n_tcr_cells": n_tcr_cells,
            "cluster_donor": cluster_donor,
            "published_label": published_label,
            "assigned_label": row.get("assigned_label", ""),
            "assigned_reference_donor": row.get("assigned_reference_donor", ""),
            "best_reference_donor": row.get("best_reference_donor", ""),
            "match_status": row.get("match_status", ""),
            "exact_allele_matches": row.get("exact_allele_matches", ""),
            "compared_alleles": row.get("compared_alleles", ""),
            "overlap_fraction": row.get("overlap_fraction", ""),
            "shared_gene_count": row.get("shared_gene_count", ""),
            "matched_genes": row.get("matched_genes", ""),
            "matched_alleles": row.get("matched_alleles", ""),
            "tied_reference_donors": row.get("tied_reference_donors", ""),
            "notes": row.get("notes", ""),
            "output_dir": str(donor_dir),
            "barcodes_path": str(donor_dir / "barcodes.txt"),
            "bam_path": str(donor_dir / "sample_alignments.bam"),
            "fastq_path": str(donor_dir / "reads.extracted.fq.gz"),
            "genotype_json_path": str(donor_dir / "hla.genotype.json"),
            "tcr_contigs_path": str(tcr_path),
        }

        alleles = summary.donors[cluster_donor]
        for gene in summary.genes:
            manifest_row[f"{gene}_1"], manifest_row[f"{gene}_2"] = alleles[gene]

        manifest_rows.append(manifest_row)

    manifest_fieldnames = [
        "pool",
        "donor",
        "n_cells",
        "n_tcr_cells",
        "cluster_donor",
        "published_label",
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
        "output_dir",
        "barcodes_path",
        "bam_path",
        "fastq_path",
        "genotype_json_path",
        "tcr_contigs_path",
    ] + allele_columns

    manifest_path = Path(args.manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fieldnames, delimiter="\t")
        writer.writeheader()
        for row in manifest_rows:
            writer.writerow(row)

    # Patient folders from an earlier run whose label changed (e.g. donor0 -> a patient
    # name after adding an HLA table). Only folders this script made are removed.
    own_files = set(expected_artifacts("x", ".", ".", ".", ".")) | {"tcr_filtered_contig_annotations.csv"}
    for folder in output_dir.iterdir():
        if folder.is_dir() and not folder.is_symlink() and folder.name not in used_labels:
            entries = list(folder.iterdir())
            if entries and all(entry.name in own_files for entry in entries):
                for entry in entries:
                    entry.unlink()
                folder.rmdir()
                print(f"Removed outdated patient folder {folder}")

    label_by_cluster = {row["cluster_donor"]: row["donor"] for row in manifest_rows}
    write_cell_table(args.cell_table, args.pool, args.clusters, label_by_cluster)

    print(f"Wrote demultiplexed manifest to {manifest_path}")
    print(f"Wrote cell-to-donor table to {args.cell_table}")


if __name__ == "__main__":
    main()
