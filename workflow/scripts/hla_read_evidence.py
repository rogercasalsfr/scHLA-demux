#!/usr/bin/env python3

"""HLA evidence from the reads behind each cluster's arcasHLA call.

arcasHLA's final call can go wrong on 10x data: in its EM step a null allele
(e.g. A*24:608N) can soak up a gene's reads, so the call comes out homozygous
and the second allele is lost, while the reads clearly show both. This script
goes back to the reads. arcasHLA saves, per cluster, every set of alleles a read
fits together with the number of reads (<cluster>_filtered.alignment.p). For
each cluster and gene, leaving null alleles out:

- the allele pair that explains the most reads gives the HLA type from reads,
  at one field (A*02 / A*24): finer fields are rarely decided by 10x reads.
  A second allele that adds fewer than 2% of the reads is not called (homozygous;
  such reads are usually ambient RNA from the other patients);
- support of each called allele: the reads lost (in % of the gene's reads) when
  it is swapped for the best other allele. Small = an alternative fits almost
  as well (e.g. A*24 vs A*23); large = clear call;
- reads that fit only one allele group (e.g. only A*02 alleles);
- fit of each patient of the HLA table: reads explained by the best pair their
  typing allows, relative to the best pair overall. 1.000 = their typing
  explains the cluster's reads as well as any HLA type can.

Outputs:
  --summary-csv   HLA type from reads per cluster, in the layout of an HLA table
  --evidence-csv  per cluster and gene: reads, call, support, closest alternative
  --fit-csv       per cluster, patient and gene: how well the patient's typing fits

Needs numpy: run it with the Python of the arcasHLA environment.
"""

import argparse
import csv
import json
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np

import hla_table


# Genes with fewer reads than this in a cluster are not used.
MIN_GENE_READS = 100
# Alleles kept per one-field group when searching for the best pair (by reads).
ALLELES_PER_GROUP = 60
# A second allele adding fewer reads than this (% of the gene's reads) is not called.
MIN_SECOND_ALLELE_GAIN = 2.0
# Groups listed in the specific-reads column.
SPECIFIC_GROUPS_SHOWN = 4

EVIDENCE_COLUMNS = [
    "cluster_donor",
    "gene",
    "hla_reads",
    "hla_from_reads",
    "explained_pct",
    "allele_support_pct",
    "closest_alternative",
    "best_alleles",
    "specific_reads",
]

FIT_COLUMNS = [
    "cluster_donor",
    "reference_donor",
    "gene",
    "typing",
    "typing_explained_pct",
    "best_explained_pct",
    "fit",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="HLA type from arcasHLA's read classes, and fit of each patient's typing."
    )
    parser.add_argument("--genotype-dir", required=True,
                        help="Folder with the clusters' <cluster>_filtered.alignment.p files.")
    parser.add_argument("--clusters", required=True,
                        help="Comma-separated cluster names (donor0,donor1,...).")
    parser.add_argument("--reference-json", required=True,
                        help="arcasHLA's dat/ref/hla.p.json (allele names).")
    parser.add_argument("--genes", default="A,B,C,DRB1", help="Genes to report.")
    parser.add_argument("--reference-csv", help="Patient HLA table (optional).")
    parser.add_argument("--pool", help="Only score the patients of this pool.")
    parser.add_argument("--matching-genes", default="A,B,C",
                        help="Genes on which patients are scored.")
    parser.add_argument("--two-field-resolution", dest="two_field_resolution",
                        action="store_true",
                        help="Use the patients' typing at two fields at most.")
    parser.add_argument("--no-two-field-resolution", dest="two_field_resolution",
                        action="store_false")
    parser.set_defaults(two_field_resolution=True)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--evidence-csv", required=True)
    parser.add_argument("--fit-csv", required=True)
    return parser.parse_args()


def write_csv(path, fieldnames, rows):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def comma_list(value):
    return [item.strip() for item in str(value).split(",") if item.strip()]


def load_allele_names(path):
    """{allele index: [allele names]} (names sharing one sequence share an index)."""
    with open(path) as handle:
        _, (_, allele_idx, _, _) = json.load(handle)
    return json.loads(allele_idx)


def load_read_classes(path):
    """{gene: [(allele indices, read count)]}; {} when arcasHLA made no alignment."""
    path = Path(path)
    if not path.exists() or not path.stat().st_size:
        return {}
    with path.open("rb") as handle:
        alignment = pickle.load(handle)
    return alignment[1]


def allele_fields(allele):
    """'A*02:01:01:01' -> ['02', '01', '01', '01'] (an expression letter is dropped)."""
    fields = allele.split("*", 1)[1].split(":")
    fields[-1] = fields[-1].rstrip("NLSCAQ") or fields[-1]
    return fields


class GeneReads:
    """One gene's read classes in one cluster."""

    def __init__(self, classes, names):
        self.names = names
        self.counts = np.array([count for _, count in classes], dtype=np.float32)
        self.total = float(self.counts.sum())
        self.rows = defaultdict(list)
        self.reads = defaultdict(float)
        self.specific = defaultdict(float)
        for row, (members, count) in enumerate(classes):
            for index in members:
                self.rows[index].append(row)
                self.reads[index] += count
            groups = {self.group(index) for index in members}
            if len(groups) == 1:
                self.specific[groups.pop()] += count

        # null alleles cannot be the expressed allele, and they are what misleads the EM
        self.usable = [index for index in self.reads
                       if not all(name.endswith("N") for name in names[index])]
        by_group = defaultdict(list)
        for index in self.usable:
            by_group[self.group(index)].append(index)
        self.candidates = [
            index for members in by_group.values()
            for index in sorted(members, key=lambda i: -self.reads[i])[:ALLELES_PER_GROUP]
        ]
        self._matching = {}

    def group(self, index):
        return self.names[index][0].split(":")[0]

    def pct(self, reads):
        return 100 * reads / self.total

    def pair_reads(self, first, second):
        """Reads explained by each pair (first[i], second[j]): reads fitting either allele."""
        def columns(indices):
            matrix = np.zeros((len(self.counts), len(indices)), dtype=np.float32)
            for column, index in enumerate(indices):
                matrix[self.rows[index], column] = 1.0
            return matrix

        m1, m2 = columns(first), columns(second)
        missed_by_first = ((1.0 - m1) * self.counts[:, None]).T
        return (self.counts @ m1)[:, None] + missed_by_first @ m2

    def matching(self, typed):
        """Usable alleles that agree with a typed allele at its resolution (most reads first)."""
        if typed not in self._matching:
            fields = allele_fields(typed)
            found = [
                index for index in self.usable
                if any(allele_fields(name)[:len(fields)] == fields for name in self.names[index])
            ]
            self._matching[typed] = sorted(found, key=lambda i: -self.reads[i])[:ALLELES_PER_GROUP]
        return self._matching[typed]

    def best(self, first, second):
        """Most reads explained by a pair from first x second (either list may be empty)."""
        if first and second:
            return float(self.pair_reads(first, second).max())
        return max((self.reads[index] for index in first or second), default=0.0)


def gene_call(gene_reads):
    """HLA type from reads for one gene, with the support of each called allele."""
    cand = gene_reads.candidates
    explained = gene_reads.pair_reads(cand, cand)
    positions = defaultdict(list)
    for position, index in enumerate(cand):
        positions[gene_reads.group(index)].append(position)
    groups = sorted(positions)

    # best allele pair for every pair of groups (a group with itself = two alleles of it)
    by_groups = {}
    for a, first in enumerate(groups):
        for second in groups[a:]:
            block = explained[np.ix_(positions[first], positions[second])]
            i, j = np.unravel_index(np.argmax(block), block.shape)
            by_groups[first, second] = (
                float(block[i, j]), cand[positions[first][i]], cand[positions[second][j]]
            )

    (first, second), (pair_reads, allele_1, allele_2) = max(
        by_groups.items(), key=lambda item: item[1][0]
    )
    top = max(cand, key=lambda index: gene_reads.reads[index])
    if gene_reads.pct(pair_reads - gene_reads.reads[top]) < MIN_SECOND_ALLELE_GAIN:
        first = second = gene_reads.group(top)
        allele_1 = allele_2 = top
        pair_reads = gene_reads.reads[top]

    # support: best pair with the allele minus best pair without it (keeping the
    # other called allele, if any)
    support = {}
    for group, other in ((first, second), (second, first)):
        def best_pair(keep):
            return max((v[0] for k, v in by_groups.items()
                        if keep(group in k) and (first == second or other in k)), default=0.0)
        support[group] = gene_reads.pct(best_pair(bool) - best_pair(lambda has: not has))

    alternative = max(
        ((k, v[0]) for k, v in by_groups.items()
         if k != (first, second) and (first != second or first not in k)),
        key=lambda item: item[1], default=None,
    )
    specific = sorted(gene_reads.specific.items(), key=lambda item: -item[1])
    return {
        "calls": [first, second],
        "hla_from_reads": f"{first} / {second}",
        "explained_pct": f"{gene_reads.pct(pair_reads):.1f}",
        "allele_support_pct": " | ".join(f"{group} {value:.1f}" for group, value in support.items()),
        "closest_alternative": (
            f"{alternative[0][0]} / {alternative[0][1]} ({gene_reads.pct(alternative[1]):.1f}%)"
            if alternative else ""
        ),
        "best_alleles": " / ".join(gene_reads.names[a][0] for a in (allele_1, allele_2)),
        "specific_reads": " | ".join(
            f"{group} {reads:.0f} ({gene_reads.pct(reads):.1f}%)"
            for group, reads in specific[:SPECIFIC_GROUPS_SHOWN]
        ),
        "best_reads": max(pair_reads, float(explained.max())),
    }


def typing_of(profile, gene, two_field_resolution):
    alleles = [allele for allele in profile.get(gene, []) if allele]
    if two_field_resolution:
        alleles = [f"{gene}*{':'.join(allele_fields(allele)[:2])}" for allele in alleles]
    return alleles


def typing_fit(gene_reads, typed, best_reads):
    """(reads explained by the best pair the typing allows, best pair overall)."""
    if len(typed) == 1:
        # one allele typed: the other one is free
        allowed = (gene_reads.matching(typed[0]), gene_reads.candidates)
    else:
        allowed = (gene_reads.matching(typed[0]), gene_reads.matching(typed[1]))
    reads = gene_reads.best(*allowed)
    return reads, max(best_reads, reads)


def main():
    args = parse_args()
    names = load_allele_names(args.reference_json)
    genes = [gene.upper() for gene in comma_list(args.genes)]
    matching_genes = [gene.upper() for gene in comma_list(args.matching_genes)]

    patients = {}
    if args.reference_csv:
        table = hla_table.read_hla_table(args.reference_csv)
        chosen = table.donors_in_pool(args.pool) if args.pool else list(table.donors)
        patients = {name: table.donors[name] for name in chosen}

    summary, evidence, fits = {}, [], []
    for cluster in comma_list(args.clusters):
        classes = load_read_classes(Path(args.genotype_dir) / f"{cluster}_filtered.alignment.p")
        summary[cluster] = {}
        for gene in genes:
            gene_classes = classes.get(gene, [])
            gene_reads = GeneReads(gene_classes, names) if gene_classes else None
            row = {"cluster_donor": cluster, "gene": gene,
                   "hla_reads": f"{gene_reads.total:.0f}" if gene_reads else "0"}
            if not gene_reads or gene_reads.total < MIN_GENE_READS or not gene_reads.candidates:
                summary[cluster][gene] = ["", ""]
                evidence.append(row)
                continue

            call = gene_call(gene_reads)
            summary[cluster][gene] = call["calls"]
            row.update({key: call[key] for key in EVIDENCE_COLUMNS if key in call})
            evidence.append(row)

            if gene not in matching_genes:
                continue
            for patient, profile in patients.items():
                typed = typing_of(profile, gene, args.two_field_resolution)
                if not typed:
                    continue
                reads, best = typing_fit(gene_reads, typed, call["best_reads"])
                fits.append({
                    "cluster_donor": cluster,
                    "reference_donor": patient,
                    "gene": gene,
                    "typing": " / ".join(typed),
                    "typing_explained_pct": f"{gene_reads.pct(reads):.2f}",
                    "best_explained_pct": f"{gene_reads.pct(best):.2f}",
                    "fit": f"{reads / best:.3f}" if best else "0.000",
                })
        print(f"{cluster}: " + "; ".join(
            f"{gene} {' / '.join(summary[cluster][gene]) or 'no call'}" for gene in genes
        ))

    hla_table.write_hla_table(args.summary_csv, genes, summary)
    write_csv(args.evidence_csv, EVIDENCE_COLUMNS, evidence)
    write_csv(args.fit_csv, FIT_COLUMNS, fits)
    print(f"Wrote {args.summary_csv}, {args.evidence_csv} and {args.fit_csv}")


if __name__ == "__main__":
    main()
