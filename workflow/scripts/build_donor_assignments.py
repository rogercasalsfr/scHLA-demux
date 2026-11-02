#!/usr/bin/env python3

"""Name each souporcell cluster after a patient of the donor HLA table.

Each cluster is scored against each patient in one of two ways:

- reads (with --read-fit-csv, the default in the pipeline): the fit of the
  patient's typing to the cluster's HLA reads, from hla_read_evidence.py. 1.000
  = the typing explains the reads as well as any HLA type can. The fit is
  averaged over the genes the patient is typed for.
- alleles: the cluster's HLA call is compared allele by allele with the
  patient's typing, at the resolution both sides have (A*03 in the table
  matches A*03:01 from arcasHLA).

Clusters are then named one-to-one: every patient names at most one cluster, and
the naming with the best total score wins. So a patient without HLA typing (or a
cluster without HLA reads) gets the one cluster left over when all the others
are named, and is reported as matched_by_elimination. With reads, a cluster is
not named after a patient whose fit is below --min-read-fit, and naming_margin
says how much lower the mean fit of the best other naming is.

match_status:
  matched                 named after the patient, from the HLA comparison
  matched_by_elimination  the only patient left for this cluster (no HLA to compare)
  manual                  named in the manual assignments file
  ambiguous               two or more patients fit (about) equally well: keeps donorN
  unmatched               no patient of the table fits: keeps donorN
  no_reference_annotation no HLA table given: keeps donorN
"""

import argparse
import csv
import itertools
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import chemistry
import hla_table


# Above this many patients the exact one-to-one search gets slow; each cluster is
# then named after its own best patient (ties and patients claimed twice are
# reported as ambiguous).
MAX_EXACT_PATIENTS = 12
# Score of leaving a cluster without a patient (allele matching). A pair with
# nothing to compare scores 0, so it is preferred; a pair whose alleles all
# differ scores lower.
UNASSIGNED = -1
# Read fits are scored in thousandths: 1.000 -> 0, 0.950 -> -50.
READ_FIT_SCALE = 1000


ASSIGNMENT_COLUMNS = [
    "cluster_donor",
    "assigned_label",
    "assigned_reference_donor",
    "best_reference_donor",
    "match_status",
    "read_fit",
    "naming_margin",
    "read_fit_by_gene",
    "exact_allele_matches",
    "compared_alleles",
    "overlap_fraction",
    "shared_gene_count",
    "matched_genes",
    "matched_alleles",
    "tied_reference_donors",
    "runner_up",
    "notes",
]

SCORE_COLUMNS = [
    "cluster_donor",
    "reference_donor",
    "read_fit",
    "read_fit_by_gene",
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
    parser.add_argument(
        "--pool",
        help="Only use the patients of this pool (when the HLA table has a pool column).",
    )
    parser.add_argument(
        "--manual-csv",
        help="Optional CSV (cluster,donor) naming clusters by hand, e.g. donor3,Patient_07.",
    )
    parser.add_argument(
        "--read-fit-csv",
        help="Fit of each patient's typing to each cluster's HLA reads "
             "(from hla_read_evidence.py). Given: clusters are named from it.",
    )
    parser.add_argument(
        "--min-read-fit",
        type=float,
        default=0.90,
        help="With --read-fit-csv: lowest mean fit for naming a cluster after a patient.",
    )
    parser.add_argument(
        "--min-margin",
        type=float,
        default=0.03,
        help="With --read-fit-csv: a cluster whose best other naming is closer than "
             "this (in mean fit) is reported as ambiguous.",
    )
    parser.add_argument(
        "--cellranger-outs",
        help="Cell Ranger outs folder of the pool: with 3' data, a warning is added to the "
             "clusters that could not be named.",
    )
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


def allele_fields(allele):
    return allele.split("*", 1)[1].split(":") if "*" in allele else [allele]


def alleles_match(first, second):
    """Same allele at the resolution both have: A*03 matches A*03:01, not A*30."""
    first, second = allele_fields(first), allele_fields(second)
    common = min(len(first), len(second))
    return first[:common] == second[:common]


def gene_overlap(cluster_alleles, reference_alleles):
    """Largest set of cluster alleles that pair up with different reference alleles."""
    small, large = sorted([sorted(cluster_alleles), sorted(reference_alleles)], key=len)
    best = []
    for order in itertools.permutations(large, len(small)):
        pairs = [(a, b) for a, b in zip(small, order) if alleles_match(a, b)]
        if len(pairs) > len(best):
            best = pairs
    cluster_side = set(cluster_alleles)
    return sorted(a if a in cluster_side else b for a, b in best)


def read_fits(path, reference_profiles, genes):
    """{(cluster, patient): (mean fit, per-gene text)} over the genes the patient is typed for."""
    by_pair = defaultdict(dict)
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            by_pair[row["cluster_donor"], row["reference_donor"]][row["gene"]] = float(row["fit"])

    fits = {}
    for (cluster, patient), by_gene in by_pair.items():
        typed = [gene for gene in genes
                 if gene in by_gene and reference_profiles.get(patient, {}).get(gene)]
        if typed:
            fits[cluster, patient] = (
                sum(by_gene[gene] for gene in typed) / len(typed),
                " | ".join(f"{gene} {by_gene[gene]:.3f}" for gene in typed),
            )
    return fits


def score_pair(cluster_name, cluster_profile, reference_name, reference_profile, genes,
               read_fit=None):
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
        overlap = gene_overlap(cluster_alleles, reference_alleles)
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
        "read_fit": f"{read_fit[0]:.3f}" if read_fit else "",
        "read_fit_by_gene": read_fit[1] if read_fit else "",
        "exact_allele_matches": exact_allele_matches,
        "compared_alleles": compared_alleles,
        "overlap_fraction": f"{overlap_fraction:.3f}",
        "shared_gene_count": shared_gene_count,
        "matched_genes": " | ".join(matched_genes),
        "matched_alleles": " | ".join(matched_alleles),
    }


def pair_score(row, by_reads):
    """Reads: 0 for a perfect fit, lower for worse ones. Alleles: shared alleles
    count most, differing alleles count a little against."""
    if by_reads:
        if not row["read_fit"]:
            return 0
        return round(READ_FIT_SCALE * (float(row["read_fit"]) - 1))
    matches = int(row["exact_allele_matches"])
    differing = int(row["compared_alleles"]) - matches
    return 2 * (100 * matches - differing)


def has_comparison(row, by_reads):
    return bool(row["read_fit"]) if by_reads else int(row["compared_alleles"]) > 0


def default_assignment(cluster_name, status, note=""):
    return {
        "cluster_donor": cluster_name,
        "assigned_label": cluster_name,
        "assigned_reference_donor": "",
        "best_reference_donor": "",
        "match_status": status,
        "read_fit": "",
        "naming_margin": "",
        "read_fit_by_gene": "",
        "exact_allele_matches": 0,
        "compared_alleles": 0,
        "overlap_fraction": "0.000",
        "shared_gene_count": 0,
        "matched_genes": "",
        "matched_alleles": "",
        "tied_reference_donors": "",
        "runner_up": "",
        "notes": note,
    }


def best_naming(clusters, patients, score, allowed, unassigned, must_name=None):
    """Best one-to-one naming: (total score, {cluster: patient or None}).
    must_name: a cluster that cannot be left without a patient."""

    @lru_cache(maxsize=None)
    def solve(index, used):
        if index == len(clusters):
            return 0, ()
        cluster = clusters[index]
        total, rest = solve(index + 1, used)
        best = (total + unassigned if cluster != must_name else float("-inf"), (None,) + rest)
        for j, patient in enumerate(patients):
            if used >> j & 1 or not allowed(cluster, patient):
                continue
            total, rest = solve(index + 1, used | 1 << j)
            if total + score[cluster, patient] > best[0]:
                best = (total + score[cluster, patient], (patient,) + rest)
        return best

    total, choice = solve(0, 0)
    return total, dict(zip(clusters, choice))


def greedy_naming(clusters, patients, score, allowed, unassigned):
    """Each cluster takes its own best patient (for tables too large for best_naming)."""
    choice = {}
    for cluster in clusters:
        options = [p for p in patients if allowed(cluster, p) and score[cluster, p] > unassigned]
        choice[cluster] = max(options, key=lambda p: score[cluster, p], default=None)
    return None, choice


def fill_from_score(assignment, row):
    for key in ("read_fit", "read_fit_by_gene", "exact_allele_matches", "compared_alleles",
                "overlap_fraction", "shared_gene_count", "matched_genes", "matched_alleles"):
        assignment[key] = row[key]


def resolve_assignments(score_rows, clusters, patients, manual, by_reads,
                        min_read_fit=0.90, min_margin=0.03):
    by_pair = {(row["cluster_donor"], row["reference_donor"]): row for row in score_rows}
    score = {pair: pair_score(row, by_reads) for pair, row in by_pair.items()}
    manual_patients = set(manual.values())
    if by_reads:
        unassigned = round(READ_FIT_SCALE * (min_read_fit - 1))
        # another naming closer than this to the best one counts as a tie
        tolerance = round(READ_FIT_SCALE * min_margin)
    else:
        unassigned, tolerance = UNASSIGNED, 1

    def allowed(cluster, patient):
        if cluster in manual:
            return patient == manual[cluster]
        return patient not in manual_patients

    exact = len(patients) <= MAX_EXACT_PATIENTS
    solver = best_naming if exact else greedy_naming
    total, choice = solver(clusters, patients, score, allowed, unassigned)

    def other_namings(cluster, patient):
        """{other patient: how much lower the best total is when this cluster takes them}."""
        others = [p for p in patients if p != patient and allowed(cluster, p)]
        if exact:
            return {
                other: total - best_naming(
                    clusters, patients, score,
                    lambda c, p: allowed(c, p) and (c != cluster or p == other), unassigned,
                    must_name=cluster,
                )[0]
                for other in others
            }
        return {other: score[cluster, patient] - score[cluster, other] for other in others}

    def runner_up_text(cluster, other):
        row = by_pair[cluster, other]
        if by_reads:
            return f"{other} (fit {row['read_fit'] or 'n/a'})"
        return f"{other} ({row['exact_allele_matches']}/{row['compared_alleles']})"

    assignments = []
    for cluster in clusters:
        # patients with HLA to compare first (an untyped patient is never the "best" match)
        row_scores = sorted(
            (p for p in patients if (cluster, p) in by_pair),
            key=lambda p: (not has_comparison(by_pair[cluster, p], by_reads), -score[cluster, p], p),
        )
        best = row_scores[0] if row_scores else ""
        patient = choice[cluster]
        runner_up = next((p for p in row_scores if p != patient), "")

        if patient is None:
            assignment = default_assignment(
                cluster,
                status="unmatched",
                note=(f"No patient's HLA typing fits this cluster's reads well enough "
                      f"(best mean fit below {min_read_fit:.2f}), or there are more clusters "
                      "than patients." if by_reads else
                      "No patient of the HLA table fits this cluster (all compared alleles differ, "
                      "or there are more clusters than patients)."),
            )
            assignment["best_reference_donor"] = best
            assignment["runner_up"] = runner_up_text(cluster, runner_up) if runner_up else ""
            if best and by_reads:
                assignment["read_fit"] = by_pair[cluster, best]["read_fit"]
                assignment["read_fit_by_gene"] = by_pair[cluster, best]["read_fit_by_gene"]
            assignments.append(assignment)
            continue

        row = by_pair[cluster, patient]
        assignment = default_assignment(cluster, status="", note="")
        fill_from_score(assignment, row)
        assignment["best_reference_donor"] = best
        assignment["runner_up"] = runner_up_text(cluster, runner_up) if runner_up else ""

        tied = []
        if cluster not in manual:
            gaps = other_namings(cluster, patient)
            tied = [other for other, gap in gaps.items() if gap < tolerance]
            # the same patient on a cluster left without a name fits about as well: e.g. one
            # untyped patient and two clusters without a match, where either could be theirs
            tied += [
                f"{other} (left without a name)" for other in clusters
                if other != cluster and choice[other] is None and allowed(other, patient)
                and score[other, patient] > unassigned
                and score[cluster, patient] - score[other, patient] < tolerance
            ]
            if not exact:
                claimed = [c for c in clusters if c != cluster and choice[c] == patient]
                if claimed:
                    tied.append(f"(also best for {', '.join(claimed)})")
            # margin against typed patients only: giving the cluster to an untyped
            # patient just costs the threshold for the cluster left without one
            typed_gaps = [gap for other, gap in gaps.items()
                          if has_comparison(by_pair[cluster, other], by_reads)]
            if by_reads and typed_gaps:
                assignment["naming_margin"] = f"{min(typed_gaps) / READ_FIT_SCALE:.3f}"

        if cluster in manual:
            status, note = "manual", "Named in the manual assignments file."
        elif tied:
            status = "ambiguous"
            note = ("Another naming fits the reads almost as well (margin below "
                    f"{min_margin:.2f})." if by_reads else
                    "Several patients fit this cluster equally well.")
            if not has_comparison(row, by_reads):
                note = (f"{patient} has no HLA typing to compare, and more than one cluster is "
                        "left without a match: it could be any of them.")
            assignment["tied_reference_donors"] = " | ".join([patient] + tied)
        elif not has_comparison(row, by_reads):
            status = "matched_by_elimination"
            note = (
                f"No HLA to compare ({patient} has no typing for the compared genes, "
                f"or this cluster has no HLA {'reads' if by_reads else 'call'}); {patient} "
                "is the only patient left once the other clusters are named."
            )
        else:
            status, note = "matched", ""

        assignment["match_status"] = status
        assignment["notes"] = note
        if status in ("matched", "matched_by_elimination", "manual"):
            assignment["assigned_label"] = patient
            assignment["assigned_reference_donor"] = patient
        assignments.append(assignment)

    unused = [p for p in patients if p not in {a["assigned_reference_donor"] for a in assignments}]
    if unused:
        print(f"Patients not given a cluster: {', '.join(unused)}")
    return assignments


def read_manual(path, clusters):
    """{cluster: patient} from a cluster,donor CSV (cluster as donor3 or 3)."""
    _, rows = hla_table.read_table_rows(path)
    manual = {}
    for line_number, row in enumerate(rows, start=2):
        row = {key.lower(): value for key, value in row.items()}
        cluster = row.get("cluster") or row.get("cluster_donor") or row.get("souporcell_cluster", "")
        patient = row.get("donor") or row.get("patient", "")
        if cluster.isdigit():
            cluster = f"donor{cluster}"
        if cluster not in clusters:
            raise SystemExit(f"{path} line {line_number}: unknown cluster {cluster!r} "
                             f"(clusters: {', '.join(clusters)})")
        if not patient:
            raise SystemExit(f"{path} line {line_number}: the donor is empty")
        if patient in manual.values():
            raise SystemExit(f"{path} line {line_number}: {patient} is given to two clusters")
        manual[cluster] = patient
    return manual


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
    clusters = sorted(cluster_profiles, key=lambda name: (len(name), name))
    manual = read_manual(args.manual_csv, clusters) if args.manual_csv else {}

    if not args.reference_csv and not manual:
        assignments = [
            default_assignment(
                cluster_name,
                status="no_reference_annotation",
                note="No reference donor HLA annotation table was provided.",
            )
            for cluster_name in clusters
        ]
    else:
        if args.reference_csv:
            reference = read_table(args.reference_csv, "donor HLA table")
            patients = reference.donors_in_pool(args.pool) if args.pool else list(reference.donors)
            shared_genes = [gene for gene in cluster_genes if gene in reference.genes]
            reference_profiles = {
                name: build_profile(reference.donors[name], shared_genes, args.two_field_resolution)
                for name in patients
            }
        else:
            patients, shared_genes, reference_profiles = [], [], {}

        for patient in manual.values():
            if patient not in reference_profiles:
                if args.reference_csv:
                    raise SystemExit(
                        f"{args.manual_csv}: {patient} is not a patient of "
                        f"{'pool ' + args.pool + ' in ' if args.pool else ''}{args.reference_csv}"
                    )
                reference_profiles[patient] = {}
                patients.append(patient)

        if args.reference_csv and not shared_genes and not manual:
            assignments = [
                default_assignment(
                    cluster_name,
                    status="no_shared_genes",
                    note=(
                        "The reference donor HLA table does not share any requested "
                        "HLA gene columns with the genotype summary."
                    ),
                )
                for cluster_name in clusters
            ]
        else:
            patients = sorted(patients)
            by_reads = bool(args.read_fit_csv)
            fits = read_fits(args.read_fit_csv, reference_profiles, shared_genes) if by_reads else {}
            for cluster_name in clusters:
                for reference_name in patients:
                    score_rows.append(
                        score_pair(
                            cluster_name,
                            cluster_profiles[cluster_name],
                            reference_name,
                            reference_profiles[reference_name],
                            shared_genes,
                            fits.get((cluster_name, reference_name)),
                        )
                    )

            assignments = resolve_assignments(
                score_rows, clusters, patients, manual, by_reads,
                args.min_read_fit, args.min_margin,
            )

    if args.cellranger_outs:
        library = chemistry.gene_expression_chemistry(args.cellranger_outs)
        if chemistry.is_five_prime(library) is False:
            print(f"WARNING: Cell Ranger chemistry is {library}. {chemistry.THREE_PRIME_NOTE}")
            for assignment in assignments:
                if assignment["match_status"] in ("unmatched", "ambiguous"):
                    assignment["notes"] = (f"{assignment['notes']} {library} library: "
                                           f"{chemistry.THREE_PRIME_NOTE}").strip()

    write_csv(args.assignment_csv, ASSIGNMENT_COLUMNS, assignments)
    write_csv(args.score_csv, SCORE_COLUMNS, score_rows)

    print(f"Wrote donor assignments to {args.assignment_csv}")
    print(f"Wrote pairwise match scores to {args.score_csv}")


if __name__ == "__main__":
    main()
