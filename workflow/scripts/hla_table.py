#!/usr/bin/env python3
"""Read and write donor HLA tables.

An HLA table has one row per donor and two columns per HLA gene, one allele
per cell:

    donor,A_1,A_2,B_1,B_2,C_1,C_2
    Patient_01,A*02:01,A*24:02,B*07:02,B*44:03,C*05:01,C*07:02
    Patient_02,A*11:01,A*11:01,B*15:01,B*40:01,,

- A homozygous gene has the same allele in both columns.
- An empty cell (or NA) means the allele is unknown.
- "A*02:01", "02:01", "HLA-A*02:01" and "A*02:01:01" are all accepted.
- Comma, semicolon and tab separated files are accepted.

The older one-column-per-gene layout (A_alleles = "A*02:01 | A*24:02") is
still accepted when reading. Tables are always written in the two-column
layout, so a table produced by the workflow can be reused as an input.
"""

import csv
import re
from pathlib import Path


PREFERRED_GENE_ORDER = ["A", "B", "C", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5",
                        "DQA1", "DQB1", "DPA1", "DPB1"]
DONOR_COLUMN_NAMES = {"donor", "donor_id", "patient", "patient_id"}
MISSING_VALUES = {"", "na", "nan", "none", "-", "?"}

TWO_COLUMN = re.compile(r"^(?:HLA-)?(?P<gene>[A-Z0-9]+)_(?P<copy>[12])$")
ONE_COLUMN = re.compile(r"^(?:HLA-)?(?P<gene>[A-Z0-9]+)_ALLELES$")
ALLELE_PATTERN = re.compile(r"^[A-Z0-9]+\*\d{2,4}(?::\d{2,4}){0,3}[A-Z]?$")


class HLATableError(ValueError):
    pass


class HLATable:
    def __init__(self, path, genes, donors, ignored_columns):
        self.path = str(path)
        self.genes = genes
        # donor name -> {gene: [allele_1, allele_2]} ("" when unknown)
        self.donors = donors
        self.ignored_columns = ignored_columns


def detect_delimiter(header_line):
    counts = {delimiter: header_line.count(delimiter) for delimiter in ",;\t"}
    delimiter = max(counts, key=counts.get)
    return delimiter if counts[delimiter] else ","


def read_table_rows(path):
    """Read a CSV/TSV file into (fieldnames, rows) with stripped values."""
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        header_line = handle.readline()
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=detect_delimiter(header_line))
        if reader.fieldnames is None:
            raise HLATableError(f"{path}: the file is empty")

        fieldnames = [name.strip() for name in reader.fieldnames]
        rows = []
        for raw_row in reader:
            row = {}
            for name, value in zip(fieldnames, (raw_row.get(key) for key in reader.fieldnames)):
                row[name] = value.strip() if isinstance(value, str) else ""
            if any(row.values()):
                rows.append(row)

    return fieldnames, rows


def ordered_genes(genes):
    genes = set(genes)
    ordered = [gene for gene in PREFERRED_GENE_ORDER if gene in genes]
    return ordered + sorted(gene for gene in genes if gene not in PREFERRED_GENE_ORDER)


def normalize_allele(gene, value):
    """Return 'GENE*xx:yy...' for an allele, '' for a missing value."""
    text = str(value or "").strip()
    if text.lower() in MISSING_VALUES:
        return ""

    text = text.upper().replace(" ", "")
    if text.startswith("HLA-"):
        text = text[4:]
    if "*" not in text:
        text = f"{gene}*{text}"

    allele_gene, fields = text.split("*", 1)
    if allele_gene != gene:
        raise HLATableError(
            f"{value!r} is in a {gene} column but is an allele of {allele_gene}"
        )

    # "2:1" -> "02:01" (spreadsheets often drop leading zeros)
    fields = ":".join(
        field.zfill(2) if field.isdigit() and len(field) == 1 else field
        for field in fields.split(":")
    )
    text = f"{gene}*{fields}"

    if not ALLELE_PATTERN.match(text):
        raise HLATableError(
            f"{value!r} does not look like an HLA allele (expected e.g. {gene}*02:01)"
        )
    return text


def split_alleles(value):
    text = str(value or "").strip()
    if text.lower() in MISSING_VALUES:
        return []
    return [item for item in re.split(r"\s*(?:\||;|,|/)\s*", text) if item]


def read_hla_table(path):
    """Read and validate an HLA table. Raises HLATableError listing every problem."""
    path = Path(path)
    if not path.is_file():
        raise HLATableError(f"{path}: file not found")

    fieldnames, rows = read_table_rows(path)
    problems = []

    donor_column = next(
        (name for name in fieldnames if name.lower() in DONOR_COLUMN_NAMES), None
    )
    if donor_column is None:
        raise HLATableError(
            f"{path}: needs a 'donor' column (found columns: {', '.join(fieldnames)})"
        )

    two_column = {}  # gene -> {"1": column, "2": column}
    one_column = {}  # gene -> column
    ignored_columns = []
    for name in fieldnames:
        if name == donor_column:
            continue
        key = name.upper()
        match = TWO_COLUMN.match(key)
        if match:
            two_column.setdefault(match.group("gene"), {})[match.group("copy")] = name
            continue
        match = ONE_COLUMN.match(key)
        if match and key != "ALL_ALLELES":
            one_column[match.group("gene")] = name
            continue
        ignored_columns.append(name)

    for gene, copies in two_column.items():
        if set(copies) != {"1", "2"}:
            missing = "2" if "1" in copies else "1"
            problems.append(f"gene {gene} has no {gene}_{missing} column")
        if gene in one_column:
            problems.append(f"gene {gene} is given both as {gene}_1/{gene}_2 and as {one_column[gene]}")

    genes = ordered_genes(set(two_column) | set(one_column))
    if not genes:
        raise HLATableError(
            f"{path}: no HLA columns found. Use two columns per gene, "
            "e.g. donor,A_1,A_2,B_1,B_2,C_1,C_2"
        )
    if not rows:
        problems.append("the table has no donor rows")

    donors = {}
    for line_number, row in enumerate(rows, start=2):
        donor = row.get(donor_column, "")
        if not donor:
            problems.append(f"line {line_number}: the donor name is empty")
            continue
        if donor in donors:
            problems.append(f"line {line_number}: donor {donor!r} appears more than once")
            continue

        profile = {}
        for gene in genes:
            if gene in two_column:
                values = [row.get(two_column[gene].get(copy, ""), "") for copy in ("1", "2")]
            else:
                values = split_alleles(row.get(one_column[gene], ""))
                if len(values) > 2:
                    problems.append(
                        f"line {line_number} ({donor}): {gene} has more than two alleles"
                    )
                    values = values[:2]
                if len(values) == 1:
                    values = values * 2
                values += [""] * (2 - len(values))

            alleles = []
            for value in values:
                try:
                    alleles.append(normalize_allele(gene, value))
                except HLATableError as error:
                    problems.append(f"line {line_number} ({donor}): {error}")
                    alleles.append("")
            profile[gene] = alleles
        donors[donor] = profile

    if problems:
        raise HLATableError(
            f"{path}:\n" + "\n".join(f"  - {problem}" for problem in problems)
        )

    return HLATable(path, genes, donors, ignored_columns)


def hla_columns(genes):
    return [f"{gene}_{copy}" for gene in genes for copy in ("1", "2")]


def write_hla_table(path, genes, donors):
    """Write {donor: {gene: [allele_1, allele_2]}} in the two-column layout."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["donor"] + hla_columns(genes))
        for donor, profile in donors.items():
            row = [donor]
            for gene in genes:
                alleles = list(profile.get(gene, [])) + ["", ""]
                row.extend(alleles[:2])
            writer.writerow(row)
