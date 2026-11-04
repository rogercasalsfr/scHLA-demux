#!/usr/bin/env python3
"""Read and write donor HLA tables.

An HLA table has one row per donor and two columns per HLA gene, one allele
per cell:

    donor,A_1,A_2,B_1,B_2,C_1,C_2
    Patient_01,A*02:01,A*24:02,B*07:02,B*44:03,C*05:01,C*07:02
    Patient_02,A*11:01,A*11:01,B*15:01,B*40:01,,

- A homozygous gene has the same allele in both columns.
- An empty cell (or NA) means the allele is unknown. A donor with no alleles at
  all (not typed) is allowed.
- "A*02:01", "02:01", "02 : 01", "HLA-A*02:01" and "A*02:01:01" are all accepted, and so
  is low-resolution typing: "A*02", "A2" or just "2" in an A column.
- Gene columns can be named A_1/A_2 or A1/A2, with or without "HLA-" or "HLA_" before.
- An optional "pool" column says which pool each donor is in, so one table can
  serve several pools.
- Comma, semicolon and tab separated files are accepted, and Excel workbooks: one sheet
  of a workbook is given as "typing.xlsx#Sheet1" (see choose_sheet).

The older one-column-per-gene layout (A_alleles = "A*02:01 | A*24:02") is
still accepted when reading. Tables are always written in the two-column
layout, so a table produced by the workflow can be reused as an input.
"""

import csv
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree


PREFERRED_GENE_ORDER = ["A", "B", "C", "DRB1", "DRB2", "DRB3", "DRB4", "DRB5",
                        "DQA1", "DQB1", "DPA1", "DPB1"]
# In order of preference, compared in lower case with spaces and hyphens as "_" ("Patient ID").
DONOR_COLUMN_NAMES = ["donor", "donor_id", "patient", "patient_id", "participant",
                      "participant_id", "subject", "subject_id", "sample", "sample_id"]
POOL_COLUMN_NAMES = {"pool"}
MISSING_VALUES = {"", "na", "nan", "none", "-", "?"}

KNOWN_GENES = set(PREFERRED_GENE_ORDER) | {"E", "F", "G"}

HLA_PREFIX = re.compile(r"^HLA[-_. ]?(?=[A-Z])")  # HLA-A1, HLA_A_1 -> A1, A_1
TWO_COLUMN = re.compile(r"^(?P<gene>[A-Z0-9]+)[_.\- ](?P<copy>[12])$")
ONE_COLUMN = re.compile(r"^(?P<gene>[A-Z0-9]+)_ALLELES$")
ALLELE_PATTERN = re.compile(r"^[A-Z0-9]+\*\d{2,4}(?::\d{2,4}){0,3}[A-Z]?$")
# P and G groups (B*39:01P, A*02:01:01G) name a group of alleles after its first allele.
# The HLA reference has no group names, so the letter is dropped.
GROUP_SUFFIX = re.compile(r"(?<=\d)[PG]$")
# Column names that were not understood but look like an HLA column (a typo such as A3).
LOOKS_LIKE_HLA = re.compile(
    r"^(?:HLA[-_ ]?)?(?:A|B|C|E|F|G|DRB\d?|DQA\d?|DQB\d?|DPA\d?|DPB\d?)[-_. ]?\d{0,2}(?:_?ALLELES?)?$"
)


class HLATableError(ValueError):
    pass


class HLATable:
    def __init__(self, path, genes, donors, ignored_columns, pools=None, group_alleles=None):
        self.path = str(path)
        self.genes = genes
        # donor name -> {gene: [allele_1, allele_2]} ("" when unknown)
        self.donors = donors
        self.ignored_columns = ignored_columns
        # donor name -> pool, when the table has a pool column ("" = any pool)
        self.pools = pools or {}
        # P/G group as written -> allele used ("B*39:01P" -> "B*39:01")
        self.group_alleles = group_alleles or {}

    @property
    def suspicious_columns(self):
        """Ignored columns whose name looks like a mistyped HLA column."""
        return [name for name in self.ignored_columns if LOOKS_LIKE_HLA.match(name.upper())]

    def donors_in_pool(self, pool):
        """Donors of this pool: all of them when the table has no pool column."""
        if not self.pools:
            return list(self.donors)
        return [donor for donor in self.donors if self.pools.get(donor, "") in ("", pool)]


def detect_delimiter(header_line):
    counts = {delimiter: header_line.count(delimiter) for delimiter in ",;\t"}
    delimiter = max(counts, key=counts.get)
    return delimiter if counts[delimiter] else ","


XLSX_MAIN = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
XLSX_REL_ID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
# Excel's built-in time formats (h:mm, h:mm:ss, [h]:mm:ss, ...).
XLSX_TIME_FORMATS = {18, 19, 20, 21, 45, 46, 47}


def split_sheet(spec):
    """'typing.xlsx#Pool1' -> ('typing.xlsx', 'Pool1'); anything else -> (spec, None)."""
    match = re.match(r"^(?P<path>.+\.xlsx)#(?P<sheet>.+)$", str(spec), re.IGNORECASE)
    if match and not Path(spec).exists():
        return match.group("path"), match.group("sheet")
    return str(spec), None


def xlsx_sheets(archive):
    """[(sheet name, its file in the archive)], in the workbook's order."""
    sheets = ElementTree.fromstring(archive.read("xl/workbook.xml")).iter(f"{XLSX_MAIN}sheet")
    relations = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {rel.get("Id"): rel.get("Target").lstrip("/") for rel in relations}
    found = []
    for sheet in sheets:
        target = targets[sheet.get(XLSX_REL_ID)]
        found.append((sheet.get("name"), target if target.startswith("xl/") else f"xl/{target}"))
    return found


def sheet_names(path):
    with zipfile.ZipFile(path) as archive:
        return [name for name, _ in xlsx_sheets(archive)]


def xlsx_time_styles(archive):
    """The cell styles that show a number as a time. Excel turns 02:01 typed in a cell into
    the time 2:01 (a fraction of a day); such cells are read as the h:mm text Excel shows."""
    if "xl/styles.xml" not in archive.namelist():
        return set()
    styles = ElementTree.fromstring(archive.read("xl/styles.xml"))
    time_formats = set(XLSX_TIME_FORMATS)
    for number_format in styles.iter(f"{XLSX_MAIN}numFmt"):
        code = re.sub(r'"[^"]*"|\[(?![hms]+\])[^\]]*\]', "", number_format.get("formatCode", "").lower())
        if "h" in code and not re.search(r"[yd]", code):
            time_formats.add(int(number_format.get("numFmtId")))
    cell_formats = styles.find(f"{XLSX_MAIN}cellXfs")
    return {
        index for index, cell_format in enumerate(cell_formats if cell_formats is not None else [])
        if int(cell_format.get("numFmtId", 0)) in time_formats
    }


def read_xlsx_rows(path, sheet=None):
    """Read one sheet (the first by default) of an .xlsx file into (fieldnames, rows),
    standard library only."""
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            for item in ElementTree.fromstring(archive.read("xl/sharedStrings.xml")).iter(f"{XLSX_MAIN}si"):
                shared.append("".join(text.text or "" for text in item.iter(f"{XLSX_MAIN}t")))

        sheets = xlsx_sheets(archive)
        if sheet is None:
            sheet, target = sheets[0]
        else:
            chosen = [entry for entry in sheets if entry[0] == sheet] or [
                entry for entry in sheets if entry[0].strip().lower() == sheet.strip().lower()
            ]
            if not chosen:
                raise HLATableError(
                    f"{path}: no sheet named {sheet!r} (sheets: {', '.join(name for name, _ in sheets)})"
                )
            sheet, target = chosen[0]
        time_styles = xlsx_time_styles(archive)
        content = ElementTree.fromstring(archive.read(target))

        table = []
        for row in content.iter(f"{XLSX_MAIN}row"):
            values = {}
            for cell in row.iter(f"{XLSX_MAIN}c"):
                letters = re.match(r"[A-Z]+", cell.get("r", "")).group(0)
                column = 0
                for letter in letters:
                    column = column * 26 + ord(letter) - 64
                kind, value = cell.get("t"), cell.find(f"{XLSX_MAIN}v")
                if kind == "s":
                    text = shared[int(value.text)]
                elif kind == "inlineStr":
                    text = "".join(t.text or "" for t in cell.iter(f"{XLSX_MAIN}t"))
                else:
                    text = value.text if value is not None and value.text else ""
                    if text and kind in (None, "n") and int(cell.get("s", 0)) in time_styles:
                        seconds = round(float(text) * 86400)
                        text = f"{seconds // 3600}:{seconds % 3600 // 60:02d}"
                        text += f":{seconds % 60:02d}" if seconds % 60 else ""
                    text = re.sub(r"^(-?\d+)\.0$", r"\1", text)  # 1234567.0 -> 1234567
                values[column - 1] = text.strip()
            if any(values.values()):
                table.append([values.get(index, "") for index in range(max(values) + 1)])

    if not table:
        raise HLATableError(f"{path}: sheet {sheet!r} is empty")
    fieldnames = [name.strip() for name in table[0]]
    rows = [dict(zip(fieldnames, values + [""] * (len(fieldnames) - len(values)))) for values in table[1:]]
    return fieldnames, rows


def read_table_rows(path, sheet=None):
    """Read a CSV/TSV (or one sheet of an .xlsx) file into (fieldnames, rows) with stripped values."""
    if str(path).lower().endswith(".xlsx"):
        return read_xlsx_rows(path, sheet)
    if sheet is not None:
        raise HLATableError(f"{path}: only an Excel .xlsx file has sheets (#{sheet})")

    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        header_line = handle.readline()
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=detect_delimiter(header_line))
        if reader.fieldnames is None:
            raise HLATableError(f"{path}: the file is empty")

        fieldnames = [name.strip() for name in reader.fieldnames]
        rows = []
        for line_number, raw_row in enumerate(reader, start=2):
            # values beyond the header: a value holding the separator (e.g. a;b in a file
            # separated by ;) that is not in quotes would otherwise be cut silently
            if any(str(extra).strip() for extra in raw_row.get(None) or []):
                raise HLATableError(
                    f"{path} line {line_number}: more values than columns; put a value that "
                    f"contains the separator {reader.reader.dialect.delimiter!r} in quotes"
                )
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

    text = re.sub(r"\s+", "", text.upper())  # "23 : 01", also non-breaking spaces
    text = re.sub(r"^HLA[-_]", "", text)
    if "*" not in text:
        # "A2" or "Cw7" (serological style) -> "2", "7"
        for prefix in [gene, "CW"] if gene == "C" else [gene]:
            if text.startswith(prefix) and text[len(prefix):][:1].isdigit():
                text = text[len(prefix):]
                break
        text = f"{gene}*{text}"

    allele_gene, fields = text.split("*", 1)
    if allele_gene != gene:
        raise HLATableError(
            f"{value!r} is in a {gene} column but is an allele of {allele_gene}"
        )
    fields = GROUP_SUFFIX.sub("", fields)

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


def find_columns(fieldnames):
    """Sort a table's columns: (donor column, pool column, {gene: {"1": column, "2": column}},
    {gene: column} for one column per gene, ignored columns)."""
    simple = {re.sub(r"[\s-]+", "_", name.strip().lower()): name for name in reversed(fieldnames)}
    donor_column = next((simple[name] for name in DONOR_COLUMN_NAMES if name in simple), None)
    pool_column = next(
        (name for name in fieldnames if name.lower() in POOL_COLUMN_NAMES), None
    )

    two_column = {}  # gene -> {"1": column, "2": column}
    one_column = {}  # gene -> column
    ignored_columns = []
    for name in fieldnames:
        if name in (donor_column, pool_column):
            continue
        key = HLA_PREFIX.sub("", name.upper())
        match = TWO_COLUMN.match(key)
        if match:
            two_column.setdefault(match.group("gene"), {})[match.group("copy")] = name
            continue
        if key not in KNOWN_GENES and key[:-1] in KNOWN_GENES and key[-1:] in ("1", "2"):
            two_column.setdefault(key[:-1], {})[key[-1]] = name  # A1, A2, DRB11
            continue
        match = ONE_COLUMN.match(key)
        if match and key != "ALL_ALLELES":
            one_column[match.group("gene")] = name
            continue
        if key in KNOWN_GENES:  # column "A" holding "A*02:01, A*24:02"
            one_column[key] = name
            continue
        ignored_columns.append(name)
    return donor_column, pool_column, two_column, one_column, ignored_columns


def choose_sheet(path, pool):
    """The sheet to read for this pool from a workbook with several sheets: the sheet named
    like the pool, else the only sheet that holds an HLA table. None for a one-sheet workbook
    or a text file (the whole file is read). Raises HLATableError when it cannot tell."""
    if not str(path).lower().endswith(".xlsx") or not Path(path).is_file():
        return None
    names = sheet_names(path)
    if len(names) == 1:
        return None
    named = [name for name in names if name.strip().lower() == pool.lower()]
    if named:
        return named[0]
    tables = []
    for name in names:
        try:
            donor_column, _, two_column, one_column, _ = find_columns(read_xlsx_rows(path, name)[0])
        except HLATableError:  # an empty sheet
            continue
        if donor_column and (two_column or one_column):
            tables.append(name)
    if len(tables) == 1:
        return tables[0]
    raise HLATableError(
        f"{path} has {len(names)} sheets ({', '.join(names)}), "
        f"{len(tables) or 'none'} of them with an HLA table, and none named {pool}: "
        f"say which sheet to use, e.g. {Path(path).name}#{(tables or names)[0]} in hla_csv"
    )


def read_hla_table(spec):
    """Read and validate an HLA table: a file, or one sheet of a workbook ("typing.xlsx#Pool1").
    Raises HLATableError listing every problem."""
    path, sheet = split_sheet(spec)
    path = Path(path)
    if not path.is_file():
        raise HLATableError(f"{path}: file not found")

    fieldnames, rows = read_table_rows(path, sheet)
    if sheet is not None:
        path = f"{path} (sheet {sheet})"
    problems = []

    donor_column, pool_column, two_column, one_column, ignored_columns = find_columns(fieldnames)
    if donor_column is None:
        raise HLATableError(
            f"{path}: needs a 'donor' column (found columns: {', '.join(fieldnames)})"
        )

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
    pools = {}
    group_alleles = {}
    for line_number, row in enumerate(rows, start=2):
        donor = row.get(donor_column, "")
        if not donor:
            problems.append(f"line {line_number}: the donor name is empty")
            continue
        pool = row.get(pool_column, "") if pool_column else ""

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
                    allele = normalize_allele(gene, value)
                except HLATableError as error:
                    problems.append(f"line {line_number} ({donor}): {error}")
                    allele = ""
                group = GROUP_SUFFIX.search(str(value or "").strip().upper())
                if allele and group:
                    group_alleles[allele + group.group()] = allele
                alleles.append(allele)
            profile[gene] = alleles

        if donor in donors:
            # the same row twice is harmless; different values are a mistake
            if donors[donor] != profile or pools.get(donor, "") != pool:
                problems.append(
                    f"line {line_number}: donor {donor!r} appears more than once, "
                    "with different values"
                )
            continue
        donors[donor] = profile
        pools[donor] = pool

    if problems:
        raise HLATableError(
            f"{path}:\n" + "\n".join(f"  - {problem}" for problem in problems)
        )

    return HLATable(path, genes, donors, ignored_columns, pools if pool_column else None,
                    group_alleles)


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
