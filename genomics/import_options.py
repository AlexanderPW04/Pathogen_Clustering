"""How to read one uploaded table: choices the user confirms on the review page.

Sprint 5 Phase 4A. Headers alone cannot say whether frequencies are
fractions or percentages, whether positions are 1-based or 0-based, which
reference a table without a reference column used, or what a depth column
counts. Nothing here is decided automatically: the review page asks, the
answer is recorded in the import provenance, and the original value is kept
in a ``Source_<field>_Original`` column whenever a value is changed.

Options (per file, all optional):

- ``frequency_scale``: "fraction" (0-1, the default the loaders validate) or
  "percent" (0-100, divided by 100 here with exact decimal arithmetic, so the
  significant digits used by the Phase 1 count-reconstruction rule are kept).
- ``position_base``: "1" (default) or "0" (each position + 1).
- ``declared_reference``: the reference accession for a table that has no
  reference column, recorded as stated by the user and later checked base by
  base against an uploaded reference FASTA.
- ``depth_definition``: what the depth column counts, as documented by the
  source; recorded for provenance only (it changes no calculation).
- ``base_count_minimum``: minimum reads for a non-reference base when variant
  records are derived from A/C/G/T counts (see `genomics.table_layouts`).
- ``layout_confirmed``: the user confirmed a wide-matrix reading.
- ``header_row``: 1-based row number of the column headers.
- ``sample_from_filename``: the file holds one sample, named by the file
  (for per-sample variant tables without a sample column).
"""

from decimal import Decimal, InvalidOperation

import pandas as pd

FRACTION, PERCENT = "fraction", "percent"
FREQUENCY_SCALES = {FRACTION: "Fractions from 0 to 1", PERCENT: "Percentages from 0 to 100"}
POSITION_BASES = {"1": "1-based (the first base is position 1)", "0": "0-based (the first base is position 0)"}
DEPTH_DEFINITIONS = {
    "": "Not stated",
    "total": "All reads covering the position (total depth)",
    "ref_alt": "Reference plus alternative reads only",
    "filtered": "Reads passing the source's quality filters",
}
DOCUMENTED_DEPTH = {
    "Fontenele et al. (2021) Supplementary Table 2": "TOTAL_DP: total depth, as defined by iVar",
    "iVar variants table": "TOTAL_DP: total depth at the position (iVar documentation)",
    "VCF": "DP / AD fields as declared in the VCF header",
}
OPTION_NAMES = ["frequency_scale", "position_base", "declared_reference", "depth_definition",
                "base_count_minimum", "layout_confirmed", "header_row", "sample_from_filename"]


def read_file_options(form, index):
    """The interpretation choices for file `index` from the review form (validated, never guessed)."""

    prefix = f"opt_{index}_"
    options = {}
    scale = form.get(prefix + "frequency_scale", "")
    if scale in FREQUENCY_SCALES:
        options["frequency_scale"] = scale
    base = form.get(prefix + "position_base", "")
    if base in POSITION_BASES:
        options["position_base"] = base
    reference = str(form.get(prefix + "declared_reference", "")).strip()[:80]
    if reference:
        options["declared_reference"] = reference
    depth = form.get(prefix + "depth_definition", "")
    if depth in DEPTH_DEFINITIONS and depth:
        options["depth_definition"] = depth
    minimum = str(form.get(prefix + "base_count_minimum", "")).strip()
    if minimum.isdigit() and 1 <= int(minimum) <= 1000:
        options["base_count_minimum"] = int(minimum)
    if form.get(prefix + "layout_confirmed") in {"on", "1", "yes", True}:
        options["layout_confirmed"] = True
    if form.get(prefix + "sample_from_filename") in {"on", "1", "yes", True}:
        options["sample_from_filename"] = True
    header = str(form.get(prefix + "header_row", "")).strip()
    if header.isdigit() and 1 <= int(header) <= 50:
        options["header_row"] = int(header)
    return options


def shift_decimal(text, places):
    """`text` divided by 10**places, written exactly (no float rounding); None if not a number."""

    try:
        value = Decimal(str(text).strip())
    except InvalidOperation:
        return None
    if not value.is_finite():
        return None
    result = value.scaleb(-places)
    written = format(result, "f")
    if "." in written:
        written = written.rstrip("0").rstrip(".") or "0"
    return written


def numbers(series):
    raw = series.fillna("").astype(str).str.strip()
    return raw, pd.to_numeric(raw.replace("", pd.NA), errors="coerce")


def frequency_problems(table, column, options):
    """Problems and warnings for the frequency column under the chosen scale."""

    problems, warnings = [], []
    raw, values = numbers(table[column])
    values = values.dropna()
    if not len(values):
        return problems, warnings
    scale = options.get("frequency_scale")
    if scale == PERCENT:
        if values.lt(0).any() or values.gt(100).any():
            problems.append(f"Column '{column}' was declared as percentages, but some values are below 0 or "
                            "above 100.")
        elif values.max() <= 1:
            warnings.append(f"Column '{column}' was declared as percentages, but every value is 1 or less "
                            f"(largest {values.max():g}). Check that these really are percentages (0-100).")
    return problems, warnings


def reference_check(table, columns, reference):
    """Share of Ref bases matching the reference FASTA when positions are read as 1- or 0-based.

    `reference` is {"id", "sequence"}. Returns {"1": share, "0": share,
    "checked": n} or None when nothing can be checked.
    """

    if not reference or not {"Position", "Ref"} <= set(columns):
        return None
    sequence = reference["sequence"]
    positions = pd.to_numeric(table[columns["Position"]], errors="coerce")
    refs = table[columns["Ref"]].astype(str).str.strip().str.upper()
    usable = positions.notna() & positions.ge(0) & refs.isin(list("ACGT"))
    if "Reference_Genome" in columns:
        usable &= table[columns["Reference_Genome"]].astype(str).str.strip().eq(reference["id"])
    if not usable.any():
        return None
    positions, refs = positions[usable].astype(int).to_numpy(), refs[usable].to_numpy()
    result = {"checked": int(len(positions))}
    for base, shift in (("1", -1), ("0", 0)):
        index = positions + shift
        inside = (index >= 0) & (index < len(sequence))
        matches = sum(1 for i, r in zip(index[inside], refs[inside]) if sequence[i] == r)
        result[base] = matches / len(positions)
    return result


def numbering_warning(check, options, name):
    """A warning when the Ref bases fit the other position numbering much better."""

    if not check or check["checked"] < 5:
        return None
    chosen = options.get("position_base", "1")
    other = "0" if chosen == "1" else "1"
    if check[other] >= 0.95 and check[chosen] < 0.5:
        return (f"{name}: only {check[chosen]:.0%} of {check['checked']} reference bases match the reference "
                f"FASTA with {POSITION_BASES[chosen].split(' (')[0]} positions, but {check[other]:.0%} match with "
                f"{POSITION_BASES[other].split(' (')[0]} positions. Check the source's documentation and choose "
                "'Position numbering' accordingly; positions are never shifted without that choice.")
    return None


def interpret_table(table, options, has_reference_column):
    """Apply the confirmed options to a mapped (standard-name) table; returns (table, provenance rows).

    Original values are kept in Source_<field>_Original columns, so every
    change can be traced and reversed.
    """

    table = table.copy()
    rows = []
    if options.get("frequency_scale") == PERCENT and "Alt_Frequency" in table:
        original = table["Alt_Frequency"].fillna("").astype(str).str.strip()
        converted = [("" if text == "" else (shift_decimal(text, 2) or text)) for text in original]
        table["Source_Alt_Frequency_Original"] = original
        table["Alt_Frequency"] = converted
        rows.append(("Alt_Frequency", "Percent values divided by 100 (exact decimal); originals in "
                                      "Source_Alt_Frequency_Original"))
    if options.get("position_base") == "0" and "Position" in table:
        original = table["Position"].fillna("").astype(str).str.strip()
        shifted = pd.to_numeric(original, errors="coerce")
        table["Source_Position_Original"] = original
        table["Position"] = [str(int(value) + 1) if pd.notna(value) and float(value).is_integer() else text
                             for value, text in zip(shifted, original)]
        rows.append(("Position", "0-based positions converted to 1-based (+1); originals in "
                                 "Source_Position_Original"))
    if not has_reference_column and options.get("declared_reference"):
        table["Reference_Genome"] = options["declared_reference"]
        table["Source_Reference_Statement"] = "Stated by the user on the review page"
        rows.append(("Reference_Genome", f"'{options['declared_reference']}' stated by the user "
                                         "(the table has no reference column)"))
    if "Read_Depth" not in table:
        table["Read_Depth"] = ""
        rows.append(("Read_Depth", "No depth column: left blank (unknown), never filled in"))
    definition = options.get("depth_definition")
    if definition:
        table["Source_Depth_Definition"] = DEPTH_DEFINITIONS[definition]
        rows.append(("Read_Depth", f"Depth definition stated by the user: {DEPTH_DEFINITIONS[definition]}"))
    return table, rows
