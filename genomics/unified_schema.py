"""Shared sample and variant tables used by both kinds of genomic evidence.

Individual sequences and pooled samples are different measurements:

- An individual (consensus) sequence gives one base call per position.
- A pooled sample gives an allele *frequency* estimated from many reads.

Both are stored in the same table layout so they can be compared at the same
reference positions, but the `Observation_Type` column keeps them distinct.

Missing information is never filled in:

- A blank `Alt_Frequency` means "not measured", which is different from 0.
- A blank `Read_Depth` means "not reported", which is different from 0 reads.
- A variant that is not listed for a sample is not proof that the sample has
  the reference base; that depends on coverage at that position.
"""

import pandas as pd


INDIVIDUAL = "Individual sequence"
POOLED = "Pooled sample"
SAMPLE_TYPES = [INDIVIDUAL, POOLED]

CONSENSUS_CALL = "Consensus base call"
FREQUENCY_CALL = "Allele frequency"

NUCLEOTIDES = ["A", "C", "G", "T"]

# Largest absolute difference between a reported pooled frequency and
# alternate reads / depth that is treated as rounding or a slightly different
# denominator (audit S3-06). Larger differences are flagged per record.
FREQUENCY_TOLERANCE = 0.01

SAMPLE_COLUMNS = [
    "Sample_ID",
    "Sample_Type",
    "Collection_Date",
    "Location_Name",
    "Latitude",
    "Longitude",
    "Coordinate_Source",
    "Geographic_Precision",
    "Reference_ID",
    "Genomic_Region",
    "Data_Source",
    "Coverage_Information",
    "Genome_Coverage_Percent",
    "Mean_Depth",
    "Quality_Status",
    "Variant_Profile_Status",
]

VARIANT_COLUMNS = [
    "Variant_Key",
    "Reference_ID",
    "Position",
    "Ref",
    "Alt",
    "Sample_ID",
    "Sample_Type",
    "Observation_Type",
    "Alt_Frequency",
    "Read_Depth",
    "Alt_Depth",
    "Callable_Status",
    "Quality_Status",
    "Collection_Date",
    "Location_Name",
    "Data_Source",
]

# Words used for one sample at one reference position.
VARIANT_OBSERVED = "Variant observed"
REFERENCE_OBSERVED = "Reference base called"  # individual sequences only
COVERED_NOT_DETECTED = "Adequate depth; variant not reported"  # pooled only
LOW_COVERAGE = "Low coverage"
NOT_COVERED = "Not covered or not callable"
COVERAGE_UNKNOWN = "Coverage unknown"
EXCLUDED_BY_QC = "Reported but excluded by QC"  # pooled only


def variant_key(reference_id, position, ref, alt):
    """Build the stable identifier Reference_ID:Position:Ref>Alt."""

    return f"{str(reference_id).strip()}:{int(position)}:{ref}>{alt}"


def make_variant_table(rows):
    """Create a variant table with the shared columns and consistent types.

    Extra columns, such as annotations published with a source table, are
    kept after the shared columns.
    """

    table = pd.DataFrame(list(rows))
    for column in VARIANT_COLUMNS:
        if column not in table.columns:
            table[column] = pd.NA
    extra_columns = [name for name in table.columns if name not in VARIANT_COLUMNS]
    table = table[VARIANT_COLUMNS + extra_columns]
    table["Position"] = pd.array(table["Position"], dtype="Int64")
    table["Alt_Frequency"] = pd.to_numeric(table["Alt_Frequency"], errors="coerce")
    table["Read_Depth"] = pd.array(
        pd.to_numeric(table["Read_Depth"], errors="coerce"), dtype="Int64"
    )
    table["Alt_Depth"] = pd.array(
        pd.to_numeric(table["Alt_Depth"], errors="coerce"), dtype="Int64"
    )
    table["Collection_Date"] = pd.to_datetime(table["Collection_Date"], errors="coerce")
    return table


def make_sample_table(rows):
    """Create a sample table with the shared columns."""

    table = pd.DataFrame(rows, columns=SAMPLE_COLUMNS)
    table["Collection_Date"] = pd.to_datetime(table["Collection_Date"], errors="coerce")
    for column in ["Latitude", "Longitude", "Genome_Coverage_Percent", "Mean_Depth"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    return table


def validate_variant_table(table):
    """Return a list of problems; an empty list means the table is usable."""

    errors = []
    missing = [column for column in VARIANT_COLUMNS if column not in table.columns]
    if missing:
        return ["Variant table is missing column(s): " + ", ".join(missing)]

    for row_number, row in enumerate(table.itertuples(index=False), start=1):
        label = f"Variant row {row_number}"
        if pd.isna(row.Position) or int(row.Position) < 1:
            errors.append(f"{label}: positions must be 1-based positive integers.")
            continue
        if row.Ref not in NUCLEOTIDES or row.Alt not in NUCLEOTIDES:
            errors.append(f"{label}: Ref and Alt must each be one of A, C, G, T.")
        elif row.Ref == row.Alt:
            errors.append(f"{label}: Ref and Alt are identical.")
        if row.Sample_Type not in SAMPLE_TYPES:
            errors.append(f"{label}: unknown Sample_Type {row.Sample_Type!r}.")
        if pd.notna(row.Alt_Frequency) and not 0 <= row.Alt_Frequency <= 1:
            errors.append(f"{label}: Alt_Frequency must be a fraction from 0 to 1.")
        for name in ["Read_Depth", "Alt_Depth"]:
            value = getattr(row, name)
            if pd.notna(value) and int(value) < 0:
                errors.append(f"{label}: {name} cannot be negative.")
        if (pd.notna(row.Read_Depth) and pd.notna(row.Alt_Depth)
                and int(row.Alt_Depth) > int(row.Read_Depth)):
            errors.append(f"{label}: Alt_Depth is larger than Read_Depth.")
        expected_key = variant_key(row.Reference_ID, row.Position, row.Ref, row.Alt)
        if row.Variant_Key != expected_key:
            errors.append(f"{label}: Variant_Key should be {expected_key}.")
        if row.Sample_Type == INDIVIDUAL and pd.notna(row.Alt_Frequency):
            errors.append(
                f"{label}: individual consensus calls do not have an allele frequency."
            )

    return errors


def depth_status(depth, minimum_depth):
    """Describe one pooled position from its read depth without guessing."""

    if depth is None or pd.isna(depth):
        return COVERAGE_UNKNOWN
    if depth <= 0:
        return NOT_COVERED
    if depth < minimum_depth:
        return LOW_COVERAGE
    return COVERED_NOT_DETECTED
