"""Load, match, clean, and validate FASTA sequences and sample metadata."""

import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from genomic_regions import fasta_region, metadata_region, regions_conflict


# Required metadata after NCBI column names are standardized. Sequence
# comparison needs an identifier that matches the FASTA and a pathogen to group
# sequences for MAFFT; nothing else is required.
REQUIRED_METADATA_COLUMNS = [
    "Accession",
    "Pathogen",
]

# Context fields. When a column is absent or a value is blank, the sample is
# still compared genetically; only the temporal or geographic screening that
# needs the value is skipped for it. Nothing is filled in.
CONTEXT_METADATA_COLUMNS = ["Collection_Date", "Location_Name", "Host_Source"]

# Coordinates are useful, but they are not required.
OPTIONAL_METADATA_COLUMNS = ["Latitude", "Longitude"]

# Complete metadata layout used in summaries and downstream analysis.
METADATA_COLUMNS = [
    "Accession",
    "Pathogen",
    "Gene_Segment",
    "Collection_Date",
    "Location_Name",
    "Latitude",
    "Longitude",
    "Host_Source",
]

# Common column names in an NCBI Virus Results Table.
METADATA_ALIASES = {
    "accession": "Accession",
    "pathogen": "Pathogen",
    "organism_name": "Pathogen",
    "gene_segment": "Gene_Segment",
    "segment": "Segment",
    "gene_name": "Gene_Name",
    "gene": "Gene",
    "region": "Region",
    "collection_date": "Collection_Date",
    "location_name": "Location_Name",
    "geo_location": "Location_Name",
    "latitude": "Latitude",
    "longitude": "Longitude",
    "host_source": "Host_Source",
    "host": "Host_Source",
    # Completeness columns are read by the complete-genome safeguard, so they
    # must be recognized however their headers are spaced or capitalized.
    "nuc_completeness": "Nuc_Completeness",
    "completeness": "Completeness",
}

REGION_COLUMNS = ["Gene_Segment", "Segment", "Gene_Name", "Gene", "Region"]

# INSDC/NCBI placeholder words that mean "no value" (plus common variants).
MISSING_VALUE_WORDS = {
    "missing", "not collected", "not provided", "not applicable", "restricted access",
    "unknown", "na", "n/a", "none", "null", "-", "?",
}


# Exact DataFrame layout expected by the later analysis stages.
INTERNAL_COLUMNS = [
    "Sample_ID",
    "Pathogen",
    "Gene_Segment",
    "Gene_Segment_Source",
    "Sequence",
    "Collection_Date",
    "Location_Name",
    "Latitude",
    "Longitude",
    "Host_Source",
    "Accession",
    # Geographic level when known from country/state columns (see
    # genomics/location_hierarchy.py); country- and state-level places are
    # too coarse for the kilometer cutoff.
    "Location_Level",
    "Location_Name_Source",
]

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}
METADATA_EXTENSIONS = {".xlsx", ".csv"}
# IUPAC ambiguity codes (for example R = A or G) are accepted, but every stage
# treats them like N: they are never counted as a difference or a variant call.
AMBIGUOUS_NUCLEOTIDES = set("RYSWKMBDHV")
ALLOWED_SEQUENCE_CHARACTERS = set("ACGTN-") | AMBIGUOUS_NUCLEOTIDES


def is_missing_value(value):
    """True for blanks and placeholder words such as 'not provided'."""

    if value is None or pd.isna(value):
        return True
    text = " ".join(str(value).split()).casefold()
    return text == "" or text in MISSING_VALUE_WORDS


def blank_mask(series):
    """Return True for missing values and strings containing only whitespace."""

    text_values = series.astype("string")
    empty_text = text_values.str.strip().eq("").fillna(True)
    return series.isna() | empty_text


def metadata_rows(mask):
    """Convert a Boolean mask to row numbers from a metadata file."""

    rows = []
    for position, is_bad in enumerate(mask):
        if is_bad:
            # Add 2 because pandas starts at 0 and the file has a header row.
            rows.append(str(position + 2))
    return ", ".join(rows)


def clean_sequence(sequence):
    """Uppercase a sequence and remove whitespace without regular expressions."""

    pieces = str(sequence).split()
    return "".join(pieces).upper()


def normalize_accession(accession):
    """Remove one final numeric accession version for matching."""

    if pd.isna(accession):
        return ""

    accession = str(accession).strip()
    base, period, version = accession.rpartition(".")

    if period and base and version.isdigit():
        return base
    return accession


def standardize_metadata_columns(data):
    """Clean header spacing/case and recognize NCBI and region aliases."""

    columns_to_rename = {}
    standardized_names = set()

    for original_name in data.columns:
        spaced_name = str(original_name).replace("_", " ").replace("-", " ")
        normalized = "_".join(spaced_name.split()).casefold()
        standard_name = METADATA_ALIASES.get(normalized, original_name)
        if standard_name in standardized_names:
            raise ValueError(f"Duplicate metadata columns for {standard_name}.")
        standardized_names.add(standard_name)
        columns_to_rename[original_name] = standard_name

    return data.rename(columns=columns_to_rename)


def load_metadata(file_path):
    """Load metadata from CSV or the first worksheet of an XLSX file."""

    file_path = Path(file_path)
    extension = file_path.suffix.lower()

    if extension == ".xlsx":
        return pd.read_excel(
            file_path, engine="openpyxl", dtype="string", keep_default_na=False
        )
    if extension == ".csv":
        return pd.read_csv(file_path, dtype="string", keep_default_na=False)

    raise ValueError("Metadata must be an .xlsx or .csv file.")


def partial_date_mask(values):
    """True for dates given only as a year or a year and month.

    Such dates (for example "2020" or "2020-06") are kept as text in the
    source but are not turned into a day, because pandas would otherwise
    silently read "2020" as 1 January 2020.
    """

    text = values.astype("string").str.strip()
    year_only = text.str.fullmatch(r"\d{4}")
    year_month = text.str.fullmatch(r"\d{4}[-/]\d{1,2}")
    month_year = text.str.fullmatch(r"[A-Za-z]{3,9}[- ]\d{4}")
    return (year_only | year_month | month_year).fillna(False).astype(bool)


def ambiguous_day_month(values):
    """Mark numeric a/b/yyyy dates when the column also shows day-first evidence.

    pandas reads 03/04/2024 as 4 March (month first) but 13/04/2024 as 13 April
    (day first). When both kinds appear, the whole column is ambiguous, so every
    a/b/yyyy value with both parts 12 or less is reported.
    """

    text = values.astype("string").str.strip()
    parts = text.str.extract(r"^(\d{1,2})[/.-](\d{1,2})[/.-](\d{4})$")
    first = pd.to_numeric(parts[0], errors="coerce")
    second = pd.to_numeric(parts[1], errors="coerce")
    day_first = (first > 12).any()
    month_first = (second > 12).any()
    both_small = first.le(12) & second.le(12)
    if not day_first:
        # Month-first only (US style) is read consistently by pandas.
        return pd.Series(False, index=values.index)
    flagged = both_small | (first.notna() if month_first else False)
    return pd.Series(flagged, index=values.index).fillna(False).astype(bool)


def validate_metadata(data, warnings=None):
    """Clean and validate metadata without deleting any rows.

    `warnings`, when a list is given, receives notes about missing optional
    context (dates, locations, hosts) that limit some analyses but do not stop
    the genetic comparison.
    """

    validated = data.copy()
    validated = standardize_metadata_columns(validated)
    # Combine Geo_Location / Country / state columns into one location per row
    # (imported here because that module imports from this one).
    from genomics.location_hierarchy import apply_location_hierarchy
    validated = apply_location_hierarchy(validated)
    errors = []
    if warnings is None:
        warnings = []

    # Check the required headers.
    missing_columns = []
    for column in REQUIRED_METADATA_COLUMNS:
        if column not in validated.columns:
            missing_columns.append(column)

    if missing_columns:
        message = "Missing required metadata column(s): "
        errors.append(message + ", ".join(missing_columns))

    absent_columns = set()
    for column in CONTEXT_METADATA_COLUMNS:
        if column not in validated.columns:
            absent_columns.add(column)
            validated[column] = pd.Series(pd.NA, index=validated.index, dtype="string")
            warnings.append(
                f"No {column} column: analyses that need it are skipped; genetic "
                "comparison still runs."
            )

    # Placeholder words such as "not provided" or "missing" (INSDC/NCBI style)
    # mean "no value" in optional context columns; they are not places or dates.
    for column in CONTEXT_METADATA_COLUMNS + OPTIONAL_METADATA_COLUMNS:
        if column not in validated.columns or column in absent_columns:
            continue
        placeholders = validated[column].map(
            lambda value: is_missing_value(value) and pd.notna(value) and str(value).strip() != ""
        ).astype(bool)
        if placeholders.any():
            warnings.append(
                f"{column} uses a 'missing value' word (for example 'not provided') at "
                f"metadata row(s): {metadata_rows(placeholders)}; treated as missing."
            )
            validated.loc[placeholders, column] = pd.NA

    # Choose an explicit region from any supported metadata column, by priority.
    region_columns = [column for column in REGION_COLUMNS if column in validated]
    region_labels = []
    region_sources = []
    region_details = []

    for row_number, (_, row) in enumerate(validated.iterrows(), start=2):
        chosen = None
        candidates = []
        for column in region_columns:
            value = row[column]
            if pd.isna(value) or not str(value).strip():
                continue
            candidate = metadata_region(value, column)
            for earlier in candidates:
                if regions_conflict(earlier, candidate):
                    errors.append(
                        f"Conflicting metadata region values at row {row_number}: "
                        f"{earlier['label']} versus {candidate['label']}."
                    )
            candidates.append(candidate)
            if chosen is None:
                chosen = candidate
        region_labels.append(chosen["label"] if chosen else pd.NA)
        region_sources.append("Metadata" if chosen else pd.NA)
        region_details.append(candidates)

    validated["Gene_Segment"] = pd.array(region_labels, dtype="string")
    validated["Gene_Segment_Source"] = pd.array(region_sources, dtype="string")
    validated["_Region_Details"] = region_details

    # Count missing values for the printed summary.
    missing_counts = {}
    for column in METADATA_COLUMNS:
        if column in validated.columns:
            missing_counts[column] = int(blank_mask(validated[column]).sum())

    # Clean and validate the text fields.
    text_columns = [
        "Accession",
        "Pathogen",
        "Location_Name",
        "Host_Source",
    ]

    for column in text_columns:
        if column not in validated.columns:
            continue

        validated[column] = validated[column].astype("string").str.strip()
        missing_values = blank_mask(validated[column])

        if missing_values.any():
            row_text = metadata_rows(missing_values)
            if column in REQUIRED_METADATA_COLUMNS:
                errors.append(f"Missing {column} at metadata row(s): {row_text}")
            elif not missing_values.all():
                warnings.append(f"Missing {column} at metadata row(s): {row_text}")
            elif column not in absent_columns:
                warnings.append(f"{column} is blank in every metadata row.")
            validated.loc[missing_values, column] = pd.NA

    # Detect duplicate metadata accessions after removing numeric versions.
    if "Accession" in validated.columns:
        original_accessions = validated["Accession"]
        normalized_accessions = original_accessions.map(normalize_accession)

        duplicate_mask = normalized_accessions.duplicated(keep=False)
        duplicate_mask = duplicate_mask & normalized_accessions.ne("")
        duplicate_values = normalized_accessions[duplicate_mask].drop_duplicates()

        for accession in duplicate_values:
            matching_rows = normalized_accessions.eq(accession).fillna(False)
            row_text = metadata_rows(matching_rows)
            original_values = original_accessions[matching_rows].astype(str).tolist()
            original_text = ", ".join(original_values)

            errors.append(
                "Duplicate metadata accession after version normalization "
                f"{accession!r} at row(s): {row_text} "
                f"(original values: {original_text})"
            )

    # Convert collection dates.
    if "Collection_Date" in validated.columns:
        original_dates = validated["Collection_Date"].copy()
        missing_dates = blank_mask(original_dates)
        partial_dates = ~missing_dates & partial_date_mask(original_dates)
        ambiguous_dates = ambiguous_day_month(original_dates.where(~missing_dates))
        if ambiguous_dates.any():
            errors.append(
                "Collection_Date mixes day-first and month-first numeric dates (for example "
                "03/04/2024 and 13/04/2024) at metadata row(s): "
                f"{metadata_rows(ambiguous_dates)}. Use YYYY-MM-DD so days apart are not "
                "calculated with the wrong convention."
            )
        converted_dates = pd.to_datetime(
            original_dates.where(~partial_dates),
            errors="coerce",
            format="mixed",
        )
        invalid_dates = ~missing_dates & ~partial_dates & converted_dates.isna()
        validated["Collection_Date"] = converted_dates

        if missing_dates.any() and not missing_dates.all():
            row_text = metadata_rows(missing_dates)
            warnings.append(
                f"Missing Collection_Date at metadata row(s): {row_text}; those samples "
                "are compared genetically but not by date."
            )
        elif missing_dates.all() and len(missing_dates) and "Collection_Date" not in absent_columns:
            warnings.append("No collection dates were given; temporal screening is skipped.")
        if partial_dates.any():
            row_text = metadata_rows(partial_dates)
            warnings.append(
                f"Collection_Date gives only a year or month at metadata row(s): {row_text}; "
                "no day is assumed, so those samples are not used for day-level comparisons."
            )

        if invalid_dates.any():
            row_text = metadata_rows(invalid_dates)
            errors.append(f"Invalid Collection_Date at metadata row(s): {row_text}")

    # Convert coordinates only when their columns are available.
    coordinate_rules = {
        "Latitude": (-90, 90),
        "Longitude": (-180, 180),
    }

    for column, limits in coordinate_rules.items():
        if column not in validated.columns:
            continue

        minimum, maximum = limits
        original_values = validated[column].copy()
        missing_values = blank_mask(original_values)
        numeric_values = pd.to_numeric(original_values, errors="coerce")

        nonnumeric_values = ~missing_values & numeric_values.isna()
        out_of_range = numeric_values.notna() & ~numeric_values.between(
            minimum,
            maximum,
        )
        validated[column] = numeric_values

        if nonnumeric_values.any():
            row_text = metadata_rows(nonnumeric_values)
            errors.append(f"Non-numeric {column} at metadata row(s): {row_text}")

        if out_of_range.any():
            row_text = metadata_rows(out_of_range)
            errors.append(
                f"{column} outside {minimum} to {maximum} "
                f"at metadata row(s): {row_text}"
            )

    return validated, errors, missing_counts


def load_and_validate_metadata(file_path, warnings=None):
    """Load and validate one metadata file."""

    metadata = load_metadata(file_path)
    return validate_metadata(metadata, warnings)


def load_fasta(file_path):
    """Read FASTA record IDs and nucleotide sequences with Biopython."""

    rows = []

    # utf-8-sig removes a byte-order mark that some Windows editors add.
    with Path(file_path).open(encoding="utf-8-sig") as fasta_file:
        records = SeqIO.parse(fasta_file, "fasta")

        for record_number, record in enumerate(records, start=1):
            row = {
                "Accession": str(record.id).strip(),
                "Normalized_Accession": normalize_accession(record.id),
                "Description": str(record.description),
                "Sequence": clean_sequence(record.seq),
                "Record_Number": record_number,
            }
            rows.append(row)

    columns = [
        "Accession",
        "Normalized_Accession",
        "Description",
        "Sequence",
        "Record_Number",
    ]
    return pd.DataFrame(rows, columns=columns)


def match_fasta_to_metadata(
    metadata, fasta_records, default_gene_segment=None, confirm_same_region=False,
    complete_genomes=False,
):
    """Match FASTA records to metadata accessions and create analysis columns."""

    errors = []
    records = fasta_records.copy()

    if records.empty:
        errors.append("The FASTA file contains no sequence records.")

    # Check for blank FASTA identifiers.
    missing_accessions = blank_mask(records["Accession"])
    if missing_accessions.any():
        record_numbers = records.loc[missing_accessions, "Record_Number"]
        record_text = ", ".join(record_numbers.astype(str))
        errors.append("Missing FASTA accession in record(s): " + record_text)

    # Check for duplicate FASTA identifiers after version normalization.
    normalized_accessions = records["Normalized_Accession"]
    duplicate_mask = normalized_accessions.duplicated(keep=False)
    duplicate_mask = duplicate_mask & normalized_accessions.ne("")
    duplicate_values = normalized_accessions[duplicate_mask].drop_duplicates()

    for accession in duplicate_values:
        matching_records = normalized_accessions.eq(accession)
        record_numbers = records.loc[matching_records, "Record_Number"]
        original_values = records.loc[matching_records, "Accession"]

        record_text = ", ".join(record_numbers.astype(str))
        original_text = ", ".join(original_values.astype(str))

        errors.append(
            "Duplicate FASTA accession after version normalization "
            f"{accession!r} in record(s): {record_text} "
            f"(original values: {original_text})"
        )

    # Validate every FASTA sequence.
    for _, record in records.iterrows():
        accession = record["Accession"]
        sequence = record["Sequence"]
        record_number = record["Record_Number"]

        if sequence == "":
            errors.append(
                f"Missing sequence for FASTA accession {accession!r} "
                f"in record {record_number}."
            )
            continue

        invalid_characters = set(sequence) - ALLOWED_SEQUENCE_CHARACTERS
        if invalid_characters:
            character_text = "".join(sorted(invalid_characters))
            errors.append(
                f"Invalid sequence character(s) {character_text!r} for FASTA "
                f"accession {accession!r} in record {record_number}. "
                "Allowed characters are A, C, G, T, N, -, and the IUPAC "
                "ambiguity codes R, Y, S, W, K, M, B, D, H, and V."
            )

    # Make dictionaries that connect normalized IDs to their original IDs.
    metadata_accessions = {}
    if "Accession" in metadata.columns:
        for value in metadata["Accession"]:
            if pd.isna(value) or str(value) == "":
                continue
            normalized = normalize_accession(value)
            metadata_accessions[normalized] = str(value)

    fasta_accessions = {}
    for _, record in records.iterrows():
        normalized = str(record["Normalized_Accession"])
        if normalized != "":
            fasta_accessions[normalized] = str(record["Accession"])

    # Report IDs found in only one of the two input files.
    metadata_without_fasta = []
    for accession in metadata_accessions:
        if accession not in fasta_accessions:
            metadata_without_fasta.append(metadata_accessions[accession])
    metadata_without_fasta.sort()

    if metadata_without_fasta:
        errors.append(
            "Metadata accession(s) without a matching FASTA sequence: "
            + ", ".join(metadata_without_fasta)
        )

    fasta_without_metadata = []
    for accession in fasta_accessions:
        if accession not in metadata_accessions:
            fasta_without_metadata.append(fasta_accessions[accession])
    fasta_without_metadata.sort()

    if fasta_without_metadata:
        errors.append(
            "FASTA accession(s) without matching metadata: "
            + ", ".join(fasta_without_metadata)
        )

    # Make lookup tables used to add FASTA information to metadata rows.
    sequence_by_accession = {}
    fasta_id_by_accession = {}
    description_by_accession = {}

    for _, record in records.iterrows():
        accession = str(record["Normalized_Accession"])

        if accession and accession not in sequence_by_accession:
            sequence_by_accession[accession] = str(record["Sequence"])
            fasta_id_by_accession[accession] = str(record["Accession"])
            description_by_accession[accession] = str(record["Description"])

    matched = metadata.copy()

    if "Accession" in matched.columns:
        normalized_metadata = matched["Accession"].map(normalize_accession)
        matched["Sample_ID"] = normalized_metadata.map(fasta_id_by_accession)
        matched["Sequence"] = normalized_metadata.map(sequence_by_accession)
    else:
        matched["Sample_ID"] = pd.Series(pd.NA, index=matched.index, dtype="string")
        matched["Sequence"] = pd.Series(pd.NA, index=matched.index, dtype="string")

    for column in INTERNAL_COLUMNS:
        if column not in matched.columns:
            matched[column] = pd.NA

    # Resolve each region before MAFFT is allowed to group any sequences.
    pathogen_count = 0
    if "Pathogen" in matched:
        pathogen_count = matched["Pathogen"].dropna().nunique()
    if confirm_same_region and pathogen_count != 1:
        errors.append("Same-region confirmation requires one pathogen only.")
    default_info = None
    if complete_genomes:
        if pathogen_count != 1:
            errors.append(
                "Complete-genome mode requires all samples to have one pathogen."
            )
        default_info = metadata_region("Complete genome")
    elif default_gene_segment and str(default_gene_segment).strip():
        default_info = metadata_region(default_gene_segment)

    for index, sample in matched.iterrows():
        accession = sample.get("Accession")
        accession_key = normalize_accession(accession)
        description = description_by_accession.get(accession_key, "")
        header_info = fasta_region(description)

        if complete_genomes:
            for metadata_info in sample.get("_Region_Details", []):
                if not metadata_info["whole"]:
                    errors.append(
                        f"Complete-genome mode conflicts with metadata for "
                        f"{accession}: {metadata_info['label']}."
                    )
            if header_info["label"] and not header_info["whole"]:
                errors.append(
                    f"Complete-genome mode conflicts with the FASTA header for "
                    f"{accession}: {header_info['label']}."
                )
            for column in ("Nuc_Completeness", "Completeness"):
                value = sample.get(column)
                if pd.notna(value) and any(
                    word in str(value).casefold() for word in ("partial", "incomplete")
                ):
                    errors.append(
                        f"Complete-genome mode conflicts with {column}="
                        f"{value} for {accession}."
                    )
            annotation = description.partition(" ")[2].casefold()
            if not header_info["whole"] and (
                "partial" in annotation
                or "incomplete" in annotation
                or "complete cds" in annotation
            ):
                errors.append(
                    f"Complete-genome mode conflicts with the FASTA description "
                    f"for {accession}: it describes a partial sequence or CDS."
                )

        if header_info["whole"] and (header_info["segment"] or header_info["genes"]):
            errors.append(f"Contradictory FASTA region annotation for {accession}.")

        existing = sample.get("Gene_Segment")
        if (pd.notna(existing)
                and str(existing).strip().casefold() == "unspecified_region"):
            if header_info["label"]:
                matched.at[index, "Gene_Segment"] = header_info["label"]
                matched.at[index, "Gene_Segment_Source"] = "FASTA header"
            elif confirm_same_region and pathogen_count == 1:
                matched.at[index, "Gene_Segment"] = "Unspecified_Region"
                matched.at[index, "Gene_Segment_Source"] = (
                    "User-confirmed unspecified region"
                )
            else:
                errors.append(
                    "Unspecified_Region requires --confirm-same-region and "
                    "one pathogen; it cannot be accepted from metadata alone."
                )
            continue

        if pd.notna(existing) and str(existing).strip():
            metadata_infos = sample.get("_Region_Details", [])
            if not metadata_infos:
                metadata_infos = [metadata_region(existing)]
            for metadata_info in metadata_infos:
                if header_info["label"] and regions_conflict(metadata_info, header_info):
                    errors.append(
                        f"Metadata and FASTA header disagree on the genomic region "
                        f"for {accession}: {metadata_info['label']} versus "
                        f"{header_info['label']}."
                    )
                    break
            matched.at[index, "Gene_Segment_Source"] = "Metadata"
        elif header_info["label"]:
            matched.at[index, "Gene_Segment"] = header_info["label"]
            matched.at[index, "Gene_Segment_Source"] = "FASTA header"
        elif default_info:
            if default_info["label"].casefold() == "unspecified_region":
                if confirm_same_region and pathogen_count == 1:
                    matched.at[index, "Gene_Segment"] = "Unspecified_Region"
                    matched.at[index, "Gene_Segment_Source"] = (
                        "User-confirmed unspecified region"
                    )
                else:
                    errors.append(
                        "Unspecified_Region requires --confirm-same-region and "
                        "one pathogen; it cannot be used as an unconfirmed default."
                    )
            else:
                matched.at[index, "Gene_Segment"] = default_info["label"]
                matched.at[index, "Gene_Segment_Source"] = (
                    "User-confirmed complete genome" if complete_genomes
                    else "User override"
                )
        elif confirm_same_region and pathogen_count == 1:
            matched.at[index, "Gene_Segment"] = "Unspecified_Region"
            matched.at[index, "Gene_Segment_Source"] = (
                "User-confirmed unspecified region"
            )
        else:
            errors.append(
                f"Cannot determine genomic region for {accession}. Add a "
                "metadata region or explicit FASTA header annotation. "
                "Select Complete genomes only if all sequences are truly "
                "complete genomes of one pathogen."
            )

    # Even if metadata uses a vague common label, FASTA evidence must not
    # reveal different segments or genes inside one MAFFT group.
    for group_names, group in matched.groupby(["Pathogen", "Gene_Segment"]):
        header_infos = []
        for accession in group["Accession"]:
            description = description_by_accession.get(normalize_accession(accession), "")
            header_infos.append(fasta_region(description))
        segments = {info["segment"] for info in header_infos if info["segment"]}
        genes = [info["genes"] for info in header_infos if info["genes"]]
        whole_count = sum(info["whole"] for info in header_infos)
        different_genes = bool(genes) and not set.intersection(*genes)
        if (len(segments) > 1 or different_genes
                or (whole_count and whole_count != len(group)
                    and (segments or genes))):
            errors.append(
                f"FASTA headers show different genomic regions within "
                f"{group_names[0]} / {group_names[1]}; MAFFT will not align "
                "this mixed group."
            )

    # Use the predictable column order required by later stages.
    matched = matched[INTERNAL_COLUMNS]
    matched.attrs["fasta_record_count"] = len(records)
    return matched, errors


def version_mismatches(metadata, fasta_records):
    """List 'metadata vs FASTA' pairs whose numeric versions differ (X.2 vs X.3)."""

    if "Accession" not in metadata.columns:
        return []
    fasta_by_base = {normalize_accession(value): str(value) for value in fasta_records["Accession"]}
    mismatched = []
    for value in metadata["Accession"].dropna():
        original = str(value).strip()
        fasta_value = fasta_by_base.get(normalize_accession(original))
        if (fasta_value and original != fasta_value
                and normalize_accession(original) != original
                and normalize_accession(fasta_value) != fasta_value):
            mismatched.append(f"{original} vs {fasta_value}")
    return mismatched


def load_and_validate_inputs(
    fasta_path, metadata_path, default_gene_segment=None, confirm_same_region=False,
    complete_genomes=False, warnings=None,
):
    """Load both inputs and return the DataFrame used by analysis stages.

    Notes about missing optional context are appended to `warnings`.
    """

    metadata, metadata_errors, missing_counts = load_and_validate_metadata(
        metadata_path, warnings
    )
    fasta_records = load_fasta(fasta_path)
    matched, matching_errors = match_fasta_to_metadata(
        metadata, fasta_records, default_gene_segment, confirm_same_region,
        complete_genomes,
    )

    missing_counts["Sample_ID"] = int(blank_mask(matched["Sample_ID"]).sum())
    missing_counts["Sequence"] = int(blank_mask(matched["Sequence"]).sum())
    missing_counts["Gene_Segment"] = int(blank_mask(matched["Gene_Segment"]).sum())

    if warnings is not None:
        mismatched = version_mismatches(metadata, fasta_records)
        if mismatched:
            warnings.append(
                "Metadata and FASTA give different accession versions for: "
                + ", ".join(mismatched[:10])
                + ". They were matched by base accession; check that the metadata "
                "describes the same sequence version."
            )

    all_errors = metadata_errors + matching_errors
    return matched, all_errors, missing_counts


def category_summary(data, column):
    """Create a compact value-and-count summary for one text column."""

    if column not in data.columns:
        return "column not available"

    available_values = ~blank_mask(data[column])
    values = data.loc[available_values, column].astype("string").str.strip()
    counts = values.value_counts().sort_index()

    if counts.empty:
        return "none"

    entries = []
    for value, count in counts.items():
        entries.append(f"{value} ({count})")

    if len(entries) > 10:
        hidden_count = len(entries) - 10
        entries = entries[:10]
        entries.append(f"... and {hidden_count} more")

    return ", ".join(entries)


def print_summary(data, errors, missing_counts):
    """Print a compact overview and every validation problem."""

    print("\nDATASET SUMMARY")
    print(f"FASTA records: {data.attrs.get('fasta_record_count', 'unknown')}")
    print(f"Metadata rows: {len(data)}")
    print(f"Matched sample rows: {data['Sequence'].notna().sum()}")
    print(f"Pathogens represented: {category_summary(data, 'Pathogen')}")
    print(f"Gene/segment groups: {category_summary(data, 'Gene_Segment')}")
    print(f"Region label sources: {category_summary(data, 'Gene_Segment_Source')}")
    for _, sample in data.head(10).iterrows():
        print(
            f"  {sample['Sample_ID']}: {sample['Gene_Segment']} "
            f"({sample['Gene_Segment_Source']})"
        )
    if len(data) > 10:
        print(f"  ... and {len(data) - 10} more samples")

    valid_dates = data["Collection_Date"].dropna()
    if valid_dates.empty:
        date_range = "no valid dates"
    else:
        first_date = valid_dates.min().date().isoformat()
        last_date = valid_dates.max().date().isoformat()
        date_range = f"{first_date} to {last_date}"

    print(f"Date range: {date_range}")
    print(f"Hosts: {category_summary(data, 'Host_Source')}")
    print(f"Locations: {category_summary(data, 'Location_Name')}")

    print("Missing values:")
    summary_columns = METADATA_COLUMNS + ["Sequence"]
    for column in summary_columns:
        if column in missing_counts:
            print(f"  {column}: {missing_counts[column]}")
        elif column in METADATA_COLUMNS:
            print(f"  {column}: column absent")

    print(f"Validation errors: {len(errors)}")
    if errors:
        for error in errors:
            print(f"  - {error}")
    else:
        print("  None")


def parse_arguments():
    """Read FASTA and metadata paths from the command line."""

    parser = argparse.ArgumentParser(
        description="Validate and match FASTA sequences with XLSX or CSV metadata."
    )
    parser.add_argument("fasta_file", type=Path)
    parser.add_argument("metadata_file", type=Path)
    parser.add_argument("--default-gene-segment", help="Region for unannotated samples")
    parser.add_argument(
        "--confirm-same-region", action="store_true",
        help="Confirm all sequences represent the same region (one pathogen only)",
    )
    return parser.parse_args()


def main():
    """Run the dual-input validator and return a process status code."""

    arguments = parse_arguments()
    fasta_path = arguments.fasta_file
    metadata_path = arguments.metadata_file

    if fasta_path.suffix.lower() not in FASTA_EXTENSIONS:
        print("ERROR: FASTA must use .fasta, .fa, or .fna.", file=sys.stderr)
        return 2

    if metadata_path.suffix.lower() not in METADATA_EXTENSIONS:
        print("ERROR: Metadata must use .xlsx or .csv.", file=sys.stderr)
        return 2

    if not fasta_path.is_file():
        print(f"ERROR: FASTA file not found: {fasta_path}", file=sys.stderr)
        return 2

    if not metadata_path.is_file():
        print(f"ERROR: Metadata file not found: {metadata_path}", file=sys.stderr)
        return 2

    try:
        data, errors, missing_counts = load_and_validate_inputs(
            fasta_path,
            metadata_path,
            arguments.default_gene_segment,
            arguments.confirm_same_region,
        )
    except Exception as error:
        print(f"ERROR: Could not read the input files: {error}", file=sys.stderr)
        return 2

    print(f"FASTA file: {fasta_path.resolve()}")
    print(f"Metadata file: {metadata_path.resolve()}")
    print_summary(data, errors, missing_counts)

    if errors:
        print("\nRESULT: FAILED validation. No metadata rows were deleted.")
        return 1

    print("\nRESULT: PASSED validation. No metadata rows were deleted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
