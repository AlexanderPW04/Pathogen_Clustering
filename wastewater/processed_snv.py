"""Load processed SNV calls and summarize what was observed in wastewater."""

from pathlib import Path

import pandas as pd

from genomics.collection_dates import (
    AUTOMATIC,
    INVALID,
    NOT_SUPPLIED,
    date_warnings,
    parse_collection_dates,
)
from genomics.column_detection import explain_snv_columns
from genomics.excel_tables import read_sheet
from genomics.file_formats import VCF_FORMAT, vcf_notes, vcf_to_snv_table
from genomics.unified_schema import FREQUENCY_TOLERANCE


# Supplementary Table 2 in Fontenele et al. (2021) uses these exact headers.
PAPER_COLUMNS = {
    "SAMPLE": "Sample_ID",
    "COLLECTION_DATE": "Collection_Date",
    "POSITION": "Position",
    "REF": "Ref",
    "ALT": "Alt",
    "ALT_FREQ": "Alt_Frequency",
    "TOTAL_DP": "Read_Depth",
}
REQUIRED_COLUMNS = [
    "Sample_ID", "Reference_Genome", "Position", "Ref", "Alt",
    "Alt_Frequency", "Read_Depth",
]
OUTPUT_COLUMNS = REQUIRED_COLUMNS + ["Collection_Date", "Location_Name"]

# Annotation columns published with Supplementary Table 2. They are kept
# exactly as published under a "Published_" name so they are never confused
# with values calculated by this application.
PUBLISHED_ANNOTATIONS = {
    "REF_CODON": "Published_Ref_Codon",
    "REF_AA": "Published_Ref_AA",
    "ALT_CODON": "Published_Alt_Codon",
    "ALT_AA": "Published_Alt_AA",
    "BIN": "Published_Bin",
    "GLOBAL_FREQUENCY": "Published_Global_Frequency",
    "USA_FREQUENCY": "Published_USA_Frequency",
    "Sy/Ns": "Published_Syn_NonSyn",
    "Annotation": "Published_Annotation",
    "Novel SNVs based on Global frequency": "Published_Novel_Global",
    "Identified in clinical derived-genomes after June 17th":
        "Published_In_Clinical_After_2020_06_17",
}

# Columns written by `ivar variants`; a SAMPLE (or Sample_ID) column must be
# added when several samples are combined into one file.
IVAR_COLUMNS = {
    "SAMPLE": "Sample_ID",
    "REGION": "Reference_Genome",
    "POS": "Position",
    "REF": "Ref",
    "ALT": "Alt",
    "ALT_FREQ": "Alt_Frequency",
    "TOTAL_DP": "Read_Depth",
    "ALT_DP": "Alt_Depth",
    "REF_DP": "Ref_Depth",
    "PVAL": "Source_Fisher_P",
    "PASS": "Source_Caller_Pass",
    "ALT_QUAL": "Source_Alt_Mean_Quality",
    "GFF_FEATURE": "Source_Gene_Feature",
    "REF_CODON": "Source_Ref_Codon",
    "REF_AA": "Source_Ref_AA",
    "ALT_CODON": "Source_Alt_Codon",
    "ALT_AA": "Source_Alt_AA",
}
IVAR_REQUIRED = {"REGION", "POS", "REF", "ALT", "ALT_FREQ", "TOTAL_DP"}
# Header spellings (letters and digits only, lower case) accepted as the
# sample column of a combined iVar table.
SAMPLE_HEADERS = {"sample", "sampleid", "samplename"}

PAPER_FORMAT = "Fontenele et al. (2021) Supplementary Table 2"
IVAR_FORMAT = "iVar variants TSV"
STANDARD_FORMAT = "Standardized processed SNV table"
VCF_SNV_FORMAT = "VCF (single-nucleotide entries converted to the standardized table)"
COMBINED_FORMAT = "Combined standardized table (some records converted from VCF)"


def read_table(path, sheet=None):
    """Read a CSV, TSV, Excel, or VCF file without losing numeric-looking sample IDs."""

    path = Path(path)
    if path.suffix.lower() == ".vcf":
        data, details = vcf_to_snv_table(path, path.stem)
        data = data.astype(str)
        data.attrs["vcf_notes"] = vcf_notes(details)
        return data
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    if path.suffix.lower() == ".tsv":
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if path.suffix.lower() == ".xlsx":
        # The paper's workbook has a descriptive first row, then its headers.
        data = read_sheet(path, sheet)
        description = data.attrs.get("preamble", "")
        data.attrs["fontenele_2021_table2"] = (
            data.attrs.get("header_row") == 1
            and "Supplementary Table 2" in description
            and "52 wastewater samples" in description
            and "SARS-CoV-2" in description
        )
        return data
    raise ValueError("Use a processed SNV .xlsx, .csv, .tsv, or SNV .vcf file.")


def significant_digits(text):
    """Significant digits written in a number such as '0.0485437' (6) or '1' (1); None if not a number."""

    text = str(text).strip().lower().lstrip("+-")
    mantissa = text.split("e")[0]
    if not mantissa or not all(c.isdigit() or c == "." for c in mantissa) or mantissa.count(".") > 1:
        return None
    digits = mantissa.replace(".", "").lstrip("0")
    return len(digits) if digits else 1


def load_processed_snvs(path, metadata_path=None, date_order=AUTOMATIC):
    """Return standardized SNV rows and visible data-quality warnings.

    Collection dates are read by `genomics.collection_dates`: ambiguous
    day/month values are kept as text and not used unless `date_order` is
    "day_first" or "month_first". The original text and a status are kept in
    Collection_Date_Original and Collection_Date_Status.
    """

    data = read_table(path)
    conversion_notes = list(data.attrs.get("vcf_notes", []))
    if conversion_notes and data.empty:
        raise ValueError("The VCF has no single-nucleotide variant entries that can be analyzed. "
                         + " ".join(conversion_notes))
    data.columns = [str(name).strip() for name in data.columns]
    if IVAR_REQUIRED.issubset(data.columns) and not {"SAMPLE", "Sample_ID"} & set(data.columns):
        sample_columns = [name for name in data.columns
                          if "".join(c for c in name.casefold() if c.isalnum()) in SAMPLE_HEADERS]
        if len(sample_columns) == 1:
            data = data.rename(columns={sample_columns[0]: "SAMPLE"})
    published = (
        data.attrs.get("fontenele_2021_table2", False)
        and set(PAPER_COLUMNS).issubset(data.columns)
        and {"STATE", "LOCATION"}.issubset(data.columns)
    )
    is_ivar = (
        not published
        and IVAR_REQUIRED.issubset(data.columns)
        and ("SAMPLE" in data.columns or "Sample_ID" in data.columns)
    )
    warnings = conversion_notes

    if published:
        input_format = PAPER_FORMAT
        data = data.rename(columns=PAPER_COLUMNS)
        data = data.rename(columns=PUBLISHED_ANNOTATIONS)
        data["Reference_Genome"] = "NC_045512.2"
        data["Location_Name"] = data["STATE"].astype(str).str.strip() + ": " + data["LOCATION"].astype(str).str.strip()
        date_info = parse_collection_dates(data["Collection_Date"], fixed_format="%d_%b_%Y")
    elif is_ivar:
        input_format = IVAR_FORMAT
        data = data.rename(columns=IVAR_COLUMNS)
        # iVar also reports insertions (+bases) and deletions (-bases). They are
        # kept out of the SNV analysis and counted visibly instead.
        indels = data["Alt"].astype(str).str.strip().str[:1].isin(["+", "-"])
        if indels.any():
            warnings.append(
                f"{int(indels.sum())} iVar insertion/deletion record(s) were set aside; "
                "only single-nucleotide variants are analyzed."
            )
            data = data.loc[~indels].copy()
        for column in ["Collection_Date", "Location_Name"]:
            if column not in data:
                data[column] = ""
        date_info = parse_collection_dates(data["Collection_Date"], date_order)
    else:
        input_format = STANDARD_FORMAT
        formats = set(data.get("Source_Input_Format", pd.Series(dtype=str)).fillna("").astype(str))
        if formats == {VCF_FORMAT}:
            input_format = VCF_SNV_FORMAT
        elif VCF_FORMAT in formats:
            input_format = COMBINED_FORMAT
        missing = sorted(set(REQUIRED_COLUMNS) - set(data.columns))
        if missing:
            raise ValueError(explain_snv_columns(data.columns)
                             or "Missing processed SNV column(s): " + ", ".join(missing))
        if "Collection_Date" not in data:
            data["Collection_Date"] = ""
        if "Location_Name" not in data:
            data["Location_Name"] = ""
        date_info = parse_collection_dates(data["Collection_Date"], date_order)
    invalid_source_dates = date_info["Status"].eq(INVALID)

    if data.empty:
        raise ValueError("The processed SNV table has no records.")

    for column in ["Sample_ID", "Reference_Genome", "Ref", "Alt", "Location_Name"]:
        data[column] = data[column].fillna("").astype(str).str.strip()

    data["Ref"] = data["Ref"].str.upper()
    data["Alt"] = data["Alt"].str.upper()
    errors = []
    # Row numbers as the user sees them in the file (headers occupy the rows above).
    first_row = data.attrs.get("header_row", 0) + 2

    def where(mask, values):
        """The first offending cell, so a large table can actually be corrected."""

        position = int(mask.to_numpy().argmax())
        return f" (first at row {position + first_row}: '{str(values.iloc[position])[:30]}')"

    if metadata_path:
        metadata = read_table(metadata_path)
        metadata.columns = [str(name).strip() for name in metadata.columns]
        if "Sample_ID" not in metadata:
            raise ValueError("Metadata needs a Sample_ID column for joining.")
        metadata["Sample_ID"] = metadata["Sample_ID"].fillna("").astype(str).str.strip()
        if metadata["Sample_ID"].eq("").any() or metadata["Sample_ID"].duplicated().any():
            raise ValueError("Metadata Sample_ID values must be present and unique.")
        unknown = sorted(set(data["Sample_ID"]) - set(metadata["Sample_ID"]))
        if unknown:
            raise ValueError("SNV sample(s) missing from metadata: " + ", ".join(unknown[:10]))
        unused = set(metadata["Sample_ID"]) - set(data["Sample_ID"])
        if unused:
            warnings.append(f"{len(unused)} metadata sample(s) have no SNV records; no mutation absence is inferred.")
        metadata = metadata.set_index("Sample_ID")
        for column in ["Collection_Date", "Location_Name"]:
            if column not in metadata:
                continue
            matched = data["Sample_ID"].map(metadata[column]).fillna("").astype(str).str.strip()
            if column == "Collection_Date":
                matched_info = parse_collection_dates(matched, date_order)
                bad = matched_info["Status"].eq(INVALID)
                if bad.any():
                    example = matched[bad.to_numpy()].iloc[0]
                    errors.append(f"{bad.sum()} record(s) get a metadata collection date that cannot be read "
                                  f"as a date (for example '{example[:30]}'; a full date such as 2020-05-19 "
                                  "is needed, and a missing year is never guessed).")
                dates = date_info["Date"]
                conflict = dates.notna() & matched_info["Date"].notna() & dates.ne(matched_info["Date"])
                if conflict.any():
                    errors.append(f"{conflict.sum()} SNV date(s) conflict with metadata.")
                # The SNV table's own value wins; metadata fills only blank SNV dates.
                use_metadata = date_info["Status"].eq(NOT_SUPPLIED) & matched_info["Status"].ne(NOT_SUPPLIED)
                date_info.loc[use_metadata] = matched_info.loc[use_metadata]
            else:
                conflict = data[column].ne("") & matched.ne("") & data[column].str.casefold().ne(matched.str.casefold())
                if conflict.any():
                    errors.append(f"{conflict.sum()} SNV location(s) conflict with metadata.")
                data[column] = data[column].where(data[column].ne(""), matched)

    # Optional Specimen_ID (Sprint 5 Phase 2) and Library_ID (Phase 3): samples that
    # share a Specimen_ID and have distinct Library_IDs are documented technical
    # replicates (separately prepared libraries of one specimen). Metadata fills
    # blank values; each sample must have one value of each.
    for name in ["Specimen_ID", "Library_ID"]:
        if metadata_path and name in metadata:
            from_metadata = data["Sample_ID"].map(metadata[name]).fillna("").astype(str).str.strip()
            if name in data:
                own = data[name].fillna("").astype(str).str.strip()
                if (own.ne("") & from_metadata.ne("") & own.ne(from_metadata)).any():
                    errors.append(f"Some SNV {name} values conflict with metadata.")
                data[name] = own.where(own.ne(""), from_metadata)
            else:
                data[name] = from_metadata
        if name in data:
            data[name] = data[name].fillna("").astype(str).str.strip()
            several = data.groupby("Sample_ID")[name].nunique().gt(1)
            if several.any():
                errors.append(f"Each sample needs one {name}; these have several: "
                              + ", ".join(several[several].index[:10]) + ".")

    if invalid_source_dates.any():
        errors.append(f"{invalid_source_dates.sum()} collection date(s) are invalid"
                      + where(invalid_source_dates, data["Collection_Date"].fillna("").astype(str)) + ".")

    for column in ["Sample_ID", "Reference_Genome", "Ref", "Alt"]:
        count = int(data[column].eq("").sum())
        if count:
            errors.append(f"{count} record(s) have no {column}.")

    for column in ["Ref", "Alt"]:
        invalid = ~data[column].isin(["A", "C", "G", "T"])
        if invalid.any():
            errors.append(f"{invalid.sum()} record(s) have invalid {column}; expected one A, C, G, or T"
                          + where(invalid, data[column]) + ".")
    same_base = data["Ref"].eq(data["Alt"])
    if same_base.any():
        errors.append(f"{same_base.sum()} record(s) have identical reference and alternative bases.")

    positions = pd.to_numeric(data["Position"], errors="coerce")
    invalid_position = positions.isna() | positions.le(0) | positions.ne(positions.round()) | positions.isin([float("inf"), float("-inf")])
    if invalid_position.any():
        errors.append(f"{invalid_position.sum()} position(s) are not positive integers"
                      + where(invalid_position, data["Position"]) + ".")
    frequencies = pd.to_numeric(data["Alt_Frequency"], errors="coerce")
    invalid_frequency = frequencies.isna() | ~frequencies.between(0, 1)
    if invalid_frequency.any():
        errors.append(f"{invalid_frequency.sum()} alternative allele frequency value(s) are outside 0 to 1 "
                      "or missing" + where(invalid_frequency, data["Alt_Frequency"]) + ".")
    raw_depths = data["Read_Depth"].fillna("").astype(str).str.strip()
    depths = pd.to_numeric(raw_depths.replace("", pd.NA), errors="coerce")
    invalid_depth_text = raw_depths.ne("") & depths.isna()
    if invalid_depth_text.any():
        errors.append(f"{invalid_depth_text.sum()} read depth value(s) are not numeric"
                      + where(invalid_depth_text, raw_depths) + ".")
    bad_depth = depths.notna() & (depths.lt(1) | depths.ne(depths.round()) | depths.isin([float("inf"), float("-inf")]))
    if bad_depth.any():
        errors.append(f"{bad_depth.sum()} read depth value(s) are not positive integers.")
    depth_missing = int(depths.isna().sum())
    if depth_missing:
        warnings.append(f"{depth_missing} record(s) have no usable read depth.")

    alt_depths = None
    if "Alt_Depth" in data:
        raw_alt = data["Alt_Depth"].fillna("").astype(str).str.strip()
        alt_depths = pd.to_numeric(raw_alt.replace("", pd.NA), errors="coerce")
        bad_alt = (raw_alt.ne("") & alt_depths.isna()) | (
            alt_depths.notna() & (alt_depths.lt(0) | alt_depths.ne(alt_depths.round()))
        )
        if bad_alt.any():
            errors.append(f"{bad_alt.sum()} alternate read count value(s) are not whole numbers.")
        elif (alt_depths.notna() & depths.notna() & alt_depths.gt(depths)).any():
            errors.append("Some alternate read counts are larger than the total read depth.")

    if errors:
        raise ValueError("Processed SNV validation failed: " + " ".join(errors))

    # A row reporting 0% (or 0 alternate reads) is not a detected variant. It is
    # set aside rather than analyzed as an observation, and it is not used as
    # coverage evidence either (coverage comes only from a per-position table).
    zero = frequencies.eq(0)
    if alt_depths is not None:
        zero = zero | alt_depths.eq(0).fillna(False)
    if zero.all():
        raise ValueError("Every record reports a 0% frequency or 0 alternate reads, so the "
                         "table contains no detected variants.")
    if zero.any():
        warnings.append(f"{int(zero.sum())} record(s) report a 0% frequency or 0 alternate reads; "
                        "they were set aside because they are not variant calls.")
        data, frequencies, depths = data.loc[~zero].copy(), frequencies[~zero], depths[~zero]
        positions, date_info = positions[~zero], date_info.loc[~zero]
        if alt_depths is not None:
            alt_depths = alt_depths[~zero]

    # IMP-9 / audit S3-06: reported frequencies should agree with the read counts
    # given. Both values are kept; each record's difference is in the variant
    # table (Frequency_Check), and no value is overwritten.
    if alt_depths is not None:
        implied = alt_depths / depths
        disagree = implied.notna() & (implied - frequencies).abs().gt(FREQUENCY_TOLERANCE)
        if disagree.any():
            warnings.append(f"{int(disagree.sum())} record(s) report a frequency that differs by more "
                            f"than {FREQUENCY_TOLERANCE} from alternate reads / depth. Both values are kept; "
                            "Sprint 3 analyses use the reported frequency, and the depth-aware analysis "
                            "does not calculate an uncertainty interval for these records.")

    data["Position"] = positions.astype(int)
    data["Alt_Frequency"] = frequencies.astype(float)
    data["Read_Depth"] = depths.astype("Int64")
    data["Collection_Date"] = date_info["Date"]
    data["Collection_Date_Original"] = date_info["Original"]
    data["Collection_Date_Status"] = date_info["Status"]
    data["Frequency_Significant_Digits"] = pd.array(
        [significant_digits(text) for text in data["Alt_Frequency"].astype(str)], dtype="Int64")
    for column in ["Collection_Date", "Location_Name"]:
        values = data[column].replace("", pd.NA) if column == "Location_Name" else data[column]
        conflicting_samples = values.groupby(data["Sample_ID"]).nunique(dropna=True).gt(1)
        if conflicting_samples.any():
            raise ValueError(f"Conflicting {column} values within sample(s): " + ", ".join(conflicting_samples[conflicting_samples].index[:10]))
        data[column] = values.groupby(data["Sample_ID"]).transform("first")
    data["Location_Name"] = data["Location_Name"].fillna("")

    keys = ["Sample_ID", "Reference_Genome", "Position", "Ref", "Alt"]
    duplicate_count = int(data.duplicated(keys).sum())
    if duplicate_count:
        conflicting = int((data.groupby(keys)["Alt_Frequency"].nunique() > 1).sum())
        warnings.append(f"{duplicate_count} repeated sample/reference/position/allele record(s) found. "
                        "A pooled sample has one frequency per allele and position, so each is counted "
                        "once (the first row); all rows stay in the downloads with a Duplicate_Status"
                        + (f" ({conflicting} report a different frequency and are flagged as conflicting)."
                           if conflicting else "."))
    warnings += date_warnings(date_info, data["Sample_ID"])
    missing_dates = int(data["Collection_Date_Status"].eq(NOT_SUPPLIED).sum())
    missing_locations = int(data["Location_Name"].eq("").sum())
    if missing_dates:
        warnings.append(f"{missing_dates} record(s) lack a collection date.")
    if missing_locations:
        warnings.append(f"{missing_locations} record(s) lack a location.")

    # Keep source-provided extras: alternate read counts and any column whose
    # name starts with Published_ (paper annotations) or Source_ (caller output).
    extra_columns = ["Collection_Date_Original", "Collection_Date_Status", "Frequency_Significant_Digits"]
    if alt_depths is not None:
        data["Alt_Depth"] = alt_depths.astype("Int64")
        extra_columns.append("Alt_Depth")
    if "Ref_Depth" in data:
        raw_ref = data["Ref_Depth"].fillna("").astype(str).str.strip()
        ref_depths = pd.to_numeric(raw_ref.replace("", pd.NA), errors="coerce")
        bad_ref = (raw_ref.ne("") & ref_depths.isna()) | (
            ref_depths.notna() & (ref_depths.lt(0) | ref_depths.ne(ref_depths.round())))
        if bad_ref.any() or (ref_depths.notna() & depths.notna() & ref_depths.gt(depths)).any():
            raise ValueError("Processed SNV validation failed: reference read counts (Ref_Depth) must be "
                             "whole numbers no larger than the total read depth.")
        data["Ref_Depth"] = ref_depths.astype("Int64")
        extra_columns.append("Ref_Depth")
    for name in ["Specimen_ID", "Library_ID"]:
        if name in data:
            extra_columns.append(name)
    for name in data.columns:
        if name.startswith(("Published_", "Source_")) and name not in extra_columns:
            extra_columns.append(name)

    result = data[OUTPUT_COLUMNS + extra_columns].copy()
    result.attrs["input_format"] = input_format
    return result, warnings


def summarize_snvs(data):
    """Count observed SNV calls; never treat an unlisted call as absent."""

    variant = ["Reference_Genome", "Position", "Ref", "Alt"]
    top = (
        data.groupby(variant, as_index=False)
        .agg(Samples_Observed=("Sample_ID", "nunique"), Records=("Sample_ID", "size"))
        .sort_values(["Samples_Observed", "Records", "Position"], ascending=[False, False, True])
        .head(10)
    )
    sample_table = (
        data.groupby("Sample_ID", as_index=False)
        .agg(
            Collection_Date=("Collection_Date", "first"),
            Location_Name=("Location_Name", "first"),
            SNV_Records=("Position", "size"),
        )
        .sort_values("Sample_ID")
    )
    observed = data.drop_duplicates(["Sample_ID"] + variant)
    sample_table["Distinct_SNVs"] = sample_table["Sample_ID"].map(observed.groupby("Sample_ID").size())
    sample_table["Collection_Date"] = sample_table["Collection_Date"].dt.strftime("%Y-%m-%d").fillna("")
    dates = data["Collection_Date"].dropna()
    depths = data["Read_Depth"].dropna()
    return {
        "samples": int(data["Sample_ID"].nunique()),
        "records": len(data),
        "unique_snvs": data[variant].drop_duplicates().shape[0],
        "date_start": dates.min().date().isoformat() if not dates.empty else None,
        "date_end": dates.max().date().isoformat() if not dates.empty else None,
        "locations": sorted(name for name in data["Location_Name"].unique() if name),
        "frequency_min": float(data["Alt_Frequency"].min()),
        "frequency_median": float(data["Alt_Frequency"].median()),
        "frequency_max": float(data["Alt_Frequency"].max()),
        "depth_count": len(depths),
        "depth_min": int(depths.min()) if not depths.empty else None,
        "depth_median": float(depths.median()) if not depths.empty else None,
        "depth_max": int(depths.max()) if not depths.empty else None,
        "top_snvs": top.to_dict("records"),
        "sample_table": sample_table.to_dict("records"),
    }
