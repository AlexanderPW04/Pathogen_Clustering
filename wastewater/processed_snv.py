"""Load processed SNV calls and summarize what was observed in wastewater."""

from pathlib import Path

import pandas as pd


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


def read_table(path):
    """Read a CSV or Excel table without losing numeric-looking sample IDs."""

    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, dtype=str, keep_default_na=False)
    if path.suffix.lower() == ".xlsx":
        # The paper's workbook has a descriptive first row, then its headers.
        first_rows = pd.read_excel(path, header=None, nrows=2, dtype=str)
        header_row = 1 if len(first_rows) > 1 and first_rows.shape[1] > 0 and str(first_rows.iloc[1, 0]).strip() == "SAMPLE" else 0
        data = pd.read_excel(path, header=header_row, dtype=str, keep_default_na=False)
        description = str(first_rows.iloc[0, 4]) if first_rows.shape[1] > 4 else ""
        data.attrs["fontenele_2021_table2"] = (
            header_row == 1
            and "Supplementary Table 2" in description
            and "52 wastewater samples" in description
            and "SARS-CoV-2" in description
        )
        return data
    raise ValueError("Use a processed SNV .xlsx or .csv file.")


def load_processed_snvs(path, metadata_path=None):
    """Return standardized SNV rows and visible data-quality warnings."""

    data = read_table(path)
    data.columns = [str(name).strip() for name in data.columns]
    published = (
        data.attrs.get("fontenele_2021_table2", False)
        and set(PAPER_COLUMNS).issubset(data.columns)
        and {"STATE", "LOCATION"}.issubset(data.columns)
    )
    warnings = []

    if published:
        data = data.rename(columns=PAPER_COLUMNS)
        data["Reference_Genome"] = "NC_045512.2"
        data["Location_Name"] = data["STATE"].astype(str).str.strip() + ": " + data["LOCATION"].astype(str).str.strip()
        dates = pd.to_datetime(data["Collection_Date"], format="%d_%b_%Y", errors="coerce")
    else:
        missing = sorted(set(REQUIRED_COLUMNS) - set(data.columns))
        if missing:
            raise ValueError("Missing processed SNV column(s): " + ", ".join(missing))
        if "Collection_Date" not in data:
            data["Collection_Date"] = ""
        if "Location_Name" not in data:
            data["Location_Name"] = ""
        dates = pd.to_datetime(data["Collection_Date"].replace("", pd.NA), errors="coerce")
    invalid_source_dates = data["Collection_Date"].fillna("").astype(str).str.strip().ne("") & dates.isna()

    if data.empty:
        raise ValueError("The processed SNV table has no records.")

    for column in ["Sample_ID", "Reference_Genome", "Ref", "Alt", "Location_Name"]:
        data[column] = data[column].fillna("").astype(str).str.strip()

    data["Ref"] = data["Ref"].str.upper()
    data["Alt"] = data["Alt"].str.upper()
    errors = []

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
                matched_dates = pd.to_datetime(matched.replace("", pd.NA), errors="coerce")
                bad = matched.ne("") & matched_dates.isna()
                if bad.any():
                    errors.append(f"{bad.sum()} metadata collection date(s) are invalid.")
                conflict = dates.notna() & matched_dates.notna() & dates.ne(matched_dates)
                if conflict.any():
                    errors.append(f"{conflict.sum()} SNV date(s) conflict with metadata.")
                dates = dates.fillna(matched_dates)
            else:
                conflict = data[column].ne("") & matched.ne("") & data[column].str.casefold().ne(matched.str.casefold())
                if conflict.any():
                    errors.append(f"{conflict.sum()} SNV location(s) conflict with metadata.")
                data[column] = data[column].where(data[column].ne(""), matched)

    if invalid_source_dates.any():
        errors.append(f"{invalid_source_dates.sum()} collection date(s) are invalid.")

    for column in ["Sample_ID", "Reference_Genome", "Ref", "Alt"]:
        count = int(data[column].eq("").sum())
        if count:
            errors.append(f"{count} record(s) have no {column}.")

    for column in ["Ref", "Alt"]:
        invalid = ~data[column].isin(["A", "C", "G", "T"])
        if invalid.any():
            errors.append(f"{invalid.sum()} record(s) have invalid {column}; expected one A, C, G, or T.")
    same_base = data["Ref"].eq(data["Alt"])
    if same_base.any():
        errors.append(f"{same_base.sum()} record(s) have identical reference and alternative bases.")

    positions = pd.to_numeric(data["Position"], errors="coerce")
    invalid_position = positions.isna() | positions.le(0) | positions.ne(positions.round()) | positions.isin([float("inf"), float("-inf")])
    if invalid_position.any():
        errors.append(f"{invalid_position.sum()} position(s) are not positive integers.")
    frequencies = pd.to_numeric(data["Alt_Frequency"], errors="coerce")
    invalid_frequency = frequencies.isna() | ~frequencies.between(0, 1)
    if invalid_frequency.any():
        errors.append(f"{invalid_frequency.sum()} alternative allele frequency value(s) are outside 0 to 1 or missing.")
    raw_depths = data["Read_Depth"].fillna("").astype(str).str.strip()
    depths = pd.to_numeric(raw_depths.replace("", pd.NA), errors="coerce")
    invalid_depth_text = raw_depths.ne("") & depths.isna()
    if invalid_depth_text.any():
        errors.append(f"{invalid_depth_text.sum()} read depth value(s) are not numeric.")
    bad_depth = depths.notna() & (depths.lt(1) | depths.ne(depths.round()) | depths.isin([float("inf"), float("-inf")]))
    if bad_depth.any():
        errors.append(f"{bad_depth.sum()} read depth value(s) are not positive integers.")
    depth_missing = int(depths.isna().sum())
    if depth_missing:
        warnings.append(f"{depth_missing} record(s) have no usable read depth.")

    if errors:
        raise ValueError("Processed SNV validation failed: " + " ".join(errors))

    data["Position"] = positions.astype(int)
    data["Alt_Frequency"] = frequencies.astype(float)
    data["Read_Depth"] = depths.astype("Int64")
    data["Collection_Date"] = dates
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
        warnings.append(f"{duplicate_count} repeated sample/reference/position/allele record(s) found; rows were retained.")
    missing_dates = int(data["Collection_Date"].isna().sum())
    missing_locations = int(data["Location_Name"].eq("").sum())
    if missing_dates:
        warnings.append(f"{missing_dates} record(s) lack a collection date.")
    if missing_locations:
        warnings.append(f"{missing_locations} record(s) lack a location.")

    return data[OUTPUT_COLUMNS].copy(), warnings


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
