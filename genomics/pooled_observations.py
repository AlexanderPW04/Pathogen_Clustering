"""Put processed pooled-sample results into the shared reference-based tables.

Three pooled inputs are supported:

1. A processed SNV table (read by `wastewater.processed_snv`). It lists only
   the variants that were called, so it cannot show coverage elsewhere.
2. Optional pooled-sample metadata with genome coverage and mean depth.
3. An optional per-position table with read depth, and optionally A/C/G/T
   counts, for each pooled sample. This is the only input that can show
   whether an unreported variant position was actually covered, and it is
   required for the paper's Yue-Clayton PCoA.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from genomics.unified_schema import (
    COVERAGE_UNKNOWN,
    FREQUENCY_CALL,
    NUCLEOTIDES,
    POOLED,
    depth_status,
    make_sample_table,
    make_variant_table,
    variant_key,
)


COUNT_COLUMNS = ["A_Count", "C_Count", "G_Count", "T_Count"]

POSITION_ALIASES = {
    "sample": "Sample_ID",
    "sample_id": "Sample_ID",
    "reference_id": "Reference_ID",
    "reference_genome": "Reference_ID",
    "reference": "Reference_ID",
    "chrom": "Reference_ID",
    "region": "Reference_ID",
    "position": "Position",
    "pos": "Position",
    "depth": "Depth",
    "read_depth": "Depth",
    "total_dp": "Depth",
    "a": "A_Count",
    "a_count": "A_Count",
    "c": "C_Count",
    "c_count": "C_Count",
    "g": "G_Count",
    "g_count": "G_Count",
    "t": "T_Count",
    "t_count": "T_Count",
}

METADATA_ALIASES = {
    "sample": "Sample_ID",
    "sample_id": "Sample_ID",
    "collection_date": "Collection_Date",
    "date": "Collection_Date",
    "location_name": "Location_Name",
    "location": "Location_Name",
    "site": "Location_Name",
    "latitude": "Latitude",
    "longitude": "Longitude",
    "genome_coverage_percent": "Genome_Coverage_Percent",
    "coverage_percent": "Genome_Coverage_Percent",
    "percent_genome_coverage": "Genome_Coverage_Percent",
    "mean_depth": "Mean_Depth",
    "average_depth": "Mean_Depth",
}


def standard_name(name, aliases):
    """Match a header without regard to case, spaces, hyphens, or underscores."""

    spaced = str(name).replace("_", " ").replace("-", " ")
    key = "_".join(spaced.split()).casefold()
    return aliases.get(key, str(name).strip())


def read_any_table(path):
    """Read a CSV, TSV, or XLSX file as text so nothing is guessed."""

    path = Path(path)
    extension = path.suffix.lower()
    if extension == ".csv":
        return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    if extension in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if extension == ".xlsx":
        return pd.read_excel(path, dtype=str, keep_default_na=False)
    raise ValueError("Use a .csv, .tsv, or .xlsx table.")


def load_position_table(path):
    """Load per-position depth (and optional base counts) for pooled samples."""

    data = read_any_table(path)
    data.columns = [standard_name(name, POSITION_ALIASES) for name in data.columns]
    if data.columns.duplicated().any():
        raise ValueError("Two position-table columns describe the same field.")

    required = ["Sample_ID", "Reference_ID", "Position"]
    missing = [name for name in required if name not in data.columns]
    if missing:
        raise ValueError("Position table is missing column(s): " + ", ".join(missing))

    has_counts = all(name in data.columns for name in COUNT_COLUMNS)
    if not has_counts and "Depth" not in data.columns:
        raise ValueError(
            "Position table needs a Depth column or all four A/C/G/T count columns."
        )
    if data.empty:
        raise ValueError("The position table has no rows.")

    errors = []
    for column in ["Sample_ID", "Reference_ID"]:
        data[column] = data[column].astype(str).str.strip()
        if data[column].eq("").any():
            errors.append(f"Some position-table rows have no {column}.")

    numeric_columns = ["Position"]
    if "Depth" in data.columns:
        numeric_columns.append("Depth")
    if has_counts:
        numeric_columns += COUNT_COLUMNS

    for column in numeric_columns:
        values = pd.to_numeric(data[column].replace("", pd.NA), errors="coerce")
        whole = values.notna() & values.eq(values.round())
        if not whole.all():
            errors.append(f"{column} must contain whole numbers in every row.")
        elif column == "Position" and values.lt(1).any():
            errors.append("Positions must be 1-based (1 or larger).")
        elif values.lt(0).any():
            errors.append(f"{column} cannot be negative.")
        data[column] = values

    if errors:
        raise ValueError("Position table validation failed: " + " ".join(errors))

    for column in numeric_columns:
        data[column] = data[column].astype(int)

    if has_counts:
        count_total = data[COUNT_COLUMNS].sum(axis=1)
        if "Depth" in data.columns:
            if (count_total > data["Depth"]).any():
                raise ValueError("A/C/G/T counts add up to more than Depth in some rows.")
        else:
            data["Depth"] = count_total

    keys = ["Sample_ID", "Reference_ID", "Position"]
    if data.duplicated(keys).any():
        raise ValueError("Each sample/reference/position may appear only once.")

    columns = keys + ["Depth"] + (COUNT_COLUMNS if has_counts else [])
    table = data[columns].sort_values(keys).reset_index(drop=True)
    table.attrs["has_base_counts"] = has_counts
    return table


def build_depth_lookup(position_table):
    """Store each sample's depths in an array indexed by reference position.

    Unlisted positions hold -1, which means "unknown", not zero reads.
    """

    lookup = {}
    if position_table is None:
        return lookup
    for (sample_id, reference_id), group in position_table.groupby(
        ["Sample_ID", "Reference_ID"], sort=True
    ):
        depths = np.full(int(group["Position"].max()) + 1, -1, dtype=np.int64)
        depths[group["Position"].to_numpy()] = group["Depth"].to_numpy()
        lookup[(sample_id, reference_id)] = depths
    return lookup


def depth_at(lookup, sample_id, reference_id, position):
    """Return the listed depth at one position, or None when it is unknown."""

    depths = lookup.get((sample_id, reference_id))
    if depths is None or position >= len(depths) or depths[position] < 0:
        return None
    return int(depths[position])


def pooled_position_status(lookup, sample_id, reference_id, position, minimum_depth):
    """Describe an unreported pooled position using only supplied depth."""

    depth = depth_at(lookup, sample_id, reference_id, position)
    if depth is None:
        return COVERAGE_UNKNOWN, None
    return depth_status(depth, minimum_depth), depth


def pooled_variant_records(snv_data, data_source):
    """Convert standardized processed SNV rows to the shared variant table."""

    passthrough = [
        name for name in snv_data.columns
        if name.startswith("Published_") or name.startswith("Source_")
    ]
    rows = []
    for record in snv_data.to_dict("records"):
        row = {
            "Variant_Key": variant_key(
                record["Reference_Genome"], record["Position"], record["Ref"], record["Alt"]
            ),
            "Reference_ID": record["Reference_Genome"],
            "Position": record["Position"],
            "Ref": record["Ref"],
            "Alt": record["Alt"],
            "Sample_ID": record["Sample_ID"],
            "Sample_Type": POOLED,
            "Observation_Type": FREQUENCY_CALL,
            "Alt_Frequency": record["Alt_Frequency"],
            "Read_Depth": record["Read_Depth"],
            "Alt_Depth": record.get("Alt_Depth", pd.NA),
            "Callable_Status": "Reported call",
            "Quality_Status": "Not evaluated",
            "Collection_Date": record["Collection_Date"],
            "Location_Name": record["Location_Name"],
            "Data_Source": data_source,
        }
        for name in passthrough:
            row[name] = record[name]
        rows.append(row)
    return make_variant_table(rows)


def load_pooled_metadata(path):
    """Read optional pooled-sample metadata (dates, sites, coverage, depth)."""

    data = read_any_table(path)
    data.columns = [standard_name(name, METADATA_ALIASES) for name in data.columns]
    if "Sample_ID" not in data.columns:
        raise ValueError("Pooled-sample metadata needs a Sample_ID column.")
    data["Sample_ID"] = data["Sample_ID"].astype(str).str.strip()
    if data["Sample_ID"].eq("").any() or data["Sample_ID"].duplicated().any():
        raise ValueError("Pooled-sample metadata Sample_ID values must be present and unique.")

    errors = []
    limits = {
        "Latitude": (-90, 90),
        "Longitude": (-180, 180),
        "Genome_Coverage_Percent": (0, 100),
        "Mean_Depth": (0, float("inf")),
    }
    for column, (low, high) in limits.items():
        if column not in data.columns:
            data[column] = pd.NA
            continue
        raw = data[column].astype(str).str.strip()
        values = pd.to_numeric(raw.replace("", pd.NA), errors="coerce")
        if (raw.ne("") & values.isna()).any():
            errors.append(f"{column} contains non-numeric values.")
        elif (values.notna() & ~values.between(low, high)).any():
            errors.append(f"{column} is outside {low} to {high}.")
        data[column] = values
    if errors:
        raise ValueError("Pooled-sample metadata validation failed: " + " ".join(errors))
    return data


def coverage_from_positions(position_table, sample_id, reference_id, reference_length,
                            coverage_depth):
    """Percent of reference positions with at least `coverage_depth` reads.

    Returns (percent, mean_depth, note). Unlisted positions are not assumed
    to be zero; if the table is incomplete the result is None.
    """

    rows = position_table[
        position_table["Sample_ID"].eq(sample_id)
        & position_table["Reference_ID"].eq(reference_id)
    ]
    if rows.empty or reference_length is None:
        return None, None, "No per-position depth for this reference."
    listed = set(rows["Position"])
    if len(listed) < reference_length or max(listed) > reference_length:
        return None, None, (
            f"Depth listed for {len(listed)} of {reference_length} reference positions; "
            "genome coverage was not calculated."
        )
    covered = int((rows["Depth"] >= coverage_depth).sum())
    percent = covered / reference_length * 100
    return percent, float(rows["Depth"].mean()), "Calculated from per-position depth."


def pooled_sample_records(snv_data, data_source, pooled_metadata=None,
                          position_table=None, reference_lengths=None,
                          coverage_depth=10):
    """Build one shared sample row per pooled sample."""

    reference_lengths = reference_lengths or {}
    metadata = None
    if pooled_metadata is not None:
        metadata = pooled_metadata.set_index("Sample_ID")

    rows = []
    for sample_id, group in snv_data.groupby("Sample_ID", sort=True):
        references = sorted(group["Reference_Genome"].unique())
        reference_id = references[0] if len(references) == 1 else "; ".join(references)
        row = {
            "Sample_ID": sample_id,
            "Sample_Type": POOLED,
            "Collection_Date": group["Collection_Date"].iloc[0],
            "Location_Name": group["Location_Name"].iloc[0],
            "Latitude": pd.NA,
            "Longitude": pd.NA,
            "Coordinate_Source": pd.NA,
            "Geographic_Precision": "Named sampling site (no coordinates)",
            "Reference_ID": reference_id,
            "Genomic_Region": "Reference genome positions",
            "Data_Source": data_source,
            "Coverage_Information": "Variant records only (coverage elsewhere unknown)",
            "Genome_Coverage_Percent": pd.NA,
            "Mean_Depth": pd.NA,
            "Quality_Status": "Not evaluated",
            "Variant_Profile_Status": "Reference positions reported by source",
        }
        if metadata is not None and sample_id in metadata.index:
            details = metadata.loc[sample_id]
            if pd.notna(details["Genome_Coverage_Percent"]):
                row["Genome_Coverage_Percent"] = details["Genome_Coverage_Percent"]
                row["Coverage_Information"] = "Genome coverage supplied in metadata"
            if pd.notna(details["Mean_Depth"]):
                row["Mean_Depth"] = details["Mean_Depth"]
            if pd.notna(details["Latitude"]) and pd.notna(details["Longitude"]):
                row["Latitude"] = details["Latitude"]
                row["Longitude"] = details["Longitude"]
                row["Coordinate_Source"] = "Provided"
                row["Geographic_Precision"] = "Provided coordinates (precision unverified)"
        if position_table is not None and len(references) == 1:
            percent, mean_depth, note = coverage_from_positions(
                position_table, sample_id, reference_id,
                reference_lengths.get(reference_id), coverage_depth,
            )
            if percent is not None:
                row["Genome_Coverage_Percent"] = percent
                row["Mean_Depth"] = mean_depth
                row["Coverage_Information"] = note
            elif (position_table["Sample_ID"] == sample_id).any():
                row["Coverage_Information"] = note
        rows.append(row)
    return make_sample_table(rows)


def base_frequencies(position_table):
    """Add A/C/G/T frequencies (count / sum of the four counts) to each row.

    Rows whose four counts sum to zero get blank frequencies: a frequency
    cannot be calculated without reads.
    """

    if not position_table.attrs.get("has_base_counts", False):
        raise ValueError("Base frequencies need A/C/G/T count columns.")
    table = position_table.copy()
    totals = table[COUNT_COLUMNS].sum(axis=1)
    for base, column in zip(NUCLEOTIDES, COUNT_COLUMNS):
        table[f"{base}_Frequency"] = (table[column] / totals).where(totals > 0)
    return table
