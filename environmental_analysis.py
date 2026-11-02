"""Explore already processed environmental measurements with simple tables."""

from io import BytesIO
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = ["Sample_ID", "Feature_ID"]
ALIASES = {
    "sample": "Sample_ID", "sample_id": "Sample_ID", "sample_identifier": "Sample_ID",
    "date": "Collection_Date", "collection_date": "Collection_Date",
    "collectiondate": "Collection_Date",
    "location": "Location_Name", "site": "Location_Name", "location_name": "Location_Name",
    "sampling_location": "Location_Name", "site_name": "Location_Name",
    "feature": "Feature_ID", "feature_identifier": "Feature_ID", "feature_id": "Feature_ID",
    "value": "Measurement_Value", "measurement": "Measurement_Value", "measurement_value": "Measurement_Value",
    "frequency": "Frequency", "feature_frequency": "Frequency",
    "alt_frequency": "Frequency", "alt_freq": "Frequency",
    "coverage": "Coverage", "coverage_depth": "Coverage",
    "read_depth": "Read_Depth", "depth": "Read_Depth",
    "quality": "Quality", "quality_flag": "Quality", "quality_score": "Quality",
    "quality_information": "Quality",
    "unit": "Unit", "units": "Unit",
}


def load_measurements(path, metric="auto"):
    """Read, standardize, and validate one row per sample and feature."""

    path = Path(path)
    if path.suffix.lower() == ".csv":
        data = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    elif path.suffix.lower() == ".xlsx":
        data = pd.read_excel(path, dtype=str, keep_default_na=False)
    else:
        raise ValueError("Choose a CSV or XLSX measurement file.")

    names = []
    for name in data.columns:
        key = "_".join(str(name).strip().casefold().replace("-", " ").replace("_", " ").split())
        names.append(ALIASES.get(key, str(name).strip()))
    data.columns = names
    if data.columns.duplicated().any():
        raise ValueError("Two columns map to the same field. Keep one column per field.")
    missing = [name for name in REQUIRED_COLUMNS if name not in data]
    if missing:
        raise ValueError("Missing column(s): " + ", ".join(missing))
    available_metrics = [name for name in ["Frequency", "Measurement_Value"] if name in data]
    if data.empty or len(data) > 10000:
        raise ValueError("Use a table containing 1 to 10,000 measurement rows.")

    for column in data:
        data[column] = data[column].fillna("").astype(str).str.strip()
    for column in ["Collection_Date", "Location_Name", "Quality", "Unit"]:
        if column not in data:
            data[column] = ""

    errors = []
    for column in ["Sample_ID", "Feature_ID"]:
        count = int(data[column].eq("").sum())
        if count:
            errors.append(f"{count} row(s) have no {column}.")

    # Excel date cells may be read as midnight timestamps; accept their date.
    raw_dates = data["Collection_Date"].str.removesuffix(" 00:00:00")
    dates = pd.to_datetime(raw_dates.replace("", pd.NA), format="%Y-%m-%d", errors="coerce")
    if (raw_dates.ne("") & dates.isna()).any():
        errors.append("Collection_Date contains invalid dates; use YYYY-MM-DD or Excel date cells.")
    data["Collection_Date"] = dates

    # Merge capitalization/spacing variants of the same location label.
    locations = {}
    for name in data["Location_Name"]:
        label = " ".join(name.split())
        locations.setdefault(label.casefold(), label)
    data["Location_Name"] = data["Location_Name"].map(
        lambda name: locations[" ".join(name.split()).casefold()]
    )

    numeric_columns = available_metrics + [
        name for name in ["Coverage", "Read_Depth"] if name in data
    ]
    for column in numeric_columns:
        raw = data[column]
        values = pd.to_numeric(raw.replace("", pd.NA), errors="coerce")
        invalid = raw.ne("") & (values.isna() | values.isin([float("inf"), float("-inf")]))
        if invalid.any():
            errors.append(f"{column} contains {int(invalid.sum())} non-numeric or non-finite value(s).")
        if column == "Frequency" and (values.notna() & ~values.between(0, 1)).any():
            errors.append("Frequency must be a fraction from 0 to 1; percentages are not converted automatically.")
        if column in {"Coverage", "Read_Depth"} and values.lt(0).any():
            errors.append(f"{column} cannot be negative.")
        data[column] = values.astype(float)

    for column in ["Collection_Date", "Location_Name"]:
        present = data[column].replace("", pd.NA)
        if present.groupby(data["Sample_ID"]).nunique().gt(1).any():
            errors.append(f"A Sample_ID has conflicting {column} values. Use a distinct ID for each collected sample.")
    if data.duplicated(["Sample_ID", "Feature_ID"]).any():
        errors.append("Repeated Sample_ID/Feature_ID rows found. Resolve replicates before uploading.")
    units = data["Unit"].replace("", pd.NA).groupby(data["Feature_ID"]).nunique()
    if units.gt(1).any():
        errors.append("A feature has different units. Convert to a common unit before comparing it.")

    if metric == "auto":
        metric = next(
            (name for name in available_metrics if data[name].notna().any()), None
        )
    if metric is not None and metric not in available_metrics:
        errors.append("The selected measurement column is not present in this file.")
    if errors:
        raise ValueError(" ".join(errors))

    return data.sort_values(["Sample_ID", "Feature_ID"]).reset_index(drop=True), metric


def summarize_measurements(data, metric):
    """Describe reported values without filling missing measurements with zero."""

    # A feature list can still be summarized without a numeric measurement.
    working = data.copy()
    value_column = metric or "_Unavailable_Value"
    if metric is None:
        working[value_column] = float("nan")
    samples = working.groupby("Sample_ID", as_index=False).agg(
        Collection_Date=("Collection_Date", "first"),
        Location_Name=("Location_Name", "first"),
        Features_Listed=("Feature_ID", "nunique"),
        Values_Available=(value_column, "count"),
    )
    samples["Missing_Values"] = samples["Features_Listed"] - samples["Values_Available"]
    observed = working.loc[working[value_column].notna()]
    dated = working.loc[working["Collection_Date"].notna()]
    temporal = dated.groupby(
        ["Feature_ID", "Location_Name", "Collection_Date"], as_index=False
    ).agg(
        Mean_Value=(value_column, "mean"),
        Samples_Listed=("Sample_ID", "nunique"),
        Values_Available=(value_column, "count"),
    )
    temporal = temporal.sort_values(["Feature_ID", "Location_Name", "Collection_Date"])
    groups = temporal.groupby(["Feature_ID", "Location_Name"])
    temporal["Previous_Date"] = groups["Collection_Date"].shift()
    temporal["Change_From_Previous"] = groups["Mean_Value"].diff()

    locations = working.loc[working["Location_Name"].ne("")].groupby(
        ["Feature_ID", "Location_Name"], as_index=False
    ).agg(
        Mean_Value=(value_column, "mean"), Minimum=(value_column, "min"),
        Maximum=(value_column, "max"), Samples_Listed=("Sample_ID", "nunique"),
        Values_Available=(value_column, "count"),
    )
    location_counts = working.loc[working["Location_Name"].ne("")].groupby(
        "Location_Name", as_index=False
    ).agg(Samples=("Sample_ID", "nunique"), Features_Listed=("Feature_ID", "nunique"))
    date_counts = dated.groupby("Collection_Date", as_index=False).agg(
        Samples=("Sample_ID", "nunique"), Features_Listed=("Feature_ID", "nunique")
    )
    missing = []
    for column in REQUIRED_COLUMNS + ["Collection_Date", "Location_Name"] + available_value_columns(data) + ["Quality", "Unit"]:
        count = int((data[column].isna() | data[column].eq("")).sum())
        missing.append({"Field": column, "Missing_Rows": count, "Total_Rows": len(data)})
    quality = data["Quality"].replace("", "Not provided").value_counts().rename_axis("Quality").reset_index(name="Rows")
    return {
        "counts": {
            "samples": int(data["Sample_ID"].nunique()),
            "samples_with_values": int(observed["Sample_ID"].nunique()),
            "features": int(data["Feature_ID"].nunique()),
            "rows": len(data),
            "frequency_values": int(data["Frequency"].count()) if "Frequency" in data else 0,
            "coverage_values": int(data["Coverage"].count()) if "Coverage" in data else 0,
            "read_depth_values": int(data["Read_Depth"].count()) if "Read_Depth" in data else 0,
            "missing_dates": int(data["Collection_Date"].isna().sum()),
            "missing_locations": int(data["Location_Name"].eq("").sum()),
        },
        "tables": {
            "measurements": data,
            "samples": samples,
            "temporal_changes": temporal,
            "location_summary": locations,
            "location_counts": location_counts,
            "date_counts": date_counts,
            "missing_data": pd.DataFrame(missing),
            "quality_summary": quality,
        },
    }


def available_value_columns(data):
    """List the numeric columns retained in the organized table."""

    return [
        name for name in ["Frequency", "Measurement_Value", "Coverage", "Read_Depth"]
        if name in data
    ]


def table_records(table):
    """Convert a summary table to simple values that Jinja can display."""

    display = table.copy()
    for column in display:
        if pd.api.types.is_datetime64_any_dtype(display[column]):
            display[column] = display[column].dt.strftime("%Y-%m-%d")
    display = display.astype(object).where(pd.notna(display), "")
    return display.to_dict("records")


def make_charts(data, metric, tables):
    """Draw sample counts and time/location plots for up to three features."""

    # Import plotting only when this mode is used; Sprint 2 needs no plotting.
    from matplotlib.dates import AutoDateLocator, ConciseDateFormatter
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    features = [] if metric is None else sorted(
        data.loc[data[metric].notna(), "Feature_ID"].unique()
    )[:3]
    figure = Figure(figsize=(10, 3.6 * (len(features) + 1)), layout="constrained")
    axes = figure.subplots(len(features) + 1, 2, squeeze=False)
    sample_axis, date_axis = axes[0]
    sample_counts = tables["samples"].head(20)
    sample_axis.barh(sample_counts["Sample_ID"], sample_counts["Features_Listed"], color="#1f5f8b")
    sample_axis.set_title("Features listed per sample (first 20)")
    sample_axis.set_xlabel("Number of features listed")
    sample_axis.xaxis.set_major_locator(MaxNLocator(integer=True))
    date_counts = tables["date_counts"]
    if date_counts.empty:
        date_axis.text(0.5, 0.5, "Collection dates unavailable", ha="center", transform=date_axis.transAxes)
    else:
        date_axis.bar(date_counts["Collection_Date"], date_counts["Samples"], color="#1f5f8b")
        locator = AutoDateLocator()
        date_axis.xaxis.set_major_locator(locator)
        date_axis.xaxis.set_major_formatter(ConciseDateFormatter(locator))
    date_axis.set_title("Samples by collection date")
    date_axis.set_ylabel("Number of samples")
    date_axis.yaxis.set_major_locator(MaxNLocator(integer=True))
    for index, feature in enumerate(features):
        dates = tables["temporal_changes"].loc[tables["temporal_changes"]["Feature_ID"].eq(feature)]
        locations = tables["location_summary"].loc[tables["location_summary"]["Feature_ID"].eq(feature)].head(6)
        time_axis, location_axis = axes[index + 1]
        for location, group in list(dates.groupby("Location_Name"))[:6]:
            time_axis.plot(group["Collection_Date"], group["Mean_Value"], marker="o", label=location or "Location unavailable")
        if dates.empty:
            time_axis.text(0.5, 0.5, "No dated measurements", ha="center", transform=time_axis.transAxes)
        else:
            locator = AutoDateLocator()
            time_axis.xaxis.set_major_locator(locator)
            time_axis.xaxis.set_major_formatter(ConciseDateFormatter(locator))
            time_axis.legend(fontsize=7)
        if locations.empty:
            location_axis.text(0.5, 0.5, "No location measurements", ha="center", transform=location_axis.transAxes)
        else:
            location_axis.barh(locations["Location_Name"], locations["Mean_Value"], color="#1f5f8b")
        units = data.loc[data["Feature_ID"].eq(feature), "Unit"]
        units = units.loc[units.ne("")]
        label = "Frequency (0–1)" if metric == "Frequency" else (units.iloc[0] if not units.empty else "Measurement value")
        time_axis.set_title(f"{str(feature)[:48]}: over time")
        time_axis.set_ylabel(label)
        location_axis.set_title(f"{str(feature)[:48]}: location means")
        location_axis.set_xlabel(label)
    image = BytesIO()
    figure.savefig(image, format="png", dpi=110)
    return image.getvalue()
