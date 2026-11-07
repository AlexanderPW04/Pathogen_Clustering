"""Follow pooled-sample variant frequencies across collection dates and sites.

Every pooled sample and every analyzed variant gets one row in a long
"observation matrix". A reported frequency is used as given. An unreported
variant is described by the per-position depth table when one was supplied
(adequate depth, low coverage, or not covered); otherwise it stays "Coverage
unknown". Unreported variants never receive a 0% frequency.

Trends here summarize the sample-level observations that exist. They are
descriptive and are not estimates of community prevalence.
"""

import numpy as np
import pandas as pd

from genomics.pooled_observations import pooled_position_status
from genomics.unified_schema import EXCLUDED_BY_QC, VARIANT_OBSERVED


MATRIX_COLUMNS = [
    "Variant_Key", "Reference_ID", "Position", "Ref", "Alt", "Sample_ID",
    "Collection_Date", "Location_Name", "Status", "Alt_Frequency", "Depth",
]


def observation_matrix(variants, samples, depth_lookup, minimum_depth):
    """Describe every included pooled sample at every included variant.

    The table has one row per variant and sample (295,000 rows for the
    published 5,677-variant, 52-sample workbook), so it is built column by
    column: one Python dictionary per row used most of a 512 MB server.
    """

    included = variants[variants["Analysis_Included"]]
    if included.empty:
        return pd.DataFrame(columns=MATRIX_COLUMNS)
    # A call that exists but failed QC is neither a measurement nor a non-detection.
    excluded = set(zip(variants.loc[~variants["Analysis_Included"], "Sample_ID"],
                       variants.loc[~variants["Analysis_Included"], "Variant_Key"]))

    # Repeated sample/variant rows (for example one position listed under two
    # annotations) keep their first reported frequency.
    reported = {}
    for row in included.to_dict("records"):
        key = (row["Sample_ID"], row["Variant_Key"])
        if key not in reported:
            reported[key] = (row["Alt_Frequency"], row["Read_Depth"])

    unique_variants = included.drop_duplicates("Variant_Key")[
        ["Variant_Key", "Reference_ID", "Position", "Ref", "Alt"]
    ].sort_values(["Reference_ID", "Position", "Alt"])
    sample_rows = samples[samples["Analysis_Included"]].sort_values(
        ["Collection_Date", "Sample_ID"]
    ).to_dict("records")

    columns = {name: [] for name in MATRIX_COLUMNS}
    for variant in unique_variants.to_dict("records"):
        for sample in sample_rows:
            if variant["Reference_ID"] not in str(sample["Reference_ID"]).split("; "):
                continue  # never place a sample on another reference's positions
            key = (sample["Sample_ID"], variant["Variant_Key"])
            if key in reported:
                status = VARIANT_OBSERVED
                frequency, depth = reported[key]
            elif key in excluded:
                status, frequency, depth = EXCLUDED_BY_QC, np.nan, None
            else:
                status, depth = pooled_position_status(
                    depth_lookup, sample["Sample_ID"], variant["Reference_ID"],
                    int(variant["Position"]), minimum_depth,
                )
                frequency = np.nan
            for name in ["Variant_Key", "Reference_ID", "Position", "Ref", "Alt"]:
                columns[name].append(variant[name])
            columns["Sample_ID"].append(sample["Sample_ID"])
            columns["Collection_Date"].append(sample["Collection_Date"])
            columns["Location_Name"].append(sample["Location_Name"])
            columns["Status"].append(status)
            columns["Alt_Frequency"].append(frequency)
            columns["Depth"].append(depth)
    matrix = pd.DataFrame({name: pd.Series(columns.pop(name), dtype=object).infer_objects()
                           for name in MATRIX_COLUMNS})
    matrix["Alt_Frequency"] = pd.to_numeric(matrix["Alt_Frequency"], errors="coerce")
    return matrix


def weekly_slope(dates, frequencies):
    """Least-squares change in reported frequency per 7 days (>= 3 points)."""

    if len(dates) < 3:
        return np.nan
    days = np.array([(date - dates[0]).days for date in dates], dtype=float)
    if np.ptp(days) == 0:
        return np.nan
    slope_per_day = np.polyfit(days, np.array(frequencies, dtype=float), 1)[0]
    return slope_per_day * 7


def frequency_trends(matrix):
    """Summarize each variant within each sampling location over time."""

    columns = [
        "Variant_Key", "Location_Name", "Samples", "Reported", "Adequate_Depth_Not_Reported",
        "Low_Or_No_Coverage", "Coverage_Unknown", "Excluded_By_QC", "First_Reported_Date",
        "Last_Reported_Date",
        "First_Reported_Frequency", "Last_Reported_Frequency", "Change_In_Reported_Frequency",
        "Slope_Per_Week",
    ]
    if matrix.empty:
        return pd.DataFrame(columns=columns)

    # Count each status per variant and site with one grouped sum (a Python
    # loop over every variant/site pair is far too slow for the paper's table).
    keys = ["Variant_Key", "Location_Name"]
    status = matrix["Status"].astype(str)
    flags = matrix[keys].copy()
    flags["Samples"] = 1
    flags["Reported"] = status.eq(VARIANT_OBSERVED).astype(int)
    flags["Adequate_Depth_Not_Reported"] = status.str.startswith("Adequate").astype(int)
    flags["Low_Or_No_Coverage"] = (status.str.startswith("Low")
                                   | status.str.startswith("Not covered")).astype(int)
    flags["Coverage_Unknown"] = status.eq("Coverage unknown").astype(int)
    flags["Excluded_By_QC"] = status.eq(EXCLUDED_BY_QC).astype(int)
    table = flags.groupby(keys, sort=True).sum()

    # First/last reported values in date order (ties broken by sample ID).
    reported = matrix[status.eq(VARIANT_OBSERVED) & matrix["Collection_Date"].notna()]
    order = keys + ["Collection_Date"] + (["Sample_ID"] if "Sample_ID" in reported else [])
    reported = reported.sort_values(order, kind="mergesort")
    groups = reported.groupby(keys, sort=True)
    counts = groups.size()
    table["First_Reported_Date"] = groups["Collection_Date"].first()
    table["Last_Reported_Date"] = groups["Collection_Date"].last()
    table["First_Reported_Frequency"] = groups["Alt_Frequency"].first()
    table["Last_Reported_Frequency"] = groups["Alt_Frequency"].last()
    change = table["Last_Reported_Frequency"] - table["First_Reported_Frequency"]
    table["Change_In_Reported_Frequency"] = change.where(counts.reindex(table.index).ge(2))
    table["Slope_Per_Week"] = np.nan
    for group_key in counts[counts >= 3].index:
        group = groups.get_group(group_key)
        table.loc[group_key, "Slope_Per_Week"] = weekly_slope(
            list(group["Collection_Date"]), list(group["Alt_Frequency"]))
    return table.reset_index()[columns]


def location_date_summary(variants, samples):
    """Count included pooled samples and reported variants by site and date."""

    included_samples = samples[samples["Analysis_Included"]]
    included_variants = variants[variants["Analysis_Included"]]
    rows = []
    for (location, date), group in included_samples.groupby(
        ["Location_Name", "Collection_Date"], sort=True, dropna=False
    ):
        sample_ids = set(group["Sample_ID"])
        records = included_variants[included_variants["Sample_ID"].isin(sample_ids)]
        rows.append({
            "Location_Name": location,
            "Collection_Date": date,
            "Pooled_Samples": len(sample_ids),
            "Variant_Records": len(records),
            "Distinct_Variants": records["Variant_Key"].nunique(),
            "Median_Reported_Frequency": (
                float(records["Alt_Frequency"].median()) if not records.empty else np.nan
            ),
        })
    return pd.DataFrame(rows, columns=[
        "Location_Name", "Collection_Date", "Pooled_Samples", "Variant_Records",
        "Distinct_Variants", "Median_Reported_Frequency",
    ])


def variants_per_sample(variants, samples):
    """Count SNV records and distinct SNVs per pooled sample (paper Fig 2A)."""

    rows = []
    for sample in samples.to_dict("records"):
        records = variants[variants["Sample_ID"].eq(sample["Sample_ID"])]
        kept = records[records["Analysis_Included"]] if "Analysis_Included" in records else records
        rows.append({
            "Sample_ID": sample["Sample_ID"],
            "Collection_Date": sample["Collection_Date"],
            "Location_Name": sample["Location_Name"],
            "SNV_Records": len(records),
            "Distinct_SNVs": records["Variant_Key"].nunique(),
            "Distinct_SNVs_After_QC": kept["Variant_Key"].nunique(),
            "Genome_Coverage_Percent": sample["Genome_Coverage_Percent"],
            "Mean_Depth": sample["Mean_Depth"],
            "Quality_Status": sample.get("Quality_Status", ""),
            "Analysis_Included": sample.get("Analysis_Included", True),
        })
    return pd.DataFrame(rows)


def depth_regression(per_sample, bootstrap_rounds=1000, permutations=2000, seed=0):
    """Linear regression of distinct SNVs on mean depth (paper Fig 2B).

    The 95% band is a bootstrap percentile interval of the fitted line, and
    the p-value is a two-sided permutation test of the slope. Both use a fixed
    random seed so the result is reproducible. Returns None when fewer than
    three included samples have a mean depth.
    """

    data = per_sample[per_sample["Analysis_Included"] & per_sample["Mean_Depth"].notna()]
    if len(data) < 3 or data["Mean_Depth"].nunique() < 2:
        return None

    x = data["Mean_Depth"].to_numpy(dtype=float)
    y = data["Distinct_SNVs_After_QC"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    r = float(np.corrcoef(x, y)[0, 1]) if np.std(y) > 0 else 0.0

    generator = np.random.default_rng(seed)
    grid = np.linspace(x.min(), x.max(), 50)
    fitted_lines = []
    for _ in range(bootstrap_rounds):
        picks = generator.integers(0, len(x), len(x))
        if np.ptp(x[picks]) == 0:
            continue
        b_slope, b_intercept = np.polyfit(x[picks], y[picks], 1)
        fitted_lines.append(b_slope * grid + b_intercept)
    fitted_lines = np.array(fitted_lines)

    extreme = 0
    for _ in range(permutations):
        shuffled_slope = np.polyfit(x, generator.permutation(y), 1)[0]
        # A tolerance keeps exact ties (equal up to rounding) counted as extreme.
        if abs(shuffled_slope) >= abs(slope) - 1e-12 * max(1.0, abs(slope)):
            extreme += 1

    return {
        "n": len(x),
        "slope": float(slope),
        "intercept": float(intercept),
        "r": r,
        "r_squared": r ** 2,
        "permutation_p": (extreme + 1) / (permutations + 1),
        "grid": grid,
        "lower": np.percentile(fitted_lines, 2.5, axis=0) if len(fitted_lines) else None,
        "upper": np.percentile(fitted_lines, 97.5, axis=0) if len(fitted_lines) else None,
        "small_sample": len(x) < 10,
    }
