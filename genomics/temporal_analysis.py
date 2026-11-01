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
    """Describe every included pooled sample at every included variant."""

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
    )

    rows = []
    for variant in unique_variants.to_dict("records"):
        for sample in sample_rows.to_dict("records"):
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
            rows.append({
                **variant,
                "Sample_ID": sample["Sample_ID"],
                "Collection_Date": sample["Collection_Date"],
                "Location_Name": sample["Location_Name"],
                "Status": status,
                "Alt_Frequency": frequency,
                "Depth": depth,
            })
    matrix = pd.DataFrame(rows, columns=MATRIX_COLUMNS)
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
    rows = []
    if matrix.empty:
        return pd.DataFrame(columns=columns)

    for (key, location), group in matrix.groupby(["Variant_Key", "Location_Name"], sort=True):
        group = group.sort_values("Collection_Date")
        reported = group[group["Status"].eq(VARIANT_OBSERVED) & group["Collection_Date"].notna()]
        status_text = group["Status"].astype(str)
        row = {
            "Variant_Key": key,
            "Location_Name": location,
            "Samples": len(group),
            "Reported": int(group["Status"].eq(VARIANT_OBSERVED).sum()),
            "Adequate_Depth_Not_Reported": int(status_text.str.startswith("Adequate").sum()),
            "Low_Or_No_Coverage": int(
                (status_text.str.startswith("Low") | status_text.str.startswith("Not covered")).sum()
            ),
            "Coverage_Unknown": int(status_text.eq("Coverage unknown").sum()),
            "Excluded_By_QC": int(status_text.eq(EXCLUDED_BY_QC).sum()),
            "First_Reported_Date": pd.NaT,
            "Last_Reported_Date": pd.NaT,
            "First_Reported_Frequency": np.nan,
            "Last_Reported_Frequency": np.nan,
            "Change_In_Reported_Frequency": np.nan,
            "Slope_Per_Week": np.nan,
        }
        if not reported.empty:
            dates = list(reported["Collection_Date"])
            values = list(reported["Alt_Frequency"])
            row["First_Reported_Date"] = dates[0]
            row["Last_Reported_Date"] = dates[-1]
            row["First_Reported_Frequency"] = values[0]
            row["Last_Reported_Frequency"] = values[-1]
            if len(values) >= 2:
                row["Change_In_Reported_Frequency"] = values[-1] - values[0]
            row["Slope_Per_Week"] = weekly_slope(dates, values)
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


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
        if abs(shuffled_slope) >= abs(slope):
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
