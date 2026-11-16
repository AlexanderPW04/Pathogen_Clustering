"""Follow pooled-sample variant frequencies across collection dates and sites.

Every pooled sample and every analyzed variant gets one row in a long
"observation matrix". A reported frequency is used as given. An unreported
variant is described by the per-position depth table when one was supplied
(adequate depth, low coverage, or not covered); otherwise it stays "Coverage
unknown".

Callable observations (audit S3-02). Trends use every observation whose
frequency is actually known, not only detections:

- a reported frequency (a detection);
- an unreported variant at adequate depth whose A/C/G/T counts show the
  alternative base: 0 alternative reads gives a frequency of 0 (an observed
  non-detection), and some reads give alternative reads / depth.

An unreported variant at adequate depth *without* base counts is not set to
0: the caller may have a reporting threshold, so the frequency is unknown
(below the caller's report). Low, absent, or unknown coverage and QC-excluded
calls are missing. Trends are descriptive, not estimates of community
prevalence and not tests of change.
"""

import numpy as np
import pandas as pd

from genomics.pooled_observations import COUNT_COLUMNS, pooled_position_status
from genomics.unified_schema import COVERED_NOT_DETECTED, EXCLUDED_BY_QC, NUCLEOTIDES, VARIANT_OBSERVED


MATRIX_COLUMNS = [
    "Variant_Key", "Reference_ID", "Position", "Ref", "Alt", "Sample_ID",
    "Collection_Date", "Location_Name", "Status", "Alt_Frequency", "Depth",
]
EVIDENCE_COLUMNS = ["Alt_Read_Count", "Callable_Frequency", "Callable_Evidence"]

REPORTED_EVIDENCE = "Reported frequency"
ZERO_ALT_READS = "Observed non-detection: 0 alternative reads at adequate depth (base counts)"
ALT_READS_NOT_CALLED = "Not reported; alternative reads in base counts (count-derived frequency)"
VALUE_UNKNOWN = "Not reported at adequate depth; alternative reads not supplied (value unknown, not 0)"


def add_callable_evidence(matrix, position_table=None):
    """Add Alt_Read_Count, Callable_Frequency, and Callable_Evidence to a matrix.

    Only explicit evidence gives a value: a reported frequency, or A/C/G/T
    counts at an adequately covered unreported position. Everything else is
    missing (NaN) with the reason in Callable_Evidence.
    """

    status = matrix["Status"].astype(str)
    observed = status.eq(VARIANT_OBSERVED)
    adequate = status.eq(COVERED_NOT_DETECTED)
    frequency = matrix["Alt_Frequency"].where(observed).astype(float)
    alt_count = pd.Series(np.nan, index=matrix.index)
    evidence = status.where(~observed, REPORTED_EVIDENCE).where(~adequate, VALUE_UNKNOWN)
    has_counts = (position_table is not None and position_table.attrs.get("has_base_counts", False)
                  and adequate.any())
    if has_counts:
        cells = matrix.loc[adequate, ["Sample_ID", "Reference_ID", "Position", "Alt"]]
        counts = cells.reset_index().merge(
            position_table[["Sample_ID", "Reference_ID", "Position", "Depth"] + COUNT_COLUMNS],
            on=["Sample_ID", "Reference_ID", "Position"], how="left").set_index("index")
        column_for = dict(zip(NUCLEOTIDES, COUNT_COLUMNS))
        reads = pd.Series([row[column_for[row["Alt"]]] for row in counts.to_dict("records")],
                          index=counts.index, dtype=float)
        depth = counts["Depth"].astype(float)
        usable = reads.notna() & depth.gt(0)
        alt_count.loc[usable.index[usable]] = reads[usable]
        frequency.loc[usable.index[usable]] = (reads / depth)[usable]
        zero = usable & reads.eq(0)
        evidence.loc[zero.index[zero]] = ZERO_ALT_READS
        evidence.loc[usable.index[usable & ~zero]] = ALT_READS_NOT_CALLED
    matrix["Alt_Read_Count"] = pd.array(alt_count.round(), dtype="Int64")
    matrix["Callable_Frequency"] = frequency
    matrix["Callable_Evidence"] = pd.Categorical(evidence)
    return matrix


def observation_matrix(variants, samples, depth_lookup, minimum_depth, position_table=None):
    """Describe every included pooled sample at every included variant.

    The table has one row per variant and sample (295,000 rows for the
    published 5,677-variant, 52-sample workbook), so it is built column by
    column: one Python dictionary per row used most of a 512 MB server.
    """

    included = variants[variants["Analysis_Included"]]
    if included.empty:
        return pd.DataFrame(columns=MATRIX_COLUMNS + EVIDENCE_COLUMNS)
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
    return add_callable_evidence(matrix, position_table)


def weekly_slope(dates, frequencies):
    """Least-squares change in reported frequency per 7 days (>= 3 points)."""

    if len(dates) < 3:
        return np.nan
    days = np.array([(date - dates[0]).days for date in dates], dtype=float)
    if np.ptp(days) == 0:
        return np.nan
    slope_per_day = np.polyfit(days, np.array(frequencies, dtype=float), 1)[0]
    return slope_per_day * 7


TREND_COLUMNS = [
    "Variant_Key", "Location_Name", "Samples", "Reported", "Adequate_Depth_Not_Reported",
    "Zero_Alt_Reads", "Alt_Reads_Not_Called", "Not_Reported_Value_Unknown", "Low_Or_No_Coverage",
    "Coverage_Unknown", "Excluded_By_QC", "Callable_Observations", "Callable_Dated",
    "First_Reported_Date", "Last_Reported_Date", "First_Callable_Date", "Last_Callable_Date",
    "First_Callable_Frequency", "Last_Callable_Frequency", "Change_In_Callable_Frequency",
    "Slope_Per_Week", "Complete_Series", "Trend_Basis",
]


def callable_frequencies(matrix):
    """The known frequency of each matrix row (NaN when it is missing)."""

    if "Callable_Frequency" in matrix:
        return matrix["Callable_Frequency"].astype(float)
    return matrix["Alt_Frequency"].where(matrix["Status"].astype(str).eq(VARIANT_OBSERVED)).astype(float)


def trend_basis(row):
    """Explain which observations a series' change and slope use."""

    if row["Callable_Dated"] < 2:
        return "Fewer than two dated callable observations: no change or trend is calculated."
    missing = int(row["Samples"] - row["Callable_Observations"])
    text = f"Uses {int(row['Callable_Dated'])} dated callable observation(s)"
    zeros = int(row["Zero_Alt_Reads"])
    if zeros:
        text += f", including {zeros} observed non-detection(s) (0 alternative reads)"
    text += "."
    if missing:
        text += (f" {missing} sample(s) without a known frequency (unknown, low, or no coverage, a "
                 "QC-excluded call, or no alternative-read counts) are left out, so the series may "
                 "over-represent samples where the variant was detected.")
    if row["Callable_Observations"] > row["Callable_Dated"]:
        text += f" {int(row['Callable_Observations'] - row['Callable_Dated'])} callable observation(s) have no usable date."
    if row["Callable_Dated"] == 2:
        text += " Two points give a change, not a trend."
    return text


def frequency_trends(matrix, keep_unknown_series=False):
    """Summarize each variant within each sampling location over time.

    Change and slope use callable observations (see the module docstring):
    reported frequencies plus base-count values at adequately covered
    unreported positions. Missing observations are counted, never filled in.

    A variant-site series without any known frequency (the variant was never
    reported there and no base counts exist) carries no frequency information;
    it is left out unless `keep_unknown_series` is True, and the number left
    out is in `table.attrs["series_without_known_frequency"]`. Its statuses
    stay in the observation matrix.
    """

    if matrix.empty:
        table = pd.DataFrame(columns=TREND_COLUMNS)
        table.attrs["series_without_known_frequency"] = 0
        return table

    # Count each status per variant and site with one grouped sum (a Python
    # loop over every variant/site pair is far too slow for the paper's table).
    keys = ["Variant_Key", "Location_Name"]
    status = matrix["Status"].astype(str)
    value = callable_frequencies(matrix)
    # Compared as a categorical column: converting 300,000 rows to text costs ~100 MB.
    evidence = (matrix["Callable_Evidence"] if "Callable_Evidence" in matrix
                else pd.Series("", index=matrix.index))
    # One-byte indicators: the matrix can have hundreds of thousands of rows.
    starts_adequate = status.str.startswith("Adequate").to_numpy(dtype=bool)
    known = value.notna().to_numpy()
    indicators = {
        "Samples": np.ones(len(matrix), dtype=bool),
        "Reported": status.eq(VARIANT_OBSERVED).to_numpy(dtype=bool),
        "Adequate_Depth_Not_Reported": starts_adequate,
        "Zero_Alt_Reads": evidence.eq(ZERO_ALT_READS).to_numpy(dtype=bool),
        "Alt_Reads_Not_Called": evidence.eq(ALT_READS_NOT_CALLED).to_numpy(dtype=bool),
        "Not_Reported_Value_Unknown": starts_adequate & ~known,
        "Low_Or_No_Coverage": (status.str.startswith("Low") | status.str.startswith("Not covered")).to_numpy(dtype=bool),
        "Coverage_Unknown": status.eq("Coverage unknown").to_numpy(dtype=bool),
        "Excluded_By_QC": status.eq(EXCLUDED_BY_QC).to_numpy(dtype=bool),
        "Callable_Observations": known,
        "Callable_Dated": known & matrix["Collection_Date"].notna().to_numpy(),
    }
    flags = pd.DataFrame({key: matrix[key].to_numpy() for key in keys})
    for name, values in indicators.items():
        flags[name] = values.astype(np.uint8)
    del indicators
    table = flags.groupby(keys, sort=True).sum()

    reported = matrix[status.eq(VARIANT_OBSERVED) & matrix["Collection_Date"].notna()]
    reported_groups = reported.sort_values(keys + ["Collection_Date"], kind="mergesort").groupby(keys)
    table["First_Reported_Date"] = reported_groups["Collection_Date"].first()
    table["Last_Reported_Date"] = reported_groups["Collection_Date"].last()

    # First/last callable values in date order (ties broken by sample ID).
    known = matrix.loc[value.notna() & matrix["Collection_Date"].notna(),
                       keys + ["Collection_Date"] + (["Sample_ID"] if "Sample_ID" in matrix else [])]
    known = known.assign(Value=value[known.index])
    known = known.sort_values(list(known.columns[:-1]), kind="mergesort")
    groups = known.groupby(keys, sort=True)
    counts = groups.size()
    table["First_Callable_Date"] = groups["Collection_Date"].first()
    table["Last_Callable_Date"] = groups["Collection_Date"].last()
    table["First_Callable_Frequency"] = groups["Value"].first()
    table["Last_Callable_Frequency"] = groups["Value"].last()
    change = table["Last_Callable_Frequency"] - table["First_Callable_Frequency"]
    table["Change_In_Callable_Frequency"] = change.where(counts.reindex(table.index).ge(2))
    table["Slope_Per_Week"] = np.nan
    for group_key in counts[counts >= 3].index:
        group = groups.get_group(group_key)
        table.loc[group_key, "Slope_Per_Week"] = weekly_slope(
            list(group["Collection_Date"]), list(group["Value"]))
    table["Complete_Series"] = table["Callable_Observations"].eq(table["Samples"])
    table = table.reset_index()
    without = table["Callable_Observations"].eq(0)
    if not keep_unknown_series:
        table = table[~without].reset_index(drop=True)
    table["Trend_Basis"] = [trend_basis(row) for row in table.to_dict("records")]
    table = table[TREND_COLUMNS]
    table.attrs["series_without_known_frequency"] = int(without.sum())
    return table


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
        # One value per sample and variant in the median (repeated rows are one observation).
        once = records.drop_duplicates(["Sample_ID", "Variant_Key"])
        rows.append({
            "Location_Name": location,
            "Collection_Date": date,
            "Pooled_Samples": len(sample_ids),
            "Variant_Records": len(records),
            "Distinct_Variants": records["Variant_Key"].nunique(),
            "Median_Reported_Frequency": (
                float(once["Alt_Frequency"].median()) if not once.empty else np.nan
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
