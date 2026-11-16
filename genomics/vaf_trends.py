"""Sprint 5 mutation-frequency time series with depth-aware uncertainty.

A mutation is identified by Reference_ID (the reference sequence, one per
genome segment) + position + Ref + Alt, the Variant_Key used everywhere in
the application, so variants from different references are never merged.
Each series is one mutation at one sampling site.

Every sample in a series is in one of three states:

- detected: the variant was reported (VAF from counts when they are valid,
  otherwise the reported frequency);
- supported non-detection: A/C/G/T counts show adequate depth with 0
  alternative reads (VAF 0, with an interval), or some alternative reads
  that the caller did not report (count-derived VAF);
- unknown: coverage unknown, low, or absent, a QC-excluded call, or an
  adequately covered unreported position without base counts. Unknown
  samples are counted and drawn apart, never given a value.

Summaries are descriptive. The earliest-to-latest change and the slope are
not statistical tests of change, and they are calculated only when there are
at least two (change) or three (slope) distinct dated time points.
"""

import numpy as np
import pandas as pd

from genomics.frequency_uncertainty import wilson_interval
from genomics.temporal_analysis import ALT_READS_NOT_CALLED, ZERO_ALT_READS
from genomics.unified_schema import VARIANT_OBSERVED


DETECTED = "Detected"
NON_DETECTION = "Supported non-detection (0 alternative reads)"
COUNTED_NOT_CALLED = "Not called; alternative reads in base counts"
UNKNOWN = "Unknown or uncallable"

POINT_COLUMNS = ["Variant_Key", "Location_Name", "Sample_ID", "Collection_Date", "State", "VAF",
                 "Alt_Read_Count", "Total_Depth", "Wilson_Lower", "Wilson_Upper", "Interval_Status"]
SUMMARY_COLUMNS = [
    "Variant_Key", "Location_Name", "Samples", "Dated_Samples", "Detections", "Supported_Non_Detections",
    "Counted_Not_Called", "Unknown_Or_Uncallable", "Known_Dated_Time_Points", "Points_With_Interval",
    "Earliest_Date", "Earliest_VAF", "Earliest_Lower", "Earliest_Upper", "Earliest_Depth",
    "Latest_Date", "Latest_VAF", "Latest_Lower", "Latest_Upper", "Latest_Depth",
    "Absolute_Change", "Slope_Per_Week", "Minimum_Depth", "Maximum_Depth", "Series_Status",
    "Interpretation_Note",
]
UNEVEN_DEPTH_RATIO = 5


def known_points(matrix, evidence):
    """One row per sample with a known value in a series (detections and supported non-detections)."""

    detected = matrix[matrix["Status"].eq(VARIANT_OBSERVED)][["Variant_Key", "Location_Name", "Sample_ID",
                                                     "Collection_Date"]]
    counted = evidence[evidence["Counted_Once"] & evidence["Analysis_Included"]][[
        "Sample_ID", "Variant_Key", "VAF", "Alt_Read_Count", "Total_Depth", "Wilson_Lower",
        "Wilson_Upper", "Interval_Status"]]
    detected = detected.merge(counted, on=["Sample_ID", "Variant_Key"], how="left")
    detected["State"] = DETECTED

    # Compared as a categorical column (a text copy of a large matrix costs ~100 MB).
    evidence_text = (matrix["Callable_Evidence"] if "Callable_Evidence" in matrix
                     else pd.Series("", index=matrix.index))
    from_counts = matrix[evidence_text.isin([ZERO_ALT_READS, ALT_READS_NOT_CALLED])]
    counts = from_counts[["Variant_Key", "Location_Name", "Sample_ID", "Collection_Date"]].copy()
    counts["VAF"] = from_counts["Callable_Frequency"].to_numpy(dtype=float)
    counts["Alt_Read_Count"] = from_counts["Alt_Read_Count"].to_numpy()
    counts["Total_Depth"] = pd.array(pd.to_numeric(from_counts["Depth"], errors="coerce"), dtype="Int64")
    lower, upper = wilson_interval(counts["Alt_Read_Count"].astype(float), counts["Total_Depth"].astype(float))
    counts["Wilson_Lower"], counts["Wilson_Upper"] = lower, upper
    counts["Interval_Status"] = "Calculated (base counts)"
    counts["State"] = np.where(evidence_text[from_counts.index].eq(ZERO_ALT_READS), NON_DETECTION,
                               COUNTED_NOT_CALLED)
    points = pd.concat([detected, counts], ignore_index=True)
    return points[POINT_COLUMNS].sort_values(["Variant_Key", "Location_Name", "Collection_Date",
                                              "Sample_ID"], kind="mergesort").reset_index(drop=True)


def series_note(row):
    notes = []
    if row["Known_Dated_Time_Points"] < 2:
        notes.append("A single time point is not a trend.")
    if row["Unknown_Or_Uncallable"]:
        notes.append(f"{int(row['Unknown_Or_Uncallable'])} sample(s) have no known value and are left out, "
                     "so the series may over-represent detections.")
    if row["Known_Dated_Time_Points"] >= 2 and row["Points_With_Interval"] < row["Known_Dated_Time_Points"]:
        notes.append("Some time points have no read-count interval.")
    if (pd.notna(row["Minimum_Depth"]) and row["Minimum_Depth"] > 0
            and row["Maximum_Depth"] / row["Minimum_Depth"] >= UNEVEN_DEPTH_RATIO):
        notes.append(f"Depth varies {row['Maximum_Depth'] / row['Minimum_Depth']:.0f}-fold between time "
                     "points; low-depth points are less precise (wider intervals).")
    return " ".join(notes) or "All samples in the series have a known value."


def series_summary(matrix, points):
    """Summarize each mutation at each site; series without any known value are left out.

    The earliest and latest values are the first and last dated samples (by
    date, then sample ID). The slope is an ordinary least-squares line through
    every dated known value, in frequency per 7 days, when there are at least
    three distinct dates. All of it is grouped arithmetic, so large tables stay fast.
    """

    keys = ["Variant_Key", "Location_Name"]
    if points.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)
    series = points[keys].drop_duplicates()
    in_series = matrix[keys + ["Collection_Date"]].merge(series, on=keys)
    table = in_series.groupby(keys).agg(Samples=("Collection_Date", "size"),
                                        Dated_Samples=("Collection_Date", "count"))
    states = points.assign(
        Detections=points["State"].eq(DETECTED).astype(int),
        Supported_Non_Detections=points["State"].eq(NON_DETECTION).astype(int),
        Counted_Not_Called=points["State"].eq(COUNTED_NOT_CALLED).astype(int),
        Known=1,
    ).groupby(keys)[["Detections", "Supported_Non_Detections", "Counted_Not_Called", "Known"]].sum()
    table = table.join(states)
    table["Unknown_Or_Uncallable"] = table["Samples"] - table["Known"]

    dated = points[points["Collection_Date"].notna()].copy()
    groups = dated.groupby(keys)
    table["Known_Dated_Time_Points"] = groups["Collection_Date"].nunique()
    table["Points_With_Interval"] = groups["Wilson_Lower"].count()
    depth = dated["Total_Depth"].astype(float)
    table["Minimum_Depth"] = depth.groupby([dated[key] for key in keys]).min()
    table["Maximum_Depth"] = depth.groupby([dated[key] for key in keys]).max()
    for prefix, ends in [("Earliest", groups.head(1)), ("Latest", groups.tail(1))]:
        ends = ends.set_index(keys)
        for column, source in [("Date", "Collection_Date"), ("VAF", "VAF"), ("Lower", "Wilson_Lower"),
                               ("Upper", "Wilson_Upper"), ("Depth", "Total_Depth")]:
            table[f"{prefix}_{column}"] = ends[source]
    table[["Known_Dated_Time_Points", "Points_With_Interval"]] = (
        table[["Known_Dated_Time_Points", "Points_With_Interval"]].fillna(0).astype(int))
    two_dates = table["Known_Dated_Time_Points"] >= 2
    table["Absolute_Change"] = (table["Latest_VAF"] - table["Earliest_VAF"]).where(two_dates)

    # Least-squares slope from grouped sums (identical to numpy.polyfit, degree 1).
    if not dated.empty:
        dated["X"] = (dated["Collection_Date"] - dated["Collection_Date"].min()).dt.days.astype(float)
        dated["Y"] = dated["VAF"].astype(float)
        dated["XX"], dated["XY"] = dated["X"] ** 2, dated["X"] * dated["Y"]
        sums = dated.groupby(keys)[["X", "Y", "XX", "XY"]].sum()
        n = dated.groupby(keys).size()
        spread = n * sums["XX"] - sums["X"] ** 2
        slope = (n * sums["XY"] - sums["X"] * sums["Y"]) / spread.where(spread > 0)
        table["Slope_Per_Week"] = (slope * 7).where(table["Known_Dated_Time_Points"] >= 3)
    else:
        table["Slope_Per_Week"] = np.nan
    table["Series_Status"] = np.select(
        [table["Known_Dated_Time_Points"] >= 3, table["Known_Dated_Time_Points"] == 2,
         table["Known_Dated_Time_Points"] == 1],
        ["Descriptive trend (3 or more dates)", "Change between two dates", "Single time point"],
        "No dated known value")
    table = table.reset_index()
    table["Interpretation_Note"] = [series_note(row) for row in table.to_dict("records")]
    return table[SUMMARY_COLUMNS]


def choose_series(summary, limit=6):
    """Series to draw: most dated known time points, then most intervals, then largest change."""

    eligible = summary[summary["Known_Dated_Time_Points"] >= 2].copy()
    if eligible.empty:
        return []
    eligible["Size"] = eligible["Absolute_Change"].abs().fillna(0)
    eligible = eligible.sort_values(["Known_Dated_Time_Points", "Points_With_Interval", "Size", "Variant_Key"],
                                    ascending=[False, False, False, True])
    return list(eligible[["Variant_Key", "Location_Name"]].head(limit).itertuples(index=False, name=None))


def unknown_dates(matrix, chosen):
    """Dates of samples with no known value, for the chosen series (drawn apart from values)."""

    wanted = pd.DataFrame(chosen, columns=["Variant_Key", "Location_Name"])
    rows = matrix.merge(wanted, on=["Variant_Key", "Location_Name"])
    value = rows["Callable_Frequency"] if "Callable_Frequency" in rows else pd.Series(np.nan, index=rows.index)
    return rows[value.isna() & rows["Collection_Date"].notna()][
        ["Variant_Key", "Location_Name", "Collection_Date", "Status"]]
