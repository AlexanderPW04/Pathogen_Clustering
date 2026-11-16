"""Sprint 5 Phase 2: sequencing uncertainty models for pooled variant frequencies.

Runs inside the Sprint 5 depth-aware analysis, on the validated evidence
table (`genomics/variant_evidence.py`). Three kinds of uncertainty result are
kept strictly apart:

1. Binomial baseline (Phase 1): the 95% Wilson interval, read sampling only.
   Always kept, never overwritten.
2. Empirically calibrated beta-binomial interval: only for samples that are
   one of two technical replicates of a specimen whose overdispersion rho was
   estimated from the replicates (`genomics/overdispersion.py`). The interval
   uses the upper 95% bound of rho, so the uncertainty of the estimate is
   included (simulation study C3).
3. Hypothetical sensitivity scenarios: the same beta-binomial interval with
   rho values typed by the user (or the defaults below). They show how much
   wider the uncertainty WOULD be if the libraries had that overdispersion.
   They are never called estimates or calibrated results.

An upload without technical replicates (for example the published
mmc2.xlsx) gets no rho estimate at all.
"""

import numpy as np
import pandas as pd

from genomics import beta_binomial as bb
from genomics.overdispersion import (
    ESTIMATED_ZERO,
    USABLE_STATUSES,
    calibrate_replicates,
)


DEFAULT_SENSITIVITY_RHOS = [0.002, 0.02, 0.05]
MAX_SENSITIVITY_VALUES = 5

# Per-observation labels are kept short: the observation table of a large workbook must stay under
# the 8 MB download limit. The full wording is in the page, the report, and the docstrings.
BASELINE_MODEL = "Binomial (Wilson)"
CALIBRATED_MODEL = "Beta-binomial (calibrated)"
SENSITIVITY_MODEL = "Beta-binomial (hypothetical rho)"
NOT_CALIBRATED_MODEL = "None (not calibrated)"
INTERVAL_TYPE = "Score interval for mean VAF (approximate)"

CALIBRATED = "Calibrated (technical replicates)"
NO_REPLICATES = "Not calibrated: no replicates in upload"
NO_REPLICATE_FOR_SAMPLE = "Not calibrated: sample not replicated"
CALIBRATION_FAILED = "Not calibrated: replicate estimate unusable"
NO_INTERVAL = "No interval"
READ_SAMPLING_ONLY = "Read sampling only"
WITH_LIBRARY_VARIATION = "Read sampling + library variation; one rho per library"

HYPOTHETICAL_SOURCE = "Hypothetical (user-specified sensitivity value; not estimated from the data)"
BORROWED_SOURCE = ("Borrowed from this upload's replicate calibration (largest upper bound; assumed, "
                   "not calibrated, for samples without their own replicates)")

EXPLANATION = ("Sequencing depth affects how precisely a mutation's frequency can be estimated. Additional "
               "sequencing variation may make the uncertainty larger than expected from read counts alone.")
NO_CALIBRATION_TEXT = ("Technical overdispersion could not be estimated from this dataset. The binomial "
                       "interval is shown as a baseline, and any overdispersion sensitivity analysis is "
                       "hypothetical.")
GENERAL_LIMITATION = ("Covers read sampling and, for beta-binomial results, library-to-library variation "
                      "only; not primer bias shared by every library, mapping or reference errors, "
                      "contamination, or biological differences between samples or dates.")

OVERDISPERSION_COLUMNS = [
    "Specimen_ID", "Baseline_Model", "Overdispersion_Model", "Rho", "Rho_Source", "Calibration_Status",
    "BetaBinomial_Lower", "BetaBinomial_Upper", "BetaBinomial_Width", "BetaBinomial_Interval_Type",
    "Uncertainty_Limitations",
]
SENSITIVITY_COLUMNS = [
    "Sample_ID", "Variant_Key", "VAF", "Alt_Read_Count", "Total_Depth", "Alt_Count_Source",
    "Assumed_Rho", "Rho_Source", "Sensitivity_Lower", "Sensitivity_Upper", "Width_Ratio_To_Wilson",
]
HYPOTHETICAL_LABEL = "Hypothetical; not calibrated"
BORROWED_LABEL = "Borrowed from replicates; not calibrated"
COMPARISON_COLUMNS = ["Model", "Rho", "Rho_Source", "Calibration", "Observations_With_Interval",
                      "Median_Depth", "Median_Width", "Median_Width_Ratio_To_Wilson", "Interpretation"]


def parse_rho_values(text):
    """Parse the sensitivity setting ("0.002, 0.02") into a sorted list; '' means none.

    Raises ValueError with a message a user can act on.
    """

    if text is None:
        return list(DEFAULT_SENSITIVITY_RHOS)
    items = [item.strip() for item in str(text).replace(";", ",").split(",") if item.strip()]
    if len(items) > MAX_SENSITIVITY_VALUES:
        raise ValueError(f"Enter at most {MAX_SENSITIVITY_VALUES} hypothetical overdispersion values.")
    values = []
    for item in items:
        try:
            value = bb.check_rho(item)
        except ValueError:
            raise ValueError(f"Hypothetical overdispersion values must be numbers from 0 up to (not "
                             f"including) 1; '{item}' is not.") from None
        if value > 0 and value not in values:
            values.append(value)
    return sorted(values)


def _limitations(source, calibrated):
    """Short per-observation limitation; GENERAL_LIMITATION applies to every row and is shown on the page."""

    text = WITH_LIBRARY_VARIATION if calibrated else READ_SAMPLING_ONLY
    return text + ("; reconstructed counts (approximate)" if source == "Reconstructed" else "")


def add_overdispersion_columns(evidence, calibration):
    """Append the Phase 2 columns to the evidence table (the Phase 1 columns are unchanged)."""

    table = evidence.copy()
    if "Specimen_ID" not in table:
        table["Specimen_ID"] = ""
    usable = calibration[calibration["Status"].isin(USABLE_STATUSES)] if len(calibration) else calibration
    rho_for_sample, specimen_of = {}, {}
    for row in usable.to_dict("records"):
        for sample in (row["Replicate_1"], row["Replicate_2"]):
            rho_for_sample[sample] = float(row["Rho_Used_For_Intervals"])
            specimen_of[sample] = row["Specimen_ID"]
    failed_samples = set()
    if len(calibration):
        for row in calibration[~calibration["Status"].isin(USABLE_STATUSES)].to_dict("records"):
            failed_samples.update(value for value in (row["Replicate_1"], row["Replicate_2"]) if value)

    has_wilson = table["Wilson_Lower"].notna().to_numpy()
    samples = table["Sample_ID"].astype(str).to_numpy()
    rho = np.array([rho_for_sample.get(sample, np.nan) for sample in samples])
    calibrated = has_wilson & np.isfinite(rho)
    lower, upper = bb.score_interval(np.where(calibrated, table["Alt_Read_Count"].astype("Float64").to_numpy(
        dtype=float, na_value=np.nan), np.nan), table["Total_Depth"].astype("Float64").to_numpy(
        dtype=float, na_value=np.nan), np.where(calibrated, rho, 0.0))

    status = []
    for index, sample in enumerate(samples):
        if not has_wilson[index]:
            status.append(f"{NO_INTERVAL} (see Interval_Status)")
        elif calibrated[index]:
            status.append(CALIBRATED)
        elif sample in failed_samples:
            status.append(CALIBRATION_FAILED)
        elif rho_for_sample or failed_samples:
            status.append(NO_REPLICATE_FOR_SAMPLE)
        else:
            status.append(NO_REPLICATES)
    sources = table["Alt_Count_Source"].astype(str).to_numpy()
    table["Baseline_Model"] = np.where(has_wilson, BASELINE_MODEL, "None (no valid counts)")
    table["Overdispersion_Model"] = np.where(calibrated, CALIBRATED_MODEL, NOT_CALIBRATED_MODEL)
    table["Rho"] = np.where(calibrated, rho, np.nan)
    table["Rho_Source"] = [
        (f"Upper 95% bound, replicates of {specimen_of[sample]}" if flag else "")
        for sample, flag in zip(samples, calibrated)]
    table["Calibration_Status"] = status
    table["BetaBinomial_Lower"] = np.where(calibrated, lower, np.nan)
    table["BetaBinomial_Upper"] = np.where(calibrated, upper, np.nan)
    table["BetaBinomial_Width"] = table["BetaBinomial_Upper"] - table["BetaBinomial_Lower"]
    table["BetaBinomial_Interval_Type"] = np.where(calibrated, INTERVAL_TYPE, "")
    table["Uncertainty_Limitations"] = [_limitations(source, flag) if wilson else ""
                                        for source, flag, wilson in zip(sources, calibrated, has_wilson)]
    return table


def sensitivity_scenarios(rhos, calibration):
    """(rho, source) pairs: the user's hypothetical values, plus a borrowed value when it is meaningful."""

    scenarios = [(value, HYPOTHETICAL_SOURCE) for value in rhos]
    if len(calibration):
        usable = calibration[calibration["Status"].isin(USABLE_STATUSES)]
        if len(usable):
            borrowed = float(usable["Rho_Used_For_Intervals"].max())
            if borrowed > 0 and borrowed not in rhos:
                scenarios.append((borrowed, BORROWED_SOURCE))
    return scenarios


def sensitivity_table(evidence, scenarios):
    """One row per observation with a Wilson interval and per scenario (long format)."""

    base = evidence[evidence["Wilson_Lower"].notna() & evidence["Counted_Once"]]
    frames = []
    alt = base["Alt_Read_Count"].astype("Float64").to_numpy(dtype=float, na_value=np.nan)
    depth = base["Total_Depth"].astype("Float64").to_numpy(dtype=float, na_value=np.nan)
    wilson_width = (base["Wilson_Upper"] - base["Wilson_Lower"]).to_numpy(dtype=float)
    for rho, source in scenarios:
        lower, upper = bb.score_interval(alt, depth, rho)
        frame = base[["Sample_ID", "Variant_Key", "VAF", "Alt_Read_Count", "Total_Depth",
                      "Alt_Count_Source"]].copy()
        frame["Assumed_Rho"] = rho
        frame["Rho_Source"] = HYPOTHETICAL_LABEL if source == HYPOTHETICAL_SOURCE else BORROWED_LABEL
        frame["Sensitivity_Lower"], frame["Sensitivity_Upper"] = lower, upper
        with np.errstate(divide="ignore", invalid="ignore"):
            frame["Width_Ratio_To_Wilson"] = (upper - lower) / wilson_width
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=SENSITIVITY_COLUMNS)
    return pd.concat(frames, ignore_index=True)[SENSITIVITY_COLUMNS]


def model_comparison(evidence, sensitivity, scenarios):
    """The compact model table shown first on the page."""

    counted = evidence[evidence["Counted_Once"]]
    wilson = counted[counted["Wilson_Lower"].notna()]
    wilson_width = wilson["Wilson_Upper"] - wilson["Wilson_Lower"]
    rows = [{"Model": BASELINE_MODEL, "Rho": 0.0, "Rho_Source": "Assumed by the binomial model",
             "Calibration": "Not applicable (read sampling only)", "Observations_With_Interval": len(wilson),
             "Median_Depth": float(wilson["Total_Depth"].astype(float).median()) if len(wilson) else np.nan,
             "Median_Width": float(wilson_width.median()) if len(wilson) else np.nan,
             "Median_Width_Ratio_To_Wilson": 1.0 if len(wilson) else np.nan,
             "Interpretation": "Baseline: the narrowest plausible interval; too narrow if libraries vary."}]
    calibrated = counted[counted["Calibration_Status"].eq(CALIBRATED)]
    if len(calibrated):
        ratio = calibrated["BetaBinomial_Width"] / (calibrated["Wilson_Upper"] - calibrated["Wilson_Lower"])
        rhos = calibrated["Rho"].astype(float)
        rows.append({"Model": CALIBRATED_MODEL,
                     "Rho": float(rhos.median()) if rhos.nunique() == 1 else np.nan,
                     "Rho_Source": (f"Upper 95% bounds, {rhos.min():.4g} to {rhos.max():.4g}"
                                    if rhos.nunique() > 1 else "Upper 95% bound of the replicate estimate"),
                     "Calibration": "Empirical (technical replicates)",
                     "Observations_With_Interval": len(calibrated),
                     "Median_Depth": float(calibrated["Total_Depth"].astype(float).median()),
                     "Median_Width": float(calibrated["BetaBinomial_Width"].median()),
                     "Median_Width_Ratio_To_Wilson": float(ratio.median()),
                     "Interpretation": "Supported result for the replicated samples only."})
    for rho, source in scenarios:
        rows_for = sensitivity[sensitivity["Assumed_Rho"].eq(rho)
                               & sensitivity["Rho_Source"].eq(HYPOTHETICAL_LABEL if source == HYPOTHETICAL_SOURCE
                                                              else BORROWED_LABEL)]
        widths = rows_for["Sensitivity_Upper"] - rows_for["Sensitivity_Lower"]
        rows.append({"Model": SENSITIVITY_MODEL if source == HYPOTHETICAL_SOURCE
                     else "Beta-binomial (borrowed rho)", "Rho": rho,
                     "Rho_Source": "Hypothetical value" if source == HYPOTHETICAL_SOURCE else "Borrowed (assumed)",
                     "Calibration": "Not calibrated (sensitivity scenario)",
                     "Observations_With_Interval": len(rows_for),
                     "Median_Depth": float(rows_for["Total_Depth"].astype(float).median()) if len(rows_for) else np.nan,
                     "Median_Width": float(widths.median()) if len(rows_for) else np.nan,
                     "Median_Width_Ratio_To_Wilson": (float(rows_for["Width_Ratio_To_Wilson"].median())
                                                      if len(rows_for) else np.nan),
                     "Interpretation": "What the uncertainty would be if the libraries had this rho."})
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def example_observations(evidence, scenarios, limit=6):
    """A few observations side by side: frequency, depth, basic and overdispersion-aware intervals."""

    counted = evidence[evidence["Counted_Once"] & evidence["Wilson_Lower"].notna()]
    calibrated = counted[counted["Calibration_Status"].eq(CALIBRATED)]
    chosen = calibrated if len(calibrated) else counted
    if not len(chosen):
        return []
    # Spread the examples over depth so the effect of depth is visible.
    chosen = chosen.sort_values("Total_Depth", kind="mergesort")
    picks = chosen.iloc[np.unique(np.linspace(0, len(chosen) - 1, min(limit, len(chosen))).astype(int))]
    reference = next(((rho, source) for rho, source in scenarios if source == HYPOTHETICAL_SOURCE),
                     scenarios[0] if scenarios else None)
    rows = []
    for row in picks.to_dict("records"):
        if row["Calibration_Status"] == CALIBRATED:
            aware = (f"{row['BetaBinomial_Lower']:.3f} to {row['BetaBinomial_Upper']:.3f} "
                     f"(calibrated, rho {row['Rho']:.3g})")
            calibrated_text = "Yes (technical replicates)"
        elif reference is not None:
            lower, upper = bb.score_interval([row["Alt_Read_Count"]], [row["Total_Depth"]], reference[0])
            aware = f"{lower[0]:.3f} to {upper[0]:.3f} (hypothetical rho {reference[0]:g})"
            calibrated_text = "No: " + row["Calibration_Status"].split(": ", 1)[-1]
        else:
            aware = "Not available (no calibration and no sensitivity value)"
            calibrated_text = "No: " + row["Calibration_Status"].split(": ", 1)[-1]
        rows.append({"Sample": row["Sample_ID"], "Mutation": str(row["Variant_Key"]).split(":", 1)[-1],
                     "Frequency": f"{row['VAF']:.3f}", "Depth": int(row["Total_Depth"]),
                     "Basic interval": f"{row['Wilson_Lower']:.3f} to {row['Wilson_Upper']:.3f}",
                     "Overdispersion-aware interval": aware, "Calibrated": calibrated_text})
    return rows


def run_overdispersion_analysis(evidence, position_table, rhos):
    """Return a dict with the extended evidence and the Phase 2 tables and summary."""

    calibration, sites = calibrate_replicates(evidence, position_table)
    extended = add_overdispersion_columns(evidence, calibration)
    scenarios = sensitivity_scenarios(rhos, calibration)
    sensitivity = sensitivity_table(extended, scenarios)
    comparison = model_comparison(extended, sensitivity, scenarios)
    usable = calibration[calibration["Status"].isin(USABLE_STATUSES)] if len(calibration) else calibration
    counted = extended[extended["Counted_Once"]]
    status_counts = counted["Calibration_Status"].value_counts()
    if len(usable):
        zero = int(usable["Status"].eq(ESTIMATED_ZERO).sum())
        headline = (f"Technical overdispersion was estimated for {len(usable)} specimen(s) with technical "
                    f"replicates (rho {usable['Rho_Estimate'].min():.3g} to {usable['Rho_Estimate'].max():.3g}"
                    + (f"; {zero} consistent with independent read sampling" if zero else "")
                    + f"). Calibrated beta-binomial intervals are given for {int(status_counts.get(CALIBRATED, 0))} "
                    "observation(s) of those replicated samples only; all other observations keep the binomial "
                    "baseline, and any sensitivity analysis for them is hypothetical.")
    else:
        headline = NO_CALIBRATION_TEXT
    reasons = []
    if len(calibration):
        for row in calibration[~calibration["Status"].isin(USABLE_STATUSES)].to_dict("records"):
            reasons.append(f"Specimen {row['Specimen_ID']}: {row['Status']}. {row['Status_Detail']}")
    elif "Specimen_ID" not in evidence or not evidence["Specimen_ID"].fillna("").astype(str).str.strip().ne("").any():
        reasons.append("No Specimen_ID column links technical replicates (separate libraries of one specimen, each with "
                       "its own Library_ID), "
                       "so the extra variation between libraries cannot be measured. Different dates, sites, or "
                       "mutations are never used in their place, because real biological differences would be "
                       "mistaken for sequencing noise.")
    if not rhos:
        reasons.append("No hypothetical overdispersion values were entered, so no sensitivity analysis was run.")
    summary = {
        "explanation": EXPLANATION,
        "headline": headline,
        "calibrated_specimens": len(usable),
        "specimens_with_replicates": len(calibration),
        "calibrated_observations": int(status_counts.get(CALIBRATED, 0)),
        "sensitivity_rhos": [rho for rho, source in scenarios if source == HYPOTHETICAL_SOURCE],
        "borrowed_rho": next((rho for rho, source in scenarios if source == BORROWED_SOURCE), None),
        "reasons": reasons,
        "comparison": comparison,
        "examples": example_observations(extended, scenarios),
    }
    return {"evidence": extended, "calibration": calibration, "sites": sites, "sensitivity": sensitivity,
            "comparison": comparison, "scenarios": scenarios, "summary": summary,
            "rho_for_sample": {sample: float(row["Rho_Used_For_Intervals"]) for row in usable.to_dict("records")
                               for sample in (row["Replicate_1"], row["Replicate_2"])}}
