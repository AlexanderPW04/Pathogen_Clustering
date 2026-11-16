"""Sprint 5, Phase 1: depth-aware variant analysis for pooled pathogen sequencing.

This step runs after the Sprint 3 pooled analyses and reuses their tables:
the quality-controlled pooled variant records, the per-position table, the
reference verification, and the observation matrix. It adds:

- `sprint5_vaf_observations.csv`: one record per observation with the
  reported frequency, read counts and their source, VAF, a 95% Wilson
  interval when valid counts exist, QC, and every reason a value is missing;
- `sprint5_temporal_summary.csv`: each mutation at each site over time;
- `sprint5_missing_information.csv`: what the uploaded data lack and what
  that prevents;
- `sprint5_depth_dependence.csv`: the same VAF at several depths.

Phase 2 (`genomics/overdispersion_analysis.py`) adds, without changing the
Phase 1 columns:

- extra columns in `sprint5_vaf_observations.csv` (model, rho, rho source,
  calibration status, beta-binomial interval, limitations);
- `sprint5_uncertainty_models.csv`: the compact model comparison;
- `sprint5_overdispersion_sensitivity.csv`: hypothetical-rho intervals;
- `sprint5_overdispersion_calibration.csv` and `sprint5_calibration_sites.csv`
  when technical replicates (a shared Specimen_ID) are present.

Phase 3 (`genomics/mutation_group_analysis.py`) adds exploratory groups of
mutations with similar frequency patterns, in their own tables and figures.

Nothing here identifies lineages or tests for change. It measures how
strongly the reads support each observed frequency and shows how those
frequencies change over time, with missing evidence kept missing.
"""

import pandas as pd

import numpy as np

from genomics import figures
from genomics.beta_binomial import score_interval
from genomics.frequency_uncertainty import depth_dependence_table
from genomics.mutation_group_analysis import add_grouping_outputs
from genomics.overdispersion_analysis import (
    DEFAULT_SENSITIVITY_RHOS,
    HYPOTHETICAL_SOURCE,
    run_overdispersion_analysis,
)
from genomics.variant_evidence import (
    BASE_COUNTS,
    FREQUENCY_ONLY,
    NO_DEPTH,
    RECONSTRUCTED,
    SUPPLIED,
    build_variant_evidence,
    validate_evidence,
)
from genomics.vaf_trends import choose_series, known_points, series_summary, unknown_dates


TIME_SERIES_PANELS = 6


def availability_messages(evidence, details, samples, position_table):
    """Plain sentences saying exactly what is missing and what that prevents."""

    messages = []
    sources = evidence["Alt_Count_Source"].value_counts()
    if sources.get(NO_DEPTH, 0) == len(evidence):
        messages.append("Read depth is unavailable, so frequency uncertainty cannot be calculated.")
    elif sources.get(NO_DEPTH, 0):
        messages.append(f"Read depth is unavailable for {sources[NO_DEPTH]} observation(s); no uncertainty is "
                        "calculated for them.")
    measured = sources.get(SUPPLIED, 0) + sources.get(BASE_COUNTS, 0)
    if not measured and sources.get(RECONSTRUCTED, 0):
        messages.append("The data contain reported frequencies but no verified read counts. Alternate read "
                        f"counts were reconstructed for {sources[RECONSTRUCTED]} observation(s) as reported "
                        f"frequency x total depth (each count is the only whole number that reproduces the "
                        f"frequency at the table's {details['precision']}-significant-digit precision). They "
                        "are derived, not measured, so their intervals are approximate.")
    elif not measured and sources.get(FREQUENCY_ONLY, 0):
        reason = ""
        if details.get("not_reproduced"):
            reason = (f" Counts were not reconstructed: {details['not_reproduced']} record(s) cannot be "
                      "reproduced as whole reads / total depth, so the frequency and depth may not share a "
                      "denominator.")
        messages.append("The data contain reported frequencies but no verified read counts, so frequency "
                        "uncertainty cannot be calculated." + reason)
    if not samples["Collection_Date"].notna().any():
        messages.append("Collection dates are unavailable, so temporal trends cannot be calculated.")
    status = evidence["Reference_Status"].astype(str)
    unverified = status.isin(["Supplied but unverified", "No reference provided", "Not evaluated"])
    if unverified.all():
        messages.append("Reference identity is unverified: positions and reference bases were not checked "
                        "against a supplied reference FASTA.")
    elif unverified.any():
        messages.append(f"Reference identity is unverified for {int(unverified.sum())} of {len(evidence)} "
                        "observations (their reference label has no matching supplied FASTA).")
    mismatched = status.str.contains("mismatch", case=False)
    if mismatched.any():
        messages.append(f"{int(mismatched.sum())} observations are on a reference with Ref-base mismatches; "
                        "records that disagree fail quality control (see reference_verification.csv).")
    if position_table is None or not position_table.attrs.get("has_base_counts", False):
        messages.append("No per-position A/C/G/T counts were supplied, so an unreported mutation is unknown "
                        "rather than a confirmed non-detection; time series use detections only.")
    return messages


def missing_information(evidence, samples, position_table):
    """One row per requirement: how many observations or samples lack it, and the effect."""

    sources = evidence["Alt_Count_Source"]
    status = evidence["Interval_Status"].astype(str)
    has_counts = position_table is not None and position_table.attrs.get("has_base_counts", False)
    rows = [
        ("Read depth", int(sources.eq(NO_DEPTH).sum()), "observations",
         "No read count can be established, so no uncertainty interval."),
        ("Measured read counts (supplied or base counts)",
         int((~sources.isin([SUPPLIED, BASE_COUNTS])).sum()), "observations",
         "Reconstructed counts give approximate intervals; frequency-only records get none."),
        ("Counts that agree with the reported frequency", int(status.str.contains("disagree").sum()),
         "observations", "No interval; both values are kept for review."),
        ("Consistent allele frequencies at a site (sum at most 1)", int(status.str.contains("sum above").sum()),
         "observations", "No interval for the alleles at that site."),
        ("Single record per sample and mutation", int(status.str.contains("repeated").sum()),
         "observations", "Counted once; conflicting repeats get no interval."),
        ("Passing quality control", int((~evidence["Analysis_Included"]).sum()), "observations",
         "Kept in the table with the reason; not used in time series."),
        ("Usable collection date", int(samples["Collection_Date"].isna().sum()), "pooled samples",
         "Not placed in time (missing, partial, or ambiguous dates are not guessed)."),
        ("Per-position A/C/G/T counts", 0 if has_counts else len(samples), "pooled samples",
         "Without them an unreported mutation is unknown, never a confirmed 0."),
        ("Verified reference", int((~evidence["Reference_Status"].eq("Verified compatible")).sum()),
         "observations", "Positions and reference bases are as reported (unverified)."),
    ]
    return pd.DataFrame(rows, columns=["Requirement", "Lacking", "Counted_As", "Effect"])


def run_depth_aware_analysis(result, pooled, matrix, reference, settings, add_figure, note):
    """Add the Sprint 5 tables, figures, notes, and summary to `result`."""

    variants, samples = pooled["variants"], pooled["samples"]
    evidence, details = build_variant_evidence(
        variants, pooled["position_table"], pooled.get("reference_status"), reference,
        pooled.get("input_format", ""), settings["min_variant_depth"])
    for problem in validate_evidence(evidence):
        note("depth_aware", "Internal consistency check failed: " + problem)
    # Phase 2: technical replicates share a Specimen_ID (blank when not supplied).
    # Phase 3: distinct Library_IDs document that the replicates are separate libraries.
    for name in ["Specimen_ID", "Library_ID"]:
        evidence[name] = (variants[name].fillna("").astype(str).to_numpy() if name in variants else "")
    uncertainty = run_overdispersion_analysis(evidence, pooled["position_table"],
                                              settings.get("overdispersion_rhos", DEFAULT_SENSITIVITY_RHOS))
    evidence = uncertainty["evidence"]
    # Phase 3's Library_ID goes last so every Phase 1-2 column keeps its position.
    evidence = evidence[[column for column in evidence.columns if column != "Library_ID"] + ["Library_ID"]]
    included_samples = samples[samples["Analysis_Included"]]
    messages = availability_messages(evidence, details, included_samples, pooled["position_table"])
    for message in messages:
        note("depth_aware", message)

    tables = result["tables"]
    tables["sprint5_vaf_observations.csv"] = evidence
    tables["sprint5_missing_information.csv"] = missing_information(evidence, included_samples,
                                                                     pooled["position_table"])
    demonstration = depth_dependence_table()
    tables["sprint5_depth_dependence.csv"] = demonstration
    tables["sprint5_uncertainty_models.csv"] = uncertainty["comparison"]
    if len(uncertainty["sensitivity"]):
        tables["sprint5_overdispersion_sensitivity.csv"] = uncertainty["sensitivity"]
    if len(uncertainty["calibration"]):
        tables["sprint5_overdispersion_calibration.csv"] = uncertainty["calibration"]
    if len(uncertainty["sites"]):
        tables["sprint5_calibration_sites.csv"] = uncertainty["sites"]

    points = known_points(matrix, evidence) if not matrix.empty else pd.DataFrame()
    rho_for_sample = uncertainty["rho_for_sample"]
    if rho_for_sample and not points.empty:
        # Calibrated beta-binomial bars for the replicated samples only.
        rho = points["Sample_ID"].map(rho_for_sample).astype(float).to_numpy()
        lower, upper = score_interval(np.where(np.isfinite(rho), points["Alt_Read_Count"].astype("Float64")
                                               .to_numpy(dtype=float, na_value=np.nan), np.nan),
                                      points["Total_Depth"].astype("Float64").to_numpy(dtype=float, na_value=np.nan),
                                      np.where(np.isfinite(rho), rho, 0.0))
        points["BetaBinomial_Lower"], points["BetaBinomial_Upper"] = lower, upper
    summary = series_summary(matrix, points) if not points.empty else pd.DataFrame()
    if not summary.empty:
        tables["sprint5_temporal_summary.csv"] = summary
    chosen = choose_series(summary, TIME_SERIES_PANELS) if not summary.empty else []

    add_figure(result, "depth_dependence", "Same frequency, different sequencing depth",
               "An illustration with exact counts, not uploaded data: each row is a VAF of 0.20 "
               "observed with a different number of reads. The 95% Wilson interval narrows as depth "
               "grows. It reflects read sampling only, not sequencing or amplification error.",
               lambda: figures.depth_dependence_figure(demonstration))
    counted = evidence[evidence["Counted_Once"]]
    with_interval = counted["Wilson_Lower"].notna()
    if with_interval.any():
        approximate = counted.loc[with_interval, "Interval_Status"].str.startswith("Approximate").any()
        add_figure(result, "vaf_depth", "Sequencing depth and frequency uncertainty",
                   "Each point is one observation with a read-count interval. Low-depth observations have "
                   "wide intervals and should not carry the same weight as high-depth ones; depth alone "
                   "does not show that a variant call is correct."
                   + (" Counts here are reconstructed from reported frequencies, so widths are approximate."
                      if approximate else ""),
                   lambda: figures.vaf_depth_figure(counted, settings["min_variant_depth"]))
    if chosen:
        gaps = unknown_dates(matrix, chosen)
        add_figure(result, "vaf_time_series", "Mutation frequencies over time with uncertainty",
                   f"Up to {TIME_SERIES_PANELS} mutation-site series with the most dated known values "
                   "(then the most intervals and the largest change); every series is in "
                   "sprint5_temporal_summary.csv. Bars are 95% read-sampling intervals"
                   + (" (approximate: counts reconstructed from reported frequencies)"
                      if details.get("used") and details.get("unique") else "")
                   + ". Gray crosses are samples with no known value; they are not zeros."
                   + (" Wide pale bands are calibrated beta-binomial intervals, shown only for samples with "
                      "technical replicates." if rho_for_sample else ""),
                   lambda: figures.vaf_time_series_figure(points, gaps, chosen))
    elif included_samples["Collection_Date"].notna().any():
        note("depth_aware", "No mutation has known values on two or more dates at one site, so no time "
                            "series is drawn.")

    phase2 = uncertainty["summary"]
    if with_interval.any():
        calibrated_rhos = sorted(set(rho_for_sample.values()))
        hypothetical = [rho for rho, source in uncertainty["scenarios"] if source == HYPOTHETICAL_SOURCE]
        reference_rho = (0.02 if 0.02 in hypothetical else (hypothetical[0] if hypothetical else None))
        add_figure(result, "uncertainty_models", "Binomial versus beta-binomial interval width",
                   "Left: the 95% interval width at a frequency of 0.20. The binomial (Wilson) interval keeps "
                   "narrowing with depth; with library-to-library variation (rho > 0) it cannot fall below a "
                   "floor set by rho. Right: each uploaded observation's binomial width against its "
                   "overdispersion-aware width"
                   + ((" (" + "; ".join(
                       (["blue: calibrated from technical replicates"] if calibrated_rhos else [])
                       + ([f"gray: hypothetical rho {reference_rho:g}, not an estimate"]
                          if reference_rho is not None else [])) + ")")
                      if calibrated_rhos or reference_rho is not None else "")
                   + ". Hypothetical curves are scenarios, never estimates from this upload.",
                   lambda: figures.uncertainty_models_figure(evidence, uncertainty["scenarios"], calibrated_rhos,
                                                             reference_rho))
    if len(uncertainty["sites"]) and uncertainty["sites"]["Used"].any():
        add_figure(result, "replicate_concordance", "Technical replicate agreement (overdispersion calibration)",
                   "Each point is one allele measured in two separately prepared libraries of the same specimen "
                   "(sites with depth above 100 in both and mean frequency 1% to 90%). Spread around the dotted "
                   "line beyond what read counts explain is the overdispersion rho; the range is a site "
                   "bootstrap. Biological differences between specimens or dates are not part of this estimate.",
                   lambda: figures.replicate_concordance_figure(uncertainty["sites"], uncertainty["calibration"]))

    # Phase 3: exploratory mutation-frequency groups (optional; never changes the results above).
    grouping = add_grouping_outputs(result, evidence, matrix, uncertainty["calibration"], included_samples,
                                    pooled["position_table"], settings, add_figure, note)

    sources = counted["Alt_Count_Source"].value_counts()
    eligible = summary[summary["Known_Dated_Time_Points"] >= 2] if not summary.empty else summary
    result["depth_aware"] = {
        "messages": messages,
        "reconstruction": details,
        "samples": int(evidence["Sample_ID"].nunique()),
        "observations": len(evidence),
        "observations_counted": len(counted),
        "with_depth": int(counted["Total_Depth"].notna().sum()),
        "supplied_counts": int(sources.get(SUPPLIED, 0) + sources.get(BASE_COUNTS, 0)),
        "reconstructed_counts": int(sources.get(RECONSTRUCTED, 0)),
        "with_interval": int(with_interval.sum()),
        "missing_required": int((~with_interval).sum()),
        "series": len(summary),
        "eligible_series": len(eligible),
        "eligible_variants": int(eligible["Variant_Key"].nunique()) if len(eligible) else 0,
        "trend_series": int((summary["Known_Dated_Time_Points"] >= 3).sum()) if len(summary) else 0,
        "uncertainty": phase2,
        "grouping": grouping,
    }
    return result["depth_aware"]
