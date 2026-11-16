"""Sprint 5 Phase 3: mutation-frequency profiles across compatible pooled samples.

A *profile* is one mutation's known frequencies across observation units.
It is the input to the exploratory grouping in `genomics/mutation_grouping.py`.

Mutation identity
-----------------
Reference_ID + position + Ref + Alt, the application's `Variant_Key`. The
Reference_ID names one reference sequence, and a segmented genome has one
reference per segment, so `Segment` is the Reference_ID. Two mutations are
compared only when they share a Reference_ID: positions on different
references or segments are never treated as homologous.

Observation units
-----------------
One unit is one biological specimen measurement:

- a single library (one `Sample_ID`); or
- a **verified** technical-replicate pair: two samples that share a
  Specimen_ID, date, location, and reference and document distinct
  Library_IDs (`overdispersion.replicate_pairs`). The two libraries measure
  the same specimen, so their counts are combined (x1 + x2 reads out of
  n1 + n2) instead of counting the specimen twice.

Samples from the same site and date that are NOT verified replicates stay
separate units. Each is flagged `Same_Site_Date_Units` so this is visible,
and nothing is silently collapsed.

States (the matrix status, never inferred)
------------------------------------------
- Known values: a detection (reported call), a supported non-detection
  (A/C/G/T counts show 0 alternative reads at adequate depth), or
  alternative reads counted in A/C/G/T data but not called.
- Unknown: coverage unknown, low, or absent; adequate depth without base
  counts (value unknown, not 0); or a call rejected by QC. Unknown cells are
  counted per mutation and never given a value. A missing variant row is
  never a zero.

Count provenance per unit
-------------------------
- "Measured": supplied alternate counts or A/C/G/T counts.
- "Reconstructed (approximate)": frequency x depth (Phase 1 rule).
- "Frequency only": a reported VAF without usable counts. It is used by the
  ordinary (frequency-only) method, never by the depth-aware method.

Memory
------
Only known values are stored, in long form (one row per mutation and unit).
The dense array that grouping needs is built later, for the eligible
mutations only.
"""

import numpy as np
import pandas as pd

from genomics.overdispersion import USABLE_STATUSES, replicate_pairs
from genomics.temporal_analysis import ALT_READS_NOT_CALLED, VALUE_UNKNOWN, ZERO_ALT_READS
from genomics.unified_schema import (
    COVERAGE_UNKNOWN,
    COVERED_NOT_DETECTED,
    EXCLUDED_BY_QC,
    LOW_COVERAGE,
    NOT_COVERED,
    VARIANT_OBSERVED,
)
from genomics.vaf_trends import COUNTED_NOT_CALLED, DETECTED, NON_DETECTION

Z_95 = 1.959963984540054

MEASURED = "Measured"
RECONSTRUCTED = "Reconstructed (approximate)"
FREQUENCY_ONLY = "Frequency only"
SINGLE_LIBRARY = "Single library"
REPLICATE_UNIT = "Specimen: 2 verified technical-replicate libraries, counts combined"

PROFILE_COLUMNS = [
    "Variant_Key", "Reference_ID", "Segment", "Position", "Ref", "Alt", "Unit_ID", "Unit_Type",
    "Sample_IDs", "Library_IDs", "Libraries_Known", "Collection_Date", "Location_Name",
    "Same_Site_Date_Units", "State", "VAF", "Alt_Read_Count", "Total_Depth", "Depth_Square_Sum",
    "Count_Basis", "Wilson_Lower", "Wilson_Upper", "Rho", "Rho_Source", "BetaBinomial_Lower",
    "BetaBinomial_Upper", "Reference_Status",
]
ELIGIBILITY_COLUMNS = [
    "Variant_Key", "Reference_ID", "Segment", "Position", "Ref", "Alt", "Annotated_Region",
    "Reference_Status", "Samples_On_Reference", "Detected", "Supported_Non_Detection",
    "Counted_Not_Called", "Unknown_Coverage", "Low_Or_No_Coverage", "Covered_Value_Unknown",
    "Rejected_By_QC", "Known_Units", "Known_Units_With_Counts", "Measured_Units", "Reconstructed_Units",
    "Frequency_Only_Units", "Eligible_Ordinary", "Eligible_Depth_Aware", "Exclusion_Reason",
]


def effective_score_interval(vaf, effective_depth, z=Z_95):
    """Wilson-form score interval for a frequency at an effective depth (arrays; NaN where invalid)."""

    p = np.asarray(vaf, dtype=float)
    n = np.asarray(effective_depth, dtype=float)
    valid = np.isfinite(p) & np.isfinite(n) & (n > 0) & (p >= 0) & (p <= 1)
    p0, n0 = np.where(valid, p, 0.0), np.where(valid, n, 1.0)
    z2 = z * z
    denominator = 1 + z2 / n0
    centre = (p0 + z2 / (2 * n0)) / denominator
    half = (z / denominator) * np.sqrt(p0 * (1 - p0) / n0 + z2 / (4 * n0 ** 2))
    lower = np.where(p0 == 0, 0.0, np.clip(centre - half, 0, 1))
    upper = np.where(p0 == 1, 1.0, np.clip(centre + half, 0, 1))
    return np.where(valid, lower, np.nan), np.where(valid, upper, np.nan)


def variance_factor(depth, depth_square_sum, rho):
    """c such that Var(VAF) = p (1 - p) c for a unit of one or more libraries (arrays).

    For libraries of depth n_l with the same rho, combined VAF = sum x_l / n has
    Var = p (1 - p) [n + rho (sum n_l^2 - n)] / n^2. With one library this is
    p (1 - p) [1 + (n - 1) rho] / n (Phase 2); with rho = 0 it is binomial, 1 / n.
    """

    n = np.asarray(depth, dtype=float)
    q = np.asarray(depth_square_sum, dtype=float)
    r = np.asarray(rho, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (n + r * (q - n)) / n ** 2


def observation_units(evidence, samples=None):
    """Return (unit table, {Sample_ID: Unit_ID}, verified pair list).

    Only verified replicate pairs (no problem from `replicate_pairs`) are combined.
    `samples` (optional) lists every included Sample_ID, including any without records.
    """

    info = evidence.drop_duplicates("Sample_ID")[["Sample_ID", "Collection_Date", "Location_Name"]]
    if samples is not None:
        extra = samples[~samples["Sample_ID"].isin(info["Sample_ID"])][["Sample_ID", "Collection_Date",
                                                                         "Location_Name"]]
        info = pd.concat([info, extra], ignore_index=True)
    libraries = {}
    if "Library_ID" in evidence:
        libraries = (evidence.assign(Library_ID=evidence["Library_ID"].fillna("").astype(str))
                     .groupby("Sample_ID")["Library_ID"].first().to_dict())
    verified = [(specimen, list(pair)) for specimen, pair, problem in replicate_pairs(evidence) if not problem]
    unit_of = {sample: sample for sample in info["Sample_ID"]}
    unit_type = {sample: SINGLE_LIBRARY for sample in info["Sample_ID"]}
    used_names = set(unit_of)
    for specimen, pair in verified:
        name = f"{specimen} [2 libraries]"
        while name in used_names:
            name += "*"
        used_names.add(name)
        for sample in pair:
            unit_of[sample] = name
        unit_type[name] = REPLICATE_UNIT
    rows = []
    for unit, group in info.assign(Unit_ID=info["Sample_ID"].map(unit_of)).groupby("Unit_ID", sort=True):
        members = sorted(group["Sample_ID"])
        rows.append({"Unit_ID": unit, "Unit_Type": unit_type.get(unit, SINGLE_LIBRARY),
                     "Sample_IDs": "; ".join(members),
                     "Library_IDs": "; ".join(libraries.get(sample, "") or "(not documented)" for sample in members),
                     "Collection_Date": group["Collection_Date"].iloc[0],
                     "Location_Name": group["Location_Name"].iloc[0]})
    units = pd.DataFrame(rows, columns=["Unit_ID", "Unit_Type", "Sample_IDs", "Library_IDs", "Collection_Date",
                                        "Location_Name"])
    site_date = units["Location_Name"].astype(str).str.casefold() + "|" + units["Collection_Date"].astype(str)
    units["Same_Site_Date_Units"] = site_date.map(site_date.value_counts()).astype(int) - 1
    return units, unit_of, verified


def known_observations(matrix, evidence):
    """Sample-level known values: one row per mutation and sample with a detection or a counted value."""

    columns = ["Variant_Key", "Sample_ID", "State", "VAF", "Alt_Read_Count", "Total_Depth", "Count_Basis"]
    observed = matrix.loc[(matrix["Status"] == VARIANT_OBSERVED).to_numpy(), ["Variant_Key", "Sample_ID"]]
    calls = evidence[evidence["Counted_Once"] & evidence["Analysis_Included"]][[
        "Sample_ID", "Variant_Key", "VAF", "Alt_Read_Count", "Total_Depth", "Alt_Count_Source", "Interval_Status"]]
    detected = observed.merge(calls, on=["Sample_ID", "Variant_Key"], how="inner")
    usable_counts = detected["Interval_Status"].astype(str).str.startswith(("Calculated", "Approximate"))
    detected["Count_Basis"] = np.where(~usable_counts, FREQUENCY_ONLY,
                                       np.where(detected["Alt_Count_Source"].eq("Reconstructed"), RECONSTRUCTED,
                                                MEASURED))
    detected["Alt_Read_Count"] = detected["Alt_Read_Count"].astype("Float64").where(usable_counts)
    detected["Total_Depth"] = detected["Total_Depth"].astype("Float64").where(usable_counts)
    detected["State"] = DETECTED

    if "Callable_Evidence" in matrix:  # compared without a text copy of the (possibly huge) matrix column
        evidence_text = matrix["Callable_Evidence"]
        from_counts = matrix[evidence_text.isin([ZERO_ALT_READS, ALT_READS_NOT_CALLED]).to_numpy()]
    else:
        evidence_text = pd.Series(dtype=object)
        from_counts = matrix.iloc[0:0]
    counted = pd.DataFrame({
        "Variant_Key": from_counts["Variant_Key"].to_numpy(),
        "Sample_ID": from_counts["Sample_ID"].to_numpy(),
        "State": np.where((evidence_text.loc[from_counts.index] == ZERO_ALT_READS).to_numpy(), NON_DETECTION,
                          COUNTED_NOT_CALLED),
        "VAF": from_counts["Callable_Frequency"].to_numpy(dtype=float),
        "Alt_Read_Count": pd.array(from_counts["Alt_Read_Count"], dtype="Float64"),
        "Total_Depth": pd.array(pd.to_numeric(from_counts["Depth"], errors="coerce"), dtype="Float64"),
        "Count_Basis": MEASURED,
    })
    table = pd.concat([detected[columns], counted[columns]], ignore_index=True)
    table["VAF"] = table["VAF"].astype(float)
    return table


def combine_units(known, unit_of, rho_for_unit):
    """Combine the known values of each mutation within each unit (counts are summed)."""

    table = known.copy()
    table["Unit_ID"] = table["Sample_ID"].map(unit_of)
    table = table[table["Unit_ID"].notna()]
    alt = table["Alt_Read_Count"].astype("Float64").to_numpy(dtype=float, na_value=np.nan)
    depth = table["Total_Depth"].astype("Float64").to_numpy(dtype=float, na_value=np.nan)
    has_counts = np.isfinite(alt) & np.isfinite(depth) & (depth > 0)
    table = table.assign(
        _x=np.where(has_counts, alt, 0.0), _n=np.where(has_counts, depth, 0.0),
        _q=np.where(has_counts, depth ** 2, 0.0), _counted=has_counts.astype(int),
        _reconstructed=(table["Count_Basis"].eq(RECONSTRUCTED) & has_counts).astype(int),
        _detected=table["State"].eq(DETECTED).astype(int), _vaf=table["VAF"].astype(float), _rows=1)
    grouped = table.groupby(["Variant_Key", "Unit_ID"], sort=False).agg(
        x=("_x", "sum"), n=("_n", "sum"), q=("_q", "sum"), counted=("_counted", "sum"),
        reconstructed=("_reconstructed", "sum"), detected=("_detected", "sum"), mean_vaf=("_vaf", "mean"),
        libraries=("_rows", "sum")).reset_index()
    with_counts = grouped["counted"].gt(0).to_numpy()
    n = grouped["n"].to_numpy(dtype=float)
    x = grouped["x"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        vaf = np.where(with_counts, x / n, grouped["mean_vaf"].to_numpy(dtype=float))
    basis = np.where(~with_counts, FREQUENCY_ONLY, np.where(grouped["reconstructed"].gt(0), RECONSTRUCTED, MEASURED))
    state = np.where(grouped["detected"].gt(0), DETECTED,
                     np.where(with_counts & (x == 0), NON_DETECTION, COUNTED_NOT_CALLED))
    q = grouped["q"].to_numpy(dtype=float)
    binomial_factor = variance_factor(n, q, 0.0)
    lower, upper = effective_score_interval(np.where(with_counts, vaf, np.nan), 1 / binomial_factor)
    rho = grouped["Unit_ID"].map(rho_for_unit).astype(float).to_numpy()
    calibrated = with_counts & np.isfinite(rho)
    bb_lower, bb_upper = effective_score_interval(np.where(calibrated, vaf, np.nan),
                                                  1 / variance_factor(n, q, np.where(calibrated, rho, 0.0)))
    return pd.DataFrame({
        "Variant_Key": grouped["Variant_Key"], "Unit_ID": grouped["Unit_ID"], "Libraries_Known": grouped["libraries"],
        "State": state, "VAF": vaf,
        "Alt_Read_Count": pd.array(np.where(with_counts, x, np.nan), dtype="Float64").round().astype("Int64"),
        "Total_Depth": pd.array(np.where(with_counts, n, np.nan), dtype="Float64").round().astype("Int64"),
        "Depth_Square_Sum": np.where(with_counts, q, np.nan), "Count_Basis": basis,
        "Wilson_Lower": lower, "Wilson_Upper": upper,
        "Rho": np.where(calibrated, rho, np.nan),
        "Rho_Source": np.where(calibrated, "Calibrated from technical replicates (upper 95% bound)", ""),
        "BetaBinomial_Lower": bb_lower, "BetaBinomial_Upper": bb_upper,
    })


def mutation_eligibility(matrix, evidence, profiles, unit_of, minimum_shared):
    """One row per mutation with its state counts and whether it can enter each grouping method."""

    # The matrix can have ~300,000 rows (mmc2): count states with integer codes and
    # bincount instead of string copies and a wide frame (they raised peak memory ~50 MB).
    codes, keys = pd.factorize(matrix["Variant_Key"], sort=False)
    # Both keep first-appearance order, so row i of `variants` is code i.
    variants = matrix.drop_duplicates("Variant_Key")[["Variant_Key", "Reference_ID", "Position", "Ref", "Alt"]]
    variants = variants.reset_index(drop=True)
    if not np.array_equal(variants["Variant_Key"].to_numpy(), np.asarray(keys)):
        raise ValueError("Internal error: mutation order does not match its codes.")
    status = matrix["Status"]
    callable_evidence = matrix["Callable_Evidence"] if "Callable_Evidence" in matrix else None

    def count(mask):
        return np.bincount(codes[np.asarray(mask, dtype=bool)], minlength=len(keys))

    def evidence_is(value):
        return (callable_evidence == value).to_numpy() if callable_evidence is not None else np.zeros(len(codes), bool)

    covered = (status == COVERED_NOT_DETECTED).to_numpy()
    table = variants.assign(
        Samples_On_Reference=np.bincount(codes, minlength=len(keys)),
        Detected=count((status == VARIANT_OBSERVED).to_numpy()),
        Supported_Non_Detection=count(evidence_is(ZERO_ALT_READS)),
        Counted_Not_Called=count(evidence_is(ALT_READS_NOT_CALLED)),
        Unknown_Coverage=count((status == COVERAGE_UNKNOWN).to_numpy()),
        Low_Or_No_Coverage=count(status.isin([LOW_COVERAGE, NOT_COVERED]).to_numpy()),
        Covered_Value_Unknown=count(covered & evidence_is(VALUE_UNKNOWN)),
        Rejected_By_QC=count((status == EXCLUDED_BY_QC).to_numpy()),
    )
    first = evidence.drop_duplicates("Variant_Key").set_index("Variant_Key")
    table["Segment"] = table["Reference_ID"]
    table["Annotated_Region"] = table["Variant_Key"].map(first["Annotated_Region"]).fillna("")
    table["Reference_Status"] = table["Variant_Key"].map(first["Reference_Status"]).fillna("Not evaluated")
    per_variant = profiles.groupby("Variant_Key")
    table["Known_Units"] = table["Variant_Key"].map(per_variant.size()).fillna(0).astype(int)
    basis = profiles.assign(Measured=profiles["Count_Basis"].eq(MEASURED).astype(int),
                            Reconstructed=profiles["Count_Basis"].eq(RECONSTRUCTED).astype(int),
                            Frequency=profiles["Count_Basis"].eq(FREQUENCY_ONLY).astype(int)).groupby("Variant_Key")
    table["Measured_Units"] = table["Variant_Key"].map(basis["Measured"].sum()).fillna(0).astype(int)
    table["Reconstructed_Units"] = table["Variant_Key"].map(basis["Reconstructed"].sum()).fillna(0).astype(int)
    table["Frequency_Only_Units"] = table["Variant_Key"].map(basis["Frequency"].sum()).fillna(0).astype(int)
    table["Known_Units_With_Counts"] = table["Measured_Units"] + table["Reconstructed_Units"]
    mismatch = table["Reference_Status"].astype(str).str.contains("mismatch", case=False)
    reasons, ordinary, depth_aware = [], [], []
    for row, bad_reference in zip(table.to_dict("records"), mismatch):
        reason = []
        if bad_reference:
            reason.append("reference has Ref-base mismatches")
        if row["Known_Units"] < minimum_shared:
            reason.append(f"known in {row['Known_Units']} observation unit(s); at least {minimum_shared} needed")
        ok_ordinary = not reason
        if ok_ordinary and row["Known_Units_With_Counts"] < minimum_shared:
            depth_reason = (f"read counts in {row['Known_Units_With_Counts']} unit(s); depth-aware grouping needs "
                            f"at least {minimum_shared}")
        else:
            depth_reason = ""
        ordinary.append(ok_ordinary)
        depth_aware.append(ok_ordinary and not depth_reason)
        reasons.append("; ".join(reason) if reason else depth_reason)
    table["Eligible_Ordinary"] = ordinary
    table["Eligible_Depth_Aware"] = depth_aware
    table["Exclusion_Reason"] = reasons
    return table[ELIGIBILITY_COLUMNS].reset_index(drop=True)


def build_mutation_profiles(matrix, evidence, calibration=None, samples=None, minimum_shared=4):
    """Return {"profiles", "eligibility", "units", "verified_pairs"} for the pooled analysis.

    `calibration` is the Phase 2 calibration table; its usable rho (upper 95%
    bound) is attached only to the units of the replicated specimens.
    """

    units, unit_of, verified = observation_units(evidence, samples)
    rho_for_unit = {}
    if calibration is not None and len(calibration):
        usable = calibration[calibration["Status"].isin(USABLE_STATUSES)]
        for row in usable.to_dict("records"):
            unit = unit_of.get(row["Replicate_1"])
            if unit is not None and unit == unit_of.get(row["Replicate_2"]):
                rho_for_unit[unit] = float(row["Rho_Used_For_Intervals"])
    known = known_observations(matrix, evidence)
    combined = combine_units(known, unit_of, rho_for_unit)
    variants = matrix.drop_duplicates("Variant_Key")[["Variant_Key", "Reference_ID", "Position", "Ref", "Alt"]]
    profiles = combined.merge(variants, on="Variant_Key", how="left").merge(
        units.drop(columns=[]), on="Unit_ID", how="left")
    profiles["Segment"] = profiles["Reference_ID"]
    status = evidence.drop_duplicates("Variant_Key").set_index("Variant_Key")["Reference_Status"]
    profiles["Reference_Status"] = profiles["Variant_Key"].map(status).fillna("Not evaluated")
    profiles = profiles[PROFILE_COLUMNS].sort_values(["Reference_ID", "Position", "Alt", "Collection_Date", "Unit_ID"],
                                                     kind="mergesort").reset_index(drop=True)
    eligibility = mutation_eligibility(matrix, evidence, profiles, unit_of, minimum_shared)
    return {"profiles": profiles, "eligibility": eligibility, "units": units, "verified_pairs": verified,
            "rho_for_unit": rho_for_unit}
