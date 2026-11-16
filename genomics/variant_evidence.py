"""Sprint 5 variant data model: one validated record per pooled variant observation.

It is built from the Sprint 3 pooled variant table after quality control, so
nothing is re-imported or recalculated twice. Each record keeps the source
values and says where every derived value came from.

Read counts can come from four places, in this order of precedence:

1. "Supplied": an alternate read count column in the file (for example iVar
   ALT_DP). A supplied reference count (REF_DP) is kept as well.
2. "Base counts": a per-position A/C/G/T table for the same sample and
   position. The reference and alternative counts are the counts of the Ref
   and Alt bases; the denominator is that table's depth.
3. "Reconstructed": reported frequency x total depth, when no count exists.
   It is used only when every frequency-only record in the table maps to one
   unique whole number of reads that reproduces the reported frequency at the
   table's stated precision (see `reconstruction_check`). Reconstructed counts
   are derived, not measured, and are labeled so.
4. "Frequency only": no count can be established; no interval is made.

A reference read count is never derived as depth - alternative reads: the
depth can include other alleles, deletions, or low-quality bases.

The observed VAF is alternative reads / total usable depth for count-based
records (the denominator of the source frequency, which includes any other
alleles at the position). The reported frequency is always kept unchanged.
"""

import math

import numpy as np
import pandas as pd

from genomics.frequency_uncertainty import wilson_interval
from genomics.pooled_observations import COUNT_COLUMNS
from genomics.unified_schema import FREQUENCY_TOLERANCE, NUCLEOTIDES


SUPPLIED = "Supplied"
BASE_COUNTS = "Base counts"
RECONSTRUCTED = "Reconstructed"
FREQUENCY_ONLY = "Frequency only"
NO_DEPTH = "No depth"
COUNT_SOURCES = [SUPPLIED, BASE_COUNTS, RECONSTRUCTED, FREQUENCY_ONLY, NO_DEPTH]

COUNT_SOURCE_TEXT = {
    SUPPLIED: "Alternate read count supplied in the file",
    BASE_COUNTS: "Counts of the Ref and Alt bases in the per-position A/C/G/T table",
    RECONSTRUCTED: ("Reconstructed as reported frequency x total depth (derived, not measured; "
                    "unique at the table's stated precision)"),
    FREQUENCY_ONLY: "No read counts: only the reported frequency is available",
    NO_DEPTH: "No read depth was supplied, so no count can be established",
}

INTERVAL_SUPPLIED = "Calculated (supplied counts)"
INTERVAL_BASE_COUNTS = "Calculated (base counts)"
INTERVAL_RECONSTRUCTED = "Approximate (reconstructed counts)"
NO_INTERVAL_FREQUENCY_ONLY = "Not calculated: no read counts"
NO_INTERVAL_DEPTH = "Not calculated: no read depth"
NO_INTERVAL_DISAGREE = "Not calculated: reported frequency and counts disagree"
NO_INTERVAL_CONFLICT = "Not calculated: conflicting repeated records"
NO_INTERVAL_REPEAT = "Not calculated: repeated record (counted once)"
NO_INTERVAL_SITE_SUM = "Not calculated: allele frequencies at this site sum above 1"

EVIDENCE_COLUMNS = [
    "Sample_ID", "Collection_Date", "Collection_Date_Original", "Collection_Date_Status",
    "Location_Name", "Organism", "Reference_ID", "Reference_Status", "Annotated_Region",
    "Position", "Ref", "Alt", "Variant_Key",
    "Ref_Read_Count", "Ref_Count_Source", "Alt_Read_Count", "Alt_Count_Source",
    "Total_Depth", "Depth_Source", "Reported_Frequency", "Count_Derived_Frequency",
    "Frequency_Difference", "Frequency_Check", "VAF", "VAF_Source",
    "Wilson_Lower", "Wilson_Upper", "Interval_Width", "Interval_Status",
    "Coverage_Status", "QC_Status", "Analysis_Included", "Duplicate_Status", "Counted_Once",
    "Alleles_At_Site", "Validation_Warnings", "Data_Source",
]


def round_significant(value, digits):
    if value == 0:
        return 0.0
    return round(value, digits - 1 - int(math.floor(math.log10(abs(value)))))


def reproduces(alt, depth, frequency, digits):
    """True when alt / depth, written with `digits` significant digits, equals the reported frequency."""

    if alt < 1 or alt > depth:
        return False
    return abs(round_significant(alt / depth, digits) - frequency) <= 1e-12 * max(1.0, frequency)


def reconstruction_check(frequencies, depths, digits):
    """Decide whether alternate counts can be reconstructed for a whole table.

    `digits` are the significant digits written in each reported frequency.
    The table's precision is the largest number of digits written: tools such
    as iVar print 6 significant digits and drop trailing zeros ('0.14' means
    0.140000). For each record the candidate count is round(f x depth). It is
    accepted when it reproduces f at that precision and neither neighbouring
    count does (so it is unique). If any record cannot be reproduced at all,
    the frequency and depth do not share a denominator and nothing is
    reconstructed.

    Returns (counts as a float array with NaN where not reconstructed, details).
    """

    frequencies = np.asarray(frequencies, dtype=float)
    depths = np.asarray(depths, dtype=float)
    digit_values = pd.array(digits, dtype="Int64")
    counts = np.full(len(frequencies), np.nan)
    details = {"eligible": len(frequencies), "precision": None, "unique": 0, "ambiguous": 0,
               "not_reproduced": 0, "used": False}
    known = [int(value) for value in digit_values if pd.notna(value)]
    if not len(frequencies) or not known:
        return counts, details
    precision = max(known)
    details["precision"] = precision
    candidates = np.full(len(frequencies), np.nan)
    for index, (frequency, depth) in enumerate(zip(frequencies, depths)):
        alt = int(round(frequency * depth))
        if not reproduces(alt, int(depth), frequency, precision):
            details["not_reproduced"] += 1
            continue
        if reproduces(alt - 1, int(depth), frequency, precision) or reproduces(alt + 1, int(depth), frequency, precision):
            details["ambiguous"] += 1
            continue
        candidates[index] = alt
        details["unique"] += 1
    if details["not_reproduced"] == 0:
        details["used"] = True
        counts = candidates
    return counts, details


def base_count_lookup(variants, position_table):
    """Ref and Alt base counts and depth for each variant record, from the A/C/G/T table."""

    empty = pd.DataFrame(index=variants.index, columns=["Base_Ref", "Base_Alt", "Base_Depth"], dtype=float)
    if position_table is None or not position_table.attrs.get("has_base_counts", False):
        return empty
    keys = ["Sample_ID", "Reference_ID", "Position"]
    cells = variants[keys + ["Ref", "Alt"]].copy()
    cells["Position"] = cells["Position"].astype(int)
    merged = cells.reset_index().merge(position_table[keys + ["Depth"] + COUNT_COLUMNS], on=keys,
                                       how="left").set_index("index")
    column_for = dict(zip(NUCLEOTIDES, COUNT_COLUMNS))
    rows = merged.to_dict("records")
    return pd.DataFrame({
        "Base_Ref": [row[column_for[row["Ref"]]] for row in rows],
        "Base_Alt": [row[column_for[row["Alt"]]] for row in rows],
        "Base_Depth": [row["Depth"] for row in rows],
    }, index=merged.index, dtype=float)


def organism_label(reference, input_format):
    if reference is not None and reference.get("Description"):
        return f"From the reference FASTA: {reference['Description']}"
    if input_format and "Fontenele" in input_format:
        return "SARS-CoV-2 (stated by the publication)"
    return "Not supplied"


def coverage_status(depth, minimum_depth):
    if pd.isna(depth):
        return "Reported call; depth not supplied"
    if depth < minimum_depth:
        return f"Reported call; depth {int(depth)} below {minimum_depth}"
    return f"Reported call; depth at least {minimum_depth}"


def build_variant_evidence(variants, position_table=None, reference_status=None, reference=None,
                           input_format="", minimum_depth=20, tolerance=FREQUENCY_TOLERANCE):
    """Return (evidence table, reconstruction details) for the pooled variant records.

    `variants` is the pooled variant table after quality control. Missing
    values stay missing: no depth, count, date, or frequency is invented.
    """

    data = variants.reset_index(drop=True)
    count = len(data)
    reported = data["Alt_Frequency"].astype(float).to_numpy()
    depth = pd.to_numeric(data["Read_Depth"], errors="coerce").astype(float).to_numpy()
    supplied = pd.to_numeric(data["Alt_Depth"], errors="coerce").astype(float).to_numpy()
    supplied_ref = (pd.to_numeric(data["Ref_Depth"], errors="coerce").astype(float).to_numpy()
                    if "Ref_Depth" in data else np.full(count, np.nan))
    bases = base_count_lookup(data, position_table)

    alt = np.full(count, np.nan)
    ref = np.full(count, np.nan)
    total = np.full(count, np.nan)
    alt_source = np.full(count, FREQUENCY_ONLY, dtype=object)
    ref_source = np.full(count, "Not supplied (not derived as depth minus alternate reads)", dtype=object)
    depth_source = np.full(count, "Total depth reported for the call", dtype=object)

    has_supplied = np.isfinite(supplied) & np.isfinite(depth)
    alt[has_supplied], total[has_supplied] = supplied[has_supplied], depth[has_supplied]
    alt_source[has_supplied] = SUPPLIED
    has_ref = np.isfinite(supplied_ref)
    ref[has_ref] = supplied_ref[has_ref]
    ref_source[has_ref] = "Supplied in the file"

    base_alt = bases["Base_Alt"].to_numpy()
    base_depth = bases["Base_Depth"].to_numpy()
    use_bases = ~has_supplied & np.isfinite(base_alt) & np.isfinite(base_depth) & (base_depth > 0)
    alt[use_bases], total[use_bases] = base_alt[use_bases], base_depth[use_bases]
    alt_source[use_bases] = BASE_COUNTS
    depth_source[use_bases] = "Depth in the per-position table"
    base_ref = bases["Base_Ref"].to_numpy()
    ref_from_bases = ~has_ref & np.isfinite(base_ref)
    ref[ref_from_bases] = base_ref[ref_from_bases]
    ref_source[ref_from_bases] = "Count of the Ref base in the per-position table"

    no_depth = ~np.isfinite(depth) & ~use_bases
    alt_source[no_depth & ~has_supplied] = NO_DEPTH
    depth_source[~np.isfinite(depth) & ~use_bases] = "Not supplied"

    eligible = (alt_source == FREQUENCY_ONLY) & np.isfinite(depth) & (depth >= 1) & (reported > 0)
    digits = (data["Frequency_Significant_Digits"] if "Frequency_Significant_Digits" in data
              else pd.Series([pd.NA] * count))
    rebuilt, details = reconstruction_check(reported[eligible], depth[eligible],
                                            list(pd.Series(digits).iloc[np.flatnonzero(eligible)]))
    eligible_index = np.flatnonzero(eligible)
    ok = np.isfinite(rebuilt)
    alt[eligible_index[ok]] = rebuilt[ok]
    total[eligible_index[ok]] = depth[eligible_index[ok]]
    alt_source[eligible_index[ok]] = RECONSTRUCTED

    with np.errstate(divide="ignore", invalid="ignore"):
        derived = np.where(np.isfinite(alt) & (total > 0), alt / total, np.nan)
    difference = reported - derived
    checks = []
    for source, gap in zip(alt_source, difference):
        if source == RECONSTRUCTED:
            checks.append("Consistent by construction (reconstructed)")
        elif not np.isfinite(gap):
            checks.append("No counts to compare")
        elif abs(gap) <= tolerance:
            checks.append("Agrees with counts")
        else:
            checks.append(f"Differs from counts by {gap:+.4f}")
    checks = np.array(checks, dtype=object)

    # Sites with several alternative alleles in one sample: each allele keeps its
    # own binomial (allele vs everything else), but frequencies must not sum above 1.
    counted = data["Counted_Once"].to_numpy(dtype=bool) if "Counted_Once" in data else np.ones(count, bool)
    site_keys = ["Sample_ID", "Reference_ID", "Position"]
    site = data[site_keys].copy()
    site["Reported"] = np.where(counted, reported, 0.0)
    site["Counted"] = counted.astype(int)
    alleles = site.groupby(site_keys)["Counted"].transform("sum").to_numpy()
    site_sum = site.groupby(site_keys)["Reported"].transform("sum").to_numpy()
    over_one = site_sum > 1 + tolerance

    duplicate = data["Duplicate_Status"].astype(str).to_numpy() if "Duplicate_Status" in data else np.full(count, "")
    conflicting = np.array(["conflicting" in text for text in duplicate])
    status = np.empty(count, dtype=object)
    for index in range(count):
        source = alt_source[index]
        if not counted[index]:
            status[index] = NO_INTERVAL_REPEAT
        elif conflicting[index]:
            status[index] = NO_INTERVAL_CONFLICT
        elif source == NO_DEPTH:
            status[index] = NO_INTERVAL_DEPTH
        elif source == FREQUENCY_ONLY:
            status[index] = NO_INTERVAL_FREQUENCY_ONLY
        elif checks[index].startswith("Differs"):
            status[index] = NO_INTERVAL_DISAGREE
        elif over_one[index]:
            status[index] = NO_INTERVAL_SITE_SUM
        else:
            status[index] = {SUPPLIED: INTERVAL_SUPPLIED, BASE_COUNTS: INTERVAL_BASE_COUNTS,
                             RECONSTRUCTED: INTERVAL_RECONSTRUCTED}[source]
    calculable = np.array([text.startswith(("Calculated", "Approximate")) for text in status])
    lower, upper = wilson_interval(np.where(calculable, alt, np.nan), np.where(calculable, total, np.nan))

    vaf = np.where(calculable, derived, reported)
    # Without usable counts the VAF shown is the reported value; Interval_Status says why.
    vaf_source = np.where(calculable, "Alternative reads / total depth", "Reported frequency")

    references = {}
    if reference_status is not None and not reference_status.empty:
        references = dict(zip(reference_status["Reference_ID"], reference_status["Status"]))
    region = [""] * count
    for column in ["Published_Annotation", "Source_Gene_Feature"]:
        if column in data:
            region = [old or ("" if pd.isna(new) else str(new)) for old, new in zip(region, data[column])]

    warnings = []
    for index, row in enumerate(data[["Reference_ID", "Analysis_Included"]].itertuples(index=False)):
        items = []
        if references.get(row.Reference_ID, "") != "Verified compatible":
            items.append("reference not verified")
        if not row.Analysis_Included:
            items.append("excluded by QC")
        if checks[index].startswith("Differs"):
            items.append("frequency and counts disagree")
        if "Repeated" in duplicate[index]:
            items.append("repeated record")
        if over_one[index]:
            items.append("site frequencies sum above 1")
        if alt_source[index] == RECONSTRUCTED:
            items.append("counts reconstructed")
        if "Collection_Date_Status" in data and data["Collection_Date_Status"].iloc[index] not in ("Parsed", "Parsed with the selected day/month order"):
            items.append("date not usable")
        warnings.append("; ".join(items))

    evidence = pd.DataFrame({
        "Sample_ID": data["Sample_ID"],
        "Collection_Date": data["Collection_Date"],
        "Collection_Date_Original": data.get("Collection_Date_Original", pd.Series([pd.NA] * count)),
        "Collection_Date_Status": data.get("Collection_Date_Status", pd.Series([pd.NA] * count)),
        "Location_Name": data["Location_Name"],
        "Organism": organism_label(reference, input_format),
        "Reference_ID": data["Reference_ID"],
        "Reference_Status": [references.get(value, "Not evaluated") for value in data["Reference_ID"]],
        "Annotated_Region": region,
        "Position": data["Position"], "Ref": data["Ref"], "Alt": data["Alt"],
        "Variant_Key": data["Variant_Key"],
        "Ref_Read_Count": pd.array(np.where(np.isfinite(ref), ref, np.nan), dtype="Float64").round().astype("Int64"),
        "Ref_Count_Source": ref_source,
        "Alt_Read_Count": pd.array(np.where(np.isfinite(alt), alt, np.nan), dtype="Float64").round().astype("Int64"),
        "Alt_Count_Source": alt_source,
        "Total_Depth": pd.array(np.where(np.isfinite(total), total, depth), dtype="Float64").round().astype("Int64"),
        "Depth_Source": depth_source,
        "Reported_Frequency": reported,
        "Count_Derived_Frequency": derived,
        "Frequency_Difference": difference,
        "Frequency_Check": checks,
        "VAF": vaf,
        "VAF_Source": vaf_source,
        "Wilson_Lower": lower,
        "Wilson_Upper": upper,
        "Interval_Width": upper - lower,
        "Interval_Status": status,
        "Coverage_Status": [coverage_status(value, minimum_depth) for value in np.where(np.isfinite(total), total, depth)],
        "QC_Status": data["Quality_Status"],
        "Analysis_Included": data["Analysis_Included"].astype(bool),
        "Duplicate_Status": duplicate,
        "Counted_Once": counted,
        "Alleles_At_Site": alleles,
        "Validation_Warnings": warnings,
        "Data_Source": data["Data_Source"],
    }, columns=EVIDENCE_COLUMNS)
    return evidence, details


def validate_evidence(evidence):
    """Return a list of internal-consistency problems (empty when the table is valid)."""

    problems = []
    missing = [name for name in EVIDENCE_COLUMNS if name not in evidence.columns]
    if missing:
        return ["Evidence table is missing column(s): " + ", ".join(missing)]
    counts = evidence[evidence["Alt_Count_Source"].isin([SUPPLIED, BASE_COUNTS, RECONSTRUCTED])]
    alt = counts["Alt_Read_Count"].astype("Float64")
    depth = counts["Total_Depth"].astype("Float64")
    if (alt.isna() | depth.isna()).any():
        problems.append("A count-based record lacks an alternate count or depth.")
    if (alt < 0).any() or (alt > depth).any():
        problems.append("An alternate count is negative or larger than the depth.")
    if (depth < 1).any():
        problems.append("A count-based record has depth below 1.")
    no_counts = evidence[evidence["Alt_Count_Source"].isin([FREQUENCY_ONLY, NO_DEPTH])]
    if no_counts["Alt_Read_Count"].notna().any() or no_counts["Wilson_Lower"].notna().any():
        problems.append("A record without counts has a count or an interval.")
    with_interval = evidence[evidence["Wilson_Lower"].notna()]
    inside = ((with_interval["Wilson_Lower"] <= with_interval["VAF"] + 1e-12)
              & (with_interval["VAF"] <= with_interval["Wilson_Upper"] + 1e-12))
    if not inside.all():
        problems.append("An observed VAF lies outside its own interval.")
    if not evidence["Reported_Frequency"].between(0, 1).all():
        problems.append("A reported frequency is outside 0 to 1.")
    return problems
