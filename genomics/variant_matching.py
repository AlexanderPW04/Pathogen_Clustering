"""Compare reference-based variants between individual and pooled evidence.

Each variant key is placed in exactly one category:

- Shared: observed in at least one individual sequence and one pooled sample.
- Individual only: carried by individual sequences and not reported in any
  pooled sample that had adequate depth at that position.
- Pooled only: reported in pooled samples while callable individual
  sequences show the reference base.
- Unresolved: the other source has no adequate coverage information at the
  position, so "not detected" cannot be distinguished from "not covered".
- Incomparable: the pooled record uses a reference that is not verified as
  the same sequence, or its reference base does not match the reference.

"Not detected" is never treated as confirmed absence. Dates and locations
are context only; they do not show which source detected a change first or
that any transmission occurred.
"""

import numpy as np
import pandas as pd

from add_temporal_geographic import haversine_km
from genomics.pooled_observations import pooled_position_status, read_any_table, standard_name
from genomics.references import reference_base, references_compatible
from genomics.sequence_variants import CALL_REFERENCE, CALL_VARIANT
from genomics.unified_schema import COVERED_NOT_DETECTED, variant_key


SHARED = "Shared: observed in both sources"
INDIVIDUAL_ONLY = "Individual only: not detected in adequately covered pooled samples"
POOLED_ONLY = "Pooled only: reference base called in individual sequences"
UNRESOLVED_POOLED = "Unresolved: pooled coverage missing or insufficient"
UNRESOLVED_INDIVIDUAL = "Unresolved: individual sequences not callable here"
INCOMPARABLE = "Incomparable: reference not verified as the same"

CATEGORY_ORDER = [SHARED, INDIVIDUAL_ONLY, POOLED_ONLY, UNRESOLVED_POOLED,
                  UNRESOLVED_INDIVIDUAL, INCOMPARABLE]


def harmonize_pooled_variants(pooled_variants, reference, registry):
    """Split pooled records into comparable rows and incomparable rows.

    Comparable rows are relabeled with the comparison reference ID only when
    the two references were shown to be the same sequence.
    """

    comparable_rows = []
    incomparable_rows = []
    for row in pooled_variants.to_dict("records"):
        ok, reason = references_compatible(row["Reference_ID"], reference["Reference_ID"], registry)
        position = int(row["Position"])
        if ok and position > reference["Length"]:
            ok, reason = False, "Position is beyond the end of the reference."
        if ok and reference_base(reference, position) != row["Ref"]:
            ok = False
            reason = (
                f"Reference base mismatch: the record says {row['Ref']}, the reference has "
                f"{reference_base(reference, position)} at position {position}."
            )
        if ok:
            row = dict(row)
            row["Original_Reference_ID"] = row["Reference_ID"]
            row["Reference_ID"] = reference["Reference_ID"]
            row["Variant_Key"] = variant_key(reference["Reference_ID"], position, row["Ref"], row["Alt"])
            comparable_rows.append(row)
        else:
            row = dict(row)
            row["Incomparable_Reason"] = reason
            incomparable_rows.append(row)
    return pd.DataFrame(comparable_rows), pd.DataFrame(incomparable_rows)


def usable_coordinates(sample):
    """Return (lat, lon) only for samples whose location precision allows distances."""

    precision = str(sample.get("Geographic_Precision", ""))
    if "too coarse" in precision or "no coordinates" in precision:
        return None
    if pd.isna(sample.get("Latitude")) or pd.isna(sample.get("Longitude")):
        return None
    return float(sample["Latitude"]), float(sample["Longitude"])


def nearest_distance(first_samples, second_samples):
    """Smallest Haversine distance between two groups, using precise points only."""

    best = None
    for first in first_samples:
        point_1 = usable_coordinates(first)
        if point_1 is None:
            continue
        for second in second_samples:
            point_2 = usable_coordinates(second)
            if point_2 is None:
                continue
            distance = haversine_km(point_1[0], point_1[1], point_2[0], point_2[1])
            if best is None or distance < best:
                best = distance
    return best


def join_names(values):
    names = sorted({str(value) for value in values if pd.notna(value) and str(value)})
    return "; ".join(names)


def compare_evidence(individual, pooled, reference, registry, minimum_depth, date_window_days):
    """Build the cross-source comparison table.

    `individual` holds the outputs of build_individual_profiles(); `pooled`
    holds QC-checked pooled variants, samples, and the depth lookup.
    """

    individual_samples = individual["samples"]
    mapped = individual_samples[
        individual_samples["Variant_Profile_Status"].eq("Mapped to reference")
    ]
    individual_by_id = {row["Sample_ID"]: row for row in mapped.to_dict("records")}
    calls = individual["calls"]

    pooled_samples = pooled["samples"][pooled["samples"]["Analysis_Included"]]
    pooled_by_id = {row["Sample_ID"]: row for row in pooled_samples.to_dict("records")}
    # Reference problems are listed for every pooled record; only records that
    # passed quality control are compared.
    comparable, incomparable = harmonize_pooled_variants(pooled["variants"], reference, registry)
    excluded_by_key = {}
    if not comparable.empty:
        kept = comparable["Analysis_Included"].astype(bool)
        for key, sample_id in zip(comparable.loc[~kept, "Variant_Key"], comparable.loc[~kept, "Sample_ID"]):
            excluded_by_key.setdefault(key, set()).add(sample_id)
        comparable = comparable[kept]
    comparable_samples = set(comparable["Sample_ID"]) if not comparable.empty else set()
    # Pooled samples whose records are on the comparison reference.
    pooled_on_reference = [
        sample_id for sample_id in pooled_by_id
        if sample_id in comparable_samples
        or reference["Reference_ID"] in str(pooled_by_id[sample_id]["Reference_ID"]).split("; ")
    ]

    individual_variants = individual["variants"]
    keys = set(individual_variants["Variant_Key"])
    carriers_by_key = {}
    for key, group in individual_variants.groupby("Variant_Key"):
        carriers_by_key[key] = list(dict.fromkeys(group["Sample_ID"]))

    reported_by_key = {}
    original_reference_by_sample = {}
    if not comparable.empty:
        keys |= set(comparable["Variant_Key"])
        for key, group in comparable.groupby("Variant_Key"):
            reported_by_key[key] = group.drop_duplicates("Sample_ID")
        for sample_id, original in zip(comparable["Sample_ID"], comparable["Original_Reference_ID"]):
            original_reference_by_sample.setdefault(sample_id, original)

    rows = []
    for key in sorted(keys, key=lambda text: (int(text.split(":")[1]), text)):
        _, position_text, change = key.split(":")
        position = int(position_text)
        ref, alt = change.split(">")

        carrier_ids = carriers_by_key.get(key, [])
        reference_calls = 0
        other_allele = 0
        not_callable = 0
        for sample_id, call_text in calls.items():
            if sample_id in carrier_ids:
                continue
            call = call_text[position - 1]
            if call == CALL_REFERENCE:
                reference_calls += 1
            elif call == CALL_VARIANT:
                other_allele += 1
            else:
                not_callable += 1

        reported = reported_by_key.get(key)
        reported_ids = list(reported["Sample_ID"]) if reported is not None else []
        adequate = 0
        low_or_none = 0
        unknown = 0
        excluded_samples = excluded_by_key.get(key, set())
        for sample_id in pooled_on_reference:
            if sample_id in reported_ids:
                continue
            if sample_id in excluded_samples:
                unknown += 1  # a QC-excluded call is not evidence either way
                continue
            # Depth tables use the pooled sample's own reference label.
            original_reference = original_reference_by_sample.get(sample_id, reference["Reference_ID"])
            status, _ = pooled_position_status(
                pooled["depth_lookup"], sample_id, original_reference, position, minimum_depth
            )
            if status == COVERED_NOT_DETECTED:
                adequate += 1
            elif status == "Coverage unknown":
                unknown += 1
            else:
                low_or_none += 1

        if carrier_ids and reported_ids:
            category = SHARED
        elif carrier_ids:
            category = INDIVIDUAL_ONLY if adequate else UNRESOLVED_POOLED
        else:
            category = POOLED_ONLY if reference_calls else UNRESOLVED_INDIVIDUAL

        carrier_samples = [individual_by_id[sample_id] for sample_id in carrier_ids
                           if sample_id in individual_by_id]
        pooled_carriers = [pooled_by_id[sample_id] for sample_id in reported_ids
                           if sample_id in pooled_by_id]
        individual_dates = [row["Collection_Date"] for row in carrier_samples
                            if pd.notna(row["Collection_Date"])]
        pooled_dates = [row["Collection_Date"] for row in pooled_carriers
                        if pd.notna(row["Collection_Date"])]
        first_individual = min(individual_dates) if individual_dates else pd.NaT
        first_pooled = min(pooled_dates) if pooled_dates else pd.NaT
        pairs_in_window = 0
        for date_1 in individual_dates:
            for date_2 in pooled_dates:
                if abs((date_2 - date_1).days) <= date_window_days:
                    pairs_in_window += 1

        frequencies = reported["Alt_Frequency"] if reported is not None else pd.Series(dtype=float)
        distance = nearest_distance(carrier_samples, pooled_carriers) if category == SHARED else None
        rows.append({
            "Variant_Key": key,
            "Reference_ID": reference["Reference_ID"],
            "Position": position,
            "Ref": ref,
            "Alt": alt,
            "Category": category,
            "Individual_Carriers": len(carrier_ids),
            "Individual_Reference_Calls": reference_calls,
            "Individual_Other_Allele": other_allele,
            "Individual_Not_Callable": not_callable,
            "Pooled_Reported": len(reported_ids),
            "Pooled_Adequate_Depth_Not_Reported": adequate,
            "Pooled_Low_Or_No_Coverage": low_or_none,
            "Pooled_Coverage_Unknown_Or_QC_Excluded": unknown,
            "Pooled_Frequency_Median": float(frequencies.median()) if len(frequencies) else np.nan,
            "Pooled_Frequency_Max": float(frequencies.max()) if len(frequencies) else np.nan,
            "First_Individual_Date": first_individual,
            "First_Pooled_Date": first_pooled,
            "Days_Pooled_Minus_Individual": (
                (first_pooled - first_individual).days
                if pd.notna(first_individual) and pd.notna(first_pooled) else np.nan
            ),
            "Date_Pairs_Within_Window": pairs_in_window,
            "Individual_Carrier_Locations": join_names(row["Location_Name"] for row in carrier_samples),
            "Pooled_Reporting_Locations": join_names(row["Location_Name"] for row in pooled_carriers),
            "Nearest_Km_Between_Carriers": round(distance, 3) if distance is not None else np.nan,
            "Note": "",
        })

    if not incomparable.empty:
        for (sample_reference, position, ref, alt), group in incomparable.groupby(
            ["Reference_ID", "Position", "Ref", "Alt"], sort=True
        ):
            rows.append({
                "Variant_Key": variant_key(sample_reference, position, ref, alt),
                "Reference_ID": sample_reference,
                "Position": int(position),
                "Ref": ref,
                "Alt": alt,
                "Category": INCOMPARABLE,
                "Pooled_Reported": group["Sample_ID"].nunique(),
                "Note": group["Incomparable_Reason"].iloc[0],
            })

    table = pd.DataFrame(rows)
    if table.empty:
        return table
    for column in ["Individual_Carriers", "Individual_Reference_Calls", "Individual_Other_Allele",
                   "Individual_Not_Callable", "Pooled_Reported", "Pooled_Adequate_Depth_Not_Reported",
                   "Pooled_Low_Or_No_Coverage", "Pooled_Coverage_Unknown_Or_QC_Excluded",
                   "Date_Pairs_Within_Window"]:
        table[column] = pd.array(table[column], dtype="Int64")
    return table


def category_counts(comparison):
    """Count variant keys per category, in a fixed order."""

    counts = comparison["Category"].value_counts() if not comparison.empty else pd.Series(dtype=int)
    return pd.DataFrame({
        "Category": CATEGORY_ORDER,
        "Variants": [int(counts.get(category, 0)) for category in CATEGORY_ORDER],
    })


def load_signatures(path):
    """Read documented variant profiles: Signature_Name, Reference_ID, Position, Ref, Alt."""

    data = read_any_table(path)
    data.columns = [standard_name(name, {"signature_name": "Signature_Name",
                                         "signature": "Signature_Name",
                                         "reference_id": "Reference_ID",
                                         "position": "Position", "ref": "Ref", "alt": "Alt",
                                         "source": "Source"})
                    for name in data.columns]
    required = ["Signature_Name", "Reference_ID", "Position", "Ref", "Alt", "Source"]
    missing = [name for name in required if name not in data.columns]
    if missing:
        raise ValueError("Variant signature table is missing column(s): " + ", ".join(missing))
    if data["Source"].astype(str).str.strip().eq("").any():
        raise ValueError("Every signature row needs a documented Source.")
    positions = pd.to_numeric(data["Position"], errors="coerce")
    if positions.isna().any() or positions.lt(1).any():
        raise ValueError("Signature positions must be whole numbers of 1 or more.")
    data["Position"] = positions.astype(int)
    for column in ["Ref", "Alt"]:
        data[column] = data[column].astype(str).str.strip().str.upper()
    return data


def compare_signatures(signatures, reference, registry, individual, pooled, minimum_depth):
    """Compare each sample with each documented profile (not a lineage call)."""

    rows = []
    calls = individual["calls"] if individual else {}
    individual_variants = individual["variants"] if individual else pd.DataFrame(columns=["Variant_Key", "Sample_ID"])
    pooled_variants = pooled["variants"] if pooled else None

    for name, group in signatures.groupby("Signature_Name", sort=True):
        keys = []
        skipped = 0
        for row in group.to_dict("records"):
            ok, _ = references_compatible(row["Reference_ID"], reference["Reference_ID"], registry)
            if ok and row["Position"] <= reference["Length"] and reference_base(reference, row["Position"]) == row["Ref"]:
                keys.append((row["Position"], variant_key(reference["Reference_ID"], row["Position"], row["Ref"], row["Alt"])))
            else:
                skipped += 1

        for sample_id, call_text in calls.items():
            carried = set(individual_variants.loc[individual_variants["Sample_ID"].eq(sample_id), "Variant_Key"])
            assessable = sum(1 for position, _ in keys if call_text[position - 1] in (CALL_REFERENCE, CALL_VARIANT))
            present = sum(1 for _, key in keys if key in carried)
            rows.append(signature_row(name, sample_id, "Individual sequence", len(keys), skipped, present, assessable))

        if pooled_variants is not None:
            included = pooled_variants[pooled_variants["Analysis_Included"]]
            comparable, _ = harmonize_pooled_variants(included, reference, registry)
            for sample in pooled["samples"][pooled["samples"]["Analysis_Included"]].to_dict("records"):
                sample_id = sample["Sample_ID"]
                reported = set()
                if not comparable.empty:
                    reported = set(comparable.loc[comparable["Sample_ID"].eq(sample_id), "Variant_Key"])
                present = sum(1 for _, key in keys if key in reported)
                assessable = present
                for position, key in keys:
                    if key in reported:
                        continue
                    status, _ = pooled_position_status(pooled["depth_lookup"], sample_id,
                                                       str(sample["Reference_ID"]), position, minimum_depth)
                    if status == COVERED_NOT_DETECTED:
                        assessable += 1
                rows.append(signature_row(name, sample_id, "Pooled sample", len(keys), skipped, present, assessable))

    return pd.DataFrame(rows, columns=[
        "Signature_Name", "Sample_ID", "Sample_Type", "Signature_Variants", "Signature_Variants_Not_Comparable",
        "Variants_Observed", "Variants_Assessable", "Fraction_Of_Assessable_Observed", "Interpretation",
    ])


def signature_row(name, sample_id, sample_type, total, skipped, present, assessable):
    return {
        "Signature_Name": name,
        "Sample_ID": sample_id,
        "Sample_Type": sample_type,
        "Signature_Variants": total,
        "Signature_Variants_Not_Comparable": skipped,
        "Variants_Observed": present,
        "Variants_Assessable": assessable,
        "Fraction_Of_Assessable_Observed": present / assessable if assessable else np.nan,
        "Interpretation": "Profile overlap only; not a lineage assignment",
    }
