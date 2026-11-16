"""Match pooled-sample metadata to variant records, and report every mismatch.

Sprint 5 Phase 4A. Samples are joined only by identical Sample_ID text
(surrounding spaces removed). Nothing is joined because dates, sites, or
similar-looking names agree: IDs that differ only by letter case, spaces,
punctuation, or leading zeros are listed as *possible* matches, never joined.

The report covers matched and unmatched samples, metadata rows without a
Sample_ID (for example a summary row such as "average:"), duplicate IDs,
values that conflict between the variant table and the metadata file, values
that conflict within one sample, specimens versus libraries, and segments.

A conflict is never settled by picking one value silently. The user chooses,
for all listed conflicts, one of:

- "metadata": the metadata file's value is used;
- "variant_table": the variant table's value is used;
- "unresolved": both values are left out (the field becomes unknown for
  those samples) and the conflict is recorded in the outputs.

Duplicate metadata IDs with different values, and different values within
one sample, can only be kept as "unresolved" (or fixed in the file).
"""

import re

import pandas as pd

from genomics.collection_dates import parse_collection_dates

FIELDS = ["Collection_Date", "Location_Name", "Specimen_ID", "Library_ID", "Latitude", "Longitude"]
RESOLUTIONS = {
    "metadata": "Use the metadata file's value",
    "variant_table": "Use the variant table's value",
    "unresolved": "Leave unresolved (the field becomes unknown for those samples)",
}


def clean(values):
    return values.fillna("").astype(str).str.strip()


def loose(text):
    """A comparison key for *hints only*: case, spaces, punctuation, and leading zeros ignored."""

    key = re.sub(r"[^0-9a-z]", "", str(text).casefold())
    return re.sub(r"(?<![0-9])0+(?=[0-9])", "", key)


def comparable(field, values, date_order):
    """Values as compared for conflicts: parsed dates, case-folded text, or numbers."""

    if field == "Collection_Date":
        parsed = parse_collection_dates(values, date_order)
        dates = parsed["Date"].dt.strftime("%Y-%m-%d")
        return dates.where(dates.notna(), values.str.casefold())
    if field in {"Latitude", "Longitude"}:
        numbers = pd.to_numeric(values, errors="coerce")
        return numbers.round(6).astype(str).where(numbers.notna(), values)
    return values.str.casefold()


def mapped(entry):
    """The entry's table with only its matched standard fields (text), or None."""

    if entry is None or entry.get("table") is None:
        return None
    columns = {match["field"]: match["column"] for match in entry.get("matches", []) if match["status"] == "matched"}
    if "Sample_ID" not in columns:
        return None
    table = pd.DataFrame({field: clean(entry["table"][column]) for field, column in columns.items()})
    return table


def sample_values(variants, field):
    """{sample: [distinct values]} of one field in the variant tables."""

    present = variants[variants[field].ne("")]
    return present.groupby("Sample_ID")[field].agg(lambda values: list(dict.fromkeys(values))).to_dict()


def match_metadata(review, options):
    """Fill review["matching"] (report rows) and review["metadata_conflicts"]; add problems."""

    entries = review["by_role"].get("snv", [])
    metadata_entry = review["by_role"].get("pooled_metadata", [None])[0]
    variant_tables = [table for table in (mapped(entry) for entry in entries) if table is not None]
    rows, conflicts = [], []
    resolution = options.get("conflict_resolution", "")
    date_order = options.get("date_order", "auto")

    def add(check, count, detail, status="OK"):
        rows.append({"Check": check, "Count": count, "Detail": detail, "Status": status})

    variant_ids = list(dict.fromkeys(str(value) for value in review.get("pooled_ids", [])))
    variants = pd.concat(variant_tables, ignore_index=True) if variant_tables else pd.DataFrame()
    for field in FIELDS:
        if field not in variants:
            variants[field] = ""

    # Different values of one field within one sample of the variant table.
    for field in FIELDS:
        for sample, values in sample_values(variants, field).items() if len(variants) else []:
            keys = set(comparable(field, pd.Series(values), date_order))
            if len(keys) > 1:
                conflicts.append({"Sample_ID": sample, "Field": field, "Kind": "Within the variant table",
                                  "Variant_Table_Value": " | ".join(values), "Metadata_Value": "",
                                  "Resolution": "unresolved" if resolution == "unresolved" else ""})

    metadata = mapped(metadata_entry)
    review["matching_summary"] = {"variant_samples": len(variant_ids)}
    if metadata_entry is None:
        add("Pooled-sample metadata file", 0, "None uploaded: dates and sites come from the variant table "
                                              "columns only (when present).", "Not supplied")
    elif metadata is None:
        add("Pooled-sample metadata file", 0, "No Sample_ID column was mapped, so nothing can be joined.",
            "Problem")
    else:
        blank = metadata["Sample_ID"].eq("")
        if blank.any():
            raw = metadata_entry["table"][blank]
            example = next((value for value in raw.astype(str).to_numpy().ravel() if value.strip()), "")
            add("Metadata rows without a Sample_ID", int(blank.sum()),
                f"Not used for joining (for example a summary row containing '{example[:30]}').", "Warning")
            review["warnings"].append(f"{metadata_entry['info']['name']}: {int(blank.sum())} row(s) have no "
                                      f"Sample_ID (for example '{example[:30]}') and are not used.")
        usable = metadata[~blank]
        ids = usable["Sample_ID"]
        repeated = ids[ids.duplicated(keep=False)].unique()
        identical, differing = [], []
        for sample in repeated:
            group = usable[ids.eq(sample)].drop(columns=["Sample_ID"])
            if len(group.drop_duplicates()) == 1:
                identical.append(sample)
                continue
            differing.append(sample)
            for field in group.columns:
                values = list(dict.fromkeys(group[field]))
                if len(set(comparable(field, pd.Series(values), date_order))) > 1:
                    conflicts.append({"Sample_ID": sample, "Field": field, "Kind": "Duplicate metadata rows",
                                      "Variant_Table_Value": "", "Metadata_Value": " | ".join(values),
                                      "Resolution": "unresolved" if resolution == "unresolved" else ""})
        if identical:
            add("Duplicate metadata IDs with identical rows", len(identical),
                "Kept once: " + ", ".join(identical[:8]), "OK")
        if differing:
            add("Duplicate metadata IDs with different values", len(differing), ", ".join(differing[:8]),
                "Conflict")
        metadata_ids = set(ids)
        matched = [sample for sample in variant_ids if sample in metadata_ids]
        unmatched = [sample for sample in variant_ids if sample not in metadata_ids]
        unused = sorted(metadata_ids - set(variant_ids))
        add("Variant samples matched to a metadata row", len(matched),
            f"{len(matched)} of {len(variant_ids)} (identical Sample_ID text)",
            "OK" if not unmatched else "Partial")
        review["matching_summary"].update(matched=len(matched), unmatched=len(unmatched), unused=len(unused))
        if unmatched:
            add("Variant samples without a metadata row", len(unmatched), ", ".join(unmatched[:10]),
                "Accepted: values unknown" if options.get("allow_unmatched_metadata") else "Problem")
        if unused:
            add("Metadata rows without variant records", len(unused), ", ".join(unused[:10])
                + ". No mutation absence is inferred for them.", "Info")
        hints = {}
        by_key = {}
        for sample in metadata_ids:
            by_key.setdefault(loose(sample), []).append(sample)
        for sample in unmatched:
            for candidate in by_key.get(loose(sample), []):
                hints[sample] = candidate
        if hints:
            add("Possible matches that were NOT joined", len(hints),
                "; ".join(f"'{a}' ~ '{b}'" for a, b in list(hints.items())[:8])
                + ". IDs must be identical; correct the file if these are the same sample.", "Check")
        review["matching_hints"] = hints
        # Conflicts between the variant table and the metadata file.
        first_rows = usable.drop_duplicates("Sample_ID").set_index("Sample_ID")
        for field in FIELDS:
            if field not in first_rows or field not in variants or not len(variants):
                continue
            for sample, values in sample_values(variants, field).items():
                if sample not in first_rows.index or sample in differing:
                    continue
                other = first_rows.at[sample, field]
                if other == "" or len(values) != 1:
                    continue
                left = comparable(field, pd.Series([values[0]]), date_order).iloc[0]
                right = comparable(field, pd.Series([other]), date_order).iloc[0]
                if left != right:
                    conflicts.append({"Sample_ID": sample, "Field": field, "Kind": "Variant table vs metadata",
                                      "Variant_Table_Value": values[0], "Metadata_Value": other,
                                      "Resolution": resolution})
        if unmatched and not options.get("allow_unmatched_metadata"):
            review["problems"].append(
                "Pooled sample(s) missing from the pooled metadata: " + ", ".join(unmatched[:8])
                + (f", ... ({len(unmatched)} in total)" if len(unmatched) > 8 else "")
                + ". Correct the IDs, or tick 'Continue without metadata for unmatched samples' under "
                  "'Sample metadata matching' (their dates and sites stay unknown; nothing is guessed).")
    # Specimens versus libraries (technical replicates are verified later with the Phase 3 rules).
    specimen_source = metadata if metadata is not None and "Specimen_ID" in metadata else variants
    if len(specimen_source) and "Specimen_ID" in specimen_source:
        table = specimen_source[specimen_source["Specimen_ID"].ne("")]
        groups = table.groupby("Specimen_ID")["Sample_ID"].nunique()
        shared = groups[groups > 1]
        if len(shared):
            libraries = table.get("Library_ID", pd.Series([""] * len(table), index=table.index))
            documented = sum(1 for specimen in shared.index
                             if libraries[table["Specimen_ID"].eq(specimen)].ne("").all()
                             and libraries[table["Specimen_ID"].eq(specimen)].nunique()
                             == table.loc[table["Specimen_ID"].eq(specimen), "Sample_ID"].nunique())
            add("Specimens with several sequenced samples", len(shared),
                f"{documented} have a distinct documented Library_ID for every sample (replicate candidates, "
                f"verified by the Phase 3 rules during the run); {len(shared) - documented} do not and are "
                "never treated as technical replicates. Specimens and libraries are kept separate.", "Info")
    references = {}
    for entry in entries:
        for reference in entry["info"].get("reference_ids", []):
            references.setdefault(reference, entry["info"]["name"])
    if len(references) > 1:
        add("Reference sequences (genome segments) in the variant data", len(references),
            ", ".join(list(references)[:10]) + ". Each is analyzed separately; positions on different "
            "references or segments are never compared.", "Info")
    statuses = review.get("pooled_date_status", {})
    unusable = {status: count for status, count in statuses.items() if status not in {"Parsed", "Not supplied"}}
    if unusable:
        add("Collection dates that cannot be used as written", sum(unusable.values()),
            "; ".join(f"{status}: {count}" for status, count in unusable.items()), "Warning")
    if conflicts:
        kinds = pd.Series([item["Kind"] for item in conflicts]).value_counts().to_dict()
        add("Conflicting metadata values", len(conflicts), "; ".join(f"{k}: {v}" for k, v in kinds.items()),
            "Resolved: " + RESOLUTIONS[resolution] if resolution else "Problem")
        needs_unresolved = [item for item in conflicts if item["Kind"] != "Variant table vs metadata"]
        if not resolution:
            review["problems"].append(
                f"{len(conflicts)} metadata value(s) conflict (listed under 'Sample metadata matching'). Choose how "
                "to handle them: no value is chosen silently.")
        elif needs_unresolved and resolution != "unresolved":
            review["problems"].append(
                f"{len(needs_unresolved)} conflict(s) are between duplicate metadata rows or within one sample; "
                "only 'Leave unresolved' (or correcting the file) can settle them.")
    review["matching"] = rows
    review["metadata_conflicts"] = conflicts
    return rows


def resolve_variant_table(table, conflicts):
    """Apply the chosen resolutions to a prepared (standard-name) variant table; returns a copy."""

    table = table.copy()
    notes = []
    for item in conflicts:
        field, sample, resolution = item["Field"], item["Sample_ID"], item["Resolution"]
        if (field not in table or resolution not in {"metadata", "unresolved"}
                or item["Kind"] == "Duplicate metadata rows"):
            continue
        rows = table["Sample_ID"].astype(str).str.strip().eq(sample)
        table.loc[rows, field] = ""
        notes.append((sample, field, resolution))
    if notes:
        text = {}
        for sample, field, resolution in notes:
            text.setdefault(sample, []).append(f"{field}: {RESOLUTIONS[resolution].lower()}")
        table["Source_Metadata_Resolution"] = table["Sample_ID"].astype(str).str.strip().map(
            lambda sample: "; ".join(text.get(sample, [])))
    return table


def resolve_metadata_table(table, conflicts, variant_ids, allow_unmatched):
    """Prepared pooled metadata: blank IDs dropped, identical duplicates kept once, resolutions applied,
    and (when the user allowed it) one empty row per unmatched variant sample."""

    table = table.copy()
    table["Sample_ID"] = clean(table["Sample_ID"])
    table = table[table["Sample_ID"].ne("")]
    for item in conflicts:
        field, sample, resolution = item["Field"], item["Sample_ID"], item["Resolution"]
        if field not in table:
            continue
        rows = table["Sample_ID"].eq(sample)
        if item["Kind"] == "Duplicate metadata rows" or resolution in {"variant_table", "unresolved"}:
            table.loc[rows, field] = ""
    table = table.drop_duplicates()
    table = table.drop_duplicates("Sample_ID")
    status = pd.Series("Joined by Sample_ID", index=table.index)
    if allow_unmatched:
        missing = [sample for sample in variant_ids if sample not in set(table["Sample_ID"])]
        if missing:
            extra = pd.DataFrame({"Sample_ID": missing})
            table = pd.concat([table, extra], ignore_index=True).fillna("")
            status = pd.concat([status.reset_index(drop=True),
                                pd.Series(["No metadata row (values unknown)"] * len(missing))], ignore_index=True)
            table = table.reset_index(drop=True)
    table["Metadata_Status"] = status.to_numpy()
    return table
