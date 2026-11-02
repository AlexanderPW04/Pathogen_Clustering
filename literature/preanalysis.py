"""Sprint 4 before the genomic analysis: find metadata gaps, apply approved values.

When the upload-page toggle is on, enrichment runs between the dataset review
and the analysis. This module builds the sample table from the prepared
inputs (no MAFFT, no GeoNames), lists what could be improved, and turns the
approved values into analysis inputs: approved individual values are written
into the prepared metadata file; approved pooled values are passed to the
integrated analysis as overrides. Original files and values are kept.
"""

import re

import pandas as pd
from Bio import SeqIO

from add_temporal_geographic import enrich_pairs
from geocode_locations import is_state_only_location
from genomics.file_formats import read_table
from genomics.pooled_observations import load_pooled_metadata, pooled_sample_records
from literature.apply import PAIR_GEOGRAPHY, apply_approved, rescreen
from literature.geography import original_geography
from literature.workflow import INDIVIDUAL, POOLED, new_session
from validate_inputs import is_missing_value, load_metadata, normalize_accession, validate_metadata
from wastewater.processed_snv import load_processed_snvs

PRECISION_COARSE = "State-level location (too coarse for distance cutoff)"
SECTIONS_BY_SOURCE = {
    INDIVIDUAL: ["1. Overview", "2. Cross-source evidence", "4. Temporal patterns", "5. Geographic context",
                 "8. Individual sequence relationships"],
    POOLED: ["1. Overview", "2. Cross-source evidence", "4. Temporal patterns", "5. Geographic context",
             "7. Population structure (colors and sites)"],
}


def location_precision(row):
    """The precision label the Sprint 2 pipeline would give before any GeoNames lookup."""

    if is_missing_value(row.get("Location_Name")):
        return "No location given"
    if is_state_only_location(row.get("Location_Name")):
        return PRECISION_COARSE
    if pd.notna(row.get("Latitude")) and pd.notna(row.get("Longitude")):
        return "Provided coordinates (precision unverified)"
    return "Named location only (precision unverified)"


def individual_samples(paths):
    """Validated individual metadata with the pipeline's Sample_ID (the FASTA record ID)."""

    if "metadata" not in paths or "fasta" not in paths:
        return None, None
    raw = read_table(paths["metadata"])
    validated, _, _ = validate_metadata(load_metadata(paths["metadata"]), [])
    with open(paths["fasta"], encoding="utf-8-sig") as handle:
        fasta_ids = {normalize_accession(record.id): str(record.id) for record in SeqIO.parse(handle, "fasta")}
    samples = validated.copy()
    samples["Sample_ID"] = samples["Accession"].map(lambda value: fasta_ids.get(normalize_accession(value)))
    samples = samples[samples["Sample_ID"].notna()].copy()
    for column in ["Latitude", "Longitude", "Location_Name", "Collection_Date"]:
        if column not in samples:
            samples[column] = pd.NA
    samples["Latitude"] = pd.to_numeric(samples["Latitude"], errors="coerce")
    samples["Longitude"] = pd.to_numeric(samples["Longitude"], errors="coerce")
    samples["Geographic_Precision"] = samples.apply(location_precision, axis=1)
    samples["Coordinate_Source"] = samples.apply(
        lambda row: "Provided" if pd.notna(row["Latitude"]) and pd.notna(row["Longitude"]) else pd.NA, axis=1)
    return samples.reset_index(drop=True), raw


def pooled_samples(paths):
    """One row per pooled sample, as the integrated analysis builds it."""

    if "snv" not in paths:
        return None, None
    snv_data, _ = load_processed_snvs(paths["snv"], paths.get("pooled_metadata"))
    metadata = load_pooled_metadata(paths["pooled_metadata"]) if paths.get("pooled_metadata") else None
    samples = pooled_sample_records(snv_data, "pooled SNV table", metadata)
    raw = read_table(paths["pooled_metadata"]) if paths.get("pooled_metadata") else None
    return samples, raw


def metadata_gaps(samples, raw_individual=None):
    """Per-sample metadata that literature could fill or make more specific."""

    raw_dates = {}
    if raw_individual is not None and "Collection_Date" in raw_individual and "Accession" in raw_individual:
        raw_dates = {normalize_accession(key): str(value) for key, value in
                     zip(raw_individual["Accession"], raw_individual["Collection_Date"])}
    rows = []
    for _, row in samples.iterrows():
        gaps = []
        geography = original_geography(row)
        if geography["level"] == 0:
            gaps.append("No collection location")
        elif geography["level"] == 1:
            gaps.append("Country-only location")
        elif geography["level"] == 2 or geography.get("coarse"):
            gaps.append("State-only location")
        elif geography["level"] < 6:
            gaps.append("Named place without coordinates (may be ambiguous; site not specified)")
        date = pd.to_datetime(row.get("Collection_Date"), errors="coerce")
        written = raw_dates.get(normalize_accession(row.get("Accession", row["Sample_ID"])), "")
        if pd.isna(date):
            if re.fullmatch(r"\s*\d{4}(-\d{1,2})?\s*", written or ""):
                gaps.append("Partial collection date (year or month only)")
            else:
                gaps.append("No collection date")
        if row.get("Sample_Source") == INDIVIDUAL and is_missing_value(row.get("Host_Source")):
            gaps.append("No sample source")
        rows.append({"Sample_ID": row["Sample_ID"], "Sample_Source": row["Sample_Source"],
                     "Location": "" if is_missing_value(row.get("Location_Name")) else str(row.get("Location_Name")),
                     "Collection_Date": "" if pd.isna(date) else date.date().isoformat(),
                     "Gaps": gaps})
    return rows


def preview_session(paths, dataset_identifiers=(), synthetic=False, label=""):
    """An enrichment session built from prepared inputs, with its metadata gaps."""

    individual, raw_individual = individual_samples(paths)
    pooled, raw_pooled = pooled_samples(paths)
    session = new_session(individual, pooled, synthetic=synthetic, label=label,
                          individual_metadata=raw_individual, pooled_metadata=raw_pooled)
    session["stage"] = "pre-analysis"
    session["gaps"] = metadata_gaps(session["samples"], raw_individual)
    session["dataset_identifiers"] = list(dataset_identifiers)
    session["eligible"] = [row["Sample_ID"] for row in session["gaps"] if row["Gaps"]]
    return session


def approved_changes(session):
    """Derived metadata and the application log for the recorded decisions (no GeoNames here:
    individual place names are resolved by the pipeline's own GeoNames step)."""

    derived, log = apply_approved(session["samples"], session["facts"], session["decisions"], geocoder=None)
    for item in log:
        item["Outcome"] = item["Outcome"].replace(
            "Applied (place name only: GeoNames not used)",
            "Applied as a place name (GeoNames resolves it during the analysis when configured; approximate)")
    return derived, log


def write_individual_metadata(path, derived):
    """Write approved individual values into the prepared metadata file."""

    changed = derived[(derived["Sample_Source"] == INDIVIDUAL) & derived["Enrichment_Sources"].notna()]
    if changed.empty:
        return []
    table = read_table(path)
    for column in ["Location_Name", "Latitude", "Longitude", "Collection_Date"]:
        if column not in table:
            table[column] = ""
    keys = table["Accession"].map(normalize_accession)
    written = []
    for _, row in changed.iterrows():
        match = keys.eq(normalize_accession(row["Sample_ID"]))
        if not match.any():
            continue
        table.loc[match, "Location_Name"] = "" if pd.isna(row["Location_Name"]) else str(row["Location_Name"])
        # GeoNames estimates are left to the pipeline; only publication coordinates are written.
        publication = str(row["Coordinate_Source"]).startswith("Publication")
        table.loc[match, "Latitude"] = str(row["Latitude"]) if publication else ""
        table.loc[match, "Longitude"] = str(row["Longitude"]) if publication else ""
        date = pd.to_datetime(row["Collection_Date"], errors="coerce")
        if pd.notna(date):
            table.loc[match, "Collection_Date"] = date.date().isoformat()
        written.append(row["Sample_ID"])
    table.to_csv(path, index=False)
    return written


def pooled_overrides(derived):
    """{Sample_ID: values} for approved pooled changes, applied inside the integrated analysis."""

    changed = derived[(derived["Sample_Source"] == POOLED) & derived["Enrichment_Sources"].notna()]
    overrides = {}
    for _, row in changed.iterrows():
        values = {}
        for column in ["Location_Name", "Latitude", "Longitude", "Coordinate_Source", "Geographic_Precision"]:
            new, old = row[column], row[f"Original_{column}"]
            if pd.notna(new) and (pd.isna(old) or new != old):
                values[column] = new
        if pd.notna(row["Collection_Date"]) and pd.isna(row["Original_Collection_Date"]):
            values["Collection_Date"] = pd.Timestamp(row["Collection_Date"])
        if values:
            overrides[str(row["Sample_ID"])] = values
    return overrides


def before_after(session, derived, run_samples, pairs, thresholds):
    """Sprint 2 screening with the original values versus the values the run used."""

    if run_samples is None or pairs is None or pairs.empty or thresholds is None:
        return None, {}
    after = run_samples.copy()
    before = run_samples.copy()
    original = derived.set_index("Sample_ID")
    for index, row in before.iterrows():
        sample_id = str(row["Sample_ID"])
        if sample_id not in original.index or pd.isna(original.at[sample_id, "Enrichment_Sources"]):
            continue
        for column in ["Location_Name", "Latitude", "Longitude", "Collection_Date", "Geographic_Precision",
                       "Coordinate_Source"]:
            before.at[index, column] = original.at[sample_id, f"Original_{column}"]
    after["Enrichment_Sources"] = after["Sample_ID"].astype(str).map(original["Enrichment_Sources"])
    genetic = pairs[[column for column in pairs.columns if column not in PAIR_GEOGRAPHY]]
    if genetic.empty:
        return None, {}
    return rescreen(before, after, genetic, thresholds)


def contribution(session, derived, log, sections_changed):
    """Counts for the results overview: did Sprint 4 improve the dataset?"""

    publications = session["publications"]
    decisions = session["decisions"]
    linked = publications[publications["Eligible_For_Enrichment"]] if publications is not None else None
    return {
        "Samples eligible for enrichment": len(session.get("eligible", [])),
        "Publications linked": 0 if linked is None else int(linked["Article_ID"].nunique()),
        "Samples with supporting publication evidence": 0 if linked is None else int(linked["Sample_ID"].nunique()),
        "Proposed metadata improvements": sum(fact["Approvable"] for fact in session["facts"]),
        "Approved": sum(choice == "approve" for choice in decisions.values()),
        "Rejected or left unresolved": sum(fact["Approvable"] and decisions.get(fact["Fact_ID"]) != "approve"
                                           for fact in session["facts"]),
        "Values changed": len([item for item in log if item["Outcome"].startswith("Applied")]),
        "Sections that used enriched metadata": ", ".join(sections_changed) or "None",
    }


def sections_using(derived):
    changed = derived[derived["Enrichment_Sources"].notna()]
    sections = []
    for source in [INDIVIDUAL, POOLED]:
        if (changed["Sample_Source"] == source).any():
            sections += [name for name in SECTIONS_BY_SOURCE[source] if name not in sections]
    return sorted(sections, key=lambda name: int(name.split(".")[0]))
