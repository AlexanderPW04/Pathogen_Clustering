"""Canonical sequencing-data tables, independent of the uploaded file's column names.

Sprint 5 Phase 4A. Every supported upload (the published Fontenele table, iVar,
VCF, standard long tables, wide sample-by-mutation matrices, and base-count
tables) is first rewritten by `genomics.input_inspection.prepare_inputs` into
the standardized columns that `wastewater.processed_snv` validates. The tables
below expose the result of that import as three linked structures:

A. ``canonical_samples.csv``: one row per sample (sequencing library) with
   its specimen, library, date, site, coordinates, organism, references
   (one per genome segment), platform, and where each value came from.
B. ``canonical_variant_observations.csv``: one row per reported SNV
   observation. Reported and calculated allele frequencies are separate
   columns; read counts are never invented. A frequency that the importer
   calculated from base counts is *not* a reported frequency.
C. ``canonical_coverage.csv``: per-position A/C/G/T counts and callable depth
   from a per-position table, kept apart from the sparse SNV observations.
   Without such a table it is empty, and an unreported variant stays
   "coverage unknown", never 0.

``import_provenance.csv`` records, for every file, which column supplied each
standard field and every transformation applied (percent to fraction, 0-based
to 1-based positions, stated reference, wide-to-long, base-count derivation).

These functions only select and rename values that the existing analyses
already computed; they change no calculation.
"""

import numpy as np
import pandas as pd

from genomics.unified_schema import POOLED

SAMPLE_COLUMNS = [
    "Sample_ID", "Specimen_ID", "Library_ID", "Sample_Type", "Collection_Date",
    "Collection_Date_Original", "Collection_Date_Status", "Collection_Site", "Latitude", "Longitude",
    "Geographic_Precision", "Organism", "Reference_IDs", "Genome_Segments", "Sequencing_Platform",
    "Genome_Coverage_Percent", "Mean_Depth", "Coverage_Information", "Reported_Observations",
    "Sample_QC", "Source_Dataset", "Metadata_Source",
]
OBSERVATION_COLUMNS = [
    "Sample_ID", "Reference_ID", "Genome_Position", "Ref", "Alt", "Variant_Key",
    "Ref_Read_Count", "Alt_Read_Count", "Total_Depth", "Reported_VAF", "Calculated_VAF",
    "VAF_Used", "Count_Provenance", "Ref_Count_Provenance", "Depth_Provenance",
    "Frequency_Provenance", "Measurement_Class", "Frequency_Check", "QC_Status", "Analysis_Included",
    "Duplicate_Status", "Reference_Status", "Source_Dataset",
]
COVERAGE_COLUMNS = [
    "Sample_ID", "Reference_ID", "Genome_Position", "A_Count", "C_Count", "G_Count", "T_Count",
    "Total_Callable_Depth", "Depth_Provenance", "Coverage_Status", "Minimum_Depth",
]
PROVENANCE_COLUMNS = ["File", "Role", "Layout", "Item", "Source", "Transformation", "Decided_By"]

# The importer writes this column when it calculated a frequency itself.
CALCULATED_BY_IMPORTER = "Calculated by the importer"

SCHEMA = {
    "canonical_samples.csv": {
        "Sample_ID": "Identifier of one sequenced sample (one library), as written in the upload",
        "Specimen_ID": "Biological specimen; several libraries of one specimen are technical replicates only "
                       "with distinct documented Library_IDs",
        "Library_ID": "Separately prepared sequencing library",
        "Collection_Date": "Parsed date (YYYY-MM-DD); blank when missing or ambiguous",
        "Collection_Date_Original": "Date text exactly as uploaded",
        "Collection_Site": "Most specific collection location supplied",
        "Latitude/Longitude": "Supplied coordinates only (never a centroid)",
        "Organism": "From the reference FASTA or the stated reference; not inferred from sequences",
        "Reference_IDs": "Reference accession(s) of the sample's observations",
        "Genome_Segments": "One reference per genome segment; segments are never compared with each other",
        "Sequencing_Platform": "From metadata when supplied",
        "Metadata_Source": "Where dates and sites came from",
    },
    "canonical_variant_observations.csv": {
        "Genome_Position": "1-based position on Reference_ID",
        "Reported_VAF": "Alternative-allele frequency reported by the source; blank when the source reported none",
        "Calculated_VAF": "Alternative reads / total depth, from measured or reconstructed counts",
        "VAF_Used": "Frequency used by the Sprint 5 uncertainty analyses",
        "Count_Provenance": "Measured, reconstructed (approximate), or frequency only",
        "Measurement_Class": "Which measurements the record has (counts, depth, frequency)",
    },
    "canonical_coverage.csv": {
        "Total_Callable_Depth": "A+C+G+T reads, or the supplied depth",
        "Coverage_Status": "Callable (depth at least the minimum), low coverage, or not covered",
    },
}


def text(series):
    return series.fillna("").astype(str).str.strip()


def measurement_class(reported, alt, depth):
    """Which measurements a record carries; nothing missing is counted as present."""

    has_freq, has_alt, has_depth = np.isfinite(reported), np.isfinite(alt), np.isfinite(depth)
    labels = np.full(len(reported), "Frequency only (no depth)", dtype=object)
    labels[has_freq & has_depth] = "Frequency and depth (counts not supplied)"
    labels[has_alt & has_depth] = "Read counts and depth"
    labels[has_alt & has_depth & ~has_freq] = "Read counts and depth (no reported frequency)"
    return labels


def canonical_observations(result):
    """Table B from the Sprint 5 evidence table and the QC'd pooled variants (same row order)."""

    pooled = result.get("pooled")
    evidence = result["tables"].get("sprint5_vaf_observations.csv")
    if pooled is None or evidence is None:
        return pd.DataFrame(columns=OBSERVATION_COLUMNS)
    variants = pooled["variants"].reset_index(drop=True)
    evidence = evidence.reset_index(drop=True)
    same_rows = (len(variants) == len(evidence)
                 and text(variants["Sample_ID"]).equals(text(evidence["Sample_ID"]))
                 and text(variants["Variant_Key"]).equals(text(evidence["Variant_Key"])))
    if not same_rows:
        raise ValueError("The evidence table no longer matches the pooled variant records.")
    source = variants.get("Source_Frequency_Provenance", pd.Series([""] * len(variants)))
    calculated = text(source).str.startswith(CALCULATED_BY_IMPORTER).to_numpy()
    reported = pd.to_numeric(evidence["Reported_Frequency"], errors="coerce").to_numpy(dtype=float)
    reported = np.where(calculated, np.nan, reported)
    alt = pd.to_numeric(evidence["Alt_Read_Count"], errors="coerce").to_numpy(dtype=float)
    depth = pd.to_numeric(evidence["Total_Depth"], errors="coerce").to_numpy(dtype=float)
    frequency_source = np.where(calculated, text(source).to_numpy(), "Reported by the source")
    table = pd.DataFrame({
        "Sample_ID": evidence["Sample_ID"],
        "Reference_ID": evidence["Reference_ID"],
        "Genome_Position": evidence["Position"],
        "Ref": evidence["Ref"],
        "Alt": evidence["Alt"],
        "Variant_Key": evidence["Variant_Key"],
        "Ref_Read_Count": evidence["Ref_Read_Count"],
        "Alt_Read_Count": evidence["Alt_Read_Count"],
        "Total_Depth": evidence["Total_Depth"],
        "Reported_VAF": reported,
        "Calculated_VAF": evidence["Count_Derived_Frequency"],
        "VAF_Used": evidence["VAF"],
        "Count_Provenance": evidence["Alt_Count_Source"],
        "Ref_Count_Provenance": evidence["Ref_Count_Source"],
        "Depth_Provenance": evidence["Depth_Source"],
        "Frequency_Provenance": frequency_source,
        "Measurement_Class": measurement_class(reported, alt, depth),
        "Frequency_Check": evidence["Frequency_Check"],
        "QC_Status": evidence["QC_Status"],
        "Analysis_Included": evidence["Analysis_Included"],
        "Duplicate_Status": evidence["Duplicate_Status"],
        "Reference_Status": evidence["Reference_Status"],
        "Source_Dataset": (variants["Source_File"].to_numpy() if "Source_File" in variants
                           else evidence["Data_Source"]),
    })
    # Original values that an import transformation changed (kept for traceability).
    for name in variants.columns:
        if name.startswith("Source_") and name.endswith("_Original"):
            table[name] = variants[name].to_numpy()
    return table


def canonical_samples(result, observations=None):
    """Table A: one row per pooled sample, values exactly as the analyses used them."""

    pooled = result.get("pooled")
    if pooled is None:
        return pd.DataFrame(columns=SAMPLE_COLUMNS)
    samples = pooled["samples"].reset_index(drop=True)
    variants = pooled["variants"]
    evidence = result["tables"].get("sprint5_vaf_observations.csv")
    by_sample = variants.groupby("Sample_ID")

    def first(table, column):
        if table is None or column not in table:
            return {}
        values = text(table[column])
        return values[values.ne("")].groupby(table.loc[values.ne(""), "Sample_ID"]).first().to_dict()

    def joined(column):
        return by_sample[column].agg(lambda values: "; ".join(sorted(set(text(values)) - {""}))).to_dict()

    ids = samples["Sample_ID"]
    references = joined("Reference_ID")
    segments = {sample: (f"{len(names.split('; '))} reference sequences (segments): {names}"
                         if names and len(names.split("; ")) > 1 else "Single reference sequence")
                for sample, names in references.items()}
    original = first(variants, "Collection_Date_Original")
    status = first(variants, "Collection_Date_Status")
    organism = first(evidence, "Organism")
    specimen = first(evidence, "Specimen_ID")
    library = first(evidence, "Library_ID")
    platform = first(pooled.get("metadata"), "Sequencing_Platform")
    counts = by_sample.size().to_dict()
    metadata = ("Pooled metadata file joined by Sample_ID" if pooled.get("metadata_supplied")
                else "Columns of the variant table")
    table = pd.DataFrame({
        "Sample_ID": ids,
        "Specimen_ID": ids.map(specimen).fillna(""),
        "Library_ID": ids.map(library).fillna(""),
        "Sample_Type": samples.get("Sample_Type", POOLED),
        "Collection_Date": samples["Collection_Date"],
        "Collection_Date_Original": ids.map(original).fillna(""),
        "Collection_Date_Status": ids.map(status).fillna(""),
        "Collection_Site": samples["Location_Name"],
        "Latitude": samples["Latitude"],
        "Longitude": samples["Longitude"],
        "Geographic_Precision": samples.get("Geographic_Precision", ""),
        "Organism": ids.map(organism).fillna(""),
        "Reference_IDs": ids.map(references).fillna(""),
        "Genome_Segments": ids.map(segments).fillna(""),
        "Sequencing_Platform": ids.map(platform).fillna(""),
        "Genome_Coverage_Percent": samples["Genome_Coverage_Percent"],
        "Mean_Depth": samples["Mean_Depth"],
        "Coverage_Information": samples["Coverage_Information"],
        "Reported_Observations": ids.map(counts).fillna(0).astype(int),
        "Sample_QC": samples.get("QC_Genome_Coverage", ""),
        "Source_Dataset": samples["Data_Source"],
        "Metadata_Source": metadata,
    })
    return table[SAMPLE_COLUMNS]


def canonical_coverage(result, minimum_depth):
    """Table C from a per-position table; empty (never zero-filled) without one."""

    pooled = result.get("pooled")
    positions = None if pooled is None else pooled.get("position_table")
    if positions is None or positions.empty:
        return pd.DataFrame(columns=COVERAGE_COLUMNS)
    table = pd.DataFrame({
        "Sample_ID": positions["Sample_ID"],
        "Reference_ID": positions["Reference_ID"],
        "Genome_Position": positions["Position"],
    })
    counts = [f"{base}_Count" for base in "ACGT"]
    has_counts = all(name in positions for name in counts)
    for name in counts:
        table[name] = positions[name] if name in positions else pd.NA
    if has_counts:
        depth = positions[counts].apply(pd.to_numeric, errors="coerce").sum(axis=1, min_count=4)
        provenance = "A+C+G+T read counts"
    else:
        depth = pd.to_numeric(positions.get("Depth"), errors="coerce")
        provenance = "Supplied depth (no base counts)"
    table["Total_Callable_Depth"] = depth
    table["Depth_Provenance"] = provenance
    status = np.where(depth.isna(), "Depth not supplied",
                      np.where(depth >= minimum_depth, "Callable",
                               np.where(depth > 0, "Low coverage", "Not covered")))
    table["Coverage_Status"] = status
    table["Minimum_Depth"] = minimum_depth
    return table[COVERAGE_COLUMNS]


def add_canonical_tables(result, provenance_rows=None):
    """Add the three canonical tables and the import provenance to the result (in place)."""

    if result.get("pooled") is None:
        return result
    minimum = result["settings"].get("min_variant_depth", 20)
    observations = canonical_observations(result)
    tables = result["tables"]
    tables["canonical_samples.csv"] = canonical_samples(result, observations)
    tables["canonical_variant_observations.csv"] = observations
    tables["canonical_coverage.csv"] = canonical_coverage(result, minimum)
    rows = list(provenance_rows or [])
    tables["import_provenance.csv"] = pd.DataFrame(rows, columns=PROVENANCE_COLUMNS)
    result["canonical"] = {
        "samples": len(tables["canonical_samples.csv"]),
        "observations": len(observations),
        "coverage_records": len(tables["canonical_coverage.csv"]),
        "measurement_classes": observations["Measurement_Class"].value_counts().to_dict(),
        "reported_vaf": int(pd.notna(observations["Reported_VAF"]).sum()),
        "calculated_vaf": int(pd.to_numeric(observations["Calculated_VAF"], errors="coerce").notna().sum()),
    }
    return result
