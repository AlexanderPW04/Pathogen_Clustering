"""Recognize table columns by their headers, without guessing.

Each file role has a list of standard fields and the header spellings
accepted for each one. Headers are compared after removing case, spaces,
underscores, hyphens, dots, and brackets, so "Collection Date",
"collection_date", and "COLLECTION-DATE" all match.

When exactly one column matches a field, it is used and the evidence is
recorded. When two or more columns match the same field, nothing is chosen:
the user must pick one on the review page. A missing optional field disables
only the analyses that need it.
"""

import pandas as pd


def squash(name):
    """Lower-case a header and keep only its letters and digits."""

    return "".join(character for character in str(name).casefold() if character.isalnum())


def field(name, aliases, required=False, purpose=""):
    return {"name": name, "aliases": {squash(name)} | {squash(alias) for alias in aliases},
            "required": required, "purpose": purpose}


SAMPLE_ALIASES = ["sample", "sample_id", "sampleid", "sample name", "samplename", "sample code",
                  "sample_accession"]
DATE_ALIASES = ["collection_date", "date", "dates", "sample_date", "sampling_date", "collection date",
                "date_collected", "collected"]
LOCATION_ALIASES = ["location_name", "location", "geo_location", "geographic_location",
                    "geo_loc_name", "site", "sampling_site", "site_name", "sampling_location"]
# Broader geographic levels are separate fields, not Location_Name candidates:
# NCBI Virus tables have Geo_Location ("USA: SD"), Country ("USA"), and a
# column named "USA" with the US state ("SD"), which describe one hierarchy.
COUNTRY_ALIASES = ["country", "country_name", "nation"]
STATE_ALIASES = ["state", "province", "state_province", "admin1", "division", "usa"]
# Technical replicates (separate libraries of one specimen) share a Specimen_ID.
SPECIMEN_ALIASES = ["specimen_id", "specimen", "specimenid", "biological_sample", "replicate_of",
                    "replicate_group"]
# Phase 3: a documented library (or technical-replicate) identifier per sample.
# Without distinct Library_IDs, a shared Specimen_ID is an unverified replicate claim.
LIBRARY_ALIASES = ["library_id", "library", "library_name", "libraryid", "library_prep_id",
                   "replicate_id", "technical_replicate", "technical_replicate_id"]
LATITUDE_ALIASES = ["latitude", "lat", "coordinates_lat", "coordinate_lat", "decimal_latitude"]
LONGITUDE_ALIASES = ["longitude", "lon", "long", "lng", "coordinates_long", "coordinates_lon",
                     "coordinate_long", "decimal_longitude"]
# Phase 4A: the instrument or platform. A column named only "sequencing" is not
# an alias: in Crits-Christoph et al. (2021) Table S1 it means enriched/unenriched.
PLATFORM_ALIASES = ["sequencing_platform", "platform", "instrument", "instrument_platform",
                    "instrument_model", "sequencer", "seq_platform"]
# Audit P1-4: segment and chromosome columns name the sequence each record is on. They
# are reference columns, so records on different segments are never put on one sequence.
REFERENCE_ALIASES = ["reference_genome", "reference", "reference_id", "ref_genome", "chrom", "scaffold",
                     "contig", "region", "chromosome", "chr", "seqid", "seq_id", "seqname", "segment"]
SEQUENCE_COLUMN_NAMES = {"".join(c for c in alias if c.isalnum()) for alias in REFERENCE_ALIASES
                         if alias not in {"reference", "reference_id", "reference_genome", "ref_genome"}}
SEQUENCE_COLUMN_NAMES |= {"molecule", "replicon", "genomesegment", "segmentname", "segmentid"}
# Phase 4A. "Site" is deliberately not a position alias: it usually names a sampling site.
POSITION_ALIASES = ["position", "pos", "reference_position", "genome_position", "genomic_position",
                    "nucleotide_position"]
# "Coverage" is accepted as depth only for per-record values; genome-coverage
# percentages have their own field and names.
DEPTH_ALIASES = ["read_depth", "depth", "dp", "total_dp", "total_depth", "coverage", "position_coverage",
                 "read_coverage", "base_coverage", "basecoverage"]

# Field lists per file role. Individual metadata keeps the Sprint 2 names
# (Accession is the identifier that must match the FASTA record IDs).
ROLE_FIELDS = {
    "metadata": [
        field("Accession", ["accession", "accession_number", "accession_id", "sequence_id",
                            "isolate_id"] + SAMPLE_ALIASES,
              required=True, purpose="Identifier that matches each FASTA record"),
        # NCBI Virus exports also have a taxonomic "Species" column; it is not
        # an alias here, so Organism_Name is used without an extra question.
        field("Pathogen", ["organism", "organism_name", "virus", "virus_name", "pathogen_name"],
              purpose="Groups sequences that may be aligned together"),
        field("Collection_Date", DATE_ALIASES, purpose="Temporal screening"),
        field("Location_Name", LOCATION_ALIASES, purpose="Most specific collection location"),
        field("Country", COUNTRY_ALIASES, purpose="Country (kept separately; fills a missing location)"),
        field("State_Province", STATE_ALIASES,
              purpose="State or province (kept separately; fills a missing location)"),
        field("Latitude", LATITUDE_ALIASES, purpose="Kilometer screening"),
        field("Longitude", LONGITUDE_ALIASES, purpose="Kilometer screening"),
        field("Host_Source", ["host", "host_source", "isolation_source", "host_name"],
              purpose="Context only"),
    ],
    "pooled_metadata": [
        field("Sample_ID", SAMPLE_ALIASES, required=True, purpose="Joins metadata to variant records"),
        field("Specimen_ID", SPECIMEN_ALIASES,
              purpose="Technical replicates of one specimen (overdispersion calibration)"),
        field("Library_ID", LIBRARY_ALIASES,
              purpose="Separately prepared library of the specimen (documents technical replicates)"),
        field("Collection_Date", DATE_ALIASES, purpose="Temporal patterns"),
        field("Location_Name", LOCATION_ALIASES, purpose="Site grouping"),
        field("Latitude", LATITUDE_ALIASES, purpose="Geographic context"),
        field("Longitude", LONGITUDE_ALIASES, purpose="Geographic context"),
        field("Genome_Coverage_Percent", ["genome_coverage_percent", "coverage_percent",
                                          "percent_genome_coverage", "genome_coverage",
                                          "pct_genome_coverage", "breadth_of_coverage"],
              purpose="Paper coverage filter (>90% at 10x)"),
        field("Mean_Depth", ["mean_depth", "average_depth", "avg_depth", "mean_read_depth",
                             "mean_coverage_depth"], purpose="Depth-versus-SNV regression"),
        field("Sequencing_Platform", PLATFORM_ALIASES, purpose="Context only (kept in the canonical sample table)"),
    ],
    "snv": [
        field("Sample_ID", SAMPLE_ALIASES, required=True, purpose="Which pooled sample"),
        # Phase 4A: a table without a reference column can still be used once the user
        # states the reference (checked later against the reference FASTA).
        field("Reference_Genome", REFERENCE_ALIASES,
              purpose="Reference coordinate system (or state it on the review page)"),
        field("Position", POSITION_ALIASES, required=True,
              purpose="1-based reference position"),
        field("Ref", ["ref", "ref_base", "reference_base", "ref_allele", "reference_allele"],
              required=True, purpose="Reference base"),
        field("Alt", ["alt", "alt_base", "alternative", "alt_allele", "alternative_allele",
                      "variant_base"], required=True, purpose="Alternative base"),
        field("Alt_Frequency", ["alt_frequency", "af", "frequency", "freq", "alt_freq",
                                "allele_frequency", "variant_frequency", "vaf", "alt_allele_frequency"],
              required=True, purpose="Variant-frequency analysis"),
        # Phase 4A: depth is optional (frequency-only tables); blank is unknown, never 0.
        field("Read_Depth", DEPTH_ALIASES,
              purpose="Depth QC, count-based frequencies, and uncertainty (column may have blanks)"),
        field("Alt_Depth", ["alt_depth", "alt_dp", "alt_count", "alt_reads", "alternative_reads",
                            "alternate_reads", "alt_read_count", "alternative_read_count"],
              purpose="Alternate read counts"),
        field("Ref_Depth", ["ref_depth", "ref_dp", "ref_count", "ref_reads", "reference_reads",
                            "ref_read_count", "reference_read_count"],
              purpose="Reference read counts"),
        field("Collection_Date", DATE_ALIASES, purpose="Temporal patterns"),
        field("Location_Name", LOCATION_ALIASES, purpose="Site grouping"),
        field("Specimen_ID", SPECIMEN_ALIASES,
              purpose="Technical replicates of one specimen (overdispersion calibration)"),
        field("Library_ID", LIBRARY_ALIASES,
              purpose="Separately prepared library of the specimen (documents technical replicates)"),
        field("Latitude", LATITUDE_ALIASES, purpose="Site coordinates (moved to the sample metadata)"),
        field("Longitude", LONGITUDE_ALIASES, purpose="Site coordinates (moved to the sample metadata)"),
    ],
    "positions": [
        field("Sample_ID", SAMPLE_ALIASES, required=True, purpose="Which pooled sample"),
        field("Reference_ID", REFERENCE_ALIASES,
              purpose="Reference coordinate system (or state it on the review page)"),
        field("Position", POSITION_ALIASES, required=True,
              purpose="1-based reference position"),
        field("Ref", ["ref", "ref_base", "reference_base", "ref_allele", "reference_allele"],
              purpose="Reference base (checked against the reference FASTA)"),
        field("Depth", DEPTH_ALIASES, purpose="Coverage at unreported positions"),
        field("A_Count", ["a", "a_count", "count_a", "a_reads"], purpose="Paper PCoA"),
        field("C_Count", ["c", "c_count", "count_c", "c_reads"], purpose="Paper PCoA"),
        field("G_Count", ["g", "g_count", "count_g", "g_reads"], purpose="Paper PCoA"),
        field("T_Count", ["t", "t_count", "count_t", "t_reads"], purpose="Paper PCoA"),
    ],
    "masked_sites": [
        field("Reference_ID", REFERENCE_ALIASES, required=True, purpose="Reference coordinate system"),
        field("Position", ["position", "pos", "reference_position"], required=True,
              purpose="Masked 1-based position"),
    ],
    "signatures": [
        field("Signature_Name", ["signature_name", "signature", "profile", "profile_name"],
              required=True, purpose="Documented profile name"),
        field("Reference_ID", REFERENCE_ALIASES, required=True, purpose="Reference coordinate system"),
        field("Position", ["position", "pos"], required=True, purpose="1-based reference position"),
        field("Ref", ["ref", "ref_base", "reference_base"], required=True, purpose="Reference base"),
        field("Alt", ["alt", "alt_base", "alternative"], required=True, purpose="Alternative base"),
        field("Source", ["source", "citation", "documented_source"], required=True,
              purpose="Where the profile is documented"),
    ],
}


def detect_columns(columns, role):
    """Match the columns of one table against the fields of one role.

    Returns a list with one dictionary per field:
    {"field", "required", "purpose", "candidates", "column", "status", "evidence"}
    where status is "matched", "ambiguous", or "missing".
    """

    results = []
    for spec in ROLE_FIELDS.get(role, []):
        candidates = [column for column in columns if squash(column) in spec["aliases"]]
        if len(candidates) == 1:
            status, column = "matched", candidates[0]
            evidence = (f"exact header '{column}'" if str(column).strip() == spec["name"]
                        else f"header '{column}' is a recognized name for {spec['name']}")
        elif candidates:
            status, column = "ambiguous", None
            evidence = ("several columns could be " + spec["name"] + ": "
                        + ", ".join(f"'{name}'" for name in candidates) + "; please choose one")
        else:
            status, column, evidence = "missing", None, "no recognized column"
        results.append({"field": spec["name"], "required": spec["required"],
                        "purpose": spec["purpose"], "candidates": candidates,
                        "column": column, "status": status, "evidence": evidence})
    return results


def explain_snv_columns(columns):
    """Why a table cannot be read as processed SNV calls, in plain words; None if it can.

    Each sentence is based on what was actually detected in the headers.
    """

    found = {match["field"] for match in detect_columns(columns, "snv") if match["status"] != "missing"}
    labels = {"Sample_ID": "pooled sample ID", "Reference_Genome": "reference genome",
              "Ref": "reference base", "Alt": "alternative base",
              "Alt_Frequency": "allele frequency", "Read_Depth": "read depth"}
    missing = [f"{label} ({name})" for name, label in labels.items() if name not in found]
    if "Position" not in found:
        return ("This spreadsheet does not contain the required variant-position column"
                + (f" (also missing: {', '.join(missing)})." if missing else "."))
    if not missing:
        return None
    if {"Alt_Frequency", "Read_Depth"} - found:
        return ("This workbook contains variant positions but not the processed SNV values required "
                f"for this analysis: no column was recognized for {', '.join(missing)}. These values "
                "cannot be estimated, so they were not filled in.")
    return f"This spreadsheet is missing required SNV columns: {', '.join(missing)}."


def resolve_choices(matches, choices, chooser="the user"):
    """Apply choices ({field: column}) to ambiguous matches; `chooser` says who chose."""

    resolved = []
    for match in matches:
        match = dict(match)
        chosen = choices.get(match["field"])
        if match["status"] == "ambiguous" and chosen in match["candidates"]:
            match["status"] = "matched"
            match["column"] = chosen
            match["evidence"] = f"header '{chosen}' chosen by {chooser} from " + ", ".join(
                f"'{name}'" for name in match["candidates"])
        resolved.append(match)
    return resolved


def unresolved(matches):
    """Return human-readable problems that must be fixed before running."""

    problems = []
    for match in matches:
        if match["status"] == "ambiguous":
            problems.append(f"Choose which column is {match['field']} ("
                            + ", ".join(match["candidates"]) + ").")
        elif match["status"] == "missing" and match["required"]:
            problems.append(f"No column was recognized as the required field {match['field']}.")
    return problems


def apply_mapping(table, matches):
    """Rename matched columns to their standard names.

    Columns that were candidates for a field but not chosen are renamed
    "Ignored_<name>" so a later loader cannot pick them up a second time.
    """

    renames = {}
    for match in matches:
        if match["status"] != "matched":
            continue
        renames[match["column"]] = match["field"]
        for other in match["candidates"]:
            if other != match["column"]:
                renames.setdefault(other, f"Ignored_{other}")
    for column in table.columns:
        if column in renames:
            continue
        # A column that already has a standard name for a different field
        # would collide with a renamed one.
        if column in renames.values():
            renames[column] = f"Ignored_{column}"
    return table.rename(columns=renames)


def value_warnings(table, matches, frequency_scale=None):
    """Check units and ranges that headers alone cannot settle.

    With `frequency_scale` "percent" (the user's Phase 4A choice) the 0-1
    frequency check is replaced by the 0-100 check in genomics.import_options.
    """

    problems = []
    columns = {match["field"]: match["column"] for match in matches if match["status"] == "matched"}

    def numbers(name):
        raw = table[columns[name]].astype(str).str.strip()
        return pd.to_numeric(raw.replace("", pd.NA), errors="coerce").dropna()

    if "Genome_Coverage_Percent" in columns:
        values = numbers("Genome_Coverage_Percent")
        if len(values) and values.max() <= 1:
            problems.append(f"Column '{columns['Genome_Coverage_Percent']}' looks like a fraction "
                            "(all values are 1 or less), but genome coverage must be a percentage "
                            "from 0 to 100. Convert it before uploading; it is not rescaled "
                            "automatically.")
    if "Alt_Frequency" in columns and frequency_scale != "percent":
        values = numbers("Alt_Frequency")
        if len(values) and values.max() > 1:
            problems.append(f"Column '{columns['Alt_Frequency']}' has values above 1 (largest {values.max():g}), "
                            "so it may be a percentage; alternative-allele frequencies are read as 0 to 1. If the "
                            "source documents percentages (0-100), choose 'Percentages' under 'How to read this "
                            "file'; values are never rescaled without that choice.")
    for name, label in [("Position", "variant positions"), ("Alt_Frequency", "allele frequencies"),
                        ("Read_Depth", "read depths")]:
        if name in columns:
            raw = table[columns[name]].astype(str).str.strip()
            bad = raw.ne("") & pd.to_numeric(raw, errors="coerce").isna()
            if bad.any():
                problems.append(f"Column '{columns[name]}' should hold {label}, but {int(bad.sum())} "
                                f"value(s) are not numbers (for example '{raw[bad].iloc[0][:30]}').")
    return problems
