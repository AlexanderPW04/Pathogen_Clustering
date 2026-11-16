"""Analysis-specific metadata readiness for the dataset review step.

Which metadata a dataset needs depends on the analysis the visitor chose:
individual sequences, pooled sequencing results, or both combined. Each field
is classified as

- Required: the analysis itself cannot run without it;
- Recommended: one named analysis step is limited or unavailable without it;
- Optional: context or extra precision; nothing stops without it.

This module only reads the review (detected columns and value counts). It
never fills values, invents coordinates, or treats a missing value as zero;
the warnings name the specific calculation that a gap affects.
"""

from validate_inputs import is_missing_value

REQUIRED, RECOMMENDED, OPTIONAL = "Required", "Recommended", "Optional"
DETECTED, PARTLY, IMPRECISE, MISSING, CHOOSE = (
    "Detected", "Partly filled", "Imprecise", "Missing", "Needs your choice")
FROM_FASTA, STATED, CONFIRMED = "From FASTA descriptions", "Stated by you", "Complete genomes confirmed"
# Badge colors on the review page: good, limited, or missing.
STATUS_TONE = {DETECTED: "run", FROM_FASTA: "run", STATED: "run", CONFIRMED: "run",
               "From FASTA record IDs": "run", PARTLY: "limited", IMPRECISE: "limited",
               MISSING: "skipped", CHOOSE: "needs"}

# Fields that linked publications can sometimes fill or make more specific
# (Sprint 4). Identifiers, sequences, and variant measurements never come
# from publications.
LITERATURE_HELP = {
    "Collection_Date": "Publications may state collection dates or date ranges.",
    "Location_Name": "Publications may name a more specific collection site.",
    "Coordinates": "Publications sometimes report sampling-site coordinates.",
    "Host_Source": "Publications may describe the sample type.",
}
NO_LITERATURE_HELP = "No: this comes from the sequencing data or your files, not from publications."
BROADER_LOCATION_HELP = ("Not searched separately; a more specific collection location found in a "
                         "publication covers it.")

MODE_LABELS = {"individual": "Individual sequences", "pooled": "Pooled sequencing results",
               "combined": "Combined individual and pooled analysis"}


def spec(key, label, level, used_for):
    return {"key": key, "label": label, "level": level, "used_for": used_for}


INDIVIDUAL_FIELDS = [
    spec("Accession", "Sample identifier", REQUIRED, "Links each metadata row to its FASTA record"),
    spec("Pathogen", "Pathogen or organism", REQUIRED,
         "Sequences are aligned and compared only within one pathogen"),
    spec("Gene_Segment", "Genomic region (segment or gene)", REQUIRED,
         "Sequence comparison keeps segments or genes apart; not needed once you confirm complete genomes"),
    spec("Collection_Date", "Collection date", RECOMMENDED, "Temporal screening (days apart)"),
    spec("Location_Name", "Collection location", RECOMMENDED,
         "Geographic-distance screening (specific places are located with GeoNames)"),
    spec("Coordinates", "Latitude and longitude", OPTIONAL, "Kilometer distances without a GeoNames lookup"),
    spec("Country", "Country", OPTIONAL, "Kept separately; fills a missing location at country level"),
    spec("State_Province", "State or province", OPTIONAL,
         "Kept separately; fills a missing location at state level"),
    spec("Host_Source", "Host or sample source", OPTIONAL, "Context only; no calculation uses it"),
]
VARIANT_FIELDS = [("Sample_ID", "Pooled sample identifier"), ("Reference_Genome", "Reference genome label"),
                  ("Position", "Reference position"), ("Ref", "Reference base"), ("Alt", "Alternative base"),
                  ("Alt_Frequency", "Variant frequency")]
# Phase 4A: depth is recommended, not required. Without it reported frequencies are
# still summarized; depth QC, count-based VAF, and uncertainty are unavailable.
POOLED_FIELDS = [spec(key, label, REQUIRED, "Pooled variant QC and frequencies")
                 for key, label in VARIANT_FIELDS] + [
    spec("Read_Depth", "Read depth", RECOMMENDED,
         "Depth QC, count-based frequencies, and uncertainty intervals"),
    spec("Collection_Date", "Collection date", RECOMMENDED, "Temporal frequency trends"),
    spec("Location_Name", "Sampling site", RECOMMENDED, "Site grouping and geographic context"),
    spec("Genome_Coverage_Percent", "Genome coverage (%)", RECOMMENDED,
         "The study's >90% genome coverage at 10x filter"),
    spec("Coordinates", "Site latitude and longitude", OPTIONAL, "Map of sampling sites"),
    spec("Mean_Depth", "Mean depth", OPTIONAL, "Depth-versus-SNV regression"),
    spec("Alt_Depth", "Alternate read counts", OPTIONAL, "Shown with each variant call"),
    spec("Base_Counts", "Per-position A/C/G/T counts (separate table)", OPTIONAL,
         "The study's PCoA and coverage-aware non-detection"),
]
# Which modes each analysis of the review's plan belongs to; the others are
# pooled analyses, except the profile overlap, which belongs to every mode.
ANALYSIS_SCOPE = {
    "Sequence comparison (MAFFT, p-distance)": {"individual", "combined"},
    "Temporal screening of sequence pairs": {"individual", "combined"},
    "Geographic-distance screening": {"individual", "combined"},
    "Reference-coordinate variant profiles": {"individual", "combined"},
    "Cross-source comparison": {"combined"},
    "Documented-profile overlap": {"individual", "pooled", "combined"},
}
POOLED_SCOPE = {"pooled", "combined"}
# Required fields whose absence stops one analysis rather than the whole run.
ANALYSIS_SPECIFIC = {"Gene_Segment": "sequence comparison", "Reference_FASTA": "cross-source comparison"}
REFERENCE_FIELD = spec("Reference_FASTA", "Shared reference genome (FASTA)", REQUIRED,
                       "Connects individual and pooled variants on one coordinate system")


def analysis_mode(review, intended):
    """The chosen mode, or the one the assigned files support when the visitor chose 'auto'."""

    if intended in MODE_LABELS:
        return intended
    roles = review["by_role"]
    if "fasta" in roles and "snv" in roles:
        return "combined"
    if "snv" in roles:
        return "pooled"
    if "fasta" in roles:
        return "individual"
    return None


def first(review, role):
    entries = review["by_role"].get(role) or [None]
    return entries[0]


def match_for(entry, field):
    for match in (entry or {}).get("matches", []):
        if match["field"] == field:
            return match
    return None


def matched_column(entry, field):
    match = match_for(entry, field)
    return match["column"] if match and match["status"] == "matched" else None


def filled(entry, field):
    """Rows with a real value in the column matched to `field` (None when there is no column)."""

    column = matched_column(entry, field)
    if column is None or entry.get("table") is None:
        return None
    return int((~entry["table"][column].map(is_missing_value)).sum())


def count_status(count, total):
    if not isinstance(count, int) or count <= 0:
        return MISSING
    return DETECTED if count >= total else PARTLY


def row(item, status, detail, source):
    """One table row; `source` names the files or columns it was read from."""

    help_text = ""
    if status not in {DETECTED, FROM_FASTA, STATED, CONFIRMED, "From FASTA record IDs"}:
        help_text = LITERATURE_HELP.get(item["key"], BROADER_LOCATION_HELP if item["key"] in {
            "Country", "State_Province"} else NO_LITERATURE_HELP)
    return {"Source": source, "Level": item["level"], "Field": item["label"], "Key": item["key"],
            "Status": status, "Tone": STATUS_TONE.get(status, "limited"), "Detail": detail,
            "Used_For": item["used_for"], "Sprint4": help_text,
            "Sprint4_Possible": item["key"] in LITERATURE_HELP and bool(help_text)}


def column_detail(column, count, total, unit):
    if column is None:
        return "No column was recognized"
    if count is None:
        return f"Column '{column}'"
    return f"Column '{column}': {count} of {total} {unit}"


def individual_rows(review, options):
    source = "Individual sequences"
    total = len(review["individual_ids"])
    entry = first(review, "metadata")
    counts = next((item for item in review["completeness"] if item["Source"] == source), {})
    annotations = review.get("individual_annotations", {})
    annotated = {field: sum(1 for values in annotations.values() if field in values)
                 for field in ("Pathogen", "Collection_Date", "Location_Name")}
    location = review.get("location_summary", {})
    region = review.get("individual_region", {})
    rows = []
    for item in INDIVIDUAL_FIELDS:
        key = item["key"]
        match = match_for(entry, key)
        column = matched_column(entry, key)
        count = filled(entry, key)
        if match and match["status"] == "ambiguous":
            rows.append(row(item, CHOOSE, match["evidence"], source))
            continue
        if key == "Accession":
            if entry is None:
                status, detail = "From FASTA record IDs", "No metadata file: the FASTA record IDs are used"
            else:
                status, detail = (DETECTED if column else MISSING), column_detail(column, count, total, "rows")
        elif key == "Pathogen":
            declared = str(options.get("declared_pathogen", "")).strip()
            if column:
                status, detail = count_status(count, total), column_detail(column, count, total, "rows")
            elif total and annotated["Pathogen"] == total:
                status, detail = FROM_FASTA, f"Annotated for all {total} sequences"
            elif declared:
                status, detail = STATED, f"'{declared}' applies to every individual sequence"
            else:
                status, detail = MISSING, ("Enter it under 'Pathogen of all individual sequences' below; "
                                           "sequence comparison cannot run without it")
        elif key == "Gene_Segment":
            headers = region.get("from_headers", 0)
            if region.get("column"):
                status, detail = DETECTED, f"Column '{region['column']}'"
            elif total and headers >= total:
                status, detail = FROM_FASTA, f"Annotated for all {total} sequences"
            elif options.get("complete_genomes"):
                status, detail = CONFIRMED, "Every sequence is compared as a complete genome"
            else:
                status = PARTLY if headers else MISSING
                detail = (f"Annotated for {headers} of {total} sequences; add a segment or gene column, "
                          "or confirm complete genomes in Advanced settings")
        elif key == "Collection_Date":
            dates = counts.get("With_Collection_Date", 0)
            status = count_status(dates, total)
            detail = (column_detail(column, dates, total, "sequences have a date") if column else
                      f"{dates} of {total} sequences annotated in FASTA descriptions" if annotated[key] else
                      "No column was recognized")
        elif key == "Location_Name":
            located = counts.get("With_Location", 0)
            coarse = location.get("coarse", 0)
            status = count_status(located, total)
            where = f"Column '{column}'" if column else (
                "Built from the country/state columns" if location else
                "FASTA descriptions" if annotated[key] else "No column was recognized")
            detail = f"{where}: {located} of {total} sequences have a location" if located else where
            if located and coarse:
                status = IMPRECISE
                example = location.get("example", "")
                detail += (f"; {coarse} only at {location_level(location)} level"
                           + (f" (for example '{example}')" if example else ""))
        elif key == "Coordinates":
            pairs = counts.get("With_Coordinates", 0)
            status = count_status(pairs, total)
            detail = f"{pairs} of {total} sequences have both latitude and longitude"
        else:
            status = count_status(count, total) if column else MISSING
            detail = column_detail(column, count, total, "rows")
        rows.append(row(item, status, detail, source))
    return rows


def location_level(location):
    if location.get("country", 0) and not location.get("state", 0):
        return "country"
    if location.get("state", 0) and not location.get("country", 0):
        return "state"
    return "state or country"


def pooled_rows(review):
    source = "Pooled samples"
    entries = review["by_role"].get("snv", [])
    metadata = first(review, "pooled_metadata")
    counts = next((item for item in review["completeness"] if item["Source"] == source), {})
    total = len(review["pooled_ids"])
    rows = []
    for item in POOLED_FIELDS:
        key = item["key"]
        if item["level"] == REQUIRED or key == "Read_Depth":
            status, detail = variant_field(entries, key)
        elif key in {"Collection_Date", "Location_Name"}:
            count = counts.get("With_Collection_Date" if key == "Collection_Date" else "With_Location", 0)
            status = count_status(count, total)
            detail = f"{count} of {total} pooled samples"
        elif key == "Genome_Coverage_Percent":
            coverage = counts.get("With_Genome_Coverage", 0)
            if isinstance(coverage, str):
                status = DETECTED if coverage.startswith("From") else MISSING
                detail = coverage
            else:
                status, detail = count_status(coverage, total), f"{coverage} of {total} pooled samples"
        elif key == "Coordinates":
            pairs = counts.get("With_Coordinates", 0)
            status = count_status(pairs, total) if isinstance(pairs, int) else MISSING
            detail = f"{pairs} of {total} pooled samples" if isinstance(pairs, int) else pairs
        elif key == "Mean_Depth":
            column = matched_column(metadata, key)
            status = count_status(filled(metadata, key), total) if column else MISSING
            detail = column_detail(column, filled(metadata, key), total, "rows")
        elif key == "Alt_Depth":
            columns = [matched_column(entry, key) for entry in entries if matched_column(entry, key)]
            status = DETECTED if columns else MISSING
            detail = f"Column '{columns[0]}'" if columns else "No column was recognized"
        else:  # Base_Counts
            tables = [entry for entry in review["by_role"].get("positions", [])
                      if all(matched_column(entry, f"{base}_Count") for base in "ACGT")]
            status = DETECTED if tables else MISSING
            detail = ("Per-position table: " + ", ".join(entry["info"]["name"] for entry in tables)
                      if tables else "No per-position table with A/C/G/T counts")
        rows.append(row(item, status, detail, source))
    return rows


def variant_field(entries, key):
    """Status of one required variant-call field across the pooled variant files."""

    if not entries:
        return MISSING, "No pooled variant file"
    found, absent = [], []
    for entry in entries:
        name = entry["info"]["name"]
        if entry.get("table") is None:
            found.append(f"read from the {entry['info']['format']} format")
            continue
        match = match_for(entry, key)
        if key == "Reference_Genome" and entry.get("declared_reference") and not (
                match and match["status"] == "matched"):
            found.append(f"stated by you: '{entry['declared_reference']}'")
        elif match and match["status"] == "matched":
            found.append(f"column '{match['column']}'")
        elif match and match["status"] == "ambiguous":
            return CHOOSE, f"{name}: {match['evidence']}"
        else:
            absent.append(name)
    if absent:
        return MISSING, "Not found in " + ", ".join(absent)
    return DETECTED, "; ".join(dict.fromkeys(found)).capitalize()


def reference_row(review):
    reference = review["reference"]
    if reference:
        return row(REFERENCE_FIELD, DETECTED, f"{reference['id']} ({reference['length']} bases)", "Both sources")
    return row(REFERENCE_FIELD, MISSING, "No reference FASTA was assigned", "Both sources")


def kept_columns(review):
    """Columns of the metadata tables that no field uses; they are carried along unchanged."""

    names = []
    region = review.get("individual_region", {}).get("column")
    for role in ("metadata", "pooled_metadata"):
        entry = first(review, role)
        if entry is None or entry.get("table") is None:
            continue
        used = {match["column"] for match in entry["matches"] if match["status"] == "matched"} | {region}
        names += [str(column) for column in entry["table"].columns if column not in used]
    return names


def consequences(rows, review):
    """Specific sentences: which calculation each gap affects, and what still runs."""

    by_key = {(item["Source"], item["Key"]): item for item in rows}
    messages = []
    individual = len(review["individual_ids"])
    location = review.get("location_summary", {})

    def status(source, key):
        item = by_key.get((source, key))
        return item["Status"] if item else None

    source = "Individual sequences"
    counts = next((item for item in review["completeness"] if item["Source"] == source), {})
    located = counts.get("With_Location", 0)
    if status(source, "Gene_Segment") in {MISSING, PARTLY}:
        messages.append("The genomic region (segment or gene) is not known for every sequence. Sequence "
                        "comparison will stop at the alignment step unless you add a segment or gene column, "
                        "or confirm complete genomes under Advanced settings before running.")
    if status(source, "Location_Name") == MISSING:
        messages.append("No collection locations were found for the individual sequences. Genetic and date "
                        "comparisons can proceed, but geographic-distance screening cannot be performed.")
    elif status(source, "Location_Name") == IMPRECISE:
        level = location_level(location)
        coarse = location.get("coarse", 0)
        if coarse >= located and by_key[(source, "Coordinates")]["Status"] == MISSING:
            messages.append(f"Collection locations are only available at the {level} level. Genetic "
                            "comparisons can proceed, but precise geographic-distance screening cannot be "
                            "performed; such pairs are reported as unresolved geography.")
        else:
            messages.append(f"{coarse} of {individual} sequences have only {level}-level locations. Genetic "
                            "comparisons include them, but pairs involving them cannot pass the kilometer "
                            "cutoff and are reported as unresolved geography.")
    if status(source, "Location_Name") in {PARTLY, IMPRECISE} and located < individual:
        messages.append(f"{individual - located} of {individual} sequences have no collection location; "
                        "pairs involving them cannot be screened by distance.")
    if status(source, "Collection_Date") == MISSING:
        messages.append("No collection dates were found for the individual sequences. Genetic comparisons "
                        "can proceed, but temporal screening (days apart) cannot be performed, so no pair "
                        "can meet all three screening criteria.")
    elif status(source, "Collection_Date") == PARTLY:
        dates = counts.get("With_Collection_Date", 0)
        messages.append(f"{individual - dates} of {individual} sequences have no collection date; pairs "
                        "involving them get no days-apart value and cannot meet the temporal criterion.")

    source = "Pooled samples"
    if status(source, "Collection_Date") == MISSING:
        messages.append("No collection dates were found for the pooled samples. Variant QC and frequencies "
                        "can be analyzed, but temporal frequency trends cannot be shown.")
    if status(source, "Location_Name") == MISSING:
        messages.append("No sampling sites were found for the pooled samples. Variant frequencies can be "
                        "analyzed, but site grouping and geographic context cannot be shown.")
    if status(source, "Genome_Coverage_Percent") == MISSING:
        messages.append("Genome coverage was not supplied for the pooled samples. Variant frequencies are "
                        "analyzed, but the study's >90% coverage filter cannot be applied; those samples "
                        "are reported as coverage not evaluable.")
    if status(source, "Base_Counts") == MISSING:
        messages.append("Without per-position A/C/G/T counts, the study's PCoA cannot be reproduced (an "
                        "exploratory SNV-set ordination is shown instead), and unreported variants stay "
                        "'coverage unknown', never 0%.")
    if status("Both sources", "Reference_FASTA") == MISSING:
        messages.append("No shared reference genome was assigned. Individual and pooled analyses can each "
                        "proceed, but cross-source variant comparison cannot be performed.")
    return messages


def metadata_readiness(review, options=None):
    """Field-by-field readiness for the chosen (or detected) analysis; None without data files."""

    options = options or {}
    mode = analysis_mode(review, options.get("intended_mode"))
    if mode is None:
        return None
    rows = []
    if mode in {"individual", "combined"}:
        rows += individual_rows(review, options)
    if mode in {"pooled", "combined"}:
        rows += pooled_rows(review)
    if mode == "combined":
        rows.append(reference_row(review))
    # Missing required identifiers, pathogen, or variant columns (and any
    # column choice) are review problems, so the run waits for them; a missing
    # region or reference stops only the analysis named in ANALYSIS_SPECIFIC.
    blocking = [f"{item['Field']} ({item['Source'].lower()})" for item in rows
                if item["Status"] == CHOOSE or (item["Level"] == REQUIRED and item["Status"] == MISSING
                                                and item["Key"] not in ANALYSIS_SPECIFIC)]
    partial = [f"{item['Field']}: without it, {ANALYSIS_SPECIFIC[item['Key']]} cannot run" for item in rows
               if item["Key"] in ANALYSIS_SPECIFIC and item["Status"] in {MISSING, PARTLY}]
    analyses = [plan for plan in review["analyses"]
                if mode in ANALYSIS_SCOPE.get(plan["Analysis"], POOLED_SCOPE) or plan["Status"] != "Unavailable"]
    enrichable = list(dict.fromkeys(item["Field"].lower() for item in rows if item["Sprint4_Possible"]))
    return {"mode": mode, "mode_label": MODE_LABELS[mode], "rows": rows, "blocking": blocking,
            "partial": partial, "analyses": analyses,
            "consequences": consequences(rows, review), "kept_columns": kept_columns(review),
            "enrichable": enrichable,
            "counts": {level: sum(1 for item in rows if item["Level"] == level and item["Status"] in {
                MISSING, CHOOSE}) for level in (REQUIRED, RECOMMENDED, OPTIONAL)}}
