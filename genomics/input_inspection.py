"""Inspect uploaded files, suggest their roles, and review what can be analyzed.

The guided upload works in three steps:

1. `inspect_file()` reads each uploaded file (already decompressed) and
   records what it contains: format, columns, identifiers, references.
2. `suggest_roles()` proposes a role for each file from its content and from
   evidence shared between files (for example metadata IDs that match FASTA
   record IDs). A file is never given a role without stated evidence; when the
   evidence is not decisive the user must choose.
3. `review_dataset()` applies the chosen roles and column mappings and
   reports problems, reference compatibility, metadata completeness, and which
   analyses can run. `prepare_inputs()` then writes standardized copies that
   the existing loaders read unchanged.

Nothing here fills in a missing date, location, reference, frequency, or depth.
"""

import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from genomic_regions import fasta_region
from genomics.column_detection import (
    apply_mapping,
    detect_columns,
    resolve_choices,
    squash,
    unresolved,
    value_warnings,
)
from genomics.location_hierarchy import location_hierarchy
from geocode_locations import is_state_only_location
from genomics.file_formats import (
    ALIGNMENT_KIND,
    EXCEL_KIND,
    FASTA_KIND,
    FASTQ_KIND,
    RAW_READ_MESSAGE,
    TABLE_KIND,
    VCF_KIND,
    detect_kind,
    file_stem,
    read_table,
    split_name,
    vcf_notes,
    vcf_to_snv_table,
)
from genomics.sequence_variants import MINIMUM_CALLABLE_POSITIONS, MINIMUM_MAPPING_IDENTITY
from validate_inputs import ALLOWED_SEQUENCE_CHARACTERS, clean_sequence, normalize_accession
from wastewater.processed_snv import (
    IVAR_REQUIRED,
    SAMPLE_HEADERS,
    PAPER_COLUMNS,
    load_processed_snvs,
    read_table as read_snv_table,
)


ROLE_LABELS = {
    "fasta": "Individual sequences (FASTA)",
    "metadata": "Individual-sequence metadata",
    "reference": "Reference genome (one FASTA sequence)",
    "snv": "Pooled variant calls (SNV table, iVar, or VCF)",
    "pooled_metadata": "Pooled-sample metadata",
    "positions": "Per-position depth or A/C/G/T counts",
    "masked_sites": "Masked problematic sites (quality control)",
    "signatures": "Documented variant profiles (annotation)",
    "ignore": "Do not use this file",
}
ROLE_ORDER = ["fasta", "metadata", "reference", "snv", "pooled_metadata", "positions",
              "masked_sites", "signatures"]
MULTI_FILE_ROLES = {"fasta", "snv"}  # several files are combined into one input
FASTA_ROLES = {"fasta", "reference"}
TABLE_ROLES = {"metadata", "pooled_metadata", "positions", "masked_sites", "signatures", "snv"}

PAPER_TABLE_FORMAT = "Fontenele et al. (2021) Supplementary Table 2"
IVAR_TABLE_FORMAT = "iVar variants table"
STANDARD_TABLE_FORMAT = "Table with recognizable columns"

PAPER_REFERENCE = "NC_045512.2"
ANNOTATION_PATTERN = re.compile(r"\[([A-Za-z_ -]+)=([^\]]+)\]")
# Explicit FASTA description annotations such as [collection_date=2020-06-01].
ANNOTATION_FIELDS = {
    "collectiondate": "Collection_Date", "date": "Collection_Date",
    "country": "Location_Name", "location": "Location_Name", "geolocname": "Location_Name",
    "host": "Host_Source", "organism": "Pathogen",
}


def preview(values, limit=12):
    values = [str(value) for value in values]
    shown = ", ".join(values[:limit])
    return shown + (f", ... ({len(values)} in total)" if len(values) > limit else "")


def fasta_annotations(description):
    """Return {standard field: value} from explicit [key=value] annotations."""

    found = {}
    for key, value in ANNOTATION_PATTERN.findall(str(description)):
        name = ANNOTATION_FIELDS.get(squash(key))
        if name and value.strip():
            found[name] = value.strip()
    return found


def empty_info(name):
    extension, compressed = split_name(name)
    return {
        "name": name, "kind": None, "format": None, "extension": extension,
        "compressed": compressed, "records": 0, "columns": [], "ids": [],
        "reference_ids": [], "lengths": [], "notes": [], "problems": [], "refused": False,
        "candidates": [], "suggested_role": None, "role_evidence": "", "role_status": "",
        "auto_choices": {}, "annotations": {}, "descriptions": {}, "id_columns": {},
        "header_hints": {},
    }


def inspect_fasta(path, info):
    records = list(SeqIO.parse(str(path), "fasta"))
    info["records"] = len(records)
    if not records:
        info["problems"].append("The FASTA file contains no sequence records.")
        return
    info["ids"] = [str(record.id).strip() for record in records]
    info["lengths"] = [len(clean_sequence(record.seq)) for record in records]
    info["descriptions"] = {str(record.id).strip(): str(record.description) for record in records}
    invalid = set()
    for record in records:
        invalid |= set(clean_sequence(record.seq)) - ALLOWED_SEQUENCE_CHARACTERS
    if invalid:
        info["problems"].append("Sequence characters that are not nucleotides: "
                                + "".join(sorted(invalid)) + ". FASTQ or protein sequences "
                                "are not accepted here.")
    for record in records:
        annotations = fasta_annotations(record.description)
        if annotations:
            info["annotations"][str(record.id).strip()] = annotations
    if info["annotations"]:
        fields = sorted({name for values in info["annotations"].values() for name in values})
        info["notes"].append("FASTA descriptions contain explicit annotations for "
                             + ", ".join(fields) + ".")
    info["format"] = f"FASTA ({len(records)} record{'s' if len(records) != 1 else ''})"
    info["candidates"] = ["fasta"] if len(records) > 1 else ["reference", "fasta"]


def inspect_vcf(path, info):
    table, details = vcf_to_snv_table(path, file_stem(info["name"]))
    info["records"] = details["records"]
    info["ids"] = details["samples"]
    info["reference_ids"] = sorted(set(table["Reference_Genome"])) or sorted(
        name for name in details["contigs"] if name)
    info["contig_lengths"] = details["contigs"]
    info["notes"] += vcf_notes(details)
    info["notes"].append(f"{details['snv_rows']} single-nucleotide allele entries can be analyzed.")
    if details["snv_rows"] == 0:
        info["problems"].append("No single-nucleotide variant entry has an interpretable frequency.")
    info["format"] = "VCF " + details["meta"].get("fileformat", "").replace("VCF", "").strip("v ")
    info["candidates"] = ["snv"]


def column_values(table, column):
    values = table[column].astype(str).str.strip()
    return sorted(set(values[values.ne("")]))


def ivar_sample_column(columns):
    """The sample column of an iVar table, whatever its capitalization, or None."""

    for column in columns:
        if column in {"SAMPLE", "Sample_ID"} or squash(column) in SAMPLE_HEADERS:
            return column
    return None


def inspect_table(path, info):
    if info["kind"] == EXCEL_KIND:
        table = read_snv_table(path)  # handles the paper's descriptive first row
    else:
        table = read_table(path, info["name"])
    table.columns = [str(name).strip() for name in table.columns]
    info["records"] = len(table)
    info["columns"] = list(table.columns)
    columns = set(table.columns)

    if (table.attrs.get("fontenele_2021_table2") and set(PAPER_COLUMNS).issubset(columns)
            and {"STATE", "LOCATION"}.issubset(columns)):
        info["format"] = PAPER_TABLE_FORMAT
        info["candidates"] = ["snv"]
        info["ids"] = column_values(table, "SAMPLE")
        info["reference_ids"] = [PAPER_REFERENCE]
        info["notes"].append("Recognized the published supplementary SNV table; its reference is "
                             f"{PAPER_REFERENCE} as stated in the paper's methods.")
        return
    if IVAR_REQUIRED.issubset(columns):
        info["format"] = IVAR_TABLE_FORMAT
        info["candidates"] = ["snv"]
        info["reference_ids"] = column_values(table, "REGION")
        sample_column = ivar_sample_column(table.columns)
        if sample_column:
            info["ids"] = column_values(table, sample_column)
        else:
            info["ids"] = [file_stem(info["name"])]
            info["notes"].append("This iVar table has no SAMPLE column, so its sample ID is taken "
                                 "from the file name: " + info["ids"][0] + ".")
        return

    info["format"] = STANDARD_TABLE_FORMAT
    matches = {role: detect_columns(table.columns, role) for role in TABLE_ROLES}

    def present(role, name):
        return any(match["field"] == name and match["status"] != "missing" for match in matches[role])

    def complete(role):
        return all(match["status"] != "missing" for match in matches[role] if match["required"])

    has_alleles = present("snv", "Ref") and present("snv", "Alt")
    candidates = []
    if complete("signatures"):
        candidates.append("signatures")
    elif complete("snv"):
        candidates.append("snv")
    if not has_alleles and complete("positions") and (
            present("positions", "Depth")
            or all(present("positions", f"{base}_Count") for base in "ACGT")):
        candidates.append("positions")
    if (not has_alleles and complete("masked_sites") and not present("positions", "Sample_ID")
            and not present("positions", "Depth")):
        candidates.append("masked_sites")
    if not present("snv", "Position"):
        hints = {
            "metadata": [name for name in ("Pathogen", "Host_Source") if present("metadata", name)],
            "pooled_metadata": [name for name in ("Genome_Coverage_Percent", "Mean_Depth")
                                if present("pooled_metadata", name)],
        }
        info["header_hints"] = hints
        for role in ("metadata", "pooled_metadata"):
            if complete(role):
                candidates.append(role)
                id_field = "Accession" if role == "metadata" else "Sample_ID"
                id_match = next(match for match in matches[role] if match["field"] == id_field)
                info["id_columns"][role] = {
                    column: column_values(table, column) for column in id_match["candidates"]}
    info["candidates"] = candidates
    for role in candidates:
        if role in {"snv", "positions"}:
            reference = next((match for match in matches[role]
                              if match["field"] in {"Reference_Genome", "Reference_ID"}), None)
            if reference and reference["status"] == "matched":
                info["reference_ids"] = column_values(table, reference["column"])
            sample = next((match for match in matches[role] if match["field"] == "Sample_ID"), None)
            if sample and sample["status"] == "matched" and role == "snv":
                info["ids"] = column_values(table, sample["column"])
    if not candidates:
        info["problems"].append("The columns do not match any supported table (see the help "
                                "list of accepted column names).")


def inspect_file(path, name):
    """Describe one decompressed upload."""

    info = empty_info(name)
    info["kind"] = detect_kind(path, name)
    info["size_kb"] = round(Path(path).stat().st_size / 1024, 1)
    try:
        if info["kind"] == FASTA_KIND:
            inspect_fasta(path, info)
        elif info["kind"] == VCF_KIND:
            inspect_vcf(path, info)
        elif info["kind"] in {TABLE_KIND, EXCEL_KIND}:
            inspect_table(path, info)
        elif info["kind"] in {FASTQ_KIND, ALIGNMENT_KIND}:
            info["refused"] = True
            info["problems"].append(RAW_READ_MESSAGE)
        else:
            info["refused"] = True
            info["problems"].append("This file type is not supported. Use FASTA, CSV/TSV/XLSX, "
                                    "or VCF (each optionally .gz).")
    except Exception as error:  # any unreadable file is refused on its own, not the whole upload
        info["refused"] = True
        info["problems"].append(f"Could not read this file: {error}")
    if info["kind"] == FASTA_KIND and info["extension"] not in {".fasta", ".fa", ".fna"}:
        info["notes"].append("The content is FASTA although the file name ends with "
                             f"'{info['extension'] or 'no extension'}'.")
    return info


def normalized_ids(values):
    return {normalize_accession(value) for value in values if str(value).strip()}


def suggest_roles(inspections):
    """Set suggested_role / role_evidence / role_status on every inspection."""

    fasta_ids = set()
    variant_ids = set()
    variant_references = {}
    for info in inspections:
        if info["candidates"] == ["fasta"]:
            fasta_ids |= normalized_ids(info["ids"])
        if info["candidates"] == ["snv"] or "positions" in info["candidates"]:
            for reference in info["reference_ids"]:
                variant_references.setdefault(reference, info["name"])
        if info["candidates"] == ["snv"]:
            variant_ids |= {str(value) for value in info["ids"]}
    metadata_ids = set()
    for info in inspections:
        for columns in info["id_columns"].values():
            for values in columns.values():
                metadata_ids |= normalized_ids(values)

    for info in inspections:
        candidates = info["candidates"]
        info["role_status"] = "needs choice"
        if info["refused"] or not candidates:
            info["suggested_role"] = "ignore"
            info["role_status"] = "refused" if info["refused"] else "needs choice"
            info["role_evidence"] = "Not a supported input." if info["refused"] else (
                "No supported role matched; choose one or leave it unused.")
            continue
        if len(candidates) == 1:
            info["suggested_role"] = candidates[0]
            info["role_status"] = "clear"
            info["role_evidence"] = role_reason(info, candidates[0])
            continue
        if set(candidates) == {"reference", "fasta"}:
            record_id = info["ids"][0]
            description = info["descriptions"].get(record_id, "")
            if record_id in variant_references:
                choose(info, "reference", f"its ID {record_id} is the reference named in "
                                          f"{variant_references[record_id]}")
            elif "reference" in description.casefold():
                choose(info, "reference", "its FASTA description says 'reference'")
            elif normalize_accession(record_id) in metadata_ids:
                choose(info, "fasta", f"its ID {record_id} appears in a metadata table")
            else:
                info["suggested_role"] = None
                info["role_evidence"] = ("One sequence: it could be the reference genome or a "
                                         "single individual sequence. Please choose.")
            continue
        if set(candidates) <= {"metadata", "pooled_metadata"}:
            choose_metadata_role(info, fasta_ids, variant_ids)
            continue
        info["suggested_role"] = None
        info["role_evidence"] = "Several roles are possible: " + ", ".join(
            ROLE_LABELS[role] for role in candidates) + ". Please choose."
    return inspections


def apply_declared_role(info, role):
    """Keep the role a person chose with a role-specific upload field.

    The choice is kept when the file's content fits it. When the content
    rules it out, no role is guessed: the conflict is explained and the
    person must choose in the review step.
    """

    info["declared_role"] = role
    label = ROLE_LABELS[role]
    if info["refused"]:
        info["role_evidence"] = f"Uploaded as '{label}', but this file type is not accepted."
        return info
    if role in info["candidates"] or (role in FASTA_ROLES and info["kind"] == FASTA_KIND):
        info["suggested_role"] = role
        info["role_status"] = "declared"
        info["role_evidence"] = f"Uploaded as '{label}'; the content fits this role."
        return info
    detected = ", ".join(ROLE_LABELS[name] for name in info["candidates"]) or "none of the supported roles"
    info["suggested_role"] = None
    info["role_status"] = "conflict"
    info["role_evidence"] = (f"Uploaded as '{label}', but the content looks like: {detected}. "
                             "Choose the correct role.")
    return info


def choose(info, role, evidence):
    info["suggested_role"] = role
    info["role_status"] = "clear"
    info["role_evidence"] = f"{ROLE_LABELS[role]}: {evidence}."


def role_reason(info, role):
    if info["format"] in {PAPER_TABLE_FORMAT, IVAR_TABLE_FORMAT} or info["kind"] == VCF_KIND:
        return f"{ROLE_LABELS[role]}: recognized {info['format']}."
    if role == "fasta":
        return f"{ROLE_LABELS[role]}: {info['records']} sequence records."
    return f"{ROLE_LABELS[role]}: its column headers match this role's required fields."


def choose_metadata_role(info, fasta_ids, variant_ids):
    """Decide between individual and pooled metadata from matching IDs, then headers."""

    best = {}
    for role, id_columns in info["id_columns"].items():
        targets = fasta_ids if role == "metadata" else {normalize_accession(v) for v in variant_ids}
        scores = {column: len(normalized_ids(values) & targets) for column, values in id_columns.items()}
        if scores and max(scores.values()) > 0:
            winners = [column for column, score in scores.items() if score == max(scores.values())]
            if len(winners) == 1:
                best[role] = (winners[0], scores[winners[0]], len(id_columns[winners[0]]))
    if len(best) == 1:
        role, (column, matched, total) = next(iter(best.items()))
        other = "FASTA records" if role == "metadata" else "pooled variant samples"
        choose(info, role, f"{matched} of {total} values in column '{column}' match the {other}")
        id_field = "Accession" if role == "metadata" else "Sample_ID"
        if len(info["id_columns"][role]) > 1:
            info["auto_choices"].setdefault(role, {})[id_field] = column
        return
    hints = {role: names for role, names in info["header_hints"].items() if names}
    if len(hints) == 1:
        role, names = next(iter(hints.items()))
        choose(info, role, "it has " + ", ".join(names) + " column(s)")
        return
    info["suggested_role"] = None
    info["role_evidence"] = ("These could be individual-sequence or pooled-sample metadata; no IDs "
                             "match another uploaded file. Please choose.")


# ----------------------------------------------------------------------------
# Dataset review
# ----------------------------------------------------------------------------

def table_for(path, info):
    if info["kind"] == EXCEL_KIND:
        table = read_snv_table(path)
    else:
        table = read_table(path, info["name"])
    table.columns = [str(name).strip() for name in table.columns]
    return table


def native_format(info):
    return info["kind"] == VCF_KIND or info["format"] in {PAPER_TABLE_FORMAT, IVAR_TABLE_FORMAT}


def review_dataset(staged, inspections, roles, choices=None, options=None):
    """Summarize the dataset under the chosen roles.

    `staged` maps file index -> decompressed path; `roles` maps index -> role;
    `choices` maps index -> {field: column} for ambiguous columns; `options`
    may hold "declared_pathogen" and "complete_genomes".
    """

    choices = choices or {}
    options = options or {}
    review = {"files": [], "problems": [], "warnings": [], "by_role": {}, "mappings": {},
              "individual_ids": [], "pooled_ids": [], "reference": None, "compatibility": [],
              "completeness": [], "analyses": [], "annotation_notes": [], "ready": False}

    for index, info in enumerate(inspections):
        role = roles.get(index, info["suggested_role"]) or None
        entry = {"index": index, "info": info, "role": role, "matches": [], "problems": []}
        review["files"].append(entry)
        if info["refused"]:
            continue
        if role is None:
            entry["problems"].append("Choose a role for this file (or 'Do not use this file').")
            continue
        if role == "ignore":
            continue
        if role not in info["candidates"] and not (role in FASTA_ROLES and info["kind"] == FASTA_KIND):
            entry["problems"].append(f"The content does not fit the role '{ROLE_LABELS[role]}'.")
        review["by_role"].setdefault(role, []).append(entry)
        if role in TABLE_ROLES and info["kind"] in {TABLE_KIND, EXCEL_KIND} and not native_format(info):
            table = table_for(staged[index], info)
            matches = resolve_choices(detect_columns(table.columns, role), choices.get(index, {}))
            matches = resolve_choices(matches, info["auto_choices"].get(role, {}),
                                      "matching its values with the other uploaded files")
            entry["matches"] = matches
            entry["problems"] += unresolved(matches) + value_warnings(table, matches)
            review["mappings"][index] = matches
            entry["table"] = table

    for role, entries in review["by_role"].items():
        if len(entries) > 1 and role not in MULTI_FILE_ROLES:
            review["problems"].append(f"{len(entries)} files were assigned to '{ROLE_LABELS[role]}'; "
                                      "only one is allowed. Set the others to 'Do not use'.")
    if "reference" in review["by_role"]:
        info = review["by_role"]["reference"][0]["info"]
        if info["kind"] != FASTA_KIND:
            review["problems"].append("The reference genome must be a FASTA file.")
        elif info["records"] != 1:
            review["problems"].append("The reference FASTA must contain exactly one sequence.")
        else:
            review["reference"] = {"id": info["ids"][0], "length": info["lengths"][0]}
            sequence_bases = set()
            for record in SeqIO.parse(str(staged[review["by_role"]["reference"][0]["index"]]), "fasta"):
                sequence_bases = set(clean_sequence(record.seq))
            if sequence_bases - set("ACGTN"):
                review["problems"].append("The reference sequence may contain only A, C, G, T, and N.")
    for entry in review["files"]:
        for problem in entry["problems"]:
            review["problems"].append(f"{entry['info']['name']}: {problem}")
        for problem in entry["info"]["problems"]:
            if not entry["info"]["refused"] and entry["role"] not in {None, "ignore"}:
                review["problems"].append(f"{entry['info']['name']}: {problem}")

    individual_review(review, options)
    pooled_review(review)
    compatibility_review(review)
    review["analyses"] = analysis_plan(review, options)
    review["ready"] = not review["problems"] and any(
        role in review["by_role"] for role in ("fasta", "snv", "reference"))
    if not review["by_role"]:
        review["problems"].append("No file has been assigned a role yet.")
    return review


def mapped_values(entry, field):
    """Values of a standard field from an assigned table, or None if unavailable."""

    for match in entry.get("matches", []):
        if match["field"] == field and match["status"] == "matched":
            values = entry["table"][match["column"]].astype(str).str.strip()
            return values
    return None


def samples_with(entry, field):
    """Number of distinct samples in a variant table with a non-blank `field`."""

    samples = mapped_values(entry, "Sample_ID")
    values = mapped_values(entry, field)
    if samples is None or values is None:
        return 0
    return int(samples[values.ne("")].nunique())


def count_present(values):
    return 0 if values is None else int(values.ne("").sum())


def location_review(review, entry):
    """Explain how the geographic columns were combined; return the number of located records."""

    columns = {match["field"]: match["column"] for match in entry.get("matches", [])
               if match["status"] == "matched"}
    if not ({"Location_Name", "Country", "State_Province"} & set(columns)):
        return None
    table = apply_mapping(entry["table"], entry["matches"])
    header = columns.get("Location_Name", "Location_Name")
    result, summary = location_hierarchy(table, location_header=header)
    total = len(result)
    broader = [f"'{columns[field]}'" for field in ("Country", "State_Province") if field in columns]
    if "Location_Name" in columns and broader:
        review["annotation_notes"].append(
            f"Geographic columns recognized: '{header}' was selected for Location_Name; "
            f"{' and '.join(broader)} were kept separately as country and state information.")
    elif broader:
        review["annotation_notes"].append(
            f"No specific location column: locations are built from {' and '.join(broader)} "
            "and labeled as country- or state-level.")
    if summary["filled"] and "Location_Name" in columns:
        review["annotation_notes"].append(
            f"{summary['filled']} record(s) without a '{header}' value use the broader country/state "
            "columns instead (labeled with that level).")
    coarse = summary["country_level"] + summary["state_level"] + int(
        (result["Location_Level"].eq("") & result["Location_Name"].map(is_state_only_location)).sum())
    if coarse:
        example = result.loc[result["Location_Level"].isin(["country", "state"]) |
                             result["Location_Name"].map(is_state_only_location), "Location_Name"].iloc[0]
        review["annotation_notes"].append(
            f"{coarse} of {total} record(s) have only country- or state-level locations (for example "
            f"'{example}'). Genetic and date comparisons continue, but these records cannot pass the "
            "kilometer cutoff; such pairs are reported as unresolved geography.")
    if summary["conflicts"]:
        ids = table["Accession"] if "Accession" in table else pd.Series(table.index, index=table.index)
        conflicts = [f"{ids[index]}: {note}" for index, note in result["Location_Conflict"].items() if note]
        review["location_conflicts"] = conflicts
        review["warnings"].append(
            f"{summary['conflicts']} record(s) have contradictory location columns; please review them "
            "(listed under 'Location values to review'). They are kept out of kilometer screening.")
    return int(result["Location_Name"].ne("").sum())


def individual_review(review, options):
    fasta_entries = review["by_role"].get("fasta", [])
    if not fasta_entries:
        return
    ids = [value for entry in fasta_entries for value in entry["info"]["ids"]]
    review["individual_ids"] = ids
    annotations = {}
    descriptions = {}
    for entry in fasta_entries:
        annotations.update(entry["info"]["annotations"])
        descriptions.update(entry["info"]["descriptions"])
    review["individual_annotations"] = annotations
    metadata = review["by_role"].get("metadata", [None])[0]
    total = len(ids)

    def from_annotations(field):
        return sum(1 for values in annotations.values() if field in values)

    if metadata is None:
        review["warnings"].append("No individual-sequence metadata: sample IDs come from the FASTA "
                                  "record IDs; dates and locations are available only where the "
                                  "FASTA descriptions annotate them explicitly.")
        pathogen_known = from_annotations("Pathogen") == total
        dates = from_annotations("Collection_Date")
        locations = from_annotations("Location_Name")
        coordinates = 0
    else:
        key = mapped_values(metadata, "Accession")
        if key is not None:
            metadata_ids = normalized_ids(key)
            fasta_ids = normalized_ids(ids)
            missing_metadata = sorted(fasta_ids - metadata_ids)
            missing_fasta = sorted(metadata_ids - fasta_ids)
            if missing_metadata:
                review["problems"].append("FASTA record(s) without a metadata row: "
                                          + preview(missing_metadata, 8))
            if missing_fasta:
                review["problems"].append("Metadata row(s) without a FASTA record: "
                                          + preview(missing_fasta, 8))
        pathogen_known = any(match["field"] == "Pathogen" and match["status"] != "missing"
                             for match in metadata.get("matches", []))
        dates = count_present(mapped_values(metadata, "Collection_Date"))
        locations = count_present(mapped_values(metadata, "Location_Name"))
        if metadata.get("table") is not None:
            located = location_review(review, metadata)
            locations = located if located is not None else locations
        latitudes = mapped_values(metadata, "Latitude")
        longitudes = mapped_values(metadata, "Longitude")
        coordinates = (0 if latitudes is None or longitudes is None
                       else int((latitudes.ne("") & longitudes.ne("")).sum()))
        for field in ("Collection_Date", "Location_Name", "Host_Source", "Pathogen"):
            column_absent = mapped_values(metadata, field) is None
            if column_absent and from_annotations(field):
                review["annotation_notes"].append(
                    f"{field} for {from_annotations(field)} sequence(s) comes from explicit "
                    "FASTA description annotations.")
                if field == "Collection_Date":
                    dates = from_annotations(field)
                if field == "Location_Name":
                    locations = from_annotations(field)
                if field == "Pathogen":
                    pathogen_known = from_annotations(field) == total
    declared = str(options.get("declared_pathogen", "")).strip()
    review["pathogen_needed"] = not pathogen_known
    if not pathogen_known and not declared:
        review["problems"].append("The pathogen of the individual sequences is not stated in "
                                  "metadata or FASTA annotations. Enter it under 'Pathogen of all "
                                  "individual sequences' (applies to every sequence).")
    elif not pathogen_known:
        review["annotation_notes"].append(f"Pathogen '{declared}' was stated by the user for every "
                                          "individual sequence.")
    region_known = metadata is not None and any(
        squash(column) in {"genesegment", "segment", "genename", "gene", "region"}
        for column in metadata["info"]["columns"])
    header_regions = sum(1 for description in descriptions.values()
                         if fasta_region(description)["label"])
    if not region_known and header_regions < total and not options.get("complete_genomes"):
        review["warnings"].append(
            f"The genomic region is annotated for {header_regions} of {total} sequences (FASTA "
            "descriptions) and no metadata region column was found. The comparison will stop "
            "unless every sequence is annotated or you confirm complete genomes in Advanced settings.")
    review["completeness"].append({
        "Source": "Individual sequences", "Samples": total,
        "With_Collection_Date": dates, "With_Location": locations,
        "With_Coordinates": coordinates, "With_Genome_Coverage": "Calculated after alignment",
    })


def pooled_review(review):
    entries = review["by_role"].get("snv", [])
    if not entries:
        return
    ids = []
    for entry in entries:
        ids += [str(value) for value in entry["info"]["ids"]]
    if len(entries) > 1 and len(ids) != len(set(ids)):
        review["problems"].append("The same pooled sample ID appears in more than one variant file.")
    if any(entry["info"]["format"] == PAPER_TABLE_FORMAT for entry in entries) and len(entries) > 1:
        review["problems"].append("The published supplementary table must be analyzed on its own.")
    review["pooled_ids"] = sorted(set(ids))
    metadata = review["by_role"].get("pooled_metadata", [None])[0]
    vcf_only = all(entry["info"]["kind"] == VCF_KIND for entry in entries)
    paper = any(entry["info"]["format"] == PAPER_TABLE_FORMAT for entry in entries)
    dates = locations = coordinates = coverage = 0
    if paper:
        dates = locations = len(review["pooled_ids"])
    else:
        dates = sum(samples_with(entry, "Collection_Date") for entry in entries)
        locations = sum(samples_with(entry, "Location_Name") for entry in entries)
    if metadata is not None:
        dates = max(dates, count_present(mapped_values(metadata, "Collection_Date")))
        locations = max(locations, count_present(mapped_values(metadata, "Location_Name")))
        latitudes = mapped_values(metadata, "Latitude")
        longitudes = mapped_values(metadata, "Longitude")
        if latitudes is not None and longitudes is not None:
            coordinates = int((latitudes.ne("") & longitudes.ne("")).sum())
        coverage = count_present(mapped_values(metadata, "Genome_Coverage_Percent"))
        key = mapped_values(metadata, "Sample_ID")
        if key is not None:
            unknown = sorted(set(review["pooled_ids"]) - set(key))
            if unknown:
                review["problems"].append("Pooled sample(s) missing from the pooled metadata: "
                                          + preview(unknown, 8))
    elif not paper:
        if vcf_only:
            review["warnings"].append("VCF files carry no collection dates or sites; add pooled-sample "
                                      "metadata for temporal and geographic analyses.")
    if review["by_role"].get("positions"):
        coverage = "From the per-position table"
    elif metadata is None:
        coverage = "Not supplied"
    review["completeness"].append({
        "Source": "Pooled samples", "Samples": len(review["pooled_ids"]),
        "With_Collection_Date": dates, "With_Location": locations,
        "With_Coordinates": coordinates if metadata is not None else "Needs pooled metadata",
        "With_Genome_Coverage": coverage,
    })


def compatibility_review(review):
    """Explain how the reference coordinates of each source will be checked."""

    reference = review["reference"]
    rows = review["compatibility"]
    if reference:
        rows.append({"Source": "Reference FASTA", "Reference": reference["id"],
                     "Status": f"Loaded: {reference['length']} bases; identity is checked by "
                               "sequence (SHA-256), not by label."})
    if review["by_role"].get("fasta"):
        rows.append({"Source": "Individual sequences", "Reference": reference["id"] if reference else "-",
                     "Status": (f"Each sequence will be aligned to {reference['id']} with MAFFT and must "
                                f"reach {MINIMUM_MAPPING_IDENTITY:.0%} identity over at least "
                                f"{MINIMUM_CALLABLE_POSITIONS} callable positions."
                                if reference else
                                "No reference FASTA: genetic distances only, no reference-coordinate "
                                "variants and no cross-source comparison.")})
    for role, label in [("snv", "Pooled variant calls"), ("positions", "Per-position table")]:
        for entry in review["by_role"].get(role, []):
            names = entry["info"]["reference_ids"]
            if not names:
                status = "Reference column not identified."
            elif not reference:
                status = ("No reference FASTA: positions are analyzed within this source only; "
                          "they are not compared with individual sequences.")
            elif names == [reference["id"]]:
                status = "Same identifier as the reference FASTA; Ref bases are also checked against it."
                lengths = entry["info"].get("contig_lengths") or {}
                declared = lengths.get(reference["id"])
                if declared and str(declared).isdigit() and int(declared) != reference["length"]:
                    status = (f"VCF contig length {declared} differs from the reference FASTA "
                              f"({reference['length']}) under the same name, so they are different "
                              "sequences. Upload the matching reference before running.")
                    review["problems"].append(entry["info"]["name"] + ": " + status)
            else:
                others = [name for name in names if name != reference["id"]]
                status = ("Different reference label(s) " + ", ".join(others) + " versus "
                          + reference["id"] + ": those records stay incomparable. Labels are never "
                          "assumed to be the same genome.")
                review["warnings"].append(f"{entry['info']['name']}: {status}")
            rows.append({"Source": f"{label} ({entry['info']['name']})",
                         "Reference": ", ".join(names) or "-", "Status": status})


def analysis_plan(review, options):
    """Which analyses can run with the assigned files, and why others cannot."""

    roles = review["by_role"]
    has_fasta = "fasta" in roles
    has_reference = review["reference"] is not None
    has_snv = "snv" in roles
    completeness = {row["Source"]: row for row in review["completeness"]}
    individual = completeness.get("Individual sequences", {})
    pooled = completeness.get("Pooled samples", {})
    positions = roles.get("positions", [])
    has_counts = any(
        all(any(match["field"] == f"{base}_Count" and match["status"] == "matched"
                for match in entry.get("matches", [])) for base in "ACGT")
        for entry in positions)

    def row(name, ready, needs, reason):
        return {"Analysis": name, "Status": "Will run" if ready is True else (
            "Limited" if ready == "limited" else "Skipped"), "Needs": needs, "Reason": reason}

    plan = []
    individual_count = individual.get("Samples", 0)
    plan.append(row("Sequence comparison (MAFFT, p-distance)", has_fasta and individual_count >= 2,
                    "Two or more individual sequences with identifiers and a pathogen",
                    f"{individual_count} individual sequence(s)" if has_fasta else "No individual FASTA"))
    date_count = individual.get("With_Collection_Date", 0)
    plan.append(row("Temporal screening of sequence pairs",
                    has_fasta and isinstance(date_count, int) and date_count >= 2,
                    "Collection dates for at least two sequences",
                    f"{date_count} of {individual_count} sequences have a date" if has_fasta else "No individual FASTA"))
    geo = individual.get("With_Coordinates", 0) + individual.get("With_Location", 0) if has_fasta else 0
    plan.append(row("Geographic-distance screening", "limited" if has_fasta and geo else False,
                    "Coordinates, or specific place names that GeoNames can locate",
                    "Coarse (state-level) places never satisfy the kilometer cutoff"
                    if has_fasta and geo else "No locations for individual sequences"))
    plan.append(row("Reference-coordinate variant profiles", has_fasta and has_reference,
                    "Individual FASTA and a reference FASTA",
                    "Ready" if has_fasta and has_reference else "Needs both files"))
    plan.append(row("Pooled variant QC and frequencies", has_snv, "A pooled variant table or VCF",
                    f"{pooled.get('Samples', 0)} pooled sample(s)" if has_snv else "No pooled variant file"))
    coverage = pooled.get("With_Genome_Coverage")
    plan.append(row("Paper coverage filter (>90% at 10x)",
                    has_snv and (bool(positions) or (isinstance(coverage, int) and coverage > 0)),
                    "Genome coverage in pooled metadata or a per-position depth table",
                    "Coverage unknown samples are reported as not evaluable" if has_snv else "No pooled data"))
    pooled_dates = pooled.get("With_Collection_Date", 0)
    plan.append(row("Temporal frequency trends", has_snv and bool(pooled_dates),
                    "Collection dates for pooled samples",
                    "Dates available" if has_snv and pooled_dates else "No pooled dates"))
    plan.append(row("Paper PCoA (Yue-Clayton)", has_snv and has_counts,
                    "Per-position A/C/G/T counts for each pooled sample",
                    "A/C/G/T counts supplied" if has_counts else
                    ("An exploratory SNV-set ordination is shown instead" if has_snv else "No pooled data")))
    plan.append(row("Coverage-aware non-detection", has_snv and bool(positions),
                    "A per-position depth table",
                    "Unreported variants can be classed as covered or not" if positions else
                    "Without it, unreported variants stay 'coverage unknown' (never 0%)"))
    plan.append(row("Cross-source comparison", has_fasta and has_reference and has_snv,
                    "Individual FASTA, reference FASTA, and pooled variants on the same reference",
                    "Ready" if has_fasta and has_reference and has_snv else
                    "Needs " + ", ".join(name for name, present in [
                        ("individual sequences", has_fasta), ("a reference FASTA", has_reference),
                        ("pooled variant calls", has_snv)] if not present)))
    plan.append(row("Documented-profile overlap", "signatures" in roles and (has_snv or has_fasta),
                    "A documented profile table", "Profile overlap only, never a lineage call"
                    if "signatures" in roles else "No profile table"))
    return plan


# ----------------------------------------------------------------------------
# Standardized copies for the existing loaders
# ----------------------------------------------------------------------------

def write_fasta(entries, staged, destination):
    with Path(destination).open("w", encoding="utf-8") as target:
        for entry in entries:
            text = Path(staged[entry["index"]]).read_text(encoding="utf-8-sig", errors="replace")
            target.write(text.rstrip("\n") + "\n")
    return destination


def prepared_metadata(review, options):
    """Individual metadata with standard names, FASTA annotations, and a stated pathogen."""

    annotations = review.get("individual_annotations", {})
    entries = review["by_role"].get("metadata")
    if entries:
        table = apply_mapping(entries[0]["table"], entries[0]["matches"])
    else:
        table = pd.DataFrame({"Accession": review["individual_ids"]})
    keys = table["Accession"].map(normalize_accession) if "Accession" in table else None
    by_id = {normalize_accession(record_id): values for record_id, values in annotations.items()}
    for field in ("Pathogen", "Collection_Date", "Location_Name", "Host_Source"):
        if field in table.columns or keys is None:
            continue  # metadata values are never overwritten by FASTA annotations
        values = keys.map(lambda key: by_id.get(key, {}).get(field, ""))
        if values.ne("").any():
            table[field] = values
            table[f"{field}_Source"] = values.ne("").map({True: "FASTA description annotation",
                                                          False: ""})
    declared = str(options.get("declared_pathogen", "")).strip()
    if "Pathogen" not in table.columns and declared:
        table["Pathogen"] = declared
        table["Pathogen_Source"] = "Stated by the user for all sequences"
    return table


def combine_variant_tables(entries, staged, folder, notes=None):
    """Load each pooled variant file and write one standardized CSV.

    Each file's import warnings are added to `notes` with the file's name.
    """

    frames = []
    for number, entry in enumerate(entries):
        single = prepare_variant_file(entry, staged, folder, f"variants_part{number}")
        data, warnings = load_processed_snvs(single)
        if notes is not None:
            notes += [f"{entry['info']['name']}: {warning}" for warning in warnings]
        data["Collection_Date"] = data["Collection_Date"].dt.strftime("%Y-%m-%d")
        frames.append(data)
    combined = pd.concat(frames, ignore_index=True).astype(object)
    combined = combined.where(combined.notna(), "")  # blanks stay blank (unknown), never 0
    path = Path(folder) / "pooled_variants_combined.csv"
    combined.to_csv(path, index=False)
    return path


def prepare_variant_file(entry, staged, folder, stem):
    info = entry["info"]
    source = Path(staged[entry["index"]])
    if info["kind"] == VCF_KIND:
        # The loader uses the file name as the sample ID of a sites-only VCF,
        # so keep the original name inside a folder of its own.
        subfolder = Path(folder) / stem
        subfolder.mkdir(exist_ok=True)
        path = subfolder / f"{file_stem(info['name'])}.vcf"
        path.write_bytes(source.read_bytes())
        return path
    if info["format"] == PAPER_TABLE_FORMAT:
        path = Path(folder) / f"{stem}.xlsx"
        path.write_bytes(source.read_bytes())
        return path
    table = entry.get("table")
    if table is None:
        table = table_for(source, info)
    if info["format"] == IVAR_TABLE_FORMAT:
        if ivar_sample_column(table.columns) is None:
            table.insert(0, "SAMPLE", info["ids"][0])
        path = Path(folder) / f"{stem}.tsv"
        table.to_csv(path, sep="\t", index=False)
        return path
    path = Path(folder) / f"{stem}.csv"
    apply_mapping(table, entry["matches"]).to_csv(path, index=False)
    return path


def prepare_inputs(staged, review, folder, options=None):
    """Write the standardized inputs; return (paths by input key, display names)."""

    options = options or {}
    folder = Path(folder)
    paths = {}
    names = {}
    roles = review["by_role"]
    if "fasta" in roles:
        paths["fasta"] = write_fasta(roles["fasta"], staged, folder / "sequences.fasta")
        names["fasta"] = ", ".join(entry["info"]["name"] for entry in roles["fasta"])
        paths["metadata"] = folder / "metadata.csv"
        prepared_metadata(review, options).to_csv(paths["metadata"], index=False)
        names["metadata"] = (roles["metadata"][0]["info"]["name"] if "metadata" in roles
                             else "Built from FASTA record IDs (no metadata file)")
    if "reference" in roles:
        paths["reference"] = write_fasta(roles["reference"], staged, folder / "reference.fasta")
        names["reference"] = roles["reference"][0]["info"]["name"]
    if "snv" in roles:
        entries = roles["snv"]
        notes = []
        paths["snv"] = (prepare_variant_file(entries[0], staged, folder, "pooled_variants")
                        if len(entries) == 1 else combine_variant_tables(entries, staged, folder, notes))
        if notes:
            names["_notes"] = notes
        names["snv"] = ", ".join(entry["info"]["name"] for entry in entries)
    for role in ("pooled_metadata", "positions", "masked_sites", "signatures"):
        if role in roles:
            entry = roles[role][0]
            paths[role] = folder / f"{role}.csv"
            apply_mapping(entry["table"], entry["matches"]).to_csv(paths[role], index=False)
            names[role] = entry["info"]["name"]
    return paths, names
