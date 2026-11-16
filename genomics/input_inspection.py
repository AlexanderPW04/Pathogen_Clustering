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
import os
from functools import wraps
from shutil import copyfile
from tempfile import TemporaryDirectory
from threading import RLock
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from genomic_regions import fasta_region
from genomics.genome_annotation import (
    GENBANK_FORMAT,
    TABLE_NOT_STATED,
    detect_annotation_format,
    load_genbank,
    load_gff3,
)
from genomics.collection_dates import date_warnings, parse_collection_dates
from genomics.column_detection import (
    ROLE_FIELDS,
    SEQUENCE_COLUMN_NAMES,
    apply_mapping,
    detect_columns,
    explain_snv_columns,
    resolve_choices,
    squash,
    unresolved,
    value_warnings,
)
from genomics.excel_tables import WorkbookTooLarge, best_sheet, data_sheets, describe_sheets, read_sheet
from genomics.import_options import (
    DEPTH_DEFINITIONS,
    DOCUMENTED_DEPTH,
    frequency_problems,
    interpret_table,
    numbering_warning,
    reference_check,
)
from genomics.table_layouts import (
    BASE_COUNT_LAYOUT,
    LONG_LAYOUT,
    WIDE_LAYOUT,
    base_count_columns,
    derive_from_base_counts,
    detect_wide_matrix,
    wide_to_long,
)
from genomics.location_hierarchy import location_hierarchy
from genomics.analysis_readiness import analysis_readiness
from genomics.metadata_matching import match_metadata, resolve_metadata_table, resolve_variant_table


# Reading an upload is the most memory-hungry step before an analysis, and it is
# not covered by the one-analysis-at-a-time slot. Several visitors reviewing
# large tables at once could exhaust a small server (audit P1-1), so file
# inspection, review, and preparation run one at a time. The lock is
# re-entrant (a review may inspect a file again) and a visitor waits at most
# INSPECTION_WAIT_SECONDS before being told the server is busy.
INSPECTION_LOCK = RLock()
INSPECTION_WAIT_SECONDS = 60


class ServerBusy(ValueError):
    """Another upload is being read; nothing was changed."""


def one_at_a_time(function):
    @wraps(function)
    def guarded(*args, **kwargs):
        if not INSPECTION_LOCK.acquire(timeout=INSPECTION_WAIT_SECONDS):
            raise ServerBusy("The server is busy reading another upload. Wait a minute and try again; "
                             "your files are kept for 30 minutes.")
        try:
            return function(*args, **kwargs)
        finally:
            INSPECTION_LOCK.release()
    return guarded
from genomics.metadata_readiness import metadata_readiness
from geocode_locations import is_state_only_location
from genomics.file_formats import (
    ALIGNMENT_KIND,
    ANNOTATION_KIND,
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
    "annotation": "Genome annotation (GenBank or GFF3) for coding consequences",
    "ignore": "Do not use this file",
}
ROLE_ORDER = ["fasta", "metadata", "reference", "annotation", "snv", "pooled_metadata", "positions",
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
        "sheets": [], "sheet": None, "sheet_choice_needed": False, "preview": [], "variant_keys": 0,
        "candidates": [], "suggested_role": None, "role_evidence": "", "role_status": "",
        "auto_choices": {}, "annotations": {}, "descriptions": {}, "id_columns": {},
        "header_hints": {}, "layout": LONG_LAYOUT, "layout_counts": {}, "not_interpreted": [],
        "header_row_choice": None,
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


def inspect_annotation(path, info):
    """Describe a GenBank or GFF3 annotation: its references and coding features (not yet verified)."""

    annotation = (load_genbank(path) if detect_annotation_format(path) == GENBANK_FORMAT else load_gff3(path))
    records = annotation["records"]
    info["format"] = annotation["format"]
    info["records"] = len(records)
    info["reference_ids"] = list(records)
    info["lengths"] = [record["Length"] for record in records.values()]
    info["candidates"] = ["annotation"]
    info["unstated_code_cds"] = sum(cds["Translation_Table_Source"] == TABLE_NOT_STATED
                                    and not cds["Record_Translation"] for cds in annotation["cds"])
    if info["unstated_code_cds"]:
        info["notes"].append(f"{info['unstated_code_cds']} coding region(s) state no genetic code (transl_table) and "
                             "have no protein sequence to check one against; their amino-acid consequences are "
                             "predicted only if you confirm the genetic code under Advanced settings.")
    for reference_id, record in records.items():
        coding = sum(cds["Reference_ID"] == reference_id for cds in annotation["cds"])
        info["notes"].append(f"Describes {reference_id} ({record['Organism']}"
                             + (f", segment {record['Segment']}" if record["Segment"] else "")
                             + f", {record['Length'] or 'unknown'} bases): {coding} coding region(s)"
                             + ("" if record["Sequence"] else "; no sequence in the file, so a matching "
                                "reference FASTA is needed") + ".")


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


def choose_sheet(path, info, sheet):
    """Which worksheet to read; several data sheets are a question for the user, not a guess."""

    sheets = describe_sheets(path)
    info["sheets"] = [item["name"] for item in sheets]
    if sheet in info["sheets"]:
        return sheet
    candidates = data_sheets(sheets)
    if len(candidates) > 1:
        info["sheet_choice_needed"] = True
        info["problems"].append("Multiple worksheets were detected ("
                                + ", ".join(f"'{item['name']}'" for item in candidates)
                                + "). Select the worksheet containing the SNV observations or sample data.")
        return candidates[0]["name"]
    return best_sheet(sheets)


def read_raw_table(path, info):
    """The uploaded table as text, with the chosen worksheet and header row."""

    header = info.get("header_row_choice")
    if info["kind"] == EXCEL_KIND:
        if header:
            table = read_sheet(path, info["sheet"], header_row=header - 1)
        else:
            table = read_snv_table(path, info["sheet"])  # handles the paper's descriptive first row
    else:
        table = read_table(path, info["name"], header_row=(header - 1) if header else None)
    table.columns = [str(name).strip() for name in table.columns]
    if info.get("sample_from_filename") and not any(
            squash(column) in SAMPLE_FIELD_ALIASES for column in table.columns):
        # The user stated that this file holds one sample, named by the file.
        table.insert(0, "Sample_ID", file_stem(info["name"]))
    return table


SAMPLE_FIELD_ALIASES = next(spec["aliases"] for spec in ROLE_FIELDS["snv"] if spec["name"] == "Sample_ID")


def inspect_table(path, info, sheet=None, header_row=None, sample_from_filename=False):
    info["header_row_choice"] = header_row
    info["sample_from_filename"] = bool(sample_from_filename)
    if info["kind"] == EXCEL_KIND:
        info["sheet"] = choose_sheet(path, info, sheet)
    table = read_raw_table(path, info)
    if info["kind"] == EXCEL_KIND:
        if len(info["sheets"]) > 1:
            info["notes"].append(f"Worksheet '{info['sheet']}' was read (the workbook has "
                                 f"{len(info['sheets'])} worksheets).")
    if header_row:
        info["notes"].append(f"Column headers are read from row {header_row}, as you chose.")
    elif table.attrs.get("header_row"):
        info["notes"].append(f"Column headers were found on row {table.attrs['header_row'] + 1} "
                             "of the worksheet. If that is wrong, set 'Header row' under 'How to read this file'.")
    info["source_records"] = len(table)
    info["source_columns"] = len(table.columns)
    info["preview"] = table.head(3).iloc[:, :10].astype(str).to_dict("records")
    wide = detect_wide_matrix(table.columns)
    if wide is not None:
        table, counts = wide_to_long(table, wide)
        info["layout"] = WIDE_LAYOUT
        info["wide"] = wide
        info["layout_counts"] = counts
        info["not_interpreted"] = list(wide["not_interpreted"])
        info["notes"].append(
            f"Wide layout detected: {counts['samples']} sample rows x {counts['mutation_columns']} mutation "
            f"columns (headers such as '{next(iter(wide['mutation_columns']))}'). {counts['reported']:,} cells "
            f"hold a value and become records; {counts['not_reported']:,} blank or NA cells mean 'not "
            "reported' and become no record (never 0%). Confirm this reading under 'How to read this file'.")
        if wide["not_interpreted"]:
            info["notes"].append(
                f"{len(wide['not_interpreted'])} column(s) describe changes that are not single-nucleotide "
                "substitutions and are not analyzed: " + preview(wide["not_interpreted"], 5))
        samples = wide["other_columns"][0]
        empty = [str(value) for value in counts.get("rows_without_records", [])]
        if empty:
            info["notes"].append(f"{len(empty)} sample row(s) have no reported single-nucleotide mutation and "
                                 f"therefore no records (first column '{samples}': {preview(empty, 5)}).")
        if counts["repeated_mutations"]:
            info["notes"].append(f"{counts['repeated_mutations']} mutation column(s) repeat a nucleotide "
                                 "change already listed under another label; repeated records are flagged later.")
    info["records"] = len(table)
    info["columns"] = list(table.columns)
    columns = set(table.columns)

    if (table.attrs.get("fontenele_2021_table2") and set(PAPER_COLUMNS).issubset(columns)
            and {"STATE", "LOCATION"}.issubset(columns)):
        info["format"] = PAPER_TABLE_FORMAT
        info["candidates"] = ["snv"]
        info["ids"] = column_values(table, "SAMPLE")
        info["variant_keys"] = count_variants(table, ["POSITION", "REF", "ALT"])
        info["reference_ids"] = [PAPER_REFERENCE]
        info["notes"].append("Recognized the published supplementary SNV table; its reference is "
                             f"{PAPER_REFERENCE} as stated in the paper's methods.")
        return
    if IVAR_REQUIRED.issubset(columns):
        info["format"] = IVAR_TABLE_FORMAT
        info["candidates"] = ["snv"]
        info["variant_keys"] = count_variants(table, ["POS", "REF", "ALT"])
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
    base_counts = base_count_columns(matches["snv"] + matches["positions"]) if wide is None else None
    if complete("signatures"):
        candidates.append("signatures")
    elif complete("snv"):
        candidates.append("snv")
    elif base_counts is not None:
        candidates.append("snv")
        info["layout"] = BASE_COUNT_LAYOUT
        info["notes"].append(
            "This table has a reference base and A/C/G/T read counts but no alternative-base column. It can "
            "be used as variant records derived from the counts (each non-reference base with at least the "
            "minimum number of reads you confirm) or as per-position coverage. Tool-specific columns such as "
            "a consensus or 'variant' base are carried along, not reinterpreted.")
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
            elif reference:
                info["reference_column_missing"] = True
            sample = next((match for match in matches[role] if match["field"] == "Sample_ID"), None)
            if sample and sample["status"] == "matched" and role == "snv":
                info["ids"] = column_values(table, sample["column"])
                fields = {match["field"]: match["column"] for match in matches["snv"]
                          if match["status"] == "matched"}
                if {"Position", "Ref", "Alt"} <= set(fields):
                    info["variant_keys"] = count_variants(
                        table, [fields["Position"], fields["Ref"], fields["Alt"]])
                elif base_counts is not None:
                    # Upper bound: every listed site with a non-reference base (exact after derivation).
                    info["variant_keys"] = count_variants(table, [fields["Position"], fields["Ref"]]) * 3
    complete_but_sample = all(present("snv", name) for name in ("Position", "Ref", "Alt", "Alt_Frequency"))
    if not candidates and not present("snv", "Sample_ID") and present("snv", "Position") and (
            complete_but_sample or base_count_columns(
                [dict(match, status="matched", column=match["column"] or "Sample_ID") if match["field"] == "Sample_ID"
                 else match for match in matches["snv"] + matches["positions"]]) is not None):
        info["needs_sample_id"] = True
        info["problems"].append(
            "This table has variant positions but no sample column. If the file holds the results of one "
            "sample, tick 'Use the file name as the sample ID' under 'How to read this file' and press "
            f"Update review (the sample would be '{file_stem(info['name'])}').")
        return
    if not candidates:
        reason = explain_snv_columns(table.columns) if looks_like_variants(table.columns) else None
        if table.attrs.get("preamble") and "Supplementary Table" in table.attrs["preamble"] and reason:
            reason += " If this is a paper's supplement, another table may hold the per-sample values."
        info["problems"].append(reason or "The columns do not match any supported table. Required columns "
                                          "depend on the file: see the file guides on the upload page.")


def count_variants(table, columns):
    """Distinct position/Ref/Alt combinations (the size of the analysis)."""

    return int(table[columns].astype(str).apply(lambda column: column.str.strip().str.upper())
               .drop_duplicates().shape[0])


def looks_like_variants(columns):
    """True when headers suggest variant data, so a missing field can be named precisely."""

    found = {match["field"] for match in detect_columns(columns, "snv") if match["status"] != "missing"}
    return bool(found & {"Position", "Ref", "Alt", "Alt_Frequency"})


@one_at_a_time
def inspect_file(path, name, sheet=None, header_row=None, sample_from_filename=False):
    """Describe one decompressed upload; `sheet`, `header_row` (1-based), and
    `sample_from_filename` are the user's choices from the review page."""

    info = empty_info(name)
    info["kind"] = detect_kind(path, name)
    info["size_kb"] = round(Path(path).stat().st_size / 1024, 1)
    try:
        if info["kind"] == FASTA_KIND:
            inspect_fasta(path, info)
        elif info["kind"] == VCF_KIND:
            inspect_vcf(path, info)
        elif info["kind"] == ANNOTATION_KIND:
            inspect_annotation(path, info)
        elif info["kind"] in {TABLE_KIND, EXCEL_KIND}:
            inspect_table(path, info, sheet, header_row, sample_from_filename)
        elif info["kind"] in {FASTQ_KIND, ALIGNMENT_KIND}:
            info["refused"] = True
            info["problems"].append(RAW_READ_MESSAGE)
        else:
            info["refused"] = True
            info["problems"].append("This file type is not supported. Use FASTA, CSV/TSV/XLSX, "
                                    "VCF, or a GenBank/GFF3 annotation (each optionally .gz).")
    except WorkbookTooLarge as error:
        info["refused"] = True
        info["problems"].append(str(error))
    except MemoryError:
        info["refused"] = True
        info["problems"].append("The uploaded dataset exceeds the available memory for this server.")
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

def table_for(path, info, role=None, file_options=None):
    """The table as the analysis sees it: wide matrices and base counts are rewritten as long records."""

    table = read_raw_table(path, info)
    if info.get("layout") == WIDE_LAYOUT:
        table, _ = wide_to_long(table, info["wide"])
    elif info.get("layout") == BASE_COUNT_LAYOUT and role == "snv":
        minimum = (file_options or {}).get("base_count_minimum")
        if minimum:
            columns = base_count_columns(detect_columns(table.columns, "snv")
                                         + detect_columns(table.columns, "positions"))
            table, summary = derive_from_base_counts(table, columns, minimum)
            info["layout_counts"] = summary
            info["base_count_source_columns"] = columns
    return table


def native_format(info):
    return info["kind"] == VCF_KIND or info["format"] in {PAPER_TABLE_FORMAT, IVAR_TABLE_FORMAT}


@one_at_a_time
def review_dataset(staged, inspections, roles, choices=None, options=None):
    """Summarize the dataset under the chosen roles.

    `staged` maps file index -> decompressed path; `roles` maps index -> role;
    `choices` maps index -> {field: column} for ambiguous columns; `options`
    may hold "declared_pathogen" and "complete_genomes".
    """

    choices = choices or {}
    options = options or {}
    review = {"options": options, "files": [], "problems": [], "warnings": [], "by_role": {}, "mappings": {},
              "individual_ids": [], "pooled_ids": [], "reference": None, "compatibility": [],
              "completeness": [], "analyses": [], "annotation_notes": [], "ready": False}

    review["reference_assigned"] = any(
        roles.get(index, info["suggested_role"]) == "reference" for index, info in enumerate(inspections))
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
            file_options = options.get("files", {}).get(index, {})
            entry["options"] = file_options
            try:
                table = table_for(staged[index], info, role, file_options)
            except ValueError as problem:
                entry["problems"].append(str(problem))
                continue
            matches = resolve_choices(detect_columns(table.columns, role), choices.get(index, {}))
            matches = resolve_choices(matches, info["auto_choices"].get(role, {}),
                                      "matching its values with the other uploaded files")
            entry["matches"] = matches
            entry["problems"] += unresolved(matches) + value_warnings(
                table, matches, file_options.get("frequency_scale"))
            interpretation_review(review, entry, table, role, file_options)
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
                sequence = clean_sequence(record.seq)
                sequence_bases = set(sequence)
                review["reference_sequence"] = {"id": info["ids"][0], "sequence": sequence}
            if sequence_bases - set("ACGTN"):
                review["problems"].append("The reference sequence may contain only A, C, G, T, and N.")
            numbering_review(review)
    for entry in review["files"]:
        for problem in entry["problems"]:
            review["problems"].append(f"{entry['info']['name']}: {problem}")
        for problem in entry["info"]["problems"]:
            if not entry["info"]["refused"] and entry["role"] not in {None, "ignore"}:
                review["problems"].append(f"{entry['info']['name']}: {problem}")

    for entry in review["by_role"].get("annotation", []):
        if entry["info"].get("unstated_code_cds") and not options.get("genetic_code"):
            review["warnings"].append(
                f"{entry['info']['name']}: {entry['info']['unstated_code_cds']} coding region(s) state no genetic "
                "code. Their amino-acid consequences will be reported as not predicted unless you choose the "
                "organism's genetic code under Advanced settings ('Genetic code for coding regions that state "
                "none'). The standard code is never assumed.")
    individual_review(review, options)
    pooled_review(review, options.get("date_order", "auto"))
    compatibility_review(review)
    review["analyses"] = analysis_plan(review, options)
    review["analysis_readiness"] = analysis_readiness(review, options)
    review["readiness"] = metadata_readiness(review, options)
    review["ready"] = not review["problems"] and any(
        role in review["by_role"] for role in ("fasta", "snv", "reference"))
    if review["ready"]:
        pooled_validation_review(review, staged, options)
        review["ready"] = not review["problems"]
    if not review["by_role"]:
        review["problems"].append("No file has been assigned a role yet.")
    return review


def pooled_validation_review(review, staged, options):
    """Run the analysis's own pooled-table validation now, so Step 2 shows what Run would refuse.

    Audit P1-3: some checks (alternative reads above depth, a non-ACGT base,
    an unreadable metadata date, ...) used to run only when the analysis
    started, after the review had said "ready", and the uploads were then
    deleted. The same function the analysis uses (`load_processed_snvs`) is
    called here on the same prepared files, so the two can never disagree.
    Nothing is changed or filled in; the problem is shown and Run stays
    disabled until it is corrected.
    """

    if "snv" not in review["by_role"]:
        return
    with TemporaryDirectory() as folder:
        try:
            paths, _ = prepare_inputs(staged, review, Path(folder), options)
            if "snv" in paths:
                load_processed_snvs(paths["snv"], paths.get("pooled_metadata"),
                                    date_order=options.get("date_order", "auto"))
        except ServerBusy:
            raise
        except ValueError as problem:
            text = str(problem).replace(str(folder), "[temporary upload]")
            review["problems"].append(
                f"The analysis would refuse the pooled variant data: {text} Correct the file and upload it "
                "again, or change the roles or choices above (for example, stop using a metadata file whose "
                "values cannot be read). Nothing was changed or guessed.")


def interpretation_review(review, entry, table, role, file_options):
    """Check the choices that headers cannot settle (layout, scale, reference, depth) for one table."""

    info = entry["info"]
    name = info["name"]
    columns = {match["field"]: match["column"] for match in entry["matches"] if match["status"] == "matched"}
    entry["columns"] = columns
    if info.get("layout") == WIDE_LAYOUT and role == "snv" and not file_options.get("layout_confirmed"):
        entry["problems"].append(
            "Confirm the wide layout under 'How to read this file': each row is one sample and each "
            "mutation column holds that sample's reported frequency (blank or NA = not reported).")
    if info.get("layout") == BASE_COUNT_LAYOUT and role == "snv" and not file_options.get("base_count_minimum"):
        entry["problems"].append(
            "Enter the minimum number of reads for a non-reference base under 'How to read this file' "
            "(use the threshold documented by the source). Variant records are derived from the A/C/G/T "
            "counts only after you confirm it.")
    if role in {"snv", "positions"}:
        field = "Reference_Genome" if role == "snv" else "Reference_ID"
        declared = file_options.get("declared_reference")
        if field not in columns:
            if not declared:
                entry["problems"].append(
                    "This table has no reference-genome column. State the reference accession the source "
                    "used (for example from its methods) under 'How to read this file'.")
            elif several_sequences(table, columns):
                column, values = several_sequences(table, columns)
                entry["problems"].append(
                    f"Column '{column}' names {len(values)} different sequences (for example "
                    f"{', '.join(values[:3])}), so one stated reference would put records from different "
                    "segments or chromosomes on one sequence. Choose that column as the reference column "
                    "under 'Columns', or split the file into one file per sequence.")
            else:
                entry["declared_reference"] = declared
                info["reference_ids"] = [declared]
                review["annotation_notes"].append(
                    f"{name}: reference '{declared}' stated by you (the table has no reference column). "
                    + ("It is checked base by base against the uploaded reference FASTA when the names match."
                       if review.get("reference_assigned") else
                       "No reference FASTA was uploaded, so it stays unverified."))
        elif declared:
            review["annotation_notes"].append(f"{name}: the stated reference '{declared}' was not used because "
                                              f"the table has its own reference column '{columns[field]}'.")
    if role == "snv" and "Alt_Frequency" in columns:
        problems, warnings = frequency_problems(table, columns["Alt_Frequency"], file_options)
        entry["problems"] += problems
        review["warnings"] += [f"{name}: {text}" for text in warnings]
    if role == "snv" and info.get("layout") != BASE_COUNT_LAYOUT:
        if "Read_Depth" not in columns:
            review["warnings"].append(
                f"{name}: no read-depth column was recognized. Reported frequencies are still summarized, but "
                "depth QC, count-based frequencies, uncertainty intervals, and the depth-aware methods are "
                "unavailable for this file. Depth is never filled in.")
        elif not file_options.get("depth_definition"):
            review["annotation_notes"].append(
                f"{name}: column '{columns['Read_Depth']}' is used as the total usable depth at each position. "
                "If the source defines it differently, state the definition under 'How to read this file'.")


def several_sequences(table, columns):
    """(column, values) when an unmapped segment/chromosome-like column holds more than one value."""

    if table is None:
        return None
    used = set(columns.values())
    for column in table.columns:
        if column in used or squash(column) not in SEQUENCE_COLUMN_NAMES:
            continue
        values = sorted({str(value).strip() for value in table[column].dropna() if str(value).strip()})
        if len(values) > 1:
            return column, values
    return None


def numbering_review(review):
    """Warn when Ref bases fit the reference only under the other position numbering."""

    reference = review.get("reference_sequence")
    for role in ("snv", "positions"):
        for entry in review["by_role"].get(role, []):
            if entry.get("table") is None or "Ref" not in entry.get("columns", {}):
                continue
            columns = dict(entry["columns"])
            if "Reference_Genome" not in columns and entry.get("declared_reference") != (reference or {}).get("id"):
                continue
            check = reference_check(entry["table"], columns, reference)
            entry["numbering_check"] = check
            warning = numbering_warning(check, entry.get("options", {}), entry["info"]["name"])
            if warning:
                review["warnings"].append(warning)


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
    state_text = int((result["Location_Level"].eq("") & result["Location_Name"].map(is_state_only_location)).sum())
    coarse = summary["country_level"] + summary["state_level"] + state_text
    example = ""
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
    located = int(result["Location_Name"].ne("").sum())
    review["location_summary"] = {"located": located, "coarse": coarse, "example": example,
                                  "state": summary["state_level"] + state_text,
                                  "country": summary["country_level"], "conflicts": summary["conflicts"]}
    return located


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
    region_columns = [] if metadata is None else [
        column for column in metadata["info"]["columns"]
        if squash(column) in {"genesegment", "segment", "genename", "gene", "region"}]
    region_known = bool(region_columns)
    header_regions = sum(1 for description in descriptions.values()
                         if fasta_region(description)["label"])
    review["individual_region"] = {"column": region_columns[0] if region_columns else None,
                                   "from_headers": header_regions}
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


def max_observation_cells():
    """Largest variant-by-sample analysis allowed (0 = no limit).

    Peak memory measured on generated tables is about 205 MB plus 0.32 MB per
    1,000 cells, so the Dockerfile sets 700,000 for a 512 MB server. Local runs
    have no limit unless MAX_OBSERVATION_CELLS is set.
    """

    try:
        return max(0, int(os.getenv("MAX_OBSERVATION_CELLS", "0")))
    except ValueError:
        return 0


def check_analysis_size(review, entries):
    """Refuse an analysis that could not finish here; data are never trimmed to fit."""

    limit = max_observation_cells()
    variants = sum(entry["info"].get("variant_keys", 0) for entry in entries)
    cells = variants * len(review["pooled_ids"])
    if limit and cells > limit:
        review["problems"].append(
            f"The uploaded dataset exceeds the available memory for this server: {variants:,} variants x "
            f"{len(review['pooled_ids']):,} samples is {cells:,} combinations, and this server can process "
            f"{limit:,}. No samples or variants are removed automatically. Run the analysis on a local "
            "computer (see the README), or upload a smaller study.")


def pooled_date_review(review, entries, date_order):
    """Flag ambiguous, partial, or two-digit-year pooled collection dates at review (audit S3-05)."""

    columns = [mapped_values(entry, "Collection_Date") for entry in entries]
    metadata = review["by_role"].get("pooled_metadata", [None])[0]
    if metadata is not None:
        columns.append(mapped_values(metadata, "Collection_Date"))
    values = [column for column in columns if column is not None]
    if not values:
        return
    parsed = parse_collection_dates(pd.concat(values, ignore_index=True), date_order)
    review["pooled_date_status"] = parsed["Status"].value_counts().to_dict()
    for message in date_warnings(parsed):
        review["warnings"].append("Pooled collection dates: " + message)


def pooled_review(review, date_order="auto"):
    entries = review["by_role"].get("snv", [])
    if not entries:
        return
    pooled_date_review(review, [entry for entry in entries
                                if entry["info"]["format"] != PAPER_TABLE_FORMAT], date_order)
    ids = []
    for entry in entries:
        ids += [str(value) for value in entry["info"]["ids"]]
    if len(entries) > 1 and len(ids) != len(set(ids)):
        review["problems"].append("The same pooled sample ID appears in more than one variant file.")
    if any(entry["info"]["format"] == PAPER_TABLE_FORMAT for entry in entries) and len(entries) > 1:
        review["problems"].append("The published supplementary table must be analyzed on its own.")
    review["pooled_ids"] = sorted(set(ids))
    check_analysis_size(review, entries)
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
        # Phase 4A fix: count only metadata rows of samples that have variant records
        # (rows of other samples, or without a Sample_ID, describe nothing analyzed).
        keys = mapped_values(metadata, "Sample_ID")
        joined = (keys.isin(set(review["pooled_ids"])) if keys is not None
                  else pd.Series(False, index=metadata["table"].index))

        def matched_present(field):
            values = mapped_values(metadata, field)
            return None if values is None else values.where(joined, "")

        dates = max(dates, count_present(matched_present("Collection_Date")))
        locations = max(locations, count_present(matched_present("Location_Name")))
        latitudes = matched_present("Latitude")
        longitudes = matched_present("Longitude")
        if latitudes is not None and longitudes is not None:
            coordinates = int((latitudes.ne("") & longitudes.ne("") & joined).sum())
        coverage = count_present(matched_present("Genome_Coverage_Percent"))
    elif not paper:
        if vcf_only:
            review["warnings"].append("VCF files carry no collection dates or sites; add pooled-sample "
                                      "metadata for temporal and geographic analyses.")
    if review["by_role"].get("positions"):
        coverage = "From the per-position table"
    elif metadata is None:
        coverage = "Not supplied"
    match_metadata(review, review.get("options", {}))
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
    for entry in review["by_role"].get("annotation", []):
        rows.append({"Source": f"Genome annotation ({entry['info']['name']})",
                     "Reference": ", ".join(entry["info"]["reference_ids"]) or "-",
                     "Status": annotation_status(review, entry["info"])})


def variant_reference_ids(review):
    names = set()
    for entry in review["by_role"].get("snv", []):
        names |= set(entry["info"]["reference_ids"])
    if review["by_role"].get("fasta") and review["reference"]:
        names.add(review["reference"]["id"])
    return names


def annotation_status(review, info):
    """How the annotation will be matched; it is applied only to the exact same accession.version."""

    described = set(info["reference_ids"])
    variants = variant_reference_ids(review)
    reference = review["reference"]
    if reference and reference["id"] in described:
        lengths = dict(zip(info["reference_ids"], info["lengths"]))
        if lengths.get(reference["id"]) and lengths[reference["id"]] != reference["length"]:
            text = (f"The annotation's {reference['id']} has {lengths[reference['id']]} bases but the reference "
                    f"FASTA has {reference['length']}: different sequences, so nothing will be annotated.")
            review["warnings"].append(f"{info['name']}: {text}")
            return text
    matched = sorted(described & variants)
    if matched:
        return ("Describes " + ", ".join(matched) + ", the reference of the uploaded variants. The sequence is "
                "checked again before use, and every reference base is verified.")
    text = ("Describes " + (", ".join(sorted(described)) or "no sequence") + ", but the variants use "
            + (", ".join(sorted(variants)) or "no identified reference")
            + ". Annotations are never transferred between references, so coding consequences will not be "
              "determined.")
    review["warnings"].append(f"{info['name']}: {text}")
    return text


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
            "Limited" if ready == "limited" else "Unavailable"), "Needs": needs, "Reason": reason}

    plan = []
    individual_count = individual.get("Samples", 0)
    plan.append(row("Sequence comparison (MAFFT, p-distance)", has_fasta and individual_count >= 2,
                    "Two or more individual sequences with identifiers and a pathogen",
                    f"{individual_count} individual sequence(s)" if has_fasta else "No individual FASTA"))
    date_count = individual.get("With_Collection_Date", 0)
    dated = has_fasta and isinstance(date_count, int) and date_count >= 2
    plan.append(row("Temporal screening of sequence pairs",
                    (True if date_count >= individual_count else "limited") if dated else False,
                    "Collection dates for at least two sequences",
                    f"{date_count} of {individual_count} sequences have a date" if has_fasta else "No individual FASTA"))
    plan.append(geographic_row(review, individual, has_fasta, row))
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
    annotation = roles.get("annotation")
    matched = bool(annotation) and bool(set(annotation[0]["info"]["reference_ids"]) & variant_reference_ids(review))
    plan.append(row("Mutation annotation and protein consequences", matched and (has_snv or has_fasta),
                    "A GenBank or GFF3 annotation of the variants' reference (same accession.version)",
                    "Coding consequences will be predicted for each distinct SNV" if matched and (has_snv or has_fasta)
                    else ("The annotation describes a different reference" if annotation else
                          "Optional: without an annotation, coding consequences are not determined")))
    plan.append(row("Documented-profile overlap", "signatures" in roles and (has_snv or has_fasta),
                    "A documented profile table", "Profile overlap only, never a lineage call"
                    if "signatures" in roles else "No profile table"))
    return plan


def geographic_row(review, individual, has_fasta, row):
    """Kilometer screening needs coordinates or specific place names; coarse places never pass."""

    name = "Geographic-distance screening"
    needs = "Coordinates, or specific place names that GeoNames can locate"
    if not has_fasta:
        return row(name, False, needs, "No individual FASTA")
    total = individual.get("Samples", 0)
    coordinates = individual.get("With_Coordinates", 0)
    located = individual.get("With_Location", 0)
    coarse = review.get("location_summary", {}).get("coarse", 0)
    specific = max(located - coarse, 0)
    if not located and not coordinates:
        return row(name, False, needs, "No locations for individual sequences; genetic and date "
                                       "comparisons still run")
    if not specific and not coordinates:
        return row(name, False, needs, "Locations are only at state or country level; genetic and date "
                                       "comparisons still run and pairs are reported as unresolved geography")
    if coordinates >= total:
        return row(name, True, needs, "Coordinates supplied for every sequence")
    parts = [f"{specific} specific place name(s), located with GeoNames during the run"]
    if coordinates:
        parts.append(f"{coordinates} with coordinates")
    if coarse:
        parts.append(f"{coarse} only at state or country level (never pass the kilometer cutoff)")
    if located < total:
        parts.append(f"{total - located} without a location")
    return row(name, "limited", needs, "; ".join(parts))


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


def combine_variant_tables(entries, staged, folder, notes=None, date_order="auto"):
    """Load each pooled variant file and write one standardized CSV.

    Each file's import warnings are added to `notes` with the file's name.
    A date that could not be used (for example an ambiguous day/month) is
    written back as its original text, so the combined table flags it again.
    """

    frames = []
    for number, entry in enumerate(entries):
        single = prepare_variant_file(entry, staged, folder, f"variants_part{number}")
        data, warnings = load_processed_snvs(single, date_order=date_order)
        data["Source_File"] = entry["info"]["name"]  # which uploaded file each record came from
        if notes is not None:
            notes += [f"{entry['info']['name']}: {warning}" for warning in warnings]
        data["Collection_Date"] = (data["Collection_Date"].dt.strftime("%Y-%m-%d")
                                   .fillna(data["Collection_Date_Original"]))
        data = data.drop(columns=["Collection_Date_Original", "Collection_Date_Status",
                                  "Frequency_Significant_Digits"])
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
        table = table_for(source, info, "snv", entry.get("options"))
    if info["format"] == IVAR_TABLE_FORMAT:
        if ivar_sample_column(table.columns) is None:
            table.insert(0, "SAMPLE", info["ids"][0])
        path = Path(folder) / f"{stem}.tsv"
        table.to_csv(path, sep="\t", index=False)
        return path
    path = Path(folder) / f"{stem}.csv"
    interpreted_table(entry, "snv").to_csv(path, index=False)
    return path


def interpreted_table(entry, role):
    """The mapped table with the user's confirmed interpretation applied (originals kept)."""

    table = apply_mapping(entry["table"], entry["matches"])
    if role == "snv" and entry.get("conflicts"):
        table = resolve_variant_table(table, entry["conflicts"])
    field = "Reference_Genome" if role == "snv" else "Reference_ID"
    has_reference = field in entry.get("columns", {})
    table, rows = interpret_table(table, entry.get("options", {}), has_reference)
    if role == "positions" and "Reference_Genome" in table and not has_reference:
        table = table.rename(columns={"Reference_Genome": "Reference_ID"})
    if role != "snv" and "Read_Depth" in table and "Read_Depth" not in entry.get("columns", {}):
        table = table.drop(columns=["Read_Depth"])
        rows = [row for row in rows if row[0] != "Read_Depth"]
    entry["interpretation"] = rows
    return table


@one_at_a_time
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
        for entry in entries:
            entry["conflicts"] = review.get("metadata_conflicts", [])
        notes = []
        paths["snv"] = (prepare_variant_file(entries[0], staged, folder, "pooled_variants")
                        if len(entries) == 1 else combine_variant_tables(entries, staged, folder, notes,
                                                                         options.get("date_order", "auto")))
        if notes:
            names["_notes"] = notes
        names["snv"] = ", ".join(entry["info"]["name"] for entry in entries)
    if "annotation" in roles:
        entry = roles["annotation"][0]
        paths["annotation"] = folder / "annotation.txt"
        copyfile(staged[entry["index"]], paths["annotation"])
        names["annotation"] = entry["info"]["name"]
    for role in ("pooled_metadata", "positions", "masked_sites", "signatures"):
        if role in roles:
            entry = roles[role][0]
            paths[role] = folder / f"{role}.csv"
            if role == "positions":
                table = interpreted_table(entry, role)
            elif role == "pooled_metadata":
                table = prepared_pooled_metadata(review, entry, options)
            else:
                table = apply_mapping(entry["table"], entry["matches"])
            table.to_csv(paths[role], index=False)
            names[role] = entry["info"]["name"]
    if "snv" in roles and "pooled_metadata" not in roles:
        derived = coordinates_from_variants(review)
        if derived is not None:
            paths["pooled_metadata"] = folder / "pooled_metadata_from_variant_table.csv"
            derived.to_csv(paths["pooled_metadata"], index=False)
            names["pooled_metadata"] = "Coordinates from the variant table columns"
    names["_provenance"] = import_provenance(review)
    return paths, names


def variant_coordinates(review):
    """{Sample_ID: (latitude, longitude)} from mapped coordinate columns of the variant tables."""

    found = {}
    for entry in review["by_role"].get("snv", []):
        columns = entry.get("columns", {})
        if entry.get("table") is None or not {"Sample_ID", "Latitude", "Longitude"} <= set(columns):
            continue
        table = entry["table"]
        samples = table[columns["Sample_ID"]].astype(str).str.strip()
        latitude = table[columns["Latitude"]].astype(str).str.strip()
        longitude = table[columns["Longitude"]].astype(str).str.strip()
        for sample, lat, lon in zip(samples, latitude, longitude):
            if lat and lon and sample not in found:
                found[sample] = (lat, lon)
    return found


def coordinates_from_variants(review):
    """Sample metadata holding only the variant table's coordinates (conflicts left blank), or None."""

    coordinates = variant_coordinates(review)
    if not coordinates:
        return None
    conflicted = {(item["Sample_ID"], item["Field"]) for item in review.get("metadata_conflicts", [])}
    rows = [{"Sample_ID": sample,
             "Latitude": "" if (sample, "Latitude") in conflicted else lat,
             "Longitude": "" if (sample, "Longitude") in conflicted else lon,
             "Metadata_Status": "Coordinates copied from the variant table"}
            for sample, (lat, lon) in coordinates.items()]
    samples = {row["Sample_ID"] for row in rows}
    rows += [{"Sample_ID": sample, "Latitude": "", "Longitude": "", "Metadata_Status": "No coordinates supplied"}
             for sample in review["pooled_ids"] if sample not in samples]
    return pd.DataFrame(rows)


def prepared_pooled_metadata(review, entry, options):
    """The pooled metadata with resolutions applied; blank coordinates filled from the variant table."""

    table = resolve_metadata_table(apply_mapping(entry["table"], entry["matches"]),
                                   review.get("metadata_conflicts", []), review["pooled_ids"],
                                   options.get("allow_unmatched_metadata"))
    coordinates = variant_coordinates(review)
    if coordinates:
        conflicted = {(item["Sample_ID"], item["Field"]) for item in review.get("metadata_conflicts", [])}
        for field, position in (("Latitude", 0), ("Longitude", 1)):
            if field not in table:
                table[field] = ""
            blank = table[field].fillna("").astype(str).str.strip().eq("")
            fill = table["Sample_ID"].map(lambda sample: "" if (sample, field) in conflicted
                                          else coordinates.get(sample, ("", ""))[position])
            table.loc[blank, field] = fill[blank]
    return table


def import_provenance(review):
    """One row per file decision: layout, worksheet, header row, field columns, and transformations."""

    rows = []
    for entry in review["files"]:
        info, role = entry["info"], entry["role"]
        if info["refused"] or role in {None, "ignore"}:
            continue
        name, layout = info["name"], info.get("layout") or info.get("format") or info["kind"]

        def add(item, source, transformation="", decided="Detected"):
            rows.append({"File": name, "Role": ROLE_LABELS.get(role, role), "Layout": layout, "Item": item,
                         "Source": source, "Transformation": transformation, "Decided_By": decided})

        add("Format", info.get("format") or info["kind"])
        if info.get("sheet"):
            add("Worksheet", info["sheet"], decided="User" if len(info.get("sheets", [])) > 1 else "Detected")
        if info.get("header_row_choice"):
            add("Header row", f"row {info['header_row_choice']}", decided="User")
        if info.get("sample_from_filename"):
            add("Sample_ID", f"file name '{file_stem(name)}'", "The file holds one sample (stated by the user)",
                "User")
        if native_format(info):
            depth = DOCUMENTED_DEPTH.get(info.get("format"), DOCUMENTED_DEPTH["VCF"] if info["kind"] == VCF_KIND else "")
            if depth:
                add("Read_Depth definition", depth, decided="Documented format")
        if info.get("layout") == WIDE_LAYOUT:
            counts = info.get("layout_counts", {})
            add("Layout", f"{counts.get('samples')} sample rows x {counts.get('mutation_columns')} mutation columns",
                f"Wide to long: {counts.get('reported')} reported cells became records; {counts.get('not_reported')} "
                "blank/NA cells became no record", "User confirmed")
        if info.get("layout") == BASE_COUNT_LAYOUT and role == "snv":
            counts = info.get("layout_counts", {})
            add("Layout", f"{counts.get('rows')} base-count rows",
                f"Each non-reference base with at least {counts.get('minimum_reads')} read(s) became one record "
                f"({counts.get('records')} records); frequency calculated as reads / (A+C+G+T)", "User confirmed")
        derived_fields = {"Alt": "each non-reference base with enough reads",
                          "Alt_Frequency": "alternative reads / (A+C+G+T)", "Read_Depth": "A+C+G+T",
                          "Alt_Depth": "reads of the alternative base", "Ref_Depth": "reads of the Ref base"}
        base_columns = info.get("base_count_source_columns") if role == "snv" else None
        for match in entry.get("matches", []):
            if match["status"] != "matched":
                continue
            decided = "User" if "chosen by the user" in match["evidence"] else (
                "Matched with other files" if "chosen by matching" in match["evidence"] else "Detected")
            if base_columns and match["field"] in derived_fields:
                counts = ", ".join(f"'{base_columns[f'{base}_Count']}'" for base in "ACGT")
                add(match["field"], f"derived from base-count columns {counts}", derived_fields[match["field"]],
                    "Importer (user-confirmed minimum)")
            elif base_columns and match["field"] in base_columns:
                if not (match["field"] == "Sample_ID" and info.get("sample_from_filename")):
                    add(match["field"], f"column '{base_columns[match['field']]}'", "recognized name", "Detected")
            elif not (match["field"] == "Sample_ID" and info.get("sample_from_filename")):
                add(match["field"], f"column '{match['column']}'", match["evidence"], decided)
        for field, text in entry.get("interpretation", []):
            add(field, "Option on the review page", text, "User")
        definition = entry.get("options", {}).get("depth_definition")
        if definition:
            add("Read_Depth definition", DEPTH_DEFINITIONS[definition], decided="User")
    return rows


def mapping_rows(review):
    """Alias kept short for the view: the provenance rows of the current review."""

    return import_provenance(review)
