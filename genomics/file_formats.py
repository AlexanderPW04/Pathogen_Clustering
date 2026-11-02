"""Recognize and read the input file formats the application accepts.

Supported inputs (see the README table for full / partial support):

- sequences: .fasta, .fa, .fna, optionally gzip-compressed (.gz);
- tables: .csv, .tsv, .txt (tab-separated), .xlsx; text tables may be .gz;
- variant calls: .vcf or .vcf.gz (single-nucleotide variants only).

Raw sequencing reads (FASTQ) and alignments (SAM/BAM/CRAM) are recognized so
they can be refused with a clear message: they need the separate offline
workflow, they are never treated as processed variant tables.

Compressed files are decompressed into a plain copy with a size limit, so every
existing loader keeps reading ordinary files.
"""

import gzip
import shutil
from pathlib import Path

import pandas as pd


SEQUENCE_EXTENSIONS = {".fasta", ".fa", ".fna"}
TABLE_EXTENSIONS = {".csv", ".tsv", ".txt", ".xlsx"}
VARIANT_EXTENSIONS = {".vcf"}
RAW_READ_EXTENSIONS = {".fastq", ".fq"}
ALIGNMENT_EXTENSIONS = {".sam", ".bam", ".cram"}
SUPPORTED_EXTENSIONS = SEQUENCE_EXTENSIONS | TABLE_EXTENSIONS | VARIANT_EXTENSIONS

# A small upload can expand a lot when decompressed; refuse anything larger.
MAX_DECOMPRESSED_BYTES = 60 * 1024 * 1024

FASTA_KIND = "FASTA sequences"
TABLE_KIND = "Table"
EXCEL_KIND = "Excel workbook"
VCF_KIND = "VCF variant calls"
FASTQ_KIND = "FASTQ raw reads"
ALIGNMENT_KIND = "Read alignment (SAM/BAM/CRAM)"
UNKNOWN_KIND = "Unrecognized"

RAW_READ_MESSAGE = (
    "Raw sequencing reads cannot be analyzed as processed results. Map and call "
    "variants first with the optional offline workflow (offline_workflow/README.md), "
    "then upload the resulting variant table or VCF."
)

VCF_FORMAT = "VCF"


def split_name(filename):
    """Return (lower-case base extension, compressed?) for a file name.

    `sample.fasta.gz` -> (".fasta", True); `table.csv` -> (".csv", False).
    """

    path = Path(str(filename).strip())
    suffixes = [suffix.lower() for suffix in path.suffixes]
    compressed = bool(suffixes) and suffixes[-1] in {".gz", ".gzip"}
    if compressed:
        suffixes = suffixes[:-1]
    return (suffixes[-1] if suffixes else ""), compressed


def file_stem(filename):
    """The file name without compression and format extensions ("S1.vcf.gz" -> "S1").

    Used only as a visible fallback sample identifier for single-sample files
    that carry no sample name; it is never used to infer scientific attributes.
    """

    name = Path(str(filename).strip()).name
    extension, compressed = split_name(name)
    if compressed:
        name = name[: name.lower().rfind(".gz")]
    if extension and name.lower().endswith(extension):
        name = name[: -len(extension)]
    return "".join(character if character.isalnum() or character in "-_." else "_"
                   for character in name) or "sample"


def is_gzip(path):
    with Path(path).open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def plain_copy(path, destination):
    """Copy `path` to `destination`, decompressing gzip data with a size limit."""

    path = Path(path)
    destination = Path(destination)
    if not is_gzip(path):
        if path.resolve() != destination.resolve():
            shutil.copyfile(path, destination)
        return destination
    written = 0
    try:
        with gzip.open(path, "rb") as source, destination.open("wb") as target:
            while True:
                block = source.read(1024 * 1024)
                if not block:
                    break
                written += len(block)
                if written > MAX_DECOMPRESSED_BYTES:
                    raise ValueError(
                        f"{path.name} expands to more than "
                        f"{MAX_DECOMPRESSED_BYTES // (1024 * 1024)} MB when decompressed; "
                        "analyze it with run_integrated_analysis.py on your own computer."
                    )
                target.write(block)
    except (OSError, EOFError) as error:
        destination.unlink(missing_ok=True)
        raise ValueError(f"{path.name} is not a readable gzip file ({error}).") from error
    except ValueError:
        destination.unlink(missing_ok=True)
        raise
    return destination


def first_lines(path, count=4):
    """Return up to `count` non-empty text lines from the start of a plain file."""

    lines = []
    with Path(path).open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            if line.strip():
                lines.append(line.rstrip("\r\n"))
            if len(lines) >= count:
                break
    return lines


def detect_kind(path, filename):
    """Decide what kind of file this is from its content, using the name only as a hint.

    `path` must already be decompressed (see `plain_copy`).
    """

    extension, _ = split_name(filename)
    with Path(path).open("rb") as handle:
        start = handle.read(4)
    if start.startswith(b"PK"):
        return EXCEL_KIND if extension == ".xlsx" else UNKNOWN_KIND
    if start.startswith(b"BAM\x01") or start.startswith(b"CRAM") or extension in ALIGNMENT_EXTENSIONS:
        return ALIGNMENT_KIND
    lines = first_lines(path)
    if not lines:
        return UNKNOWN_KIND
    if lines[0].startswith("##fileformat=VCF") or (lines[0].startswith("#CHROM") and "POS" in lines[0]):
        return VCF_KIND
    if lines[0].startswith(">"):
        return FASTA_KIND
    if lines[0].startswith("@") and len(lines) >= 3 and lines[2].startswith("+"):
        return FASTQ_KIND
    if extension in RAW_READ_EXTENSIONS:
        return FASTQ_KIND
    if extension in TABLE_EXTENSIONS:
        return TABLE_KIND
    if "," in lines[0] or "\t" in lines[0]:
        return TABLE_KIND
    return UNKNOWN_KIND


def table_separator(path, filename):
    extension, _ = split_name(filename)
    if extension == ".csv":
        return ","
    if extension in {".tsv", ".txt"}:
        return "\t"
    header = first_lines(path, 1)
    header = header[0] if header else ""
    return "\t" if header.count("\t") > header.count(",") else ","


def read_table(path, filename=None, nrows=None):
    """Read a plain CSV/TSV/XLSX table as text (nothing is converted or guessed)."""

    filename = filename or Path(path).name
    if detect_kind(path, filename) == EXCEL_KIND:
        return pd.read_excel(path, dtype=str, keep_default_na=False, nrows=nrows)
    return pd.read_csv(path, sep=table_separator(path, filename), dtype=str,
                       keep_default_na=False, encoding="utf-8-sig", nrows=nrows)


def write_csv(table, path):
    table.to_csv(path, index=False)
    return Path(path)


# ----------------------------------------------------------------------------
# VCF
# ----------------------------------------------------------------------------

BASES = {"A", "C", "G", "T"}


def parse_info(text):
    info = {}
    if text in {"", "."}:
        return info
    for item in text.split(";"):
        key, _, value = item.partition("=")
        info[key] = value if _ else True
    return info


def numbers(text):
    """Split a comma-separated VCF value into floats; '.' becomes None."""

    values = []
    for item in str(text).split(","):
        item = item.strip()
        try:
            values.append(None if item in {"", "."} else float(item))
        except ValueError:
            values.append(None)
    return values


def pick(values, index):
    return values[index] if values is not None and index < len(values) else None


def allele_measurement(fields, info, alt_index, alt_count, use_info):
    """Return (frequency, alt_reads, depth, field_used) for one alternative allele.

    Only allele-specific values are used, so multi-allelic records are
    interpreted correctly or not at all. Order of preference:
    FORMAT/AD (one value per allele), FORMAT/AO+RO (FreeBayes), FORMAT/AF
    (one value per ALT), iVar's FORMAT/ALT_FREQ (single ALT only), and the
    LoFreq-style INFO/AF only for a single-sample or sites-only file (with
    several samples, INFO values describe the whole cohort, not one sample).
    """

    depth = pick(numbers(fields["DP"]), 0) if "DP" in fields else None
    if depth is None and use_info and "DP" in info:
        depth = pick(numbers(info["DP"]), 0)

    if "AD" in fields:
        counts = numbers(fields["AD"])
        if len(counts) == alt_count + 1 and None not in counts:
            alt_reads = counts[alt_index + 1]
            total = sum(counts)
            if total > 0:
                return alt_reads / total, alt_reads, depth if depth is not None else total, "FORMAT/AD"
    if "AO" in fields and "RO" in fields:
        alt_counts = numbers(fields["AO"])
        reference_count = pick(numbers(fields["RO"]), 0)
        if len(alt_counts) == alt_count and None not in alt_counts and reference_count is not None:
            total = reference_count + sum(alt_counts)
            if total > 0:
                alt_reads = alt_counts[alt_index]
                return alt_reads / total, alt_reads, depth if depth is not None else total, "FORMAT/AO,RO"
    if "AF" in fields:
        frequencies = numbers(fields["AF"])
        if len(frequencies) == alt_count and frequencies[alt_index] is not None:
            return frequencies[alt_index], None, depth, "FORMAT/AF"
    if alt_count == 1 and "ALT_FREQ" in fields:
        frequency = pick(numbers(fields["ALT_FREQ"]), 0)
        if frequency is not None:
            alt_reads = pick(numbers(fields.get("ALT_DP", ".")), 0)
            return frequency, alt_reads, depth, "FORMAT/ALT_FREQ"
    if use_info and "AF" in info and info["AF"] is not True:
        frequencies = numbers(info["AF"])
        if len(frequencies) == alt_count and frequencies[alt_index] is not None:
            return frequencies[alt_index], None, depth, "INFO/AF"
    return None, None, depth, None


def vcf_to_snv_table(path, fallback_sample_id):
    """Convert a VCF into standardized processed-SNV rows.

    Returns (table, details). `details` lists what was read and every record
    that was set aside, so nothing is dropped silently.
    """

    details = {"meta": {}, "contigs": {}, "samples": [], "records": 0, "snv_rows": 0,
                "set_aside": {}, "fields_used": set(), "sample_id_source": "VCF sample columns"}

    def set_aside(reason, count=1):
        details["set_aside"][reason] = details["set_aside"].get(reason, 0) + count

    rows = []
    sample_names = []
    header_seen = False
    with Path(path).open("r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("##"):
                key, _, value = line[2:].partition("=")
                if key == "contig" and "ID=" in value:
                    items = dict(part.split("=", 1) for part in value.strip("<>").split(",")
                                 if "=" in part)
                    details["contigs"][items.get("ID")] = items.get("length")
                elif key in {"fileformat", "source", "reference"}:
                    details["meta"][key] = value
                continue
            if line.startswith("#CHROM"):
                columns = line.lstrip("#").split("\t")
                if columns[:8] != ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
                    raise ValueError("The VCF header line does not have the standard eight columns.")
                sample_names = columns[9:] if len(columns) > 9 else []
                header_seen = True
                continue
            if not header_seen:
                raise ValueError("The VCF has records before its #CHROM header line.")
            fields = line.split("\t")
            if len(fields) < 8:
                set_aside("malformed line (fewer than 8 columns)")
                continue
            details["records"] += 1
            chrom, pos, record_id, ref, alt_text, qual, filter_text, info_text = fields[:8]
            ref = ref.upper()
            alts = [alt.upper() for alt in alt_text.split(",")]
            info = parse_info(info_text)
            if ref not in BASES:
                set_aside("not a single-nucleotide reference allele (indel, MNP, or N)")
                continue
            format_keys = fields[8].split(":") if len(fields) > 8 else []
            entries = fields[9:] if sample_names else [None]
            names = sample_names if sample_names else [fallback_sample_id]
            use_info = len(names) == 1
            for name, entry in zip(names, entries):
                values = dict(zip(format_keys, entry.split(":"))) if entry is not None else {}
                for index, alt in enumerate(alts):
                    if alt not in BASES:
                        set_aside("ALT is not a single nucleotide (indel, symbolic, or '*')")
                        continue
                    frequency, alt_reads, depth, used = allele_measurement(
                        values, info, index, len(alts), use_info)
                    if used is None:
                        set_aside("no allele-specific frequency or read counts (for example GT only)")
                        continue
                    if frequency <= 0 or (alt_reads is not None and alt_reads <= 0):
                        set_aside("sample entry has no alternate reads (not treated as a variant call)")
                        continue
                    details["fields_used"].add(used)
                    rows.append({
                        "Sample_ID": name,
                        "Reference_Genome": chrom,
                        "Position": pos,
                        "Ref": ref,
                        "Alt": alt,
                        "Alt_Frequency": round(float(frequency), 6),
                        "Read_Depth": "" if depth is None else int(depth),
                        "Alt_Depth": "" if alt_reads is None else int(alt_reads),
                        "Collection_Date": "",
                        "Location_Name": "",
                        "Source_VCF_ID": "" if record_id == "." else record_id,
                        "Source_VCF_Qual": "" if qual == "." else qual,
                        "Source_VCF_Filter": filter_text,
                        "Source_Frequency_Field": used,
                        "Source_Input_Format": VCF_FORMAT,
                    })
    if not header_seen:
        raise ValueError("This file has no #CHROM header line, so it is not a readable VCF.")
    details["samples"] = sample_names or [fallback_sample_id]
    if not sample_names:
        details["sample_id_source"] = "file name (the VCF has no sample columns)"
    details["snv_rows"] = len(rows)
    details["fields_used"] = sorted(details["fields_used"])
    columns = ["Sample_ID", "Reference_Genome", "Position", "Ref", "Alt", "Alt_Frequency",
               "Read_Depth", "Alt_Depth", "Collection_Date", "Location_Name", "Source_VCF_ID",
               "Source_VCF_Qual", "Source_VCF_Filter", "Source_Frequency_Field",
               "Source_Input_Format"]
    return pd.DataFrame(rows, columns=columns), details


def vcf_notes(details):
    """Plain-language notes about a VCF conversion, for the review page and warnings."""

    notes = []
    if details["fields_used"]:
        notes.append("Allele frequencies came from " + ", ".join(details["fields_used"]) + ".")
    for reason, count in sorted(details["set_aside"].items()):
        notes.append(f"{count} VCF allele entr{'y was' if count == 1 else 'ies were'} set aside: {reason}.")
    if details["sample_id_source"].startswith("file name"):
        notes.append("The VCF has no sample columns, so its sample ID is taken from the file name.")
    notes.append("VCF files do not carry collection dates or sites; add pooled-sample metadata for "
                 "temporal and geographic analyses. FILTER values are recorded, not applied.")
    return notes
