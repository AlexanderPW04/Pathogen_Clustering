"""Turn individual consensus sequences into reference-coordinate variant profiles.

The Sprint 2 alignment compares samples with each other, so its alignment
column numbers are not reference genome positions. Here each sample is
aligned with MAFFT to one verified reference sequence. The aligned reference
row is then walked column by column: a column where the reference has a base
advances the 1-based reference position; a column where the reference has a
gap is an insertion in the sample and has no reference position.

Each reference position receives one call code for the sample:

    R  callable base that matches the reference
    V  callable base that differs from the reference (a substitution)
    N  N or an IUPAC ambiguity code: no confident call
    D  gap inside the sequence (deletion relative to the reference)
    .  outside the sequenced span (partial sequences)

Only R and V positions are "callable". Insertions and deletions are listed
separately and are not compared across evidence sources.
"""

from pathlib import Path
from io import StringIO
from tempfile import TemporaryDirectory
import subprocess

from Bio import SeqIO
import pandas as pd

from align_sequences import run_mafft
from genomics.unified_schema import (
    CONSENSUS_CALL,
    INDIVIDUAL,
    make_sample_table,
    make_variant_table,
    variant_key,
)


CALL_REFERENCE = "R"
CALL_VARIANT = "V"
CALL_AMBIGUOUS = "N"
CALL_DELETION = "D"
CALL_OUTSIDE = "."

CALL_DESCRIPTIONS = {
    CALL_REFERENCE: "Reference base called",
    CALL_VARIANT: "Substitution called",
    CALL_AMBIGUOUS: "N or ambiguous base (not callable)",
    CALL_DELETION: "Deletion relative to reference (not callable)",
    CALL_OUTSIDE: "Outside sequenced span (not callable)",
}

# Prototype safeguards for accepting a reference mapping. They are not taken
# from the Fontenele et al. paper; they stop unrelated sequences (for example a
# different pathogen) from receiving false reference positions.
MINIMUM_CALLABLE_POSITIONS = 20
MINIMUM_MAPPING_IDENTITY = 0.90


def map_alignment_to_reference(reference_row, sample_row):
    """Convert one reference/sample alignment into reference-position calls."""

    reference_row = reference_row.upper()
    sample_row = sample_row.upper()
    if len(reference_row) != len(sample_row):
        raise ValueError("Aligned reference and sample rows have different lengths.")

    # The sequenced span runs from the first to the last aligned sample base.
    letter_columns = [index for index, base in enumerate(sample_row) if base != "-"]
    if letter_columns:
        first_column = letter_columns[0]
        last_column = letter_columns[-1]
    else:
        first_column = len(sample_row)
        last_column = -1

    calls = []
    substitutions = []
    insertions = []
    reference_position = 0
    pending_insertion = ""

    for column, (reference_base, sample_base) in enumerate(zip(reference_row, sample_row)):
        if reference_base == "-":
            # An inserted sample base has no reference coordinate.
            if sample_base != "-":
                pending_insertion += sample_base
            continue

        if pending_insertion:
            insertions.append((reference_position, pending_insertion))
            pending_insertion = ""

        reference_position += 1
        if column < first_column or column > last_column:
            calls.append(CALL_OUTSIDE)
        elif sample_base == "-":
            calls.append(CALL_DELETION)
        elif sample_base in "ACGT" and reference_base in "ACGT":
            if sample_base == reference_base:
                calls.append(CALL_REFERENCE)
            else:
                calls.append(CALL_VARIANT)
                substitutions.append((reference_position, reference_base, sample_base))
        else:
            calls.append(CALL_AMBIGUOUS)

    if pending_insertion:
        insertions.append((reference_position, pending_insertion))

    call_text = "".join(calls)
    return {
        "calls": call_text,
        "substitutions": substitutions,
        "insertions": insertions,
        "deletions": find_runs(call_text, CALL_DELETION),
    }


def find_runs(call_text, code):
    """Return (start, end) 1-based reference ranges where `code` repeats."""

    runs = []
    start = None
    for position, call in enumerate(call_text, start=1):
        if call == code and start is None:
            start = position
        elif call != code and start is not None:
            runs.append((start, position - 1))
            start = None
    if start is not None:
        runs.append((start, len(call_text)))
    return runs


def summarize_calls(call_text):
    """Count callable positions and identity to the reference."""

    reference_calls = call_text.count(CALL_REFERENCE)
    variant_calls = call_text.count(CALL_VARIANT)
    callable_positions = reference_calls + variant_calls
    identity = reference_calls / callable_positions if callable_positions else None

    span = [position for position, call in enumerate(call_text, start=1)
            if call != CALL_OUTSIDE]
    return {
        "callable_positions": callable_positions,
        "identity": identity,
        "ambiguous_positions": call_text.count(CALL_AMBIGUOUS),
        "deleted_positions": call_text.count(CALL_DELETION),
        "span_start": span[0] if span else None,
        "span_end": span[-1] if span else None,
        "callable_percent": callable_positions / len(call_text) * 100 if call_text else 0.0,
    }


def mapping_problem(summary):
    """Explain why a mapping cannot be trusted, or return None."""

    if summary["callable_positions"] < MINIMUM_CALLABLE_POSITIONS:
        return (
            f"only {summary['callable_positions']} callable positions "
            f"(at least {MINIMUM_CALLABLE_POSITIONS} required)"
        )
    if summary["identity"] < MINIMUM_MAPPING_IDENTITY:
        return (
            f"{summary['identity'] * 100:.1f}% identity to the reference "
            f"(at least {MINIMUM_MAPPING_IDENTITY * 100:.0f}% required)"
        )
    return None


def same_bases(first, second):
    """Compare ungapped sequences, treating every non-ACGT character as N."""

    def simplify(sequence):
        return "".join(base if base in "ACGT" else "N" for base in sequence.upper()
                       if base != "-")

    return simplify(first) == simplify(second)


def align_with_reference(reference, sequence, mafft_path):
    """Align one sample sequence to the reference and return both rows."""

    with TemporaryDirectory() as folder:
        input_path = Path(folder) / "reference_pair.fasta"
        input_path.write_text(
            f">REFERENCE\n{reference['Sequence']}\n>SAMPLE\n{sequence}\n",
            encoding="utf-8",
        )
        try:
            aligned_text = run_mafft(mafft_path, input_path)
        except (OSError, subprocess.CalledProcessError) as error:
            details = getattr(error, "stderr", "") or str(error)
            raise RuntimeError(f"MAFFT reference alignment failed: {details.strip()}") from error

    rows = {}
    for record in SeqIO.parse(StringIO(aligned_text), "fasta"):
        rows[str(record.id)] = "".join(str(record.seq).split()).upper()
    if set(rows) != {"REFERENCE", "SAMPLE"}:
        raise RuntimeError("MAFFT did not return the reference and sample rows.")
    if len(rows["REFERENCE"]) != len(rows["SAMPLE"]):
        raise RuntimeError("MAFFT returned unequal reference/sample row lengths.")
    if not same_bases(rows["REFERENCE"], reference["Sequence"]):
        raise RuntimeError("MAFFT changed the reference sequence.")
    if not same_bases(rows["SAMPLE"], sequence):
        raise RuntimeError("MAFFT changed the sample sequence.")
    return rows["REFERENCE"], rows["SAMPLE"]


def build_individual_profiles(data, reference, mafft_path, data_source):
    """Create shared variant and sample tables for individual sequences.

    `data` is the validated Sprint 2 table with each sample's original
    (unaligned) sequence. Samples whose mapping cannot be established keep
    their Sprint 2 genetic-distance results but get no variant profile.
    """

    variant_rows = []
    sample_rows = []
    indel_rows = []
    calls_by_sample = {}

    for sample in data.to_dict("records"):
        sample_id = str(sample["Sample_ID"])
        sequence = str(sample["Sequence"]).replace("-", "")
        reference_row, sample_row = align_with_reference(reference, sequence, mafft_path)
        mapping = map_alignment_to_reference(reference_row, sample_row)
        if len(mapping["calls"]) != reference["Length"]:
            raise RuntimeError("Reference positions were lost during alignment mapping.")

        summary = summarize_calls(mapping["calls"])
        problem = mapping_problem(summary)
        if problem is None:
            status = "Mapped to reference"
            calls_by_sample[sample_id] = mapping["calls"]
        else:
            status = f"Reference mapping not established: {problem}"

        sample_rows.append({
            "Sample_ID": sample_id,
            "Sample_Type": INDIVIDUAL,
            "Collection_Date": sample.get("Collection_Date"),
            "Location_Name": sample.get("Location_Name"),
            "Latitude": sample.get("Latitude"),
            "Longitude": sample.get("Longitude"),
            "Coordinate_Source": sample.get("Coordinate_Source"),
            "Geographic_Precision": sample.get("Geographic_Precision"),
            "Reference_ID": reference["Reference_ID"],
            "Genomic_Region": sample.get("Gene_Segment"),
            "Data_Source": data_source,
            "Coverage_Information": "Callable positions from consensus base calls",
            "Genome_Coverage_Percent": summary["callable_percent"],
            "Mean_Depth": pd.NA,
            "Quality_Status": (
                f"{summary['callable_positions']} callable positions; "
                f"{summary['ambiguous_positions']} ambiguous; "
                f"{summary['deleted_positions']} deleted"
            ),
            "Variant_Profile_Status": status,
            "Sequenced_Span": (
                f"{summary['span_start']}-{summary['span_end']}"
                if summary["span_start"] else ""
            ),
            "Reference_Identity_Percent": (
                summary["identity"] * 100 if summary["identity"] is not None else pd.NA
            ),
        })

        if problem is not None:
            continue

        for position, ref, alt in mapping["substitutions"]:
            variant_rows.append({
                "Variant_Key": variant_key(reference["Reference_ID"], position, ref, alt),
                "Reference_ID": reference["Reference_ID"],
                "Position": position,
                "Ref": ref,
                "Alt": alt,
                "Sample_ID": sample_id,
                "Sample_Type": INDIVIDUAL,
                "Observation_Type": CONSENSUS_CALL,
                "Alt_Frequency": pd.NA,
                "Read_Depth": pd.NA,
                "Alt_Depth": pd.NA,
                "Callable_Status": "Callable",
                "Quality_Status": "Consensus base call (no read-level quality supplied)",
                "Collection_Date": sample.get("Collection_Date"),
                "Location_Name": sample.get("Location_Name"),
                "Data_Source": data_source,
            })
        for after_position, bases in mapping["insertions"]:
            indel_rows.append({
                "Sample_ID": sample_id, "Event": "Insertion",
                "Reference_Positions": f"after {after_position}", "Bases": bases,
                "Note": "Not compared across sources",
            })
        for start, end in mapping["deletions"]:
            indel_rows.append({
                "Sample_ID": sample_id, "Event": "Deletion",
                "Reference_Positions": f"{start}-{end}", "Bases": "",
                "Note": "Positions are not callable for substitutions",
            })

    samples = make_sample_table(sample_rows)
    extra_columns = ["Sequenced_Span", "Reference_Identity_Percent"]
    for column in extra_columns:
        samples[column] = [row[column] for row in sample_rows]
    indels = pd.DataFrame(
        indel_rows, columns=["Sample_ID", "Event", "Reference_Positions", "Bases", "Note"]
    )
    return {
        "variants": make_variant_table(variant_rows),
        "samples": samples,
        "calls": calls_by_sample,
        "indels": indels,
    }


def individual_position_status(calls_by_sample, sample_id, position):
    """Return the call code description at one reference position."""

    calls = calls_by_sample.get(sample_id)
    if calls is None or position < 1 or position > len(calls):
        return None
    return calls[position - 1]
