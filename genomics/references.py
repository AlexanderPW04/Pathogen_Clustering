"""Load reference genomes and decide whether two references are compatible.

Genomic positions can only be compared when they refer to the same reference
sequence. Two different accession labels (for example MN908947.3 and
NC_045512.2) are accepted as the same reference only when both sequences were
supplied and are exactly identical. Labels alone are never trusted.
"""

import hashlib
from pathlib import Path

from Bio import SeqIO


REFERENCE_CHARACTERS = set("ACGTN")


def split_reference_version(reference_id):
    """Split 'NC_045512.2' into ('NC_045512', '2'); the version may be blank."""

    reference_id = str(reference_id).strip()
    base, period, version = reference_id.rpartition(".")
    if period and base and version.isdigit():
        return base, version
    return reference_id, ""


def sequence_checksum(sequence):
    """Return a SHA-256 fingerprint so identical sequences can be recognized."""

    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def load_reference_fasta(path):
    """Read exactly one reference sequence and describe it."""

    path = Path(path)
    with path.open(encoding="utf-8") as fasta_file:
        records = list(SeqIO.parse(fasta_file, "fasta"))

    if len(records) != 1:
        raise ValueError(
            f"A reference FASTA must contain exactly one sequence; found {len(records)}."
        )

    record = records[0]
    sequence = "".join(str(record.seq).split()).upper()
    if not sequence:
        raise ValueError("The reference sequence is empty.")

    invalid = set(sequence) - REFERENCE_CHARACTERS
    if invalid:
        raise ValueError(
            "The reference may contain only A, C, G, T, or N; found "
            + "".join(sorted(invalid))
            + "."
        )

    reference_id = str(record.id).strip()
    base, version = split_reference_version(reference_id)
    return {
        "Reference_ID": reference_id,
        "Accession": base,
        "Version": version,
        "Description": str(record.description),
        "Sequence": sequence,
        "Length": len(sequence),
        "Checksum": sequence_checksum(sequence),
        "Source_File": path.name,
    }


def build_reference_registry(references):
    """Index reference descriptions by their identifier."""

    registry = {}
    for reference in references:
        reference_id = reference["Reference_ID"]
        if reference_id in registry:
            if registry[reference_id]["Checksum"] != reference["Checksum"]:
                raise ValueError(
                    f"Two different sequences use the reference ID {reference_id}."
                )
            continue
        registry[reference_id] = reference
    return registry


def references_compatible(reference_id_1, reference_id_2, registry=None):
    """Return (True/False, reason) for comparing positions on two references."""

    first = str(reference_id_1).strip()
    second = str(reference_id_2).strip()
    if not first or not second:
        return False, "A reference identifier is missing."

    if first == second:
        return True, "Same reference identifier."

    registry = registry or {}
    if first in registry and second in registry:
        if registry[first]["Checksum"] == registry[second]["Checksum"]:
            return True, "Different identifiers, but the supplied sequences are identical."
        return False, "The supplied reference sequences differ."

    first_base, first_version = split_reference_version(first)
    second_base, second_version = split_reference_version(second)
    if first_base == second_base:
        return False, (
            f"Different versions of {first_base} ({first_version or 'none'} versus "
            f"{second_version or 'none'}); sequence identity was not verified."
        )
    return False, (
        "Different reference identifiers; their sequence identity was not "
        "verified with supplied reference sequences."
    )


def reference_base(reference, position):
    """Return the reference nucleotide at a 1-based position."""

    if position < 1 or position > reference["Length"]:
        raise ValueError(
            f"Position {position} is outside {reference['Reference_ID']} "
            f"(1 to {reference['Length']})."
        )
    return reference["Sequence"][position - 1]
