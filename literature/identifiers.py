"""Recognize identifiers in sample metadata that can link a sample to literature."""

import re

import pandas as pd

from validate_inputs import is_missing_value

PATTERNS = {
    "BioSample": re.compile(r"^SAM[NED][A-Z]?\d+$"),
    "BioProject": re.compile(r"^PRJ[NED][A-Z]\d+$"),
    "SRA": re.compile(r"^[SED]R[RXSP]\d{5,}$"),
    "PMCID": re.compile(r"^PMC\d+$", re.IGNORECASE),
    "PMID": re.compile(r"^\d{1,9}$"),
    "DOI": re.compile(r"^10\.\d{4,9}/\S+$"),
    # INSDC nucleotide accessions such as MN908947.3, OQ123456, NC_045512.2.
    "Nucleotide": re.compile(r"^([A-Z]{1,2}\d{5,6}|[A-Z]{2}\d{8}|[A-Z]{2}_\d{6,9})(\.\d+)?$"),
}

# Column names (lowercase, without spaces or underscores) and the identifier
# type they hold. "Publication" columns may hold PMIDs, PMCIDs, DOIs, or
# project-specific article IDs (for example the synthetic SYNPUB-001).
COLUMN_TYPES = {
    "accession": "Nucleotide", "genbankaccession": "Nucleotide", "nucleotideaccession": "Nucleotide",
    "biosample": "BioSample", "biosampleaccession": "BioSample",
    "bioproject": "BioProject", "bioprojectaccession": "BioProject",
    "sra": "SRA", "sraaccession": "SRA", "srarun": "SRA", "run": "SRA", "runaccession": "SRA",
    "pmid": "PMID", "pubmedid": "PMID", "pmcid": "PMCID", "doi": "DOI",
    "publication": "Publication", "publicationid": "Publication", "reference": "Publication",
    "article": "Publication", "articleid": "Publication", "citation": "Publication",
}
# Other names a paper might use for the same sample.
ALIAS_COLUMNS = {"isolate", "strain", "samplename", "isolatename", "siteid", "site"}


def column_key(name):
    return re.sub(r"[\s_\-]", "", str(name)).casefold()


def is_standard(kind, value):
    """True when a value has the official format for its identifier type."""

    pattern = PATTERNS.get(kind)
    return pattern is not None and pattern.match(str(value).strip()) is not None


def classify_value(value):
    """Return the identifier type of one value, or None."""

    text = str(value).strip()
    for kind in ["BioSample", "BioProject", "SRA", "PMCID", "DOI", "Nucleotide"]:
        if PATTERNS[kind].match(text):
            return kind
    return None


def sample_identifiers(row):
    """All linkable identifiers and aliases for one sample row.

    Returns (identifiers, aliases): identifiers is a list of (type, value,
    column); aliases are names that may label the sample inside article text.
    """

    identifiers = []
    aliases = {str(row["Sample_ID"]).strip()}
    for column, value in row.items():
        if value is None or (not isinstance(value, str) and pd.isna(value)) or is_missing_value(value):
            continue
        text = str(value).strip()
        if not text:
            continue
        key = column_key(column)
        declared = COLUMN_TYPES.get(key)
        if key in ALIAS_COLUMNS:
            aliases.add(text)
        if declared == "Publication":
            kind = classify_value(text) or ("PMID" if PATTERNS["PMID"].match(text) else "Publication")
        elif declared == "PMID":
            kind = "PMID" if PATTERNS["PMID"].match(text) else None
        elif declared:
            # A value in a declared column keeps that type even when it is not
            # an NCBI-format identifier (for example a project-specific
            # accession); live NCBI lookups skip such values.
            kind = classify_value(text) or declared
        elif key in {"sampleid", "sample"}:
            kind = classify_value(text)
        else:
            continue
        if kind:
            identifiers.append((kind, text, str(column)))
            if kind in {"Nucleotide", "BioSample", "SRA"}:
                aliases.add(text)
                aliases.add(text.split(".")[0])  # accession without its version
    unique = list(dict.fromkeys(identifiers))
    return unique, sorted(alias for alias in aliases if alias)


def parse_identifiers(text, limit=10):
    """Publication or study identifiers typed by the user: [(type, value)].

    PMIDs, PMCIDs, DOIs, BioProject, BioSample, SRA, and nucleotide accessions
    are recognized; anything else is kept as a project-specific publication ID.
    """

    found = []
    for value in re.split(r"[\s,;]+", str(text or "")):
        value = re.sub(r"^(pmid|pmcid|doi|bioproject):", "", value.strip().rstrip("."), flags=re.IGNORECASE)
        if not value or value.casefold() in {"pmid", "pmcid", "doi", "bioproject", "and"}:
            continue  # labels, not identifiers
        kind = classify_value(value) or ("PMID" if PATTERNS["PMID"].match(value) else "Publication")
        if (kind, value) not in found:
            found.append((kind, value))
    return found[:limit]


def mentions(text, alias):
    """True when an alias appears in text as a whole identifier, not inside a longer one."""

    pattern = r"(?<![A-Za-z0-9_.\-])" + re.escape(alias) + r"(?![A-Za-z0-9_])"
    return re.search(pattern, text, flags=re.IGNORECASE) is not None


def linkable_columns(frame):
    """Columns that hold identifiers or alternative names for samples."""

    return [column for column in frame.columns
            if column_key(column) in COLUMN_TYPES or column_key(column) in ALIAS_COLUMNS]
