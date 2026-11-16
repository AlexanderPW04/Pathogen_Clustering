"""Sprint 5, Phase 4: load a verified genome annotation for SNV consequence prediction.

Supported sources are an NCBI-style GenBank flat file (parsed with Biopython,
already a project dependency) or a GFF3 file (parsed here). Gene and coding
boundaries come only from these files; nothing is inferred from the sequence.

Coordinate conventions
----------------------
- GenBank locations (``join(266..13468,13468..21555)``) and GFF3 columns 4-5
  are 1-based and inclusive. They are stored that way in every table here.
- Biopython's ``FeatureLocation`` is 0-based with an exclusive end, so a
  part is converted as ``(int(part.start) + 1, int(part.end))``.
- Sequences are Python strings, so the base at 1-based position ``p`` is
  ``sequence[p - 1]``.
- Offsets inside a coding sequence ("coding offsets") are 0-based and count
  from the first base of the start codon on the coding strand: codon number
  ``offset // 3 + 1`` and codon position ``offset % 3 + 1``.

Each coding feature (CDS) is checked before it is used. A CDS is marked
unsupported, with the reason, rather than translated when its structure
cannot be interpreted reliably: partial ends, a non-zero starting phase,
pseudogenes, translation exceptions (for example RNA editing or
selenocysteine), parts out of genomic order (origin-spanning or
trans-spliced), overlapping parts without a documented ribosomal slippage,
inconsistent GFF3 phases, incomplete codons, internal stop codons, or a
translation that disagrees with the record's own protein sequence.

Genetic code (final-polish 1C)
------------------------------
The translation table of each CDS is never assumed silently. Its source is
recorded in ``Translation_Table_Source``:

- stated in the annotation (GenBank ``/transl_table``, GFF3 ``transl_table``);
- not stated, but the standard code (table 1) reproduces the record's own
  ``/translation`` exactly, which verifies it for that CDS;
- confirmed by the user (``genetic_code`` setting) for CDS that state none;
- otherwise not verified: the CDS is marked unsupported and no amino-acid
  consequence is predicted for it.
"""

import hashlib
import io
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqFeature import AfterPosition, BeforePosition

from genomics.references import split_reference_version


GENBANK_FORMAT = "GenBank flat file"
GFF3_FORMAT = "GFF3"
SUPPORTED = "Supported"
TABLE_STATED = "Stated in the annotation"
TABLE_NOT_STATED = "Not stated in the annotation"
TABLE_VERIFIED = "Not stated; the standard code (table 1) reproduces the record's own protein sequence"
TABLE_CONFIRMED = "Confirmed by the user (Advanced settings: genetic code)"
TABLE_UNVERIFIED_DETAIL = (
    "The annotation states no genetic code (transl_table) and has no protein sequence to check one against, "
    "so amino-acid consequences are not predicted for this coding region. If the organism's genetic code is "
    "known, choose it under Advanced settings ('Genetic code for coding regions that state none').")
UNSUPPORTED = "Unsupported"
NUCLEOTIDES = set("ACGT")
COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVN", "TGCAYRMKSWVHDBN")
SLIPPAGE = "ribosomal slippage"
# Feature types kept as named noncoding context (UTRs, structured RNA elements, genes).
CONTEXT_TYPES = {"5'UTR", "3'UTR", "five_prime_UTR", "three_prime_UTR", "stem_loop", "ncRNA",
                 "misc_RNA", "regulatory", "repeat_region", "gene"}
PEPTIDE_TYPES = {"mat_peptide", "mature_protein_region_of_CDS", "sig_peptide", "signal_peptide_region_of_CDS"}

SOURCE_COLUMNS = ["Reference_ID", "Accession", "Version", "Organism", "Segment", "Reference_Length",
                  "Sequence_Source", "Sequence_Check", "Gene", "Gene_ID", "Locus_Tag", "CDS_ID", "Product",
                  "Strand", "Parts_1_Based", "Phases", "Coding_Length", "Codons", "Translation_Table",
                  "Translation_Table_Source", "Structure", "Mature_Peptides", "Status", "Status_Detail", "Verification",
                  "Annotation_Format", "Annotation_File"]


def detect_annotation_format(path):
    """Return GENBANK_FORMAT or GFF3_FORMAT from the first non-empty line, or None."""

    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith("##gff-version") and text.split()[-1].startswith("3"):
                return GFF3_FORMAT
            if text.startswith("LOCUS "):
                return GENBANK_FORMAT
            return None
    return None


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def clean_sequence(text):
    return "".join(str(text).split()).upper().replace("U", "T")


def new_record(reference_id, organism="", segment="", length=None, sequence=None, source=""):
    base, version = split_reference_version(reference_id)
    return {"Reference_ID": reference_id, "Accession": base, "Version": version,
            "Organism": organism or "Not stated", "Segment": segment, "Length": length,
            "Sequence": sequence, "Sequence_Source": source if sequence else "",
            "Sequence_Check": "", "Compatible": True}


def new_cds(reference_id, strand, parts, **values):
    cds = {"Reference_ID": reference_id, "Strand": strand, "Parts": parts, "Phases": [],
           "Gene": "", "Gene_ID": "", "Locus_Tag": "", "CDS_ID": "", "Product": "",
           "Translation_Table": 1, "Translation_Table_Source": TABLE_NOT_STATED, "Codon_Start": 1, "Exceptions": [], "Slippage": False,
           "Record_Translation": None, "Partial": False, "Pseudo": False, "Transl_Except": False,
           "Status": SUPPORTED, "Status_Detail": "", "Verification": "", "Coding_Sequence": None,
           "Start_Codon_Verified": False, "Mature_Peptides": []}
    cds.update(values)
    return cds


def genetic_code_choices():
    """(table id, NCBI name) for every translation table Biopython knows, for the settings menu."""

    return [(table_id, table.names[0]) for table_id, table in sorted(unambiguous_dna_by_id.items())]


def parse_genetic_code(text):
    """The user's confirmed genetic code: None when not confirmed, else a known NCBI table id."""

    text = str(text or "").strip()
    if not text:
        return None
    if not text.isdigit() or int(text) not in unambiguous_dna_by_id:
        raise ValueError("Choose a genetic code from the list (an NCBI translation table), or leave it as "
                         "'Not confirmed'.")
    return int(text)


def translation_table(stated):
    """Translation_Table and its source from the annotation's own value ('' when absent)."""

    text = str(stated or "").strip()
    if text.isdigit():
        return {"Translation_Table": int(text), "Translation_Table_Source": f"{TABLE_STATED} (transl_table={text})"}
    if text:  # present but unreadable: keep it visible; structure_problem refuses unknown tables
        return {"Translation_Table": -1, "Translation_Table_Source": f"{TABLE_STATED} (unreadable value '{text[:20]}')"}
    return {"Translation_Table": 1, "Translation_Table_Source": TABLE_NOT_STATED}


def gene_id_from(xrefs):
    for value in xrefs:
        if str(value).startswith("GeneID:"):
            return str(value).split(":", 1)[1]
    return ""


# ---------------------------------------------------------------------------
# GenBank
# ---------------------------------------------------------------------------

def qualifier(feature, name, default=""):
    values = feature.qualifiers.get(name)
    return str(values[0]).strip() if values else default


def location_parts(location):
    """1-based inclusive (start, end) parts in transcription order, plus a partial flag."""

    parts, partial = [], False
    for part in location.parts:  # Biopython keeps join() parts in biological order
        if isinstance(part.start, (BeforePosition, AfterPosition)) or isinstance(
                part.end, (BeforePosition, AfterPosition)):
            partial = True
        parts.append((int(part.start) + 1, int(part.end)))
    return parts, partial


def load_genbank(path):
    """Parse every record of a GenBank file into the annotation model."""

    annotation = {"format": GENBANK_FORMAT, "source_file": Path(path).name, "records": {}, "cds": [],
                  "peptides": [], "context": [], "notes": [], "problems": []}
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if not records:
        raise ValueError("The GenBank file contains no records.")
    for record in records:
        reference_id = str(record.id).strip()
        try:
            sequence = clean_sequence(record.seq)
        except Exception:  # noqa: BLE001 - an undefined sequence (CONTIG records) has no bases
            sequence = None
        source = next((feature for feature in record.features if feature.type == "source"), None)
        organism = str(record.annotations.get("organism", "")) or (qualifier(source, "organism") if source else "")
        segment = qualifier(source, "segment") if source else ""
        if reference_id in annotation["records"]:
            raise ValueError(f"The GenBank file contains {reference_id} more than once.")
        annotation["records"][reference_id] = new_record(
            reference_id, organism, segment, len(record.seq), sequence or None, "GenBank record sequence")
        for feature in record.features:
            strand = "-" if feature.location.strand == -1 else "+"
            parts, partial = location_parts(feature.location)
            if feature.type == "CDS":
                exceptions = [str(value).strip() for value in feature.qualifiers.get("exception", [])]
                slippage = "ribosomal_slippage" in feature.qualifiers or any(
                    value.casefold() == SLIPPAGE for value in exceptions)
                xrefs = feature.qualifiers.get("db_xref", [])
                annotation["cds"].append(new_cds(
                    reference_id, strand, parts, Gene=qualifier(feature, "gene"),
                    Gene_ID=gene_id_from(xrefs), Locus_Tag=qualifier(feature, "locus_tag"),
                    CDS_ID=qualifier(feature, "protein_id"), Product=qualifier(feature, "product"),
                    **translation_table(qualifier(feature, "transl_table")),
                    Codon_Start=int(qualifier(feature, "codon_start", "1") or 1),
                    Exceptions=[value for value in exceptions if value.casefold() != SLIPPAGE],
                    Slippage=slippage, Record_Translation=qualifier(feature, "translation") or None,
                    Partial=partial or "partial" in feature.qualifiers,
                    Pseudo="pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers,
                    Transl_Except="transl_except" in feature.qualifiers))
            elif feature.type in PEPTIDE_TYPES:
                annotation["peptides"].append({"Reference_ID": reference_id, "Strand": strand, "Parts": parts,
                                               "Product": qualifier(feature, "product"),
                                               "Protein_ID": qualifier(feature, "protein_id"),
                                               "Parent": ""})
            elif feature.type in CONTEXT_TYPES:
                annotation["context"].append({"Reference_ID": reference_id, "Type": feature.type,
                                              "Name": qualifier(feature, "gene") or qualifier(feature, "product")
                                              or feature.type, "Start": min(p[0] for p in parts),
                                              "End": max(p[1] for p in parts), "Strand": strand})
    return annotation


# ---------------------------------------------------------------------------
# GFF3
# ---------------------------------------------------------------------------

def gff_attributes(text):
    values = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            values[key.strip()] = unquote(value.strip())
    return values


def load_gff3(path):
    """Parse a GFF3 file (and its optional ##FASTA section) into the annotation model."""

    annotation = {"format": GFF3_FORMAT, "source_file": Path(path).name, "records": {}, "cds": [],
                  "peptides": [], "context": [], "notes": [], "problems": []}
    regions, species, cds_lines, genes, peptides = {}, "", {}, {}, {}
    fasta_lines = None
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for number, line in enumerate(handle, start=1):
            if fasta_lines is not None:
                fasta_lines.append(line)
                continue
            text = line.rstrip("\n\r")
            if not text.strip():
                continue
            if text.startswith("##FASTA"):
                fasta_lines = []
                continue
            if text.startswith("##sequence-region"):
                fields = text.split()
                if len(fields) == 4 and fields[3].isdigit():
                    regions.setdefault(fields[1], {})["Length"] = int(fields[3])
                continue
            if text.startswith("##species"):
                species = text.split(None, 1)[1].strip() if len(text.split()) > 1 else ""
                continue
            if text.startswith("#"):
                continue
            fields = text.split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF3 line {number} does not have 9 tab-separated columns.")
            seqid, _, kind, start, end, _, strand, phase, attributes = fields
            if not (start.isdigit() and end.isdigit()) or int(start) < 1 or int(end) < int(start):
                raise ValueError(f"GFF3 line {number} has invalid coordinates {start}..{end}.")
            start, end = int(start), int(end)
            attrs = gff_attributes(attributes)
            seqid = unquote(seqid)
            if kind in {"region", "source"}:
                info = regions.setdefault(seqid, {})
                info.setdefault("Organism", attrs.get("organism", ""))
                info.setdefault("Segment", attrs.get("segment", ""))
                info.setdefault("Taxon", next((x.split(":", 1)[1] for x in attrs.get("Dbxref", "").split(",")
                                               if x.startswith("taxon:")), ""))
                info.setdefault("Length", end if start == 1 else None)
            elif kind == "CDS":
                key = (seqid, attrs.get("ID") or f"line-{number}")
                entry = cds_lines.setdefault(key, {"strand": strand, "attrs": attrs, "parts": []})
                if entry["strand"] != strand:
                    entry["mixed_strand"] = True
                entry["parts"].append((start, end, phase))
            elif kind in PEPTIDE_TYPES:  # multi-line peptides (across a frameshift) share one ID
                key = (seqid, attrs.get("ID") or f"line-{number}")
                peptide = peptides.setdefault(key, {"Reference_ID": seqid, "Strand": strand, "Parts": [],
                                                    "Product": attrs.get("product", attrs.get("Name", "")),
                                                    "Protein_ID": attrs.get("protein_id", ""),
                                                    "Parent": attrs.get("Parent", "")})
                peptide["Parts"].append((start, end))
            if kind == "gene" and attrs.get("ID"):
                genes[attrs["ID"]] = attrs
            if kind in CONTEXT_TYPES:
                annotation["context"].append({"Reference_ID": seqid, "Type": kind,
                                              "Name": attrs.get("gene") or attrs.get("Name") or attrs.get("product")
                                              or kind, "Start": start, "End": end, "Strand": strand})
    for peptide in peptides.values():
        peptide["Parts"].sort(key=lambda part: part[0], reverse=(peptide["Strand"] == "-"))
        annotation["peptides"].append(peptide)
    sequences = {}
    if fasta_lines:
        for record in SeqIO.parse(io.StringIO("".join(fasta_lines)), "fasta"):
            sequences[str(record.id).strip()] = clean_sequence(record.seq)
    taxon_note = ""
    if species:
        taxon_note = species.rsplit("id=", 1)[-1] if "id=" in species else species
    seqids = list(dict.fromkeys(list(regions) + [key[0] for key in cds_lines] + list(sequences)))
    for seqid in seqids:
        info = regions.get(seqid, {})
        organism = info.get("Organism") or (f"NCBI taxon {info.get('Taxon') or taxon_note}"
                                            if info.get("Taxon") or taxon_note else "")
        sequence = sequences.get(seqid)
        length = info.get("Length") or (len(sequence) if sequence else None)
        annotation["records"][seqid] = new_record(seqid, organism, info.get("Segment", ""), length, sequence,
                                                  "GFF3 ##FASTA section")
        if sequence and info.get("Length") and len(sequence) != info["Length"]:
            annotation["records"][seqid]["Compatible"] = False
            annotation["records"][seqid]["Sequence_Check"] = (
                f"The ##FASTA sequence has {len(sequence)} bases but ##sequence-region declares {info['Length']}.")
    for (seqid, _), entry in cds_lines.items():
        attrs, strand = entry["attrs"], entry["strand"]
        lines = sorted(entry["parts"], key=lambda item: item[0], reverse=(strand == "-"))
        parts = [(start, end) for start, end, _ in lines]
        phases = [int(phase) if phase in {"0", "1", "2"} else None for _, _, phase in lines]
        exceptions = [value.strip() for value in attrs.get("exception", "").split(",") if value.strip()]
        parent = genes.get(attrs.get("Parent", ""), {})
        xrefs = (attrs.get("Dbxref", "") + "," + parent.get("Dbxref", "")).split(",")
        cds = new_cds(seqid, strand, parts, Phases=phases,
                      Gene=attrs.get("gene") or parent.get("gene") or parent.get("Name", ""),
                      Gene_ID=gene_id_from(xrefs),
                      Locus_Tag=attrs.get("locus_tag") or parent.get("locus_tag", ""),
                      CDS_ID=attrs.get("protein_id") or attrs.get("ID", "").removeprefix("cds-"),
                      Product=attrs.get("product", ""),
                      **translation_table(attrs.get("transl_table", "")),
                      Codon_Start=(phases[0] + 1) if phases and phases[0] is not None else 1,
                      Exceptions=[value for value in exceptions if value.casefold() != SLIPPAGE],
                      Slippage=any(value.casefold() == SLIPPAGE for value in exceptions),
                      Partial=(attrs.get("partial", "").casefold() == "true" or "start_range" in attrs
                               or "end_range" in attrs),
                      Pseudo=attrs.get("pseudo", "").casefold() == "true",
                      Transl_Except="transl_except" in attrs)
        if entry.get("mixed_strand") or strand not in {"+", "-"}:
            cds["Status"], cds["Status_Detail"] = UNSUPPORTED, "The CDS parts do not share one strand."
        if None in phases:
            cds["Status"], cds["Status_Detail"] = UNSUPPORTED, "A CDS line has no valid phase (0, 1, or 2)."
        annotation["cds"].append(cds)
    return annotation


# ---------------------------------------------------------------------------
# Reference matching and coding-feature verification
# ---------------------------------------------------------------------------

def attach_reference(annotation, reference):
    """Use the uploaded reference FASTA where the annotation describes the same accession.version.

    A record that carries its own sequence must be identical to the FASTA; a
    record without a sequence (GFF3 without ##FASTA) borrows the FASTA only
    when the declared length matches. Labels alone never make two references
    compatible.
    """

    if reference is None:
        return
    record = annotation["records"].get(reference["Reference_ID"])
    if record is None:
        annotation["notes"].append(
            f"The reference FASTA {reference['Reference_ID']} is not described by the annotation file (which "
            "describes " + ", ".join(sorted(annotation["records"])) + "), so its variants are not annotated.")
        return
    if record["Sequence"]:
        checksum = hashlib.sha256(record["Sequence"].encode("ascii", "replace")).hexdigest()
        if checksum == reference["Checksum"]:
            record["Sequence_Check"] = "Identical to the uploaded reference FASTA"
        else:
            record["Compatible"] = False
            record["Sequence_Check"] = ("The annotation's sequence differs from the uploaded reference FASTA with "
                                        "the same identifier, so its coordinates are not used.")
    elif record["Length"] and record["Length"] != reference["Length"]:
        record["Compatible"] = False
        record["Sequence_Check"] = (f"The annotation declares {record['Length']} bases but the reference FASTA has "
                                    f"{reference['Length']}, so they are different sequences.")
    else:
        record["Sequence"] = reference["Sequence"]
        record["Length"] = reference["Length"]
        record["Sequence_Source"] = "Uploaded reference FASTA (same accession.version and length)"
        record["Sequence_Check"] = "Annotation has no sequence; the matching reference FASTA is used"


def coding_offsets(cds, position):
    """0-based coding offsets of a 1-based genomic position (several where parts overlap)."""

    offsets, before = [], 0
    for start, end in cds["Parts"]:
        if start <= position <= end:
            offsets.append(before + (position - start if cds["Strand"] == "+" else end - position))
        before += end - start + 1
    return offsets


def build_coding_sequence(cds, sequence):
    pieces = [sequence[start - 1:end] for start, end in cds["Parts"]]
    if cds["Strand"] == "-":
        pieces = [reverse_complement(piece) for piece in pieces]
    return "".join(pieces)


def translate_codon(codon, table):
    """One-letter amino acid, '*' for a stop, or None when the codon is not plain A/C/G/T."""

    if len(codon) != 3 or not set(codon) <= NUCLEOTIDES:
        return None
    codons = unambiguous_dna_by_id[table]
    return "*" if codon in codons.stop_codons else codons.forward_table[codon]


def structure_label(cds):
    if len(cds["Parts"]) == 1:
        return "Single interval"
    if cds["Slippage"]:
        return f"Joined ({len(cds['Parts'])} parts): programmed ribosomal frameshift"
    return f"Joined ({len(cds['Parts'])} parts)"


def mark(cds, detail):
    cds["Status"], cds["Status_Detail"] = UNSUPPORTED, detail


def structure_problem(cds, record):
    """The first reason this CDS cannot be interpreted reliably, or ''."""

    if cds["Pseudo"]:
        return "Annotated as a pseudogene; no protein is predicted."
    if cds["Partial"]:
        return "Partial coding sequence (an incomplete 5' or 3' end)."
    if cds["Codon_Start"] != 1:
        return f"Translation starts at codon position {cds['Codon_Start']} (a partial 5' end)."
    if cds["Transl_Except"]:
        return "The record lists a translation exception (/transl_except), so plain translation would be wrong."
    if cds["Exceptions"]:
        return "The record lists an exception (" + "; ".join(cds["Exceptions"]) + ") that changes how the " \
               "genome is read, so codons cannot be derived from genomic coordinates alone."
    if cds["Translation_Table"] not in unambiguous_dna_by_id:
        return f"Unknown translation table {cds['Translation_Table']}."
    for (start, end), (next_start, next_end) in zip(cds["Parts"], cds["Parts"][1:]):
        ordered = next_start > start if cds["Strand"] == "+" else next_end < end
        if not ordered:
            return "The CDS parts are not in genomic order (origin-spanning or trans-spliced); not supported."
        overlap = next_start <= end if cds["Strand"] == "+" else next_end >= start
        if overlap and not cds["Slippage"]:
            return "CDS parts overlap without a documented ribosomal slippage."
    if record is None or not record["Compatible"]:
        return "The annotation's reference sequence is not compatible with the uploaded reference."
    if not record["Sequence"]:
        return "No reference sequence is available to derive codons."
    if any(end > len(record["Sequence"]) for _, end in cds["Parts"]):
        return "The CDS extends beyond the end of the reference sequence."
    if cds["Phases"]:  # GFF3 phases must agree with the part lengths
        before = 0
        for (start, end), phase in zip(cds["Parts"], cds["Phases"]):
            if phase != (3 - before % 3) % 3:
                return "The GFF3 phase column disagrees with the lengths of the earlier CDS parts."
            before += end - start + 1
    return ""


def verify_cds(cds, record):
    """Check one CDS and set Status, Status_Detail, Verification, and its coding sequence."""

    if cds["Status"] == UNSUPPORTED:
        return
    problem = structure_problem(cds, record)
    if problem:
        mark(cds, problem)
        return
    if cds["Translation_Table_Source"] == TABLE_NOT_STATED and not cds["Record_Translation"]:
        mark(cds, TABLE_UNVERIFIED_DETAIL)
        return
    coding = build_coding_sequence(cds, record["Sequence"])
    table = cds["Translation_Table"]
    codon_table = unambiguous_dna_by_id[table]
    if len(coding) % 3:
        mark(cds, f"The coding length ({len(coding)} bases) is not a whole number of codons.")
        return
    protein = [translate_codon(coding[index:index + 3], table) for index in range(0, len(coding), 3)]
    if protein[-1] != "*":
        mark(cds, "The CDS does not end with a stop codon in the reference sequence.")
        return
    internal = [index + 1 for index, amino in enumerate(protein[:-1]) if amino == "*"]
    if internal:
        mark(cds, f"Internal stop codon(s) at codon {internal[0]} in the reference translation (for example a "
                  "frameshift or editing the annotation does not describe).")
        return
    cds["Start_Codon_Verified"] = coding[:3] in codon_table.start_codons
    ambiguous = sum(amino is None for amino in protein)
    if cds["Record_Translation"]:
        predicted = "".join(amino or "X" for amino in protein[:-1])
        expected = cds["Record_Translation"]
        if cds["Start_Codon_Verified"] and expected[:1] == "M":
            predicted = "M" + predicted[1:]  # an alternative start codon is translated as methionine
        if predicted != expected:
            mismatch = next((i + 1 for i, (a, b) in enumerate(zip(predicted, expected)) if a != b),
                            min(len(predicted), len(expected)) + 1)
            mark(cds, f"Translating the reference disagrees with the record's protein at residue {mismatch}; "
                      "the coding structure is not interpreted.")
            return
        cds["Verification"] = ("Translation of the reference coding sequence matches the record's /translation "
                               f"({len(expected)} amino acids)")
        if cds["Translation_Table_Source"] == TABLE_NOT_STATED:
            cds["Translation_Table_Source"] = TABLE_VERIFIED
    else:
        cds["Verification"] = ("Internal checks only (GFF3 has no protein sequence): whole codons, one terminal "
                               "stop, no internal stop" + (", start codon present" if cds["Start_Codon_Verified"]
                                                           else ", first codon is not a start codon"))
    if ambiguous:
        cds["Verification"] += f"; {ambiguous} codon(s) contain ambiguous reference bases"
    cds["Coding_Sequence"] = coding


def map_peptides(annotation):
    """Attach mature peptides that lie in frame inside each supported CDS (codon ranges, 1-based)."""

    for cds in annotation["cds"]:
        if cds["Status"] != SUPPORTED:
            continue
        seen = set()
        for peptide in annotation["peptides"]:
            if peptide["Reference_ID"] != cds["Reference_ID"] or peptide["Strand"] != cds["Strand"]:
                continue
            if peptide["Parent"] and cds["CDS_ID"] and not peptide["Parent"].endswith(cds["CDS_ID"]):
                continue
            span = peptide_span(cds, peptide["Parts"])
            if span is None:
                continue
            key = (peptide["Product"], span)
            if key not in seen:
                seen.add(key)
                cds["Mature_Peptides"].append({"Product": peptide["Product"], "Protein_ID": peptide["Protein_ID"],
                                               "First_Codon": span[0], "Last_Codon": span[1]})
        cds["Mature_Peptides"].sort(key=lambda item: item["First_Codon"])


def peptide_span(cds, parts):
    """(first codon, last codon) when the peptide parts map contiguously and in frame, else None."""

    expected = None
    first = None
    for start, end in parts:
        five, three = (start, end) if cds["Strand"] == "+" else (end, start)
        found = None
        for offset_five in coding_offsets(cds, five):
            for offset_three in coding_offsets(cds, three):
                if offset_three - offset_five == end - start and (expected is None or offset_five == expected):
                    found = (offset_five, offset_three)
        if found is None:
            return None
        first = found[0] if first is None else first
        expected = found[1] + 1
    last = expected - 1
    if first % 3 or (last + 1) % 3:
        return None
    return first // 3 + 1, (last + 1) // 3


def finalize_annotation(annotation, reference=None, genetic_code=None):
    """Match the uploaded reference, verify every CDS, and map mature peptides.

    `genetic_code` is the user's confirmed NCBI translation table for CDS whose annotation states none
    (None: not confirmed). A table stated in the annotation is never overridden.
    """

    attach_reference(annotation, reference)
    if genetic_code is not None:
        for cds in annotation["cds"]:
            if cds["Translation_Table_Source"] == TABLE_NOT_STATED:
                cds["Translation_Table"] = int(genetic_code)
                cds["Translation_Table_Source"] = f"{TABLE_CONFIRMED} (table {int(genetic_code)})"
    for cds in annotation["cds"]:
        verify_cds(cds, annotation["records"].get(cds["Reference_ID"]))
    unverified = sum(cds["Status_Detail"] == TABLE_UNVERIFIED_DETAIL for cds in annotation["cds"])
    if unverified:
        annotation["notes"].append(
            f"{unverified} coding region(s) state no genetic code and carry no protein sequence to verify one, so "
            "their amino-acid consequences are not predicted (the standard code is never assumed). Choose the "
            "organism's genetic code under Advanced settings if it is known.")
    confirmed = sorted({cds["Translation_Table"] for cds in annotation["cds"]
                        if cds["Translation_Table_Source"].startswith(TABLE_CONFIRMED)})
    if confirmed:
        annotation["notes"].append("Genetic code confirmed by you for coding regions that state none: table "
                                   + ", ".join(str(value) for value in confirmed) + ".")
    map_peptides(annotation)
    if not annotation["cds"]:
        annotation["notes"].append("The annotation file contains no coding (CDS) features, so every SNV is "
                                   "reported as outside annotated coding regions or not annotated.")
    return annotation


def load_annotation(path, reference=None, genetic_code=None):
    """Load, match, and verify a GenBank or GFF3 annotation file (see finalize_annotation)."""

    kind = detect_annotation_format(path)
    if kind == GENBANK_FORMAT:
        annotation = load_genbank(path)
    elif kind == GFF3_FORMAT:
        annotation = load_gff3(path)
    else:
        raise ValueError("The annotation file must be a GenBank flat file (starting with LOCUS) or GFF3 "
                         "(starting with ##gff-version 3).")
    return finalize_annotation(annotation, reference, genetic_code)


def source_table(annotation):
    """One row per coding feature: what it is, where it is, and whether it is used."""

    rows = []
    for cds in annotation["cds"]:
        record = annotation["records"].get(cds["Reference_ID"], new_record(cds["Reference_ID"]))
        coding_length = sum(end - start + 1 for start, end in cds["Parts"])
        rows.append({
            "Reference_ID": cds["Reference_ID"], "Accession": record["Accession"], "Version": record["Version"],
            "Organism": record["Organism"], "Segment": record["Segment"], "Reference_Length": record["Length"],
            "Sequence_Source": record["Sequence_Source"] or "None available",
            "Sequence_Check": record["Sequence_Check"], "Gene": cds["Gene"], "Gene_ID": cds["Gene_ID"],
            "Locus_Tag": cds["Locus_Tag"], "CDS_ID": cds["CDS_ID"], "Product": cds["Product"],
            "Strand": cds["Strand"], "Parts_1_Based": ",".join(f"{start}..{end}" for start, end in cds["Parts"]),
            "Phases": ",".join(str(value) for value in cds["Phases"]) if cds["Phases"] else
            f"codon_start={cds['Codon_Start']}", "Coding_Length": coding_length,
            "Codons": coding_length // 3 if coding_length % 3 == 0 else None,
            "Translation_Table": cds["Translation_Table"],
            "Translation_Table_Source": cds["Translation_Table_Source"], "Structure": structure_label(cds),
            "Mature_Peptides": "; ".join(f"{item['Product']} ({item['First_Codon']}-{item['Last_Codon']})"
                                         for item in cds["Mature_Peptides"]),
            "Status": cds["Status"], "Status_Detail": cds["Status_Detail"], "Verification": cds["Verification"],
            "Annotation_Format": annotation["format"], "Annotation_File": annotation["source_file"]})
    return pd.DataFrame(rows, columns=SOURCE_COLUMNS)
