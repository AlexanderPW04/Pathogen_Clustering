"""Sprint 5, Phase 4: predicted coding consequences of single-nucleotide variants.

Each distinct SNV is annotated once, against one verified reference:

1. its reference must be described by the annotation (the same
   accession.version; an uploaded FASTA must be identical to the annotation's
   sequence, see `genome_annotation.attach_reference`);
2. Ref and Alt must be single, different A/C/G/T bases;
3. the reference base at the 1-based position must equal Ref;
4. every supported CDS containing the position gives one consequence, so a
   base in overlapping reading frames keeps all of them;
5. the codon is read on the coding strand (Alt is complemented for "-" genes),
   translated with the CDS's translation table, and classified.

Categories follow the Sequence Ontology definitions:

- start lost (SO:0002012): any change to the verified start codon;
- stop gained (SO:0001587): a sense codon becomes a stop codon;
- stop lost (SO:0001578): the terminal stop codon becomes a sense codon;
- stop retained (SO:0001567): the stop codon changes but remains a stop;
- synonymous (SO:0001819): the encoded amino acid is unchanged;
- missense (SO:0001583): a different amino acid, same protein length;
- noncoding: outside every annotated CDS. It is never called synonymous.

Everything else is reported as unsupported with its reason, for example an
ambiguous base, an indel or multi-base change, a base read twice because of
programmed ribosomal frameshifting, or a CDS whose structure was not verified.

Each consequence describes one SNV alone against the reference. Pooled reads
do not show which SNVs share a viral genome, so when other observed SNVs fall
in the same codon their combined effect is not predicted (`linkage_flags`).
"""

import pandas as pd

from genomics.genome_annotation import (
    COMPLEMENT,
    NUCLEOTIDES,
    SUPPORTED,
    coding_offsets,
    translate_codon,
)


SYNONYMOUS = "Synonymous"
MISSENSE = "Missense"
STOP_GAINED = "Stop gained"
STOP_LOST = "Stop lost"
STOP_RETAINED = "Stop retained"
START_LOST = "Start lost"
NONCODING = "Noncoding"
UNSUPPORTED_CONSEQUENCE = "Unsupported or ambiguous"
NOT_ANNOTATED = "Not annotated"
OVERLAP_DIFFERS = "Differs between overlapping coding regions"
CODING_CONSEQUENCES = [SYNONYMOUS, MISSENSE, STOP_GAINED, STOP_LOST, STOP_RETAINED, START_LOST]
OTHER_CODING = [STOP_GAINED, STOP_LOST, STOP_RETAINED, START_LOST]
CATEGORY_ORDER = CODING_CONSEQUENCES + [OVERLAP_DIFFERS, NONCODING, UNSUPPORTED_CONSEQUENCE, NOT_ANNOTATED]

ANNOTATED = "Annotated"
PARTLY_ANNOTATED = "Partly annotated"
REF_VERIFIED = "Reference base verified"

CONSEQUENCE_COLUMNS = [
    "Variant_Key", "Reference_ID", "Segment", "Position", "Ref", "Alt", "Nucleotide_Change", "Reference_Check",
    "Gene", "Gene_ID", "Locus_Tag", "CDS_ID", "Product", "Strand", "Coding_Position", "Coding_Change",
    "Codon_Number", "Codon_Position", "Ref_Codon", "Alt_Codon", "Ref_Amino_Acid", "Alt_Amino_Acid",
    "Amino_Acid_Change", "Mature_Peptide", "Mature_Peptide_Position", "Consequence", "Consequence_Detail",
    "Shared_Codon_With", "Linkage_Note",
]
THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His",
    "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Stop",
}


def classify_codon_change(codon_number, ref_amino, alt_amino, start_verified):
    """The consequence category of one codon change (codon_number is 1-based)."""

    if codon_number == 1 and start_verified:
        return START_LOST
    if ref_amino == "*":
        return STOP_RETAINED if alt_amino == "*" else STOP_LOST
    if alt_amino == "*":
        return STOP_GAINED
    return SYNONYMOUS if ref_amino == alt_amino else MISSENSE


def blank_row(variant, record):
    row = dict.fromkeys(CONSEQUENCE_COLUMNS, "")
    row.update({"Variant_Key": variant["Variant_Key"], "Reference_ID": variant["Reference_ID"],
                "Segment": record["Segment"] if record else "", "Position": variant["Position"],
                "Ref": variant["Ref"], "Alt": variant["Alt"],
                "Nucleotide_Change": f"{variant['Position']}{variant['Ref']}>{variant['Alt']}"})
    for name in ["Coding_Position", "Codon_Number", "Codon_Position", "Mature_Peptide_Position"]:
        row[name] = None
    return row


def cds_row(base, cds, offset, alt):
    """One consequence row for an SNV inside one supported CDS."""

    row = dict(base)
    row.update({"Gene": cds["Gene"], "Gene_ID": cds["Gene_ID"], "Locus_Tag": cds["Locus_Tag"],
                "CDS_ID": cds["CDS_ID"], "Product": cds["Product"], "Strand": cds["Strand"]})
    coding = cds["Coding_Sequence"]
    codon_index, codon_position = divmod(offset, 3)
    ref_codon = coding[3 * codon_index:3 * codon_index + 3]
    coding_alt = alt if cds["Strand"] == "+" else alt.translate(COMPLEMENT)
    alt_codon = ref_codon[:codon_position] + coding_alt + ref_codon[codon_position + 1:]
    row.update({"Coding_Position": offset + 1, "Codon_Number": codon_index + 1, "Codon_Position": codon_position + 1,
                "Coding_Change": f"c.{offset + 1}{coding[offset]}>{coding_alt}",
                "Ref_Codon": ref_codon, "Alt_Codon": alt_codon})
    ref_amino = translate_codon(ref_codon, cds["Translation_Table"])
    alt_amino = translate_codon(alt_codon, cds["Translation_Table"])
    if ref_amino is None or alt_amino is None:
        row["Consequence"] = UNSUPPORTED_CONSEQUENCE
        row["Consequence_Detail"] = "The reference codon contains an ambiguous base, so it cannot be translated."
        return row
    if codon_index == 0 and not cds["Start_Codon_Verified"]:
        row["Consequence"] = UNSUPPORTED_CONSEQUENCE
        row["Consequence_Detail"] = ("The first codon is not a recognized start codon, so a change to it is not "
                                     "interpreted.")
        return row
    consequence = classify_codon_change(codon_index + 1, ref_amino, alt_amino, cds["Start_Codon_Verified"])
    row.update({"Ref_Amino_Acid": THREE_LETTER[ref_amino], "Alt_Amino_Acid": THREE_LETTER[alt_amino],
                "Amino_Acid_Change": f"{ref_amino}{codon_index + 1}{alt_amino}", "Consequence": consequence})
    if consequence == START_LOST:
        row["Consequence_Detail"] = ("The start codon changes; whether translation can begin elsewhere is not "
                                     "predicted.")
    for peptide in cds["Mature_Peptides"]:
        if peptide["First_Codon"] <= codon_index + 1 <= peptide["Last_Codon"]:
            row["Mature_Peptide"] = peptide["Product"]
            row["Mature_Peptide_Position"] = codon_index + 1 - peptide["First_Codon"] + 1
            break
    return row


def context_names(annotation, reference_id, position):
    names = [item["Name"] if item["Type"] == "gene" else (
        item["Type"] if item["Name"] == item["Type"] else f"{item['Type']} ({item['Name']})")
             for item in annotation["context"]
             if item["Reference_ID"] == reference_id and item["Start"] <= position <= item["End"]]
    return "; ".join(dict.fromkeys(names))


def annotate_snv(variant, annotation, by_reference):
    """All consequence rows for one distinct SNV (a dict with Variant_Key, Reference_ID, Position, Ref, Alt)."""

    record = annotation["records"].get(variant["Reference_ID"])
    base = blank_row(variant, record)
    ref, alt = str(variant["Ref"]).upper(), str(variant["Alt"]).upper()

    def single(consequence, detail, check=""):
        base.update({"Consequence": consequence, "Consequence_Detail": detail, "Reference_Check": check})
        return [base]

    if len(ref) != 1 or len(alt) != 1:
        return single(UNSUPPORTED_CONSEQUENCE, "Not a single-nucleotide substitution (an indel or multi-base "
                                               "change); not annotated in this phase.")
    if ref not in NUCLEOTIDES or alt not in NUCLEOTIDES or ref == alt:
        return single(UNSUPPORTED_CONSEQUENCE, "Ref and Alt must be different A, C, G, or T bases; ambiguous "
                                               "bases are not interpreted.")
    if record is None:
        described = ", ".join(sorted(annotation["records"])) or "no sequences"
        return single(NOT_ANNOTATED, f"The annotation describes {described}, not {variant['Reference_ID']}; "
                                     "annotations are never transferred between different references.",
                      "Reference not described by the annotation")
    if not record["Compatible"]:
        return single(NOT_ANNOTATED, record["Sequence_Check"], "Annotation and reference are incompatible")
    if not record["Sequence"]:
        return single(NOT_ANNOTATED, "Neither the annotation nor an uploaded reference FASTA supplies this "
                                     "reference's sequence, so the reference base and codons cannot be checked.",
                      "Not checked (no reference sequence)")
    position = int(variant["Position"])
    if not 1 <= position <= len(record["Sequence"]):
        return single(NOT_ANNOTATED, f"Position {position} is outside {variant['Reference_ID']} "
                                     f"(1 to {len(record['Sequence'])}).", "Position outside the reference")
    observed = record["Sequence"][position - 1]
    if observed not in NUCLEOTIDES:
        return single(UNSUPPORTED_CONSEQUENCE, f"The reference base at {position} is ambiguous ({observed}).",
                      f"Ambiguous reference base ({observed})")
    if observed != ref:
        return single(NOT_ANNOTATED, f"The reference has {observed} at {position}, not {ref}; the record is not "
                                     "annotated.", f"Mismatch: reference has {observed}")
    base["Reference_Check"] = REF_VERIFIED
    rows = []
    for cds in by_reference.get(variant["Reference_ID"], []):
        offsets = coding_offsets(cds, position)
        if not offsets:
            continue
        if cds["Status"] != SUPPORTED:
            row = dict(base)
            row.update({"Gene": cds["Gene"], "Gene_ID": cds["Gene_ID"], "Locus_Tag": cds["Locus_Tag"],
                        "CDS_ID": cds["CDS_ID"], "Product": cds["Product"], "Strand": cds["Strand"],
                        "Consequence": UNSUPPORTED_CONSEQUENCE,
                        "Consequence_Detail": "Coding structure not supported: " + cds["Status_Detail"]})
            rows.append(row)
        elif len(offsets) > 1:
            row = dict(base)
            row.update({"Gene": cds["Gene"], "Gene_ID": cds["Gene_ID"], "Locus_Tag": cds["Locus_Tag"],
                        "CDS_ID": cds["CDS_ID"], "Product": cds["Product"], "Strand": cds["Strand"],
                        "Consequence": UNSUPPORTED_CONSEQUENCE,
                        "Consequence_Detail": "This base is read twice in this protein because of programmed "
                                              "ribosomal frameshifting (coding positions "
                                              + ", ".join(str(value + 1) for value in offsets)
                                              + "); a single consequence is not assigned."})
            rows.append(row)
        else:
            rows.append(cds_row(base, cds, offsets[0], alt))
    if not rows:
        context = context_names(annotation, variant["Reference_ID"], position)
        base["Consequence"] = NONCODING
        base["Consequence_Detail"] = ("Outside every annotated coding region"
                                      + (f" (within: {context})" if context else "") + ".")
        base["Gene"] = ""
        rows.append(base)
    return rows


def linkage_flags(rows, samples_by_variant=None):
    """Mark SNVs whose codon also contains another observed SNV at a different position.

    `rows` are consequence dictionaries (changed in place). `samples_by_variant`
    optionally maps each Variant_Key to the set of samples that reported it;
    then only SNVs reported together in at least one sample are flagged.
    """

    codons = {}
    for row in rows:
        if row["CDS_ID"] and row["Codon_Number"] is not None:
            codons.setdefault((row["Reference_ID"], row["CDS_ID"], row["Codon_Number"]), []).append(row)
    for group in codons.values():
        if len({row["Position"] for row in group}) < 2:
            continue
        for row in group:
            others = list(dict.fromkeys(other["Variant_Key"] for other in group
                                        if other["Position"] != row["Position"]))
            where = ""
            if samples_by_variant is not None:
                mine = samples_by_variant.get(row["Variant_Key"], set())
                together = {other: len(mine & samples_by_variant.get(other, set())) for other in others}
                others = [other for other in others if together[other]]
                if not others:
                    continue
                where = f" (reported together in up to {max(together[other] for other in others)} sample(s))"
            row["Shared_Codon_With"] = "; ".join(others)
            row["Linkage_Note"] = (
                "Another observed SNV changes this codon" + where + ". Each change is predicted alone; pooled "
                "reads do not show whether they occur on the same genome, so a combined amino-acid change is "
                "not predicted.")
    return rows


def combined_consequence(consequences):
    """One category per SNV that keeps disagreement between overlapping coding regions visible."""

    found = list(dict.fromkeys(consequences))
    if len(found) == 1:
        return found[0]
    return OVERLAP_DIFFERS


def annotate_variants(variants, annotation, samples_by_variant=None):
    """Annotate each distinct (Reference_ID, Position, Ref, Alt) once.

    Returns (consequences: one row per SNV per overlapping CDS, summary: one
    row per SNV). `variants` needs Variant_Key, Reference_ID, Position, Ref,
    and Alt; duplicates are dropped before annotation so repeated sample
    records are never translated again.
    """

    unique = variants[["Variant_Key", "Reference_ID", "Position", "Ref", "Alt"]].drop_duplicates("Variant_Key")
    by_reference = {}
    for cds in annotation["cds"]:
        by_reference.setdefault(cds["Reference_ID"], []).append(cds)
    rows = []
    for variant in unique.to_dict("records"):
        rows.extend(annotate_snv(variant, annotation, by_reference))
    linkage_flags(rows, samples_by_variant)
    summary = summarize_variants(rows)
    table = pd.DataFrame(rows, columns=CONSEQUENCE_COLUMNS)
    for name in ["Position", "Coding_Position", "Codon_Number", "Codon_Position", "Mature_Peptide_Position"]:
        table[name] = pd.to_numeric(table[name], errors="coerce").astype("Int64")
    return table, summary


SUMMARY_COLUMNS = ["Variant_Key", "Reference_ID", "Segment", "Position", "Ref", "Alt", "Nucleotide_Change",
                   "Reference_Check", "Annotation_Status", "Genes", "Amino_Acid_Changes", "Mature_Peptides",
                   "Consequence_Category", "Consequences_By_CDS", "Coding_Regions", "Shared_Codon_With",
                   "Detail"]


def summarize_variants(rows):
    """One row per SNV (from consequence dictionaries): a combined category plus every per-CDS consequence."""

    groups = {}
    for row in rows:
        groups.setdefault(row["Variant_Key"], []).append(row)
    summary = []
    for key, group in groups.items():
        first = group[0]
        consequences = [row["Consequence"] for row in group]
        if NOT_ANNOTATED in consequences:
            status = "Not annotated"
        elif UNSUPPORTED_CONSEQUENCE in consequences:
            status = PARTLY_ANNOTATED if len(set(consequences)) > 1 else "Not annotated"
        else:
            status = ANNOTATED
        coding = [row for row in group if row["CDS_ID"]]
        summary.append({
            "Variant_Key": key, "Reference_ID": first["Reference_ID"], "Segment": first["Segment"],
            "Position": first["Position"], "Ref": first["Ref"], "Alt": first["Alt"],
            "Nucleotide_Change": first["Nucleotide_Change"], "Reference_Check": first["Reference_Check"],
            "Annotation_Status": status,
            "Genes": "; ".join(dict.fromkeys(row["Gene"] for row in coding if row["Gene"])),
            "Amino_Acid_Changes": "; ".join(dict.fromkeys(
                f"{row['Gene'] or row['CDS_ID']}:{row['Amino_Acid_Change']}" for row in coding
                if row["Amino_Acid_Change"])),
            "Mature_Peptides": "; ".join(dict.fromkeys(
                f"{row['Mature_Peptide']} {row['Mature_Peptide_Position']}" for row in coding
                if row["Mature_Peptide"])),
            "Consequence_Category": combined_consequence(consequences),
            "Consequences_By_CDS": "; ".join(f"{row['CDS_ID']}: {row['Consequence']}" for row in coding)
            if coding else consequences[0],
            "Coding_Regions": len(coding),
            "Shared_Codon_With": "; ".join(dict.fromkeys(
                item for row in group for item in row["Shared_Codon_With"].split("; ") if item)),
            "Detail": " ".join(dict.fromkeys(row["Consequence_Detail"] for row in group if row["Consequence_Detail"])),
        })
    table = pd.DataFrame(summary, columns=SUMMARY_COLUMNS)
    table["Position"] = pd.to_numeric(table["Position"], errors="coerce").astype("Int64")
    return table
