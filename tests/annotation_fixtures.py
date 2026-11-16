"""Hand-built SYNTHETIC annotation fixtures for the Phase 4 tests (not biological data).

Two invented segments, SEGA.1 (85 bases) and SEGB.1 (55 bases), are assembled
from pieces so every coordinate below can be checked by hand. Expected codons
and amino acids in the tests were derived by hand from these pieces, not by
the code under test.

SEGA.1 (1-based, inclusive):
     1-9   CCCCCCCCC                     5'UTR (noncoding context)
    10-27  ATG GAT TGG CCA TGC TAA       P1, + strand: M D W P C *
    21-32  ATG CTA AGC TGA               P3, + strand, overlaps P1 in another frame: M L S *
    28-33  GCTGAC                         spacer (28-32 are inside P3)
    34-51  reverse complement of          P2, - strand: complement(34..51)
           ATG TTT GCA TAC CGT TAG        M F A Y R *  (coding offset o <-> genome 51 - o)
    52-60  GGGNGGGGG                      stem_loop; N at 55
    61-69  reverse complement of GGCTTTTAG   P7 part 2 (- strand, spliced)
    70-75  CCCCCC                         P7 intron
    76-81  reverse complement of ATGCCA      P7 part 1: coding ATGCCA|GGCTTTTAG = M P G F *
    82-85  GGGG

SEGB.1:
     1-6   TTTTTT
     7-11  ATGGC                          P4 part 1 (+, spliced)
    12-17  GTAAGT                         P4 intron
    18-27  TNCAAATTAA                     P4 part 2 (GFF3 phase 1): coding ATG GCT NCA AAT TAA = M A X N *
    28-29  CC
    30-41  ATG TAA GCT TAA                P5: internal stop -> unsupported
    42-50  ATG AAA TAA                    P6: exception "RNA editing" -> unsupported
    51-55  CCCCC
"""

from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def rc(sequence):
    """Reverse complement written independently of the code under test."""

    return sequence[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


SEGA = ("C" * 9 + "ATGGATTGGCCATGCTAA" + "GCTGAC" + rc("ATGTTTGCATACCGTTAG") + "GGGNGGGGG"
        + rc("GGCTTTTAG") + "CCCCCC" + rc("ATGCCA") + "GGGG")
SEGB = "TTTTTT" + "ATGGC" + "GTAAGT" + "TNCAAATTAA" + "CC" + "ATGTAAGCTTAA" + "ATGAAATAA" + "CCCCC"
assert len(SEGA) == 85 and len(SEGB) == 55

GFF3_FEATURES = """##gff-version 3
##sequence-region SEGA.1 1 85
##sequence-region SEGB.1 1 55
SEGA.1\ttest\tregion\t1\t85\t.\t+\t.\tID=SEGA.1;organism=Synthetic test virus;segment=A
SEGA.1\ttest\tfive_prime_UTR\t1\t9\t.\t+\t.\tID=utr5
SEGA.1\ttest\tgene\t10\t27\t.\t+\t.\tID=gene-g1;Name=g1;gene=g1;Dbxref=GeneID:111
SEGA.1\ttest\tCDS\t10\t27\t.\t+\t0\tID=cds-P1;Parent=gene-g1;product=protein one;protein_id=P1
SEGA.1\ttest\tCDS\t21\t32\t.\t+\t0\tID=cds-P3;gene=g3;product=overlapping protein;protein_id=P3
SEGA.1\ttest\tCDS\t34\t51\t.\t-\t0\tID=cds-P2;gene=g2;product=minus-strand protein;protein_id=P2
SEGA.1\ttest\tstem_loop\t52\t60\t.\t+\t.\tID=sl1
SEGA.1\ttest\tCDS\t76\t81\t.\t-\t0\tID=cds-P7;gene=g7;product=spliced minus protein;protein_id=P7
SEGA.1\ttest\tCDS\t61\t69\t.\t-\t{p7_phase}\tID=cds-P7;gene=g7;product=spliced minus protein;protein_id=P7
SEGB.1\ttest\tregion\t1\t55\t.\t+\t.\tID=SEGB.1;organism=Synthetic test virus;segment=B
SEGB.1\ttest\tCDS\t7\t11\t.\t+\t0\tID=cds-P4;gene=g4;product=spliced protein;protein_id=P4
SEGB.1\ttest\tCDS\t18\t27\t.\t+\t{p4_phase}\tID=cds-P4;gene=g4;product=spliced protein;protein_id=P4
SEGB.1\ttest\tCDS\t30\t41\t.\t+\t0\tID=cds-P5;gene=g5;product=internal stop protein;protein_id=P5
SEGB.1\ttest\tCDS\t42\t50\t.\t+\t0\tID=cds-P6;gene=g6;product=edited protein;protein_id=P6;exception=RNA editing
"""


def write_gff3(folder, name="synthetic.gff3", with_fasta=True, p4_phase=1, p7_phase=0):
    text = GFF3_FEATURES.format(p4_phase=p4_phase, p7_phase=p7_phase)
    if with_fasta:
        text += f"##FASTA\n>SEGA.1\n{SEGA}\n>SEGB.1\n{SEGB}\n"
    path = Path(folder) / name
    path.write_text(text, encoding="utf-8")
    return path


def cds(location, gene, protein_id, translation=None, **extra):
    qualifiers = {"gene": [gene], "protein_id": [protein_id], "product": [f"{gene} protein"],
                  "codon_start": ["1"], "transl_table": ["1"]}
    if translation is not None:
        qualifiers["translation"] = [translation]
    qualifiers.update({key: [value] for key, value in extra.items()})
    return SeqFeature(location, type="CDS", qualifiers=qualifiers)


def write_genbank(folder, name="synthetic.gb", p1_translation="MDWPC", extra_records=()):
    """The same design as a GenBank file; complement(join()) is written by Biopython."""

    sega = SeqRecord(Seq(SEGA), id="SEGA.1", name="SEGA", description="Synthetic test segment A")
    sega.annotations.update(molecule_type="RNA", organism="Synthetic test virus")
    sega.features = [
        SeqFeature(FeatureLocation(0, 85, strand=1), type="source",
                   qualifiers={"organism": ["Synthetic test virus"], "segment": ["A"]}),
        SeqFeature(FeatureLocation(0, 9, strand=1), type="5'UTR"),
        cds(FeatureLocation(9, 27, strand=1), "g1", "P1", p1_translation),
        cds(FeatureLocation(20, 32, strand=1), "g3", "P3", "MLS"),
        cds(FeatureLocation(33, 51, strand=-1), "g2", "P2", "MFAYR"),
        cds(CompoundLocation([FeatureLocation(75, 81, strand=-1), FeatureLocation(60, 69, strand=-1)]),
            "g7", "P7", "MPGF"),
    ]
    segb = SeqRecord(Seq(SEGB), id="SEGB.1", name="SEGB", description="Synthetic test segment B")
    segb.annotations.update(molecule_type="RNA", organism="Synthetic test virus")
    segb.features = [
        SeqFeature(FeatureLocation(0, 55, strand=1), type="source",
                   qualifiers={"organism": ["Synthetic test virus"], "segment": ["B"]}),
        cds(CompoundLocation([FeatureLocation(6, 11, strand=1), FeatureLocation(17, 27, strand=1)]),
            "g4", "P4", "MAXN"),
        cds(FeatureLocation(29, 41, strand=1), "g5", "P5"),
        cds(FeatureLocation(41, 50, strand=1), "g6", "P6", "MK", exception="RNA editing"),
    ]
    path = Path(folder) / name
    SeqIO.write([sega, segb, *extra_records], path, "genbank")
    return path


def write_fasta(folder, name, record_id, sequence):
    path = Path(folder) / name
    path.write_text(f">{record_id}\n{sequence}\n", encoding="utf-8")
    return path
