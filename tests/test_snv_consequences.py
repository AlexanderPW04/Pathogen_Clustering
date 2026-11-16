"""Sprint 5 Phase 4: predicted SNV consequences, checked against hand-derived codons.

Every expected codon and amino acid below was worked out by hand: either from
the synthetic pieces in tests/annotation_fixtures.py, or for the real SARS-CoV-2
reference by slicing the FASTA string at coordinates computed with plain
arithmetic (CDS start + 3 x (codon - 1)), independently of the code under test.
Well-known SARS-CoV-2 examples (S D614G at A23403G, ORF1ab P4715L / nsp12 P323L
at C14408T, ORF8 Q27* at C27972T, N R203K / G204R at 28881-28883) are widely
reported names; the tests check that the arithmetic reproduces them.
"""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotation_fixtures import SEGA, write_genbank, write_gff3  # noqa: E402

from genomics.genome_annotation import load_annotation  # noqa: E402
from genomics.references import load_reference_fasta  # noqa: E402
from genomics.snv_consequences import (  # noqa: E402
    ANNOTATED,
    MISSENSE,
    NONCODING,
    NOT_ANNOTATED,
    OVERLAP_DIFFERS,
    PARTLY_ANNOTATED,
    START_LOST,
    STOP_GAINED,
    STOP_LOST,
    STOP_RETAINED,
    SYNONYMOUS,
    UNSUPPORTED_CONSEQUENCE,
    annotate_variants,
    classify_codon_change,
)

ANNOTATION = Path(__file__).resolve().parents[1] / "examples" / "annotation"


def variants(rows):
    frame = pd.DataFrame(rows, columns=["Reference_ID", "Position", "Ref", "Alt"])
    frame["Variant_Key"] = (frame["Reference_ID"] + ":" + frame["Position"].astype(str) + ":" + frame["Ref"]
                            + ">" + frame["Alt"])
    return frame


def one(table, key, cds=None):
    rows = table[table["Variant_Key"].eq(key)]
    if cds is not None:
        rows = rows[rows["CDS_ID"].eq(cds)]
    assert len(rows) == 1, (key, cds, len(rows))
    return rows.iloc[0]


class ClassificationRuleTest(unittest.TestCase):
    def test_rules(self):
        self.assertEqual(classify_codon_change(5, "D", "G", True), MISSENSE)
        self.assertEqual(classify_codon_change(5, "D", "D", True), SYNONYMOUS)
        self.assertEqual(classify_codon_change(5, "W", "*", True), STOP_GAINED)
        self.assertEqual(classify_codon_change(9, "*", "Q", True), STOP_LOST)
        self.assertEqual(classify_codon_change(9, "*", "*", True), STOP_RETAINED)
        self.assertEqual(classify_codon_change(1, "M", "V", True), START_LOST)
        self.assertEqual(classify_codon_change(1, "M", "M", True), START_LOST)  # any change to the start codon


class SyntheticConsequenceTest(unittest.TestCase):
    """Hand-derived expectations on the two synthetic segments."""

    CASES = [
        # (reference, position, ref, alt, CDS, codon number, ref codon, alt codon, change, consequence)
        ("SEGA.1", 14, "A", "G", "P1", 2, "GAT", "GGT", "D2G", MISSENSE),
        ("SEGA.1", 15, "T", "C", "P1", 2, "GAT", "GAC", "D2D", SYNONYMOUS),
        ("SEGA.1", 18, "G", "A", "P1", 3, "TGG", "TGA", "W3*", STOP_GAINED),
        ("SEGA.1", 26, "A", "C", "P1", 6, "TAA", "TCA", "*6S", STOP_LOST),
        ("SEGA.1", 27, "A", "G", "P1", 6, "TAA", "TAG", "*6*", STOP_RETAINED),
        ("SEGA.1", 10, "A", "G", "P1", 1, "ATG", "GTG", "M1V", START_LOST),
        # Minus strand: genome position = 51 - coding offset; Alt is complemented.
        ("SEGA.1", 47, "A", "G", "P2", 2, "TTT", "TCT", "F2S", MISSENSE),
        ("SEGA.1", 48, "A", "T", "P2", 2, "TTT", "ATT", "F2I", MISSENSE),
        ("SEGA.1", 46, "A", "G", "P2", 2, "TTT", "TTC", "F2F", SYNONYMOUS),
        ("SEGA.1", 40, "G", "T", "P2", 4, "TAC", "TAA", "Y4*", STOP_GAINED),
        ("SEGA.1", 35, "T", "G", "P2", 6, "TAG", "TCG", "*6S", STOP_LOST),
        ("SEGA.1", 34, "C", "T", "P2", 6, "TAG", "TAA", "*6*", STOP_RETAINED),
        ("SEGA.1", 51, "T", "C", "P2", 1, "ATG", "GTG", "M1V", START_LOST),
        # Minus-strand spliced CDS (complement(join(61..69,76..81))).
        ("SEGA.1", 76, "T", "C", "P7", 2, "CCA", "CCG", "P2P", SYNONYMOUS),
        ("SEGA.1", 69, "C", "A", "P7", 3, "GGC", "TGC", "G3C", MISSENSE),
        ("SEGA.1", 62, "T", "C", "P7", 5, "TAG", "TGG", "*5W", STOP_LOST),
        # Plus-strand spliced CDS; codon 2 spans the intron (genome 10, 11 | 18).
        ("SEGB.1", 18, "T", "A", "P4", 2, "GCT", "GCA", "A2A", SYNONYMOUS),
        ("SEGB.1", 11, "C", "A", "P4", 2, "GCT", "GAT", "A2D", MISSENSE),
        ("SEGB.1", 24, "T", "G", "P4", 4, "AAT", "AAG", "N4K", MISSENSE),
    ]

    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        path = Path(cls.folder.name)
        cls.annotations = {"GFF3": load_annotation(write_gff3(path)), "GenBank": load_annotation(write_genbank(path))}

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def annotate(self, rows, kind="GFF3", samples=None):
        return annotate_variants(variants(rows), self.annotations[kind], samples)

    def test_hand_derived_codons_both_formats(self):
        rows = [case[:4] for case in self.CASES]
        for kind in ["GFF3", "GenBank"]:
            table, _ = self.annotate(rows, kind)
            for reference, position, ref, alt, cds, codon, ref_codon, alt_codon, change, consequence in self.CASES:
                with self.subTest(kind=kind, position=position, cds=cds):
                    row = one(table, f"{reference}:{position}:{ref}>{alt}", cds)
                    self.assertEqual(row["Codon_Number"], codon)
                    self.assertEqual(row["Ref_Codon"], ref_codon)
                    self.assertEqual(row["Alt_Codon"], alt_codon)
                    self.assertEqual(row["Amino_Acid_Change"], change)
                    self.assertEqual(row["Consequence"], consequence)
                    self.assertEqual(row["Reference_Check"], "Reference base verified")

    def test_codon_position_and_coding_change(self):
        table, _ = self.annotate([("SEGA.1", 47, "A", "G"), ("SEGA.1", 46, "A", "G")])
        row = one(table, "SEGA.1:47:A>G")
        self.assertEqual((row["Coding_Position"], row["Codon_Position"], row["Strand"]), (5, 2, "-"))
        self.assertEqual(row["Coding_Change"], "c.5T>C")
        self.assertEqual(one(table, "SEGA.1:46:A>G")["Codon_Position"], 3)

    def test_noncoding_is_never_synonymous(self):
        table, summary = self.annotate([("SEGA.1", 9, "C", "T"), ("SEGA.1", 33, "C", "T"),
                                        ("SEGB.1", 14, "A", "G")])
        self.assertTrue(table["Consequence"].eq(NONCODING).all())
        self.assertTrue(summary["Consequence_Category"].eq(NONCODING).all())
        self.assertIn("five_prime_UTR", one(table, "SEGA.1:9:C>T")["Consequence_Detail"])
        self.assertEqual(one(table, "SEGB.1:14:A>G")["Amino_Acid_Change"], "")  # intron: no amino acid

    def test_overlapping_frames_keep_every_consequence(self):
        table, summary = self.annotate([("SEGA.1", 22, "T", "C"), ("SEGA.1", 25, "T", "C"),
                                        ("SEGA.1", 24, "C", "T")])
        self.assertEqual(one(table, "SEGA.1:22:T>C", "P1")["Amino_Acid_Change"], "C5R")
        self.assertEqual(one(table, "SEGA.1:22:T>C", "P3")["Consequence"], START_LOST)
        self.assertEqual(one(table, "SEGA.1:25:T>C", "P1")["Consequence"], STOP_LOST)
        self.assertEqual(one(table, "SEGA.1:25:T>C", "P3")["Amino_Acid_Change"], "L2P")
        categories = summary.set_index("Variant_Key")["Consequence_Category"]
        self.assertEqual(categories["SEGA.1:22:T>C"], OVERLAP_DIFFERS)
        self.assertEqual(categories["SEGA.1:25:T>C"], OVERLAP_DIFFERS)
        self.assertEqual(categories["SEGA.1:24:C>T"], SYNONYMOUS)  # synonymous in both frames
        self.assertEqual(summary.set_index("Variant_Key").loc["SEGA.1:24:C>T", "Coding_Regions"], 2)

    def test_same_position_on_different_segments_is_annotated_separately(self):
        table, summary = self.annotate([("SEGA.1", 24, "C", "T"), ("SEGB.1", 24, "T", "G")])
        self.assertEqual(set(table[table["Position"].eq(24)]["Reference_ID"]), {"SEGA.1", "SEGB.1"})
        self.assertEqual(summary.set_index("Variant_Key").loc["SEGB.1:24:T>G", "Segment"], "B")

    def test_reference_mismatch_and_wrong_accession_are_not_annotated(self):
        table, summary = self.annotate([("SEGA.1", 14, "G", "A"), ("SEGA.2", 14, "A", "G"),
                                        ("NC_045512.2", 14, "A", "G"), ("SEGA.1", 999, "A", "G")])
        mismatch = one(table, "SEGA.1:14:G>A")
        self.assertEqual(mismatch["Consequence"], NOT_ANNOTATED)
        self.assertEqual(mismatch["Reference_Check"], "Mismatch: reference has A")
        version = one(table, "SEGA.2:14:A>G")
        self.assertEqual(version["Consequence"], NOT_ANNOTATED)
        self.assertIn("never transferred", version["Consequence_Detail"])
        self.assertEqual(one(table, "NC_045512.2:14:A>G")["Consequence"], NOT_ANNOTATED)
        self.assertEqual(one(table, "SEGA.1:999:A>G")["Reference_Check"], "Position outside the reference")
        self.assertTrue(summary["Annotation_Status"].eq("Not annotated").all())
        self.assertTrue(table["Amino_Acid_Change"].eq("").all())

    def test_ambiguous_bases_indels_and_unsupported_structures(self):
        table, _ = self.annotate([("SEGB.1", 20, "C", "T"), ("SEGB.1", 19, "A", "G"), ("SEGA.1", 14, "A", "AT"),
                                  ("SEGA.1", 14, "A", "R"), ("SEGB.1", 35, "A", "G"), ("SEGB.1", 46, "A", "G"),
                                  ("SEGA.1", 14, "A", "A")])
        self.assertIn("ambiguous base", one(table, "SEGB.1:20:C>T")["Consequence_Detail"])
        self.assertEqual(one(table, "SEGB.1:19:A>G")["Reference_Check"], "Ambiguous reference base (N)")
        self.assertIn("indel", one(table, "SEGA.1:14:A>AT")["Consequence_Detail"])
        self.assertIn("ambiguous bases are not interpreted", one(table, "SEGA.1:14:A>R")["Consequence_Detail"])
        self.assertIn("Internal stop", one(table, "SEGB.1:35:A>G")["Consequence_Detail"])
        self.assertIn("RNA editing", one(table, "SEGB.1:46:A>G")["Consequence_Detail"])
        self.assertTrue(table["Consequence"].eq(UNSUPPORTED_CONSEQUENCE).all())
        self.assertTrue(table["Amino_Acid_Change"].eq("").all())

    def test_adjacent_snvs_in_one_codon_are_flagged_not_combined(self):
        rows = [("SEGA.1", 14, "A", "G"), ("SEGA.1", 15, "T", "C"), ("SEGA.1", 18, "G", "A")]
        together = {"SEGA.1:14:A>G": {"s1", "s2"}, "SEGA.1:15:T>C": {"s2"}, "SEGA.1:18:G>A": {"s1"}}
        table, summary = self.annotate(rows, samples=together)
        first = one(table, "SEGA.1:14:A>G")
        self.assertEqual(first["Shared_Codon_With"], "SEGA.1:15:T>C")
        self.assertIn("reported together in up to 1 sample(s)", first["Linkage_Note"])
        self.assertIn("combined amino-acid change is not predicted", first["Linkage_Note"])
        self.assertEqual(first["Amino_Acid_Change"], "D2G")  # each SNV alone: GGT, never GGC
        self.assertEqual(one(table, "SEGA.1:15:T>C")["Amino_Acid_Change"], "D2D")
        self.assertEqual(one(table, "SEGA.1:18:G>A")["Shared_Codon_With"], "")
        apart = {"SEGA.1:14:A>G": {"s1"}, "SEGA.1:15:T>C": {"s2"}, "SEGA.1:18:G>A": {"s1"}}
        table, _ = self.annotate(rows, samples=apart)
        self.assertTrue(table["Shared_Codon_With"].eq("").all())  # never reported together
        table, _ = self.annotate(rows)
        self.assertEqual(one(table, "SEGA.1:15:T>C")["Shared_Codon_With"], "SEGA.1:14:A>G")

    def test_each_distinct_snv_is_annotated_once(self):
        rows = [("SEGA.1", 14, "A", "G")] * 500 + [("SEGA.1", 15, "T", "C")] * 300
        table, summary = self.annotate(rows)
        self.assertEqual(len(table), 2)
        self.assertEqual(len(summary), 2)


class RealSarsCov2ConsequenceTest(unittest.TestCase):
    """NC_045512.2 (real RefSeq) with hand arithmetic on the FASTA string."""

    @classmethod
    def setUpClass(cls):
        reference = load_reference_fasta(ANNOTATION / "NC_045512.2.fasta")
        cls.sequence = reference["Sequence"]
        cls.annotation = load_annotation(ANNOTATION / "NC_045512.2.gb", reference)
        cls.gff = load_annotation(ANNOTATION / "NC_045512.2.gff3", reference)

    def codon(self, cds_start, number):
        start = cds_start + 3 * (number - 1)  # 1-based first base of the codon
        return self.sequence[start - 1:start + 2]

    def annotate(self, rows, annotation=None, samples=None):
        return annotate_variants(variants(rows), annotation or self.annotation, samples)

    def test_spike_d614g(self):
        self.assertEqual(self.codon(21563, 614), "GAT")  # S CDS 21563..25384; 23403 is codon position 2
        for annotation in [self.annotation, self.gff]:
            table, _ = self.annotate([("NC_045512.2", 23403, "A", "G")], annotation)
            row = one(table, "NC_045512.2:23403:A>G")
            self.assertEqual((row["Gene"], row["Codon_Number"], row["Codon_Position"]), ("S", 614, 2))
            self.assertEqual((row["Ref_Codon"], row["Alt_Codon"], row["Amino_Acid_Change"]), ("GAT", "GGT", "D614G"))
            self.assertEqual(row["Consequence"], MISSENSE)

    def test_orf1ab_p4715l_after_the_frameshift_is_nsp12_p323l(self):
        # In pp1ab, 14408 lies in part 2: coding offset 13203 + (14408 - 13468) = 14143 -> codon 4715, position 2.
        offset = 13203 + (14408 - 13468)
        self.assertEqual((offset // 3 + 1, offset % 3 + 1), (4715, 2))
        table, _ = self.annotate([("NC_045512.2", 14408, "C", "T")])
        self.assertEqual(len(table), 1)  # 14408 is beyond ORF1a (266..13483)
        row = table.iloc[0]
        self.assertEqual((row["CDS_ID"], row["Amino_Acid_Change"]), ("YP_009724389.1", "P4715L"))
        self.assertEqual(row["Ref_Codon"], self.sequence[14407 - 1:14407 + 2])  # 14407..14409 in the shifted frame
        self.assertEqual((row["Mature_Peptide"], row["Mature_Peptide_Position"]),
                         ("RNA-dependent RNA polymerase", 323))

    def test_c3037t_is_synonymous_in_orf1a_and_orf1ab(self):
        offset = 3037 - 266
        self.assertEqual((offset // 3 + 1, offset % 3 + 1), (924, 3))
        self.assertEqual(self.codon(266, 924), "TTC")
        table, summary = self.annotate([("NC_045512.2", 3037, "C", "T")])
        self.assertEqual(set(table["CDS_ID"]), {"YP_009724389.1", "YP_009725295.1"})
        self.assertTrue(table["Amino_Acid_Change"].eq("F924F").all())
        self.assertTrue(table["Mature_Peptide"].eq("nsp3").all())
        self.assertTrue(table["Mature_Peptide_Position"].eq(106).all())  # 924 - 819 + 1
        self.assertEqual(summary.iloc[0]["Consequence_Category"], SYNONYMOUS)
        self.assertEqual(summary.iloc[0]["Amino_Acid_Changes"], "ORF1ab:F924F")

    def test_5utr_c241t_is_noncoding(self):
        table, summary = self.annotate([("NC_045512.2", 241, "C", "T")])
        self.assertEqual(summary.iloc[0]["Consequence_Category"], NONCODING)
        self.assertIn("5'UTR", table.iloc[0]["Consequence_Detail"])

    def test_nucleocapsid_r203k_g204r_and_shared_codon(self):
        self.assertEqual(self.codon(28274, 203), "AGG")  # N CDS 28274..29533
        self.assertEqual(self.codon(28274, 204), "GGA")
        rows = [("NC_045512.2", 28881, "G", "A"), ("NC_045512.2", 28882, "G", "A"), ("NC_045512.2", 28883, "G", "C")]
        samples = {"NC_045512.2:28881:G>A": {"ww1"}, "NC_045512.2:28882:G>A": {"ww1"},
                   "NC_045512.2:28883:G>C": {"ww1"}}
        table, _ = self.annotate(rows, samples=samples)
        self.assertEqual(one(table, "NC_045512.2:28881:G>A")["Amino_Acid_Change"], "R203K")
        self.assertEqual(one(table, "NC_045512.2:28882:G>A")["Amino_Acid_Change"], "R203R")
        self.assertEqual(one(table, "NC_045512.2:28883:G>C")["Amino_Acid_Change"], "G204R")
        self.assertEqual(one(table, "NC_045512.2:28881:G>A")["Shared_Codon_With"], "NC_045512.2:28882:G>A")
        self.assertEqual(one(table, "NC_045512.2:28883:G>C")["Shared_Codon_With"], "")  # codon 204 alone

    def test_stop_gain_stop_loss_and_start_loss(self):
        self.assertEqual(self.codon(27894, 27), "CAA")  # ORF8 27894..28259
        self.assertEqual(self.sequence[29531 - 1:29533], "TAA")  # N stop codon
        self.assertEqual(self.sequence[21563 - 1:21565], "ATG")  # S start codon
        table, _ = self.annotate([("NC_045512.2", 27972, "C", "T"), ("NC_045512.2", 29531, "T", "C"),
                                  ("NC_045512.2", 21563, "A", "G")])
        self.assertEqual(one(table, "NC_045512.2:27972:C>T")["Amino_Acid_Change"], "Q27*")
        self.assertEqual(one(table, "NC_045512.2:27972:C>T")["Consequence"], STOP_GAINED)
        self.assertEqual(one(table, "NC_045512.2:29531:T>C")["Amino_Acid_Change"], "*420Q")
        self.assertEqual(one(table, "NC_045512.2:29531:T>C")["Consequence"], STOP_LOST)
        self.assertEqual(one(table, "NC_045512.2:21563:A>G")["Consequence"], START_LOST)

    def test_orf7a_orf7b_overlap_keeps_both(self):
        # ORF7b starts at 27756 (ATG); the same base is codon 121, position 3, of ORF7a (27394..27759).
        self.assertEqual(self.sequence[27756 - 1:27758], "ATG")
        offset = 27756 - 27394
        self.assertEqual((offset // 3 + 1, offset % 3 + 1), (121, 3))
        ref_codon = self.codon(27394, 121)
        table, summary = self.annotate([("NC_045512.2", 27756, "A", "G")])
        orf7b = table[table["Gene"].eq("ORF7b")].iloc[0]
        orf7a = table[table["Gene"].eq("ORF7a")].iloc[0]
        self.assertEqual(orf7b["Consequence"], START_LOST)
        self.assertEqual(orf7a["Ref_Codon"], ref_codon)
        self.assertEqual(orf7a["Alt_Codon"], ref_codon[:2] + "G")
        self.assertEqual(summary.iloc[0]["Consequence_Category"], OVERLAP_DIFFERS)

    def test_frameshift_base_gets_no_single_pp1ab_consequence(self):
        base = self.sequence[13468 - 1]
        alt = "A" if base != "A" else "G"
        table, summary = self.annotate([("NC_045512.2", 13468, base, alt)])
        pp1ab = one(table, f"NC_045512.2:13468:{base}>{alt}", "YP_009724389.1")
        self.assertEqual(pp1ab["Consequence"], UNSUPPORTED_CONSEQUENCE)
        self.assertIn("read twice", pp1ab["Consequence_Detail"])
        pp1a = one(table, f"NC_045512.2:13468:{base}>{alt}", "YP_009725295.1")
        self.assertNotEqual(pp1a["Consequence"], UNSUPPORTED_CONSEQUENCE)
        self.assertEqual(summary.iloc[0]["Annotation_Status"], PARTLY_ANNOTATED)

    def test_sars_cov_2_annotation_is_not_applied_to_influenza(self):
        influenza = load_annotation(ANNOTATION / "influenza_A_PR8_segments.gb")
        table, summary = self.annotate([("NC_045512.2", 23403, "A", "G")], influenza)
        self.assertEqual(summary.iloc[0]["Annotation_Status"], "Not annotated")
        self.assertEqual(table.iloc[0]["Gene"], "")


class RealInfluenzaConsequenceTest(unittest.TestCase):
    """Segmented reference (real RefSeq records): spliced M2 overlaps M1 in another frame."""

    @classmethod
    def setUpClass(cls):
        cls.annotation = load_annotation(ANNOTATION / "influenza_A_PR8_segments.gb")
        cls.segment7 = cls.annotation["records"]["NC_002016.1"]["Sequence"]

    def test_segment_7_m1_and_spliced_m2(self):
        # M1 26..784: 760 -> offset 734 -> codon 245, position 3.
        # M2 join(26..51,740..1007): 760 -> offset 26 + (760 - 740) = 46 -> codon 16, position 2.
        base = self.segment7[760 - 1]
        alt = "C" if base != "C" else "T"
        m1_codon = self.segment7[26 + 3 * 244 - 1:26 + 3 * 244 + 2]
        m2_coding = self.segment7[26 - 1:51] + self.segment7[740 - 1:1007]
        m2_codon = m2_coding[45:48]
        table, summary = annotate_variants(variants([("NC_002016.1", 760, base, alt)]), self.annotation)
        m1 = table[table["Gene"].eq("M1")].iloc[0]
        m2 = table[table["Gene"].eq("M2")].iloc[0]
        self.assertEqual((m1["Codon_Number"], m1["Codon_Position"], m1["Ref_Codon"]), (245, 3, m1_codon))
        self.assertEqual((m2["Codon_Number"], m2["Codon_Position"], m2["Ref_Codon"]), (16, 2, m2_codon))
        self.assertEqual(m2["Alt_Codon"], m2_codon[0] + alt + m2_codon[2])
        self.assertEqual(summary.iloc[0]["Segment"], "7")
        self.assertEqual(summary.iloc[0]["Annotation_Status"], ANNOTATED)

    def test_two_segments_in_one_table(self):
        segment4 = self.annotation["records"]["NC_002017.1"]["Sequence"]
        rows = [("NC_002016.1", 100, self.segment7[99], "A" if self.segment7[99] != "A" else "G"),
                ("NC_002017.1", 100, segment4[99], "A" if segment4[99] != "A" else "G")]
        table, summary = annotate_variants(variants(rows), self.annotation)
        genes = summary.set_index("Reference_ID")["Genes"]
        self.assertEqual(genes["NC_002016.1"], "M1")
        self.assertEqual(genes["NC_002017.1"], "HA")
        self.assertEqual(set(summary["Segment"]), {"7", "4"})


class SequenceSanityTest(unittest.TestCase):
    def test_fixture_pieces(self):
        self.assertEqual(SEGA[9:27], "ATGGATTGGCCATGCTAA")
        self.assertEqual(SEGA[20:32], "ATGCTAAGCTGA")  # P3 in the +1 frame of P1


if __name__ == "__main__":
    unittest.main()
