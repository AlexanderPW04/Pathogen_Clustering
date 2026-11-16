"""Sprint 5 Phase 4: loading and verifying genome annotations (GenBank and GFF3).

Real data: the official NCBI RefSeq files in examples/annotation/ (SARS-CoV-2
NC_045512.2 as GenBank and GFF3, and the eight influenza A/Puerto Rico/8/1934
segments as GenBank). Synthetic data: tests/annotation_fixtures.py.
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotation_fixtures import SEGA, rc, write_fasta, write_genbank, write_gff3  # noqa: E402

from genomics.genome_annotation import (  # noqa: E402
    GENBANK_FORMAT,
    GFF3_FORMAT,
    SUPPORTED,
    UNSUPPORTED,
    coding_offsets,
    detect_annotation_format,
    load_annotation,
    source_table,
)
from genomics.references import load_reference_fasta  # noqa: E402

ANNOTATION = Path(__file__).resolve().parents[1] / "examples" / "annotation"
SARS2_GB = ANNOTATION / "NC_045512.2.gb"
SARS2_GFF = ANNOTATION / "NC_045512.2.gff3"
SARS2_FASTA = ANNOTATION / "NC_045512.2.fasta"
INFLUENZA_GB = ANNOTATION / "influenza_A_PR8_segments.gb"


def by_id(annotation):
    return {cds["CDS_ID"]: cds for cds in annotation["cds"]}


class SyntheticAnnotationTest(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.path = Path(self.folder.name)

    def test_format_detection(self):
        self.assertEqual(detect_annotation_format(write_gff3(self.path)), GFF3_FORMAT)
        self.assertEqual(detect_annotation_format(write_genbank(self.path)), GENBANK_FORMAT)
        fasta = write_fasta(self.path, "x.fasta", "SEGA.1", SEGA)
        self.assertIsNone(detect_annotation_format(fasta))
        with self.assertRaises(ValueError):
            load_annotation(fasta)

    def test_gff3_and_genbank_give_the_same_model(self):
        gff = by_id(load_annotation(write_gff3(self.path), genetic_code=1))
        genbank = by_id(load_annotation(write_genbank(self.path)))
        for protein in ["P1", "P2", "P3", "P4", "P7"]:
            self.assertEqual(gff[protein]["Parts"], genbank[protein]["Parts"], protein)
            self.assertEqual(gff[protein]["Strand"], genbank[protein]["Strand"], protein)
            self.assertEqual(gff[protein]["Coding_Sequence"], genbank[protein]["Coding_Sequence"], protein)
            self.assertEqual(gff[protein]["Status"], SUPPORTED, protein)
            self.assertEqual(genbank[protein]["Status"], SUPPORTED, protein)
        # Hand-derived coding sequences (the pieces in annotation_fixtures.py).
        self.assertEqual(gff["P1"]["Coding_Sequence"], "ATGGATTGGCCATGCTAA")
        self.assertEqual(gff["P2"]["Coding_Sequence"], "ATGTTTGCATACCGTTAG")
        self.assertEqual(gff["P7"]["Coding_Sequence"], "ATGCCAGGCTTTTAG")
        self.assertEqual(gff["P4"]["Coding_Sequence"], "ATGGCTNCAAATTAA")
        # Transcription order: the minus-strand spliced CDS starts at the higher part.
        self.assertEqual(gff["P7"]["Parts"], [(76, 81), (61, 69)])
        self.assertIn("Translation of the reference coding sequence matches", genbank["P1"]["Verification"])
        self.assertIn("Internal checks only", gff["P1"]["Verification"])
        self.assertIn("1 codon(s) contain ambiguous reference bases", gff["P4"]["Verification"])

    def test_segments_and_organism(self):
        annotation = load_annotation(write_genbank(self.path))
        self.assertEqual(set(annotation["records"]), {"SEGA.1", "SEGB.1"})
        self.assertEqual(annotation["records"]["SEGA.1"]["Segment"], "A")
        self.assertEqual(annotation["records"]["SEGB.1"]["Segment"], "B")
        self.assertEqual(annotation["records"]["SEGB.1"]["Organism"], "Synthetic test virus")
        gff = load_annotation(write_gff3(self.path), genetic_code=1)
        self.assertEqual(gff["records"]["SEGB.1"]["Segment"], "B")
        self.assertEqual(gff["records"]["SEGA.1"]["Length"], 85)

    def test_coordinate_conversion_one_based_inclusive(self):
        p1 = by_id(load_annotation(write_genbank(self.path)))["P1"]
        self.assertEqual(p1["Parts"], [(10, 27)])  # Biopython 9..27 (0-based, end-exclusive) -> 10..27
        self.assertEqual(coding_offsets(p1, 10), [0])  # first base of the start codon
        self.assertEqual(coding_offsets(p1, 27), [17])  # last base of the stop codon
        self.assertEqual(coding_offsets(p1, 9), [])
        self.assertEqual(coding_offsets(p1, 28), [])
        p2 = by_id(load_annotation(write_genbank(self.path)))["P2"]
        self.assertEqual(coding_offsets(p2, 51), [0])  # minus strand: offset counts from the high end
        self.assertEqual(coding_offsets(p2, 34), [17])
        p4 = by_id(load_annotation(write_gff3(self.path), genetic_code=1))["P4"]
        self.assertEqual(coding_offsets(p4, 11), [4])
        self.assertEqual(coding_offsets(p4, 18), [5])  # first base after the intron continues codon 2
        self.assertEqual(coding_offsets(p4, 14), [])  # intron

    def test_unsupported_structures_are_marked_with_reasons(self):
        annotation = by_id(load_annotation(write_gff3(self.path), genetic_code=1))
        self.assertEqual(annotation["P5"]["Status"], UNSUPPORTED)
        self.assertIn("Internal stop codon", annotation["P5"]["Status_Detail"])
        self.assertEqual(annotation["P6"]["Status"], UNSUPPORTED)
        self.assertIn("RNA editing", annotation["P6"]["Status_Detail"])
        self.assertIsNone(annotation["P5"]["Coding_Sequence"])

    def test_inconsistent_gff3_phase_is_unsupported(self):
        bad = by_id(load_annotation(write_gff3(self.path, p4_phase=0), genetic_code=1))["P4"]
        self.assertEqual(bad["Status"], UNSUPPORTED)
        self.assertIn("phase", bad["Status_Detail"])
        good = by_id(load_annotation(write_gff3(self.path, p4_phase=1), genetic_code=1))["P4"]
        self.assertEqual(good["Status"], SUPPORTED)

    def test_translation_disagreeing_with_record_is_unsupported(self):
        p1 = by_id(load_annotation(write_genbank(self.path, p1_translation="MDWPA")))["P1"]
        self.assertEqual(p1["Status"], UNSUPPORTED)
        self.assertIn("residue 5", p1["Status_Detail"])

    def test_partial_out_of_order_and_overlapping_parts_are_unsupported(self):
        text = write_gff3(self.path, with_fasta=False).read_text(encoding="utf-8")
        extra = ("SEGA.1\ttest\tCDS\t10\t27\t.\t+\t0\tID=cds-PART;protein_id=PART;partial=true\n"
                 "SEGA.1\ttest\tCDS\t10\t21\t.\t+\t0\tID=cds-OVER;protein_id=OVER\n"
                 "SEGA.1\ttest\tCDS\t19\t27\t.\t+\t0\tID=cds-OVER;protein_id=OVER\n")
        path = self.path / "extra.gff3"
        path.write_text(text + extra + f"##FASTA\n>SEGA.1\n{SEGA}\n", encoding="utf-8")
        annotation = by_id(load_annotation(path, genetic_code=1))
        self.assertIn("Partial", annotation["PART"]["Status_Detail"])
        self.assertIn("overlap without a documented ribosomal slippage", annotation["OVER"]["Status_Detail"])
        from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        circular = SeqRecord(Seq(SEGA), id="CIRC.1", name="CIRC", description="origin-spanning test")
        circular.annotations["molecule_type"] = "DNA"
        circular.features = [SeqFeature(CompoundLocation([FeatureLocation(60, 85, strand=1),
                                                          FeatureLocation(0, 5, strand=1)]),
                                        type="CDS", qualifiers={"protein_id": ["WRAP"]})]
        wrapped = by_id(load_annotation(write_genbank(self.path, "wrap.gb", extra_records=[circular])))
        self.assertIn("not in genomic order", wrapped["WRAP"]["Status_Detail"])

    def test_gff3_without_sequence_uses_only_a_matching_reference(self):
        path = write_gff3(self.path, with_fasta=False)
        alone = load_annotation(path, genetic_code=1)
        self.assertIsNone(alone["records"]["SEGA.1"]["Sequence"])
        self.assertIn("No reference sequence", by_id(alone)["P1"]["Status_Detail"])
        matching = load_reference_fasta(write_fasta(self.path, "a.fasta", "SEGA.1", SEGA))
        used = load_annotation(path, matching, genetic_code=1)
        self.assertEqual(by_id(used)["P1"]["Status"], SUPPORTED)
        self.assertIn("matching reference FASTA is used", used["records"]["SEGA.1"]["Sequence_Check"])
        shorter = load_reference_fasta(write_fasta(self.path, "b.fasta", "SEGA.1", SEGA[:80]))
        refused = load_annotation(path, shorter, genetic_code=1)
        self.assertFalse(refused["records"]["SEGA.1"]["Compatible"])
        self.assertEqual(by_id(refused)["P1"]["Status"], UNSUPPORTED)

    def test_sequence_must_be_identical_to_an_uploaded_reference(self):
        changed = SEGA[:20] + ("A" if SEGA[20] != "A" else "C") + SEGA[21:]
        reference = load_reference_fasta(write_fasta(self.path, "c.fasta", "SEGA.1", changed))
        annotation = load_annotation(write_genbank(self.path), reference)
        self.assertFalse(annotation["records"]["SEGA.1"]["Compatible"])
        self.assertIn("differs from the uploaded reference FASTA", annotation["records"]["SEGA.1"]["Sequence_Check"])
        other = load_reference_fasta(write_fasta(self.path, "d.fasta", "OTHER.1", SEGA))
        unrelated = load_annotation(write_genbank(self.path), other)
        self.assertTrue(any("not described by the annotation" in text for text in unrelated["notes"]))

    def test_source_table_documents_every_cds(self):
        table = source_table(load_annotation(write_gff3(self.path), genetic_code=1))
        self.assertEqual(len(table), 7)
        self.assertEqual(table.set_index("CDS_ID").loc["P7", "Parts_1_Based"], "76..81,61..69")
        self.assertEqual(table.set_index("CDS_ID").loc["P4", "Phases"], "0,1")
        self.assertEqual(table.set_index("CDS_ID").loc["P2", "Strand"], "-")
        self.assertEqual(set(table["Status"]), {SUPPORTED, UNSUPPORTED})


class RealAnnotationTest(unittest.TestCase):
    """Official NCBI RefSeq records (real)."""

    @classmethod
    def setUpClass(cls):
        cls.reference = load_reference_fasta(SARS2_FASTA)
        cls.genbank = load_annotation(SARS2_GB, cls.reference)
        cls.gff = load_annotation(SARS2_GFF, cls.reference, genetic_code=1)

    def test_reference_identity(self):
        record = self.genbank["records"]["NC_045512.2"]
        self.assertEqual(record["Accession"], "NC_045512")
        self.assertEqual(record["Version"], "2")
        self.assertEqual(record["Length"], 29903)
        self.assertEqual(record["Organism"], "Severe acute respiratory syndrome coronavirus 2")
        self.assertEqual(record["Sequence_Check"], "Identical to the uploaded reference FASTA")

    def test_every_cds_reproduces_the_ncbi_protein(self):
        cds = self.genbank["cds"]
        self.assertEqual(len(cds), 12)
        for item in cds:
            self.assertEqual(item["Status"], SUPPORTED, item["CDS_ID"])
            self.assertTrue(item["Verification"].startswith("Translation of the reference coding sequence matches"))
        pp1ab = by_id(self.genbank)["YP_009724389.1"]
        self.assertEqual(pp1ab["Parts"], [(266, 13468), (13468, 21555)])
        self.assertTrue(pp1ab["Slippage"])
        self.assertEqual(len(pp1ab["Record_Translation"]), 7096)

    def test_gff3_matches_genbank(self):
        genbank, gff = by_id(self.genbank), by_id(self.gff)
        self.assertEqual(set(genbank), set(gff))
        for key in genbank:
            self.assertEqual(genbank[key]["Parts"], gff[key]["Parts"], key)
            self.assertEqual(genbank[key]["Coding_Sequence"], gff[key]["Coding_Sequence"], key)
            self.assertEqual(gff[key]["Status"], SUPPORTED, key)

    def test_frameshift_base_is_read_twice(self):
        pp1ab = by_id(self.genbank)["YP_009724389.1"]
        self.assertEqual(coding_offsets(pp1ab, 13468), [13202, 13203])
        self.assertEqual(coding_offsets(pp1ab, 13467), [13201])
        self.assertEqual(coding_offsets(pp1ab, 13469), [13204])

    def test_mature_peptides_match_ncbi_protein_coordinates(self):
        # NCBI's GFF3 IDs give each mature peptide's protein coordinates, e.g.
        # "id-YP_009724389.1:181..818"; they were not used to compute these spans.
        expected = {}
        for line in SARS2_GFF.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) == 9 and fields[2] == "mature_protein_region_of_CDS" and "YP_009724389.1:" in fields[8]:
                span = fields[8].split("ID=id-YP_009724389.1:", 1)[1].split(";", 1)[0]
                first, last = (int(value) for value in span.split(".."))
                expected[(first, last)] = True
        for annotation in [self.genbank, self.gff]:
            spans = {(item["First_Codon"], item["Last_Codon"])
                     for item in by_id(annotation)["YP_009724389.1"]["Mature_Peptides"]}
            self.assertEqual(spans, set(expected))
        rdrp = next(item for item in by_id(self.genbank)["YP_009724389.1"]["Mature_Peptides"]
                    if item["Product"] == "RNA-dependent RNA polymerase")
        self.assertEqual((rdrp["First_Codon"], rdrp["Last_Codon"]), (4393, 5324))  # spans the frameshift

    def test_influenza_segments_splicing_and_frameshift(self):
        annotation = load_annotation(INFLUENZA_GB)
        self.assertEqual(len(annotation["records"]), 8)
        self.assertEqual(sorted(record["Segment"] for record in annotation["records"].values()),
                         [str(number) for number in range(1, 9)])
        self.assertEqual(len(annotation["cds"]), 12)
        self.assertTrue(all(item["Status"] == SUPPORTED for item in annotation["cds"]))
        genes = {item["Gene"]: item for item in annotation["cds"]}
        self.assertEqual(genes["M2"]["Parts"], [(26, 51), (740, 1007)])  # spliced
        self.assertEqual(genes["NEP"]["Parts"], [(27, 56), (529, 864)])
        self.assertEqual(genes["PA-X"]["Parts"], [(25, 597), (599, 784)])  # +1 ribosomal frameshift
        self.assertTrue(genes["PA-X"]["Slippage"])
        self.assertEqual(genes["M2"]["Reference_ID"], "NC_002016.1")

    @unittest.skipUnless(os.environ.get("ANNOTATION_MPOX_GENBANK"),
                         "set ANNOTATION_MPOX_GENBANK to NCBI NC_063383.1 (GenBank) to check minus-strand genes")
    def test_mpox_minus_strand_genes_reproduce_ncbi_proteins(self):
        annotation = load_annotation(os.environ["ANNOTATION_MPOX_GENBANK"])
        minus = [item for item in annotation["cds"] if item["Strand"] == "-"]
        self.assertGreater(len(minus), 50)
        self.assertTrue(all(item["Status"] == SUPPORTED for item in annotation["cds"]))
        sequence = annotation["records"]["NC_063383.1"]["Sequence"]
        gene = minus[0]
        start, end = gene["Parts"][0]
        self.assertEqual(gene["Coding_Sequence"][:3], rc(sequence[end - 3:end]))


if __name__ == "__main__":
    unittest.main()
