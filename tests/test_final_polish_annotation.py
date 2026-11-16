"""Final polish, Milestone 1: gene-track labels, the influenza gene "NA", and the genetic code.

1A  The genome-position gene track uses verified CDS coordinates from an uploaded
    annotation; without one, only source labels forming one uninterrupted run are drawn
    (a label such as "UTR" that recurs between genes is left unannotated).
1B  Text such as the influenza gene "NA" stays text in every download; only empty cells
    are missing.
1C  The genetic code is never assumed: it is stated by the annotation, verified against the
    record's own protein, or confirmed by the user; otherwise no amino-acid consequence.

Expected values are worked out by hand (codon tables 1 and 4 from NCBI), independently of
the code under test. All sequences are SYNTHETIC except the bundled NCBI NC_045512.2 files.
"""

import io
import tempfile
import unittest
from pathlib import Path

import openpyxl
import pandas as pd

from compare_genetic_distances import read_pair_table
from genomics.genome_annotation import (
    SUPPORTED,
    TABLE_CONFIRMED,
    TABLE_STATED,
    TABLE_UNVERIFIED_DETAIL,
    TABLE_VERIFIED,
    UNSUPPORTED,
    load_annotation,
    parse_genetic_code,
    source_table,
)
from genomics.integrated_analysis import annotation_file_bands, annotation_ranges, gene_track
from genomics.references import load_reference_fasta
from genomics.report_exports import StoredTables, build_workbook
from genomics.snv_consequences import STOP_GAINED, SYNONYMOUS, UNSUPPORTED_CONSEQUENCE, annotate_variants

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION = ROOT / "examples" / "annotation"

# A 16-base SYNTHETIC genome with one CDS at 3..14: ATG TGG AAA TAA (Met-Trp-Lys-stop).
GENOME = "CC" + "ATGTGGAAATAA" + "GG"
# Position 8 is the third base of codon 2 (TGG). G>A makes TGA: a stop under the standard code
# (table 1) but tryptophan under the Mycoplasma/Spiroplasma code (table 4), so W2* versus W2W.
SNV = pd.DataFrame([{"Variant_Key": "CODE.1:8:G>A", "Reference_ID": "CODE.1", "Position": 8, "Ref": "G",
                     "Alt": "A"}])


def write(folder, name, text):
    path = Path(folder) / name
    path.write_text(text, encoding="utf-8")
    return path


def gff3(folder, attributes=""):
    return write(folder, "code.gff3", "##gff-version 3\n##sequence-region CODE.1 1 16\n"
                                      f"CODE.1\ttest\tCDS\t3\t14\t.\t+\t0\tID=cds-C1;gene=codeA;protein_id=C1"
                                      f"{attributes}\n##FASTA\n>CODE.1\n{GENOME}\n")


def genbank(folder, qualifiers):
    lines = "\n".join(f"                     /{text}" for text in qualifiers)
    return write(folder, "code.gb", f"""LOCUS       CODE                      16 bp    DNA     linear   SYN 01-JAN-2026
DEFINITION  Synthetic genetic-code test.
ACCESSION   CODE
VERSION     CODE.1
FEATURES             Location/Qualifiers
     source          1..16
                     /organism="Synthetic test organism"
     CDS             3..14
                     /gene="codeA"
                     /protein_id="C1.1"
{lines}
ORIGIN
        1 {GENOME.lower()}
//
""")


def consequence(annotation):
    rows, summary = annotate_variants(SNV, annotation)
    coding = rows[rows["Gene"].eq("codeA")]
    return (coding["Consequence"].iloc[0] if len(coding) else None,
            coding["Amino_Acid_Change"].iloc[0] if len(coding) else None,
            summary["Consequence_Category"].iloc[0])


class GeneticCodeTest(unittest.TestCase):
    """1C: no silent standard code."""

    def setUp(self):
        self.folder = tempfile.mkdtemp()

    def test_gff3_without_a_stated_code_predicts_no_amino_acid_change(self):
        annotation = load_annotation(gff3(self.folder))
        cds = annotation["cds"][0]
        self.assertEqual(cds["Status"], UNSUPPORTED)
        self.assertEqual(cds["Status_Detail"], TABLE_UNVERIFIED_DETAIL)
        self.assertTrue(any("state no genetic code" in text for text in annotation["notes"]))
        category, change, summary = consequence(annotation)
        self.assertFalse(change)  # no codon is translated for an unsupported CDS
        self.assertNotEqual(summary, STOP_GAINED)  # the false "W2*" of the old default is gone
        self.assertEqual(summary, UNSUPPORTED_CONSEQUENCE)

    def test_the_users_confirmed_code_is_used_and_recorded(self):
        four = load_annotation(gff3(self.folder), genetic_code=4)
        self.assertEqual(four["cds"][0]["Translation_Table"], 4)
        self.assertTrue(four["cds"][0]["Translation_Table_Source"].startswith(TABLE_CONFIRMED))
        self.assertEqual(consequence(four)[:2], (SYNONYMOUS, "W2W"))
        one = load_annotation(gff3(self.folder), genetic_code=1)
        self.assertEqual(consequence(one)[:2], (STOP_GAINED, "W2*"))
        self.assertIn("Translation_Table_Source", source_table(four).columns)

    def test_a_stated_code_is_used_and_never_overridden(self):
        stated = load_annotation(gff3(self.folder, ";transl_table=4"), genetic_code=1)
        cds = stated["cds"][0]
        self.assertEqual(cds["Translation_Table"], 4)
        self.assertEqual(cds["Translation_Table_Source"], f"{TABLE_STATED} (transl_table=4)")
        self.assertEqual(consequence(stated)[:2], (SYNONYMOUS, "W2W"))
        genbank_stated = load_annotation(genbank(self.folder, ["transl_table=4", 'translation="MWK"']))
        self.assertEqual(genbank_stated["cds"][0]["Status"], SUPPORTED)
        self.assertEqual(consequence(genbank_stated)[:2], (SYNONYMOUS, "W2W"))

    def test_genbank_protein_verifies_the_standard_code(self):
        verified = load_annotation(genbank(self.folder, ['translation="MWK"']))
        cds = verified["cds"][0]
        self.assertEqual((cds["Status"], cds["Translation_Table_Source"]), (SUPPORTED, TABLE_VERIFIED))
        self.assertEqual(consequence(verified)[:2], (STOP_GAINED, "W2*"))
        # A record whose own protein disagrees with table 1 is not translated with it.
        mismatched = load_annotation(genbank(self.folder, ['translation="MRK"']))
        self.assertEqual(mismatched["cds"][0]["Status"], UNSUPPORTED)
        self.assertFalse(consequence(mismatched)[1])

    def test_genbank_without_protein_or_code_is_not_translated(self):
        bare = load_annotation(genbank(self.folder, []))
        self.assertEqual(bare["cds"][0]["Status_Detail"], TABLE_UNVERIFIED_DETAIL)

    def test_real_ncbi_genbank_keeps_its_annotations(self):
        reference = load_reference_fasta(ANNOTATION / "NC_045512.2.fasta")
        annotation = load_annotation(ANNOTATION / "NC_045512.2.gb", reference)
        self.assertTrue(all(cds["Status"] == SUPPORTED for cds in annotation["cds"]))
        self.assertTrue(all(cds["Translation_Table_Source"] == TABLE_VERIFIED for cds in annotation["cds"]))
        d614g = pd.DataFrame([{"Variant_Key": "NC_045512.2:23403:A>G", "Reference_ID": "NC_045512.2",
                               "Position": 23403, "Ref": "A", "Alt": "G"}])
        rows, _ = annotate_variants(d614g, annotation)
        self.assertEqual(rows["Amino_Acid_Change"].iloc[0], "D614G")

    def test_settings_value_is_validated(self):
        self.assertIsNone(parse_genetic_code(""))
        self.assertEqual(parse_genetic_code("11"), 11)
        for bad in ["0", "99", "x", "-1"]:
            with self.assertRaises(ValueError):
                parse_genetic_code(bad)


class GeneTrackTest(unittest.TestCase):
    """1A: verified bands, and no 'UTR' band across genes."""

    def test_recurring_source_labels_are_not_drawn_across_genes(self):
        # Fontenele et al. (2021) label every intergenic position 'UTR' (25385 ... 29557).
        summary = pd.DataFrame({"Position": [21564, 25354, 25385, 25394, 26215, 26500, 28275, 29532, 29557],
                                "Gene_Or_Annotation": ["S", "S", "UTR", "ORF3a", "ORF3a", "UTR", "N", "N", "UTR"]})
        bands = annotation_ranges(summary)
        self.assertNotIn("UTR", [band[0] for band in bands])
        self.assertEqual(bands, [("S", 21564, 25355), ("ORF3a", 25394, 26216), ("N", 28275, 29533)])

    def test_annotation_file_bands_are_the_ncbi_cds_coordinates(self):
        reference = load_reference_fasta(ANNOTATION / "NC_045512.2.fasta")
        annotation = load_annotation(ANNOTATION / "NC_045512.2.gb", reference)
        bands = {band[0]: band[1:] for band in annotation_file_bands(annotation, "NC_045512.2")}
        # NCBI NC_045512.2 CDS (1-based, inclusive): S 21563..25384, ORF3a 25393..26220, N 28274..29533,
        # ORF1ab join(266..13468,13468..21555) merged with ORF1a 266..13483.
        self.assertEqual(bands["S"], (21563, 25385))
        self.assertEqual(bands["ORF3a"], (25393, 26221))
        self.assertEqual(bands["N"], (28274, 29534))
        self.assertEqual(bands["ORF1ab"], (266, 21556))
        self.assertNotIn("UTR", bands)
        covered = [(start, end) for start, end in bands.values()]
        self.assertFalse(any(start <= 25389 < end for start, end in covered))  # 25385-25392 unannotated
        summary = pd.DataFrame({"Position": [25390], "Gene_Or_Annotation": ["UTR"]})
        _, caption = gene_track(summary, annotation, "NC_045512.2")
        self.assertIn("coding regions (CDS) from the uploaded annotation", caption)
        self.assertEqual(annotation_file_bands(annotation, "OTHER.1"), [])


class GeneNameTextTest(unittest.TestCase):
    """1B: 'NA' is a gene name, an empty cell is missing."""

    TABLE = ("Variant_Key,Genes,Amino_Acid_Changes,Annotation_Status\n"
             "SEG6:100:A>G,NA,NA:Gly454Asp,Annotated\n"
             "SEG4:50:C>T,HA,HA:Gln46Stop,Annotated\n"
             "SEG7:8:T>G,M1; M2,M1:Leu3Arg; M2:Leu3Arg,Annotated\n"
             "SEG8:20:A>C,NS,NS:Asp2Glu,Annotated\n"
             "SEG9:5:G>T,,,Not annotated\n")

    def test_stored_tables_keep_na_and_short_symbols_and_blank_is_missing(self):
        stored = StoredTables(lambda name: self.TABLE.encode() if name == "x.csv" else None)
        table = stored.get("x.csv")
        self.assertEqual(list(table["Genes"][:4]), ["NA", "HA", "M1; M2", "NS"])
        self.assertTrue(pd.isna(table["Genes"].iloc[4]))
        self.assertTrue(pd.isna(table["Amino_Acid_Changes"].iloc[4]))
        self.assertEqual(table["Amino_Acid_Changes"].iloc[0], "NA:Gly454Asp")

    def test_workbook_keeps_na(self):
        table = StoredTables(lambda name: self.TABLE.encode()).get("x.csv")
        content, _ = build_workbook({"sprint5_annotated_mutations.csv": table},
                                    {"findings": [], "notes": {}, "readiness": [], "labels": {}})
        book = openpyxl.load_workbook(io.BytesIO(content), read_only=True)
        sheet = next(name for name in book.sheetnames if name != "Summary")
        rows = list(book[sheet].iter_rows(values_only=True))
        genes = [row[1] for row in rows[1:]]
        self.assertEqual(genes, ["NA", "HA", "M1; M2", "NS", None])

    def test_sprint2_pair_tables_keep_na_segments(self):
        path = Path(tempfile.mkdtemp()) / "pairs.csv"
        path.write_text("Sample_ID_1,Sample_ID_2,Pathogen,Gene_Segment,P_Distance\n"
                        "001,NA,Influenza A virus,NA,0.01\nA,B,Influenza A virus,,0.02\n", encoding="utf-8")
        pairs = read_pair_table(path)
        self.assertEqual(list(pairs.iloc[0, :4]), ["001", "NA", "Influenza A virus", "NA"])
        self.assertEqual(pairs["Gene_Segment"].iloc[1], "")


if __name__ == "__main__":
    unittest.main()
