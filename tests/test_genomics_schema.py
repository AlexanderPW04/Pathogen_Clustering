"""Tests for the shared reference-based sample and variant representation."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from genomics.pooled_observations import (
    base_frequencies,
    build_depth_lookup,
    coverage_from_positions,
    depth_at,
    load_position_table,
    pooled_position_status,
    pooled_sample_records,
    pooled_variant_records,
)
from genomics.references import (
    build_reference_registry,
    load_reference_fasta,
    reference_base,
    references_compatible,
)
from genomics.unified_schema import (
    COVERAGE_UNKNOWN,
    COVERED_NOT_DETECTED,
    INDIVIDUAL,
    LOW_COVERAGE,
    NOT_COVERED,
    POOLED,
    make_variant_table,
    validate_variant_table,
    variant_key,
)


def write_text(folder, name, text):
    path = Path(folder) / name
    path.write_text(text, encoding="utf-8")
    return path


def snv_rows(**changes):
    row = {
        "Sample_ID": "W1", "Reference_Genome": "REF.1", "Position": 5,
        "Ref": "A", "Alt": "G", "Alt_Frequency": 0.25, "Read_Depth": 80,
        "Collection_Date": pd.Timestamp("2026-01-05"), "Location_Name": "Site A",
    }
    row.update(changes)
    return pd.DataFrame([row])


class ReferenceTest(unittest.TestCase):
    def test_reference_is_loaded_with_one_based_positions(self):
        with tempfile.TemporaryDirectory() as folder:
            path = write_text(folder, "ref.fasta", ">REF.1 invented\nacgtn\nACGT\n")
            reference = load_reference_fasta(path)
        self.assertEqual(reference["Reference_ID"], "REF.1")
        self.assertEqual(reference["Accession"], "REF")
        self.assertEqual(reference["Version"], "1")
        self.assertEqual(reference["Length"], 9)
        self.assertEqual(reference_base(reference, 1), "A")
        self.assertEqual(reference_base(reference, 9), "T")
        with self.assertRaises(ValueError):
            reference_base(reference, 0)

    def test_reference_rejects_multiple_records_and_gaps(self):
        with tempfile.TemporaryDirectory() as folder:
            two = write_text(folder, "two.fasta", ">A\nACGT\n>B\nACGT\n")
            gapped = write_text(folder, "gap.fasta", ">A\nAC-GT\n")
            with self.assertRaisesRegex(ValueError, "exactly one"):
                load_reference_fasta(two)
            with self.assertRaisesRegex(ValueError, "only A, C, G, T"):
                load_reference_fasta(gapped)

    def test_aliases_need_identical_supplied_sequences(self):
        with tempfile.TemporaryDirectory() as folder:
            first = load_reference_fasta(write_text(folder, "a.fasta", ">MN908947.3\nACGTACGT\n"))
            same = load_reference_fasta(write_text(folder, "b.fasta", ">NC_045512.2\nACGTACGT\n"))
            other = load_reference_fasta(write_text(folder, "c.fasta", ">OTHER.1\nACGTACGA\n"))
        registry = build_reference_registry([first, same, other])

        ok, reason = references_compatible("MN908947.3", "NC_045512.2", registry)
        self.assertTrue(ok)
        self.assertIn("identical", reason)

        ok, _ = references_compatible("MN908947.3", "OTHER.1", registry)
        self.assertFalse(ok)

        # Without supplied sequences, labels alone are not accepted as aliases.
        ok, reason = references_compatible("MN908947.3", "NC_045512.2")
        self.assertFalse(ok)
        self.assertIn("not", reason)

        ok, reason = references_compatible("NC_045512.1", "NC_045512.2")
        self.assertFalse(ok)
        self.assertIn("versions", reason)

        self.assertTrue(references_compatible("REF.1", "REF.1")[0])
        self.assertFalse(references_compatible("", "REF.1")[0])


class VariantTableTest(unittest.TestCase):
    def test_variant_key_and_validation(self):
        self.assertEqual(variant_key("REF.1", 23403, "A", "G"), "REF.1:23403:A>G")
        good = make_variant_table([{
            "Variant_Key": "REF.1:5:A>G", "Reference_ID": "REF.1", "Position": 5,
            "Ref": "A", "Alt": "G", "Sample_ID": "W1", "Sample_Type": POOLED,
            "Alt_Frequency": 0.0, "Read_Depth": 30, "Alt_Depth": 0,
            "Published_Annotation": "ORF1ab",
        }])
        self.assertEqual(validate_variant_table(good), [])
        self.assertIn("Published_Annotation", good.columns)
        self.assertEqual(good.loc[0, "Alt_Frequency"], 0.0)

        bad = make_variant_table([
            {"Variant_Key": "REF.1:0:A>G", "Reference_ID": "REF.1", "Position": 0,
             "Ref": "A", "Alt": "G", "Sample_ID": "W1", "Sample_Type": POOLED},
            {"Variant_Key": "wrong", "Reference_ID": "REF.1", "Position": 3,
             "Ref": "A", "Alt": "A", "Sample_ID": "S1", "Sample_Type": INDIVIDUAL,
             "Alt_Frequency": 1.0},
            {"Variant_Key": "REF.1:4:C>T", "Reference_ID": "REF.1", "Position": 4,
             "Ref": "C", "Alt": "T", "Sample_ID": "W1", "Sample_Type": POOLED,
             "Alt_Frequency": 1.5, "Read_Depth": 10, "Alt_Depth": 12},
        ])
        errors = " ".join(validate_variant_table(bad))
        self.assertIn("1-based", errors)
        self.assertIn("identical", errors)
        self.assertIn("Variant_Key should be", errors)
        self.assertIn("do not have an allele frequency", errors)
        self.assertIn("fraction from 0 to 1", errors)
        self.assertIn("larger than Read_Depth", errors)


class PooledObservationTest(unittest.TestCase):
    def test_processed_snvs_become_shared_records_without_filling_gaps(self):
        data = pd.concat([
            snv_rows(),
            snv_rows(Sample_ID="W2", Alt_Frequency=0.0, Read_Depth=pd.NA,
                     Published_Annotation="ORF1ab"),
        ], ignore_index=True)
        data["Read_Depth"] = data["Read_Depth"].astype("Int64")
        variants = pooled_variant_records(data, "unit test")
        self.assertEqual(validate_variant_table(variants), [])
        self.assertEqual(variants.loc[0, "Variant_Key"], "REF.1:5:A>G")
        self.assertEqual(variants.loc[1, "Alt_Frequency"], 0.0)  # supplied zero stays zero
        self.assertTrue(pd.isna(variants.loc[1, "Read_Depth"]))  # unknown stays unknown
        self.assertTrue(pd.isna(variants.loc[0, "Alt_Depth"]))
        self.assertEqual(variants.loc[1, "Published_Annotation"], "ORF1ab")

        samples = pooled_sample_records(data, "unit test")
        self.assertEqual(list(samples["Sample_ID"]), ["W1", "W2"])
        self.assertTrue(samples["Genome_Coverage_Percent"].isna().all())
        self.assertIn("coverage elsewhere unknown", samples.loc[0, "Coverage_Information"])

    def test_position_table_separates_zero_low_unknown_and_covered(self):
        text = (
            "Sample_ID,Reference_ID,Position,A,C,G,T\n"
            "W1,REF.1,1,30,0,0,0\n"
            "W1,REF.1,2,0,0,0,0\n"
            "W1,REF.1,3,0,5,0,0\n"
            "W1,REF.1,4,0,0,15,15\n"
        )
        with tempfile.TemporaryDirectory() as folder:
            table = load_position_table(write_text(folder, "counts.csv", text))
        self.assertTrue(table.attrs["has_base_counts"])
        self.assertEqual(list(table["Depth"]), [30, 0, 5, 30])

        lookup = build_depth_lookup(table)
        self.assertEqual(depth_at(lookup, "W1", "REF.1", 2), 0)
        self.assertIsNone(depth_at(lookup, "W1", "REF.1", 9))
        self.assertIsNone(depth_at(lookup, "W9", "REF.1", 1))

        self.assertEqual(pooled_position_status(lookup, "W1", "REF.1", 1, 20)[0],
                         COVERED_NOT_DETECTED)
        self.assertEqual(pooled_position_status(lookup, "W1", "REF.1", 2, 20)[0], NOT_COVERED)
        self.assertEqual(pooled_position_status(lookup, "W1", "REF.1", 3, 20)[0], LOW_COVERAGE)
        self.assertEqual(pooled_position_status(lookup, "W1", "REF.1", 7, 20)[0],
                         COVERAGE_UNKNOWN)

        frequencies = base_frequencies(table)
        self.assertAlmostEqual(frequencies.loc[3, "G_Frequency"], 0.5)
        self.assertTrue(pd.isna(frequencies.loc[1, "A_Frequency"]))  # zero depth: no frequency

        percent, mean_depth, _ = coverage_from_positions(table, "W1", "REF.1", 4, 10)
        self.assertAlmostEqual(percent, 50.0)
        self.assertAlmostEqual(mean_depth, 16.25)
        percent, _, note = coverage_from_positions(table, "W1", "REF.1", 10, 10)
        self.assertIsNone(percent)  # incomplete tables are not treated as zero depth
        self.assertIn("4 of 10", note)

    def test_position_table_rejects_bad_values(self):
        with tempfile.TemporaryDirectory() as folder:
            no_depth = write_text(folder, "a.csv", "Sample_ID,Reference_ID,Position\nW1,R,1\n")
            negative = write_text(folder, "b.csv", "Sample_ID,Reference_ID,Position,Depth\nW1,R,1,-2\n")
            zero_pos = write_text(folder, "c.csv", "Sample_ID,Reference_ID,Position,Depth\nW1,R,0,5\n")
            duplicate = write_text(
                folder, "d.csv", "Sample_ID,Reference_ID,Position,Depth\nW1,R,1,5\nW1,R,1,6\n"
            )
            with self.assertRaisesRegex(ValueError, "Depth column"):
                load_position_table(no_depth)
            with self.assertRaisesRegex(ValueError, "negative"):
                load_position_table(negative)
            with self.assertRaisesRegex(ValueError, "1-based"):
                load_position_table(zero_pos)
            with self.assertRaisesRegex(ValueError, "only once"):
                load_position_table(duplicate)


if __name__ == "__main__":
    unittest.main()
