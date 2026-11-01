"""Tests for reference-coordinate variant profiles from individual sequences."""

import shutil
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from genomics.references import sequence_checksum
from genomics.sequence_variants import (
    CALL_AMBIGUOUS,
    CALL_DELETION,
    CALL_OUTSIDE,
    CALL_REFERENCE,
    CALL_VARIANT,
    build_individual_profiles,
    individual_position_status,
    map_alignment_to_reference,
    summarize_calls,
)
from genomics.unified_schema import validate_variant_table


REF = "ATGCGTACGTTAGCCGATCGATTACGGCTAGCTAGGCTAA"  # 40 invented bases
REFERENCE = {
    "Reference_ID": "SYNREF.1",
    "Sequence": REF,
    "Length": len(REF),
    "Checksum": sequence_checksum(REF),
}


def change(sequence, position, base):
    """Replace a 1-based position."""

    return sequence[: position - 1] + base + sequence[position:]


class AlignmentMappingTest(unittest.TestCase):
    def test_substitutions_use_one_based_reference_positions(self):
        sample = change(change(REF, 1, "C"), 40, "G")
        mapping = map_alignment_to_reference(REF, sample)
        self.assertEqual(mapping["substitutions"], [(1, "A", "C"), (40, "A", "G")])
        self.assertEqual(mapping["calls"][0], CALL_VARIANT)
        self.assertEqual(mapping["calls"][1], CALL_REFERENCE)
        self.assertEqual(len(mapping["calls"]), 40)

    def test_partial_sequence_is_not_assumed_to_start_at_position_one(self):
        sample_row = "-----" + REF[5:35] + "-----"
        mapping = map_alignment_to_reference(REF, sample_row)
        self.assertEqual(mapping["calls"], CALL_OUTSIDE * 5 + CALL_REFERENCE * 30 + CALL_OUTSIDE * 5)
        summary = summarize_calls(mapping["calls"])
        self.assertEqual((summary["span_start"], summary["span_end"]), (6, 35))
        self.assertEqual(summary["callable_positions"], 30)
        self.assertEqual(mapping["substitutions"], [])

    def test_deletions_and_insertions_are_kept_apart_from_substitutions(self):
        deleted = REF[:20] + "---" + REF[23:]
        mapping = map_alignment_to_reference(REF, deleted)
        self.assertEqual(mapping["calls"][20:23], CALL_DELETION * 3)
        self.assertEqual(mapping["deletions"], [(21, 23)])
        self.assertEqual(mapping["substitutions"], [])

        reference_row = REF[:20] + "--" + REF[20:]
        sample_row = REF[:20] + "GG" + REF[20:]
        mapping = map_alignment_to_reference(reference_row, sample_row)
        self.assertEqual(len(mapping["calls"]), 40)  # alignment columns are not positions
        self.assertEqual(mapping["calls"], CALL_REFERENCE * 40)
        self.assertEqual(mapping["insertions"], [(20, "GG")])

    def test_n_and_iupac_codes_are_not_callable(self):
        sample = change(change(REF, 5, "N"), 6, "R")
        mapping = map_alignment_to_reference(REF, sample)
        self.assertEqual(mapping["calls"][4:6], CALL_AMBIGUOUS * 2)
        self.assertEqual(mapping["substitutions"], [])
        self.assertEqual(summarize_calls(mapping["calls"])["ambiguous_positions"], 2)


def fake_reference_mafft(mafft_path, input_path):
    """Stand-in for MAFFT: pad a substring of REF to its reference position."""

    lines = Path(input_path).read_text(encoding="utf-8").split()
    sample = lines[3]
    if len(sample) == len(REF):
        sample_row = sample
    elif sample in REF:
        start = REF.index(sample)
        sample_row = "-" * start + sample + "-" * (len(REF) - start - len(sample))
    else:
        sample_row = (sample + "-" * len(REF))[: len(REF)]
    return f">REFERENCE\n{REF}\n>SAMPLE\n{sample_row}\n"


class IndividualProfileTest(unittest.TestCase):
    def sample_table(self):
        return pd.DataFrame({
            "Sample_ID": ["FULL", "PARTIAL", "AMBIGUOUS", "UNRELATED"],
            "Sequence": [
                change(REF, 10, "A"),
                REF[5:35],
                change(change(REF, 10, "A"), 12, "N"),
                "T" * 40,
            ],
            "Collection_Date": pd.to_datetime(["2026-03-01"] * 4),
            "Location_Name": ["Tempe, Arizona"] * 4,
            "Latitude": [33.4] * 4,
            "Longitude": [-111.9] * 4,
            "Coordinate_Source": ["Provided"] * 4,
            "Geographic_Precision": ["Provided coordinates (precision unverified)"] * 4,
            "Gene_Segment": ["Complete genome"] * 4,
        })

    def test_profiles_are_built_only_for_trustworthy_mappings(self):
        with patch("genomics.sequence_variants.run_mafft", side_effect=fake_reference_mafft):
            result = build_individual_profiles(self.sample_table(), REFERENCE, "mafft", "test")

        variants = result["variants"]
        self.assertEqual(validate_variant_table(variants), [])
        self.assertEqual(sorted(variants["Variant_Key"].unique()), ["SYNREF.1:10:T>A"])
        self.assertEqual(sorted(variants["Sample_ID"]), ["AMBIGUOUS", "FULL"])
        self.assertTrue(variants["Alt_Frequency"].isna().all())

        samples = result["samples"].set_index("Sample_ID")
        self.assertEqual(samples.loc["PARTIAL", "Sequenced_Span"], "6-35")
        self.assertIn("not established", samples.loc["UNRELATED", "Variant_Profile_Status"])
        self.assertNotIn("UNRELATED", result["calls"])

        # The partial sequence did not cover position 10, so it is unknown there.
        self.assertEqual(individual_position_status(result["calls"], "PARTIAL", 2), CALL_OUTSIDE)
        self.assertEqual(individual_position_status(result["calls"], "PARTIAL", 10), CALL_REFERENCE)
        self.assertEqual(individual_position_status(result["calls"], "AMBIGUOUS", 12), CALL_AMBIGUOUS)
        self.assertIsNone(individual_position_status(result["calls"], "UNRELATED", 10))


@unittest.skipUnless(shutil.which("mafft"), "Real MAFFT is not installed on PATH.")
class RealMafftReferenceAlignmentTest(unittest.TestCase):
    def test_real_mafft_maps_a_partial_sequence_with_a_substitution(self):
        sample = change(REF, 20, "T" if REF[19] != "T" else "A")[3:37]
        data = pd.DataFrame({
            "Sample_ID": ["REAL"], "Sequence": [sample],
            "Collection_Date": [pd.Timestamp("2026-03-01")], "Location_Name": ["X"],
            "Gene_Segment": ["Complete genome"],
        })
        result = build_individual_profiles(data, REFERENCE, shutil.which("mafft"), "real")
        self.assertEqual(result["samples"].loc[0, "Sequenced_Span"], "4-37")
        self.assertEqual(list(result["variants"]["Position"]), [20])


if __name__ == "__main__":
    unittest.main()
