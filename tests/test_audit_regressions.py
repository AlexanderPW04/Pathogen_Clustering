"""Regression tests for defects found in the pre-Sprint 4 audit (SPRINT2_3_FULL_AUDIT.md).

Each expected value below was worked out by hand from the inputs, not copied
from the function under test. Issue IDs refer to the audit's issue register.
MAFFT is replaced by a pass-through (equal-length inputs) and GeoNames is off.
"""

import gzip
import io
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from add_temporal_geographic import enrich_pairs, haversine_km  # noqa: E402
from align_sequences import run_mafft  # noqa: E402
from compare_genetic_distances import compare_sequences  # noqa: E402
from genomics.file_formats import plain_copy, vcf_to_snv_table  # noqa: E402
from genomics.input_inspection import inspect_file, review_dataset, suggest_roles  # noqa: E402
from genomics.quality_control import apply_pooled_quality_control  # noqa: E402
from genomics.temporal_analysis import depth_regression  # noqa: E402
from genomics.variant_matching import usable_coordinates  # noqa: E402
from geocode_locations import is_state_only_location  # noqa: E402
from run_pipeline import output_paths, run_pipeline  # noqa: E402
from validate_inputs import validate_metadata  # noqa: E402
from wastewater.processed_snv import COMBINED_FORMAT, load_processed_snvs  # noqa: E402

SEQUENCE = "ACGTACGTACGTACGTACGTACGTACGT"  # 28 bases


def passthrough(command, **kwargs):
    class Result:
        stdout = Path(command[-1]).read_text()
    return Result()


def run_two_samples(ids, sequences, locations, dates=("2024-01-01", "2024-01-03"), extra=None):
    """Run the Sprint 2 pipeline on two samples; return (status, cluster rows)."""

    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        fasta = folder / "s.fasta"
        fasta.write_text("".join(f">{i}\n{s}\n" for i, s in zip(ids, sequences)))
        metadata = folder / "m.csv"
        pd.DataFrame({"Accession": ids, "Pathogen": "V", "Gene_Segment": "Segment 4",
                      "Collection_Date": list(dates), "Location_Name": locations,
                      "Host_Source": "Human", **(extra or {})}).to_csv(metadata, index=False)
        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}), \
                patch("align_sequences.shutil.which", return_value="mafft"), \
                patch("align_sequences.subprocess.run", side_effect=passthrough), \
                redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            status = run_pipeline(fasta, metadata, 0.05, 14, 50)
        if status != 0:
            return status, None
        return status, pd.read_csv(output_paths(metadata)[3], converters={"Sample_ID": str})


class Sprint2RegressionTest(unittest.TestCase):
    def test_s2_1_numeric_looking_ids_survive_the_csv_round_trip(self):
        status, clusters = run_two_samples(["001", "NA"], [SEQUENCE, SEQUENCE],
                                           ["Tempe, AZ", "Tempe, AZ"])
        self.assertEqual(status, 0)
        self.assertEqual(sorted(clusters["Sample_ID"]), ["001", "NA"])

    def test_s2_2_coarse_or_placeholder_locations_never_form_strict_clusters(self):
        for location in ["USA: California", "USA", "not provided", "USA: AZ"]:
            status, clusters = run_two_samples(["A", "B"], [SEQUENCE, SEQUENCE], [location, location])
            self.assertEqual(status, 0, location)
            self.assertTrue(clusters.empty, location)
        self.assertTrue(is_state_only_location("United States: New Mexico"))
        self.assertFalse(is_state_only_location("USA: AZ, Tempe"))
        self.assertFalse(is_state_only_location("Tempe, Arizona"))

    def test_s2_3_one_shared_called_base_is_not_genetic_relatedness(self):
        coordinates = {"Latitude": [33.0, 33.0], "Longitude": [-111.0, -111.0]}
        status, clusters = run_two_samples(["A", "B"], ["A" + "N" * 27, SEQUENCE],
                                           ["Site 1", "Site 1"], extra=coordinates)
        self.assertEqual(status, 0)
        self.assertTrue(clusters.empty)
        # 20 or more comparable positions still qualify.
        status, clusters = run_two_samples(["A", "B"], [SEQUENCE, SEQUENCE], ["Site 1", "Site 1"],
                                           extra=coordinates)
        self.assertEqual(len(clusters), 2)

    def test_s2_5_placeholder_words_are_missing_values_not_errors(self):
        warnings = []
        data, errors, _ = validate_metadata(pd.DataFrame({
            "Accession": ["A", "B"], "Pathogen": ["V", "V"],
            "Collection_Date": ["2024-01-01", "missing"], "Location_Name": ["Tempe", "not provided"],
            "Latitude": ["33.4", "NA"], "Longitude": ["-111.9", "NA"]}, dtype="string"), warnings)
        self.assertEqual(errors, [])
        self.assertTrue(pd.isna(data["Collection_Date"].iloc[1]))
        self.assertTrue(pd.isna(data["Location_Name"].iloc[1]))
        self.assertTrue(any("missing value" in warning for warning in warnings))

    def test_s2_6_mixed_day_first_and_month_first_dates_are_rejected(self):
        frame = {"Accession": ["A", "B"], "Pathogen": ["V", "V"]}
        _, errors, _ = validate_metadata(pd.DataFrame(
            {**frame, "Collection_Date": ["03/04/2024", "13/04/2024"]}, dtype="string"), [])
        self.assertTrue(any("day-first and month-first" in error for error in errors))
        # US month-first dates alone are read consistently and stay accepted.
        data, errors, _ = validate_metadata(pd.DataFrame(
            {**frame, "Collection_Date": ["03/04/2024", "03/15/2024"]}, dtype="string"), [])
        self.assertEqual(errors, [])
        self.assertEqual(data["Collection_Date"].iloc[0], pd.Timestamp("2024-03-04"))

    def test_s2_10_days_apart_uses_calendar_days_in_either_order(self):
        metadata = pd.DataFrame({"Sample_ID": ["A", "B"],
                                 "Collection_Date": pd.to_datetime(["2024-01-01 23:00", "2024-01-02 01:00"]),
                                 "Location_Name": [pd.NA, pd.NA]})
        for first, second in [("A", "B"), ("B", "A")]:
            pairs = pd.DataFrame({"Sample_ID_1": [first], "Sample_ID_2": [second]})
            self.assertEqual(enrich_pairs(metadata, pairs)["Days_Apart"].iloc[0], 1)

    def test_s2_9_batch_file_and_input_in_folders_with_spaces(self):
        folder = Path(tempfile.mkdtemp(prefix="space test "))
        self.addCleanup(shutil.rmtree, folder, True)
        (folder / "my mafft").mkdir()
        batch = folder / "my mafft" / "mafft.bat"
        batch.write_text('@echo off\r\ntype "%~3"\r\n')
        (folder / "in put").mkdir()
        fasta = folder / "in put" / "group.fasta"
        fasta.write_text(">A\nACGT\n")
        if os.name != "nt":
            self.skipTest("cmd.exe batch files are Windows-only")
        self.assertEqual(run_mafft(str(batch), fasta).strip(), ">A\nACGT")

    def test_known_answer_p_distance_and_haversine(self):
        # Gaps and N are not comparable: positions 1-4 compare, 1 differs.
        result = compare_sequences("ACGTN-", "ACGAAA")
        self.assertEqual((result["Differences"], result["Comparable_Positions"]), (1, 4))
        self.assertAlmostEqual(result["P_Distance"], 0.25)
        # One degree of latitude on a 6371 km sphere = 6371 * pi / 180 km.
        self.assertAlmostEqual(haversine_km(0, 0, 1, 0), 6371 * np.pi / 180, places=6)


class ImportRegressionTest(unittest.TestCase):
    def setUp(self):
        self.folder = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.folder, True)

    def write(self, name, text):
        path = self.folder / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_imp_1_lowercase_ivar_sample_column_keeps_samples_apart(self):
        path = self.write("ivar_all.tsv", "sample\tREGION\tPOS\tREF\tALT\tALT_FREQ\tTOTAL_DP\n"
                                          "WW1\tR.1\t5\tA\tG\t0.5\t100\nWW2\tR.1\t5\tA\tG\t0.1\t80\n")
        data, _ = load_processed_snvs(path)
        self.assertEqual(sorted(data["Sample_ID"]), ["WW1", "WW2"])
        info = inspect_file(path, "ivar_all.tsv")
        self.assertEqual(info["ids"], ["WW1", "WW2"])

    def test_imp_2_zero_frequency_rows_are_not_detections(self):
        path = self.write("calls.csv", "Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth\n"
                                       "W1,R,5,A,G,0,500\nW1,R,9,C,T,0.4,500\n")
        data, warnings = load_processed_snvs(path)
        self.assertEqual(list(data["Position"]), [9])
        self.assertTrue(any("0% frequency" in warning for warning in warnings))

    def test_imp_4_corrupt_gzip_is_a_readable_error_with_the_user_file_name(self):
        raw = self.folder / "raw_0"
        raw.write_bytes(gzip.compress(b"a,b\n1,2\n")[:12] + b"\x00garbage" * 4)
        with self.assertRaisesRegex(ValueError, "pooled.csv.gz is not a readable gzip file"):
            plain_copy(raw, self.folder / "plain.csv", display_name="pooled.csv.gz")

    def test_imp_5_7_vcf_lines_without_sample_fields_and_genotype_info_af(self):
        header = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
        table, details = vcf_to_snv_table(self.write("a.vcf", header + "R\t5\t.\tA\tG\t.\tPASS\tAF=0.5\n"), "a")
        self.assertTrue(table.empty)
        self.assertIn("record has no FORMAT/sample fields", details["set_aside"])
        single = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                  "R\t5\t.\tA\tG\t.\tPASS\tAF=0.5\tGT:AD\t0/1:.,.\n")
        table, _ = vcf_to_snv_table(self.write("b.vcf", single), "b")
        self.assertTrue(table.empty)  # INFO/AF of a sample-column VCF is not a read frequency

    def test_imp_6_depth_is_the_allele_count_total_when_frequency_uses_it(self):
        vcf = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
               "R\t5\t.\tA\tG\t.\tPASS\t.\tAD:DP\t0,30:25\n")
        table, _ = vcf_to_snv_table(self.write("c.vcf", vcf), "c")
        row = table.iloc[0]
        self.assertEqual((row["Alt_Frequency"], row["Alt_Depth"], row["Read_Depth"]), (1.0, 30, 30))

    def test_app_4_infinite_depth_is_ignored(self):
        vcf = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
               "R\t5\t.\tA\tG\t.\tPASS\tAF=0.2;DP=1e400\n")
        table, _ = vcf_to_snv_table(self.write("d.vcf", vcf), "d")
        self.assertEqual(table["Read_Depth"].iloc[0], "")

    def test_imp_8_mixed_combined_table_is_not_labeled_vcf(self):
        path = self.write("mixed.csv", "Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,"
                                       "Read_Depth,Source_Input_Format\nW1,R,5,A,G,0.2,50,VCF\n"
                                       "W2,R,6,A,C,0.3,40,\n")
        data, _ = load_processed_snvs(path)
        self.assertEqual(data.attrs["input_format"], COMBINED_FORMAT)

    def test_app_2_non_fasta_file_set_as_reference_is_a_problem_not_a_crash(self):
        path = self.write("one.csv", "Sample_ID,Value\nA,1\n")
        info = inspect_file(path, "one.csv")
        suggest_roles([info])
        review = review_dataset({0: path}, [info], {0: "reference"})
        self.assertIn("The reference genome must be a FASTA file.", review["problems"])

    def test_imp_3_contig_length_mismatch_blocks_the_run(self):
        reference = self.write("ref.fasta", ">R.1 reference\n" + "ACGT" * 10 + "\n")
        vcf = self.write("p.vcf", "##fileformat=VCFv4.2\n##contig=<ID=R.1,length=999>\n"
                                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nR.1\t5\t.\tA\tG\t.\tPASS\tAF=0.2;DP=50\n")
        infos = [inspect_file(reference, "ref.fasta"), inspect_file(vcf, "p.vcf")]
        suggest_roles(infos)
        review = review_dataset({0: reference, 1: vcf}, infos, {})
        self.assertFalse(review["ready"])
        self.assertTrue(any("different sequences" in problem for problem in review["problems"]))


class ScienceRegressionTest(unittest.TestCase):
    def test_sci_3_state_level_coordinates_give_no_kilometre_distance(self):
        self.assertIsNone(usable_coordinates({"Location_Name": "USA: AZ", "Latitude": 33.4,
                                              "Longitude": -111.9,
                                              "Geographic_Precision": "Provided coordinates (precision unverified)"}))
        self.assertEqual(usable_coordinates({"Location_Name": "Tempe", "Latitude": 33.4,
                                             "Longitude": -111.9, "Geographic_Precision": "Provided"}),
                         (33.4, -111.9))

    def test_sci_4_masked_list_for_another_reference_is_not_a_pass(self):
        variants = pd.DataFrame({"Variant_Key": ["R.1:5:A>G"], "Reference_ID": ["R.1"], "Position": [5],
                                 "Ref": ["A"], "Alt": ["G"], "Sample_ID": ["W1"], "Read_Depth": [50]})
        samples = pd.DataFrame({"Sample_ID": ["W1"], "Genome_Coverage_Percent": [95.0]})
        checked, _, _ = apply_pooled_quality_control(variants, samples, masked_sites={("OTHER", 5)})
        self.assertTrue(checked["QC_Masked_Site"].iloc[0].startswith("Not evaluable"))

    def test_sci_6_permutation_p_counts_tied_slopes(self):
        # Every permutation of y gives |slope| >= the observed |slope| (hand check
        # over all 6 orders), so the permutation p-value must be 1.
        per_sample = pd.DataFrame({"Mean_Depth": [0.1, 0.2, 0.3],
                                   "Distinct_SNVs_After_QC": [0.1, 0.7, 0.3],
                                   "Analysis_Included": [True, True, True]})
        fit = depth_regression(per_sample, bootstrap_rounds=10, permutations=200)
        self.assertAlmostEqual(fit["permutation_p"], 1.0)
        # Least squares by hand: slope = Sxy / Sxx = 0.02 / 0.02 = 1.0, intercept = 0.166...
        self.assertAlmostEqual(fit["slope"], 1.0)
        self.assertAlmostEqual(fit["intercept"], 1 / 6)


if __name__ == "__main__":
    unittest.main()
