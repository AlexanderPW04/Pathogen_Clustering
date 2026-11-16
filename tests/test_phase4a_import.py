"""Sprint 5 Phase 4A: generalized import, metadata matching, readiness, and canonical tables.

All data here are small SYNTHETIC tables written by the tests; the real
published datasets are checked in test_phase4a_published.py (optional).
"""

import os
import shutil
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from genomics.analysis_readiness import AVAILABLE, EXPERIMENTAL, PARTIAL, UNAVAILABLE
from genomics.column_detection import detect_columns
from genomics.import_options import interpret_table, read_file_options, shift_decimal
from genomics.mapping_profiles import apply_profile, build_profile, read_profile
from genomics.metadata_matching import loose
from genomics.table_layouts import (
    BASE_COUNT_LAYOUT,
    WIDE_LAYOUT,
    derive_from_base_counts,
    detect_wide_matrix,
    parse_mutation_header,
    wide_to_long,
)
from validation.check_published_datasets import run_import

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_SEQUENCE = "ACGTACGTAC" * 12  # 120 bases, SYNTHETIC
SAMPLES = ["S1", "S2", "S3"]
# (sample, position, ref, alt, alt reads, depth); every Ref matches REFERENCE_SEQUENCE.
RECORDS = [("S1", 2, "C", "T", 10, 100), ("S1", 7, "G", "A", 30, 120), ("S2", 2, "C", "T", 25, 50),
           ("S2", 17, "G", "T", 8, 80), ("S3", 7, "G", "A", 40, 160), ("S3", 23, "G", "C", 12, 60)]
DATES = {"S1": "2026-01-05", "S2": "2026-01-12", "S3": "2026-01-19"}
SITES = {"S1": "Site North", "S2": "Site North", "S3": "Site South"}


def long_table(**renames):
    rows = [{"Sample_ID": s, "Reference_Genome": "SYN.1", "Position": p, "Ref": r, "Alt": a,
             "Alt_Frequency": repr(n / d), "Read_Depth": d, "Alt_Depth": n, "Collection_Date": DATES[s],
             "Location_Name": SITES[s]} for s, p, r, a, n, d in RECORDS]
    table = pd.DataFrame(rows)
    return table.rename(columns=renames)


class Workspace(unittest.TestCase):
    def setUp(self):
        self.folder = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.folder, ignore_errors=True)
        self.reference = self.write_text("synthetic_reference.fasta", ">SYN.1 SYNTHETIC reference\n"
                                         + REFERENCE_SEQUENCE + "\n")

    def write_text(self, name, text):
        path = self.folder / name
        path.write_text(text, encoding="utf-8")
        return path

    def write(self, name, table, sep=","):
        path = self.folder / name
        if name.endswith(".xlsx"):
            table.to_excel(path, index=False)
        else:
            table.to_csv(path, index=False, sep=sep)
        return path

    def run_files(self, files, options=None, settings=None):
        review, result, _ = run_import(files, options, settings)
        return review, result

    def observations(self, result):
        table = result["tables"]["canonical_variant_observations.csv"]
        return table.sort_values(["Sample_ID", "Genome_Position", "Alt"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Layout readers (unit level)
# ---------------------------------------------------------------------------

class LayoutTest(unittest.TestCase):
    def test_mutation_headers_are_parsed_only_when_they_state_one_snv(self):
        self.assertEqual(parse_mutation_header("C10029T"), {"Ref": "C", "Position": 10029, "Alt": "T"})
        self.assertEqual(parse_mutation_header("ORF1ab:T3255I::C10029T")["Position"], 10029)
        self.assertEqual(parse_mutation_header("10029C>T")["Alt"], "T")
        self.assertEqual(parse_mutation_header("g.10029C>T")["Ref"], "C")
        for header in ["Site", "position", "ORF1ab:L5261-::T16047C,T16047G", "del21765", "C0T", "AC10029TG"]:
            self.assertIsNone(parse_mutation_header(header), header)

    def test_wide_matrix_needs_several_mutation_columns_and_other_columns(self):
        self.assertIsNone(detect_wide_matrix(["Sample", "C10T"]))
        self.assertIsNone(detect_wide_matrix(["C10T", "G20A", "T30C"]))  # no sample column
        layout = detect_wide_matrix(["Sample", "Date", "C10T", "G20A", "T30C", "M:L93-::C26801G,C26801T"])
        self.assertEqual(len(layout["mutation_columns"]), 3)
        self.assertEqual(layout["not_interpreted"], ["M:L93-::C26801G,C26801T"])

    def test_na_cells_become_no_record_never_zero(self):
        table = pd.DataFrame({"Sample": ["A", "B"], "C10T": ["0.5", "NA"], "G20A": ["", "0.25"],
                              "T30C": ["nan", "N/A"]})
        long, counts = wide_to_long(table, detect_wide_matrix(table.columns))
        self.assertEqual(len(long), 2)
        self.assertEqual(counts["not_reported"], 4)
        self.assertEqual(sorted(long["Alt_Frequency"]), ["0.25", "0.5"])
        self.assertFalse((long["Alt_Frequency"] == "0").any())

    def test_base_count_derivation_uses_measured_counts_and_the_confirmed_minimum(self):
        table = pd.DataFrame({"sample": ["S1", "S1"], "genome_position": ["2", "7"], "ref_base": ["C", "G"],
                              "A": ["1", "30"], "C": ["90", "0"], "G": ["0", "88"], "T": ["9", "2"],
                              "position_coverage": ["100", "120"], "var_freq": ["0.9", "0.25"]})
        columns = {"Sample_ID": "sample", "Position": "genome_position", "Ref": "ref_base",
                   "A_Count": "A", "C_Count": "C", "G_Count": "G", "T_Count": "T"}
        long, summary = derive_from_base_counts(table, columns, 2)
        self.assertEqual(list(zip(long["Position"], long["Alt"], long["Alt_Depth"])),
                         [("2", "T", 9), ("7", "A", 30), ("7", "T", 2)])
        self.assertEqual(list(long["Read_Depth"]), [100, 120, 120])
        self.assertEqual(list(long["Ref_Depth"]), [90, 88, 88])
        self.assertIn("Calculated by the importer", long["Source_Frequency_Provenance"].iloc[0])
        # Columns that could be mistaken for the derived values are kept under Source_ names.
        self.assertIn("Source_position_coverage", long.columns)
        self.assertEqual(summary["depth_checks"]["position_coverage"], "equals A+C+G+T on every derived record")
        self.assertIn("var_freq", long.columns)  # tool-specific columns are carried, never reinterpreted
        with self.assertRaises(ValueError):
            derive_from_base_counts(table.assign(A=["x", "1"]), columns, 1)


class OptionTest(unittest.TestCase):
    def test_percent_conversion_is_exact_and_keeps_the_original(self):
        self.assertEqual(shift_decimal("12.5", 2), "0.125")
        self.assertEqual(shift_decimal("100", 2), "1")
        self.assertEqual(shift_decimal("4.85437", 2), "0.0485437")  # significant digits kept
        self.assertIsNone(shift_decimal("abc", 2))
        table = pd.DataFrame({"Alt_Frequency": ["12.5", "4.85437"], "Position": ["9", "0"]})
        converted, rows = interpret_table(table, {"frequency_scale": "percent", "position_base": "0"}, True)
        self.assertEqual(list(converted["Alt_Frequency"]), ["0.125", "0.0485437"])
        self.assertEqual(list(converted["Source_Alt_Frequency_Original"]), ["12.5", "4.85437"])
        self.assertEqual(list(converted["Position"]), ["10", "1"])
        self.assertEqual(list(converted["Source_Position_Original"]), ["9", "0"])
        self.assertEqual(converted["Read_Depth"].tolist(), ["", ""])  # never filled in
        self.assertEqual(len(rows), 3)

    def test_form_options_are_validated(self):
        form = {"opt_0_frequency_scale": "percent", "opt_0_position_base": "7", "opt_0_base_count_minimum": "0",
                "opt_0_header_row": "3", "opt_0_declared_reference": " NC_045512.2 ", "opt_0_layout_confirmed": "on"}
        options = read_file_options(form, 0)
        self.assertEqual(options, {"frequency_scale": "percent", "header_row": 3, "declared_reference": "NC_045512.2",
                                   "layout_confirmed": True})

    def test_site_is_not_a_position_and_sequencing_is_not_a_platform(self):
        matches = {m["field"]: m for m in detect_columns(["Sample", "Site", "Ref", "Alt", "AF", "sequencing"], "snv")}
        self.assertEqual(matches["Position"]["status"], "missing")
        meta = {m["field"]: m for m in detect_columns(["Sample", "sequencing", "platform"], "pooled_metadata")}
        self.assertEqual(meta["Sequencing_Platform"]["column"], "platform")

    def test_loose_keys_are_hints_only(self):
        self.assertEqual(loose("WW-001"), loose("ww_1"))
        self.assertNotEqual(loose("WW-001"), loose("WW-0011"))


# ---------------------------------------------------------------------------
# End to end through the website's import functions
# ---------------------------------------------------------------------------

class EquivalentLayoutsTest(Workspace):
    """Scientifically equivalent files with different layouts give the same downstream values."""

    VALUE_COLUMNS = ["Sample_ID", "Genome_Position", "Ref", "Alt", "Alt_Read_Count", "Total_Depth",
                     "Calculated_VAF", "VAF_Used", "Count_Provenance"]

    def baseline(self):
        review, result = self.run_files([{"path": self.write("standard.csv", long_table()), "role": "snv"},
                                         {"path": self.reference, "role": "reference"}])
        self.assertTrue(review["ready"], review["problems"])
        return result

    def test_renamed_columns_tsv_and_excel_give_identical_vaf(self):
        expected = self.observations(self.baseline())
        renamed = long_table(Sample_ID="SampleName", Position="POS", Ref="REF", Alt="ALT", Alt_Frequency="AF",
                             Read_Depth="DP", Alt_Depth="ALT_DP", Reference_Genome="CHROM",
                             Collection_Date="Sampling_Date")
        for name in ["renamed.tsv", "renamed.xlsx"]:
            path = self.write(name, renamed, sep="\t")
            review, result = self.run_files([{"path": path, "role": "snv"}, {"path": self.reference, "role": "reference"}])
            self.assertTrue(review["ready"], review["problems"])
            got = self.observations(result)
            pd.testing.assert_frame_equal(got[self.VALUE_COLUMNS], expected[self.VALUE_COLUMNS], check_dtype=False)

    def test_percent_frequencies_need_a_choice_then_match_fractions(self):
        expected = self.observations(self.baseline())
        table = long_table()
        table["Alt_Frequency"] = [shift_decimal(text, -2) for text in table["Alt_Frequency"]]
        path = self.write("percent.csv", table)
        review, _ = self.run_files([{"path": path, "role": "snv"}])
        self.assertFalse(review["ready"])
        self.assertTrue(any("may be a percentage" in text for text in review["problems"]))
        review, result = self.run_files([{"path": path, "role": "snv", "options": {"frequency_scale": "percent"}},
                                         {"path": self.reference, "role": "reference"}])
        self.assertTrue(review["ready"], review["problems"])
        got = self.observations(result)
        np.testing.assert_allclose(got["Reported_VAF"], expected["Reported_VAF"], rtol=0, atol=1e-15)
        pd.testing.assert_frame_equal(got[self.VALUE_COLUMNS], expected[self.VALUE_COLUMNS], check_dtype=False)
        self.assertIn("Source_Alt_Frequency_Original", got.columns)
        bad = table.assign(Alt_Frequency=["150"] + list(table["Alt_Frequency"][1:]))
        review, _ = self.run_files([{"path": self.write("bad.csv", bad), "role": "snv",
                                     "options": {"frequency_scale": "percent"}}])
        self.assertTrue(any("below 0 or above 100" in text for text in review["problems"]))

    def test_wide_matrix_equals_the_long_table_for_reported_frequencies(self):
        expected = self.observations(self.baseline())
        wide = {"samplename": SAMPLES, "dates": [DATES[s] for s in SAMPLES],
                "location_name": [SITES[s] for s in SAMPLES]}
        for s, p, r, a, n, d in RECORDS:
            wide.setdefault(f"GENE:X1Y::{r}{p}{a}", ["NA"] * 3)[SAMPLES.index(s)] = repr(n / d)
        path = self.write("wide.xlsx", pd.DataFrame(wide))
        files = [{"path": path, "role": "snv", "options": {"declared_reference": "SYN.1"}},
                 {"path": self.reference, "role": "reference"}]
        review, _ = self.run_files(files)
        self.assertFalse(review["ready"])  # the wide reading must be confirmed first
        files[0]["options"]["layout_confirmed"] = True
        review, result = self.run_files(files)
        self.assertTrue(review["ready"], review["problems"])
        self.assertEqual(review["files"][0]["info"]["layout"], WIDE_LAYOUT)
        got = self.observations(result)
        self.assertEqual(len(got), len(RECORDS))
        np.testing.assert_allclose(got["Reported_VAF"], expected["Reported_VAF"])
        self.assertTrue(got["Total_Depth"].isna().all())  # a wide matrix carries no depth; none invented
        self.assertTrue((got["Measurement_Class"] == "Frequency only (no depth)").all())
        samples = result["tables"]["canonical_samples.csv"].set_index("Sample_ID")
        self.assertEqual(str(samples.loc["S3", "Collection_Date"])[:10], "2026-01-19")
        readiness = {row["Analysis"]: row["Status"] for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Wilson uncertainty intervals"], UNAVAILABLE)
        self.assertEqual(readiness["Reported VAF summaries"], AVAILABLE)
        self.assertIn("Read_Depth", {row["Item"] for row in result["tables"]["import_provenance.csv"].to_dict("records")})

    def test_base_count_table_gives_measured_counts(self):
        rows = []
        for s, p, r, a, n, d in RECORDS:
            counts = {base: 0 for base in "ACGT"}
            counts[r], counts[a] = d - n, n
            rows.append({"sample": s, "genome_position": p, "ref_base": r, **counts})
        path = self.write("counts.csv", pd.DataFrame(rows))
        files = [{"path": path, "role": "snv", "options": {"declared_reference": "SYN.1"}},
                 {"path": self.reference, "role": "reference"}]
        review, _ = self.run_files(files)
        self.assertFalse(review["ready"])
        self.assertTrue(any("minimum number of reads" in text for text in review["problems"]))
        files[0]["options"]["base_count_minimum"] = 1
        review, result = self.run_files(files)
        self.assertTrue(review["ready"], review["problems"])
        self.assertEqual(review["files"][0]["info"]["layout"], BASE_COUNT_LAYOUT)
        got = self.observations(result)
        expected = self.observations(self.baseline())
        pd.testing.assert_frame_equal(got[["Sample_ID", "Genome_Position", "Alt", "Alt_Read_Count", "Total_Depth"]],
                                      expected[["Sample_ID", "Genome_Position", "Alt", "Alt_Read_Count",
                                                "Total_Depth"]], check_dtype=False)
        self.assertTrue(got["Reported_VAF"].isna().all())  # calculated, never shown as reported
        readiness = {row["Analysis"]: row["Status"] for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Reported VAF summaries"], PARTIAL)
        self.assertEqual(readiness["Count-based VAF calculations"], AVAILABLE)
        self.assertEqual(readiness["Full-count Yue-Clayton dissimilarity"], UNAVAILABLE)

    def test_zero_based_positions_warn_then_convert(self):
        table = long_table()
        table["Position"] = table["Position"] - 1
        path = self.write("zero.csv", table)
        review, _ = self.run_files([{"path": path, "role": "snv"}, {"path": self.reference, "role": "reference"}])
        self.assertTrue(any("match with 0-based positions" in text for text in review["warnings"]))
        review, result = self.run_files([{"path": path, "role": "snv", "options": {"position_base": "0"}},
                                         {"path": self.reference, "role": "reference"}])
        got = self.observations(result)
        self.assertEqual(sorted(got["Genome_Position"]), sorted(p for _, p, *_ in RECORDS))
        self.assertTrue(result["tables"]["reference_verification.csv"]["Status"].eq("Verified compatible").all())

    def test_header_row_choice_reads_a_table_below_notes(self):
        path = self.folder / "notes_above.csv"
        path.write_text("Table S9. Variant calls\nGenerated by a pipeline\n" + long_table().to_csv(index=False),
                        encoding="utf-8")
        review, _ = self.run_files([{"path": path, "role": "snv"}])
        self.assertFalse(review["ready"])
        review, result = self.run_files([{"path": path, "role": "snv", "options": {"header_row": 3}},
                                         {"path": self.reference, "role": "reference"}])
        self.assertTrue(review["ready"], review["problems"])
        self.assertEqual(len(self.observations(result)), len(RECORDS))


class ReferenceAndCoverageTest(Workspace):
    def test_missing_reference_column_needs_a_statement(self):
        path = self.write("noref.csv", long_table().drop(columns=["Reference_Genome"]))
        review, _ = self.run_files([{"path": path, "role": "snv"}])
        self.assertTrue(any("no reference-genome column" in text for text in review["problems"]))
        review, result = self.run_files([{"path": path, "role": "snv", "options": {"declared_reference": "OTHER.9"}},
                                         {"path": self.reference, "role": "reference"}])
        status = result["tables"]["reference_verification.csv"].set_index("Reference_ID")["Status"]
        self.assertNotEqual(status["OTHER.9"], "Verified compatible")  # a different label is never assumed equal

    def test_segmented_genome_keeps_segments_apart(self):
        table = long_table()
        table.loc[table["Position"].eq(23), "Reference_Genome"] = "SYN_SEG2.1"
        review, result = self.run_files([{"path": self.write("segments.csv", table), "role": "snv"}])
        checks = {row["Check"]: row for row in review["matching"]}
        self.assertEqual(checks["Reference sequences (genome segments) in the variant data"]["Count"], 2)
        samples = result["tables"]["canonical_samples.csv"].set_index("Sample_ID")
        self.assertIn("2 reference sequences (segments)", samples.loc["S3", "Genome_Segments"])

    def test_sparse_table_has_no_coverage_and_never_zero_fills(self):
        review, result = self.run_files([{"path": self.write("sparse.csv", long_table()), "role": "snv"}])
        self.assertEqual(len(result["tables"]["canonical_coverage.csv"]), 0)
        readiness = {row["Analysis"]: row["Status"] for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Coverage-aware non-detection"], UNAVAILABLE)
        matrix = result["tables"]["pooled_observation_matrix.csv"]
        unreported = matrix[matrix["Status"].ne("Variant observed")]
        self.assertTrue(unreported["Alt_Frequency"].isna().all())

    def test_full_base_count_table_fills_coverage_records(self):
        rows = [{"Sample_ID": s, "Reference_ID": "SYN.1", "Position": i + 1,
                 **{f"{b}_Count": (50 if REFERENCE_SEQUENCE[i] == b else 0) for b in "ACGT"}}
                for s in SAMPLES for i in range(len(REFERENCE_SEQUENCE))]
        review, result = self.run_files([{"path": self.write("sparse.csv", long_table()), "role": "snv"},
                                         {"path": self.write("counts.csv", pd.DataFrame(rows)), "role": "positions"},
                                         {"path": self.reference, "role": "reference"}])
        coverage = result["tables"]["canonical_coverage.csv"]
        self.assertEqual(len(coverage), 3 * len(REFERENCE_SEQUENCE))
        self.assertTrue(coverage["Coverage_Status"].eq("Callable").all())
        readiness = {row["Analysis"]: row["Status"] for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Full-count Yue-Clayton dissimilarity"], AVAILABLE)
        self.assertEqual(readiness["Coverage-aware non-detection"], AVAILABLE)


class MetadataMatchingTest(Workspace):
    def split_files(self, metadata):
        variants = self.write("variants.csv", long_table().drop(columns=["Collection_Date", "Location_Name"]))
        return [{"path": variants, "role": "snv"}, {"path": self.write("meta.csv", metadata), "role": "pooled_metadata"},
                {"path": self.reference, "role": "reference"}]

    def metadata(self):
        return pd.DataFrame({"Sample_ID": SAMPLES, "Collection_Date": [DATES[s] for s in SAMPLES],
                             "Location_Name": [SITES[s] for s in SAMPLES]})

    def test_separate_metadata_gives_the_same_results_as_a_combined_table(self):
        _, combined = self.run_files([{"path": self.write("combined.csv", long_table()), "role": "snv"},
                                      {"path": self.reference, "role": "reference"}])
        review, separate = self.run_files(self.split_files(self.metadata()))
        self.assertTrue(review["ready"], review["problems"])
        for name in ["sprint5_vaf_observations.csv", "temporal_trends.csv", "pooled_samples_with_qc.csv"]:
            left, right = combined["tables"][name], separate["tables"][name]
            shared = [column for column in left.columns if column not in {"Data_Source"}]
            pd.testing.assert_frame_equal(left[shared].reset_index(drop=True), right[shared].reset_index(drop=True),
                                          check_dtype=False)

    def test_unmatched_samples_block_until_accepted_and_stay_unknown(self):
        metadata = self.metadata().iloc[:2]
        review, _ = self.run_files(self.split_files(metadata))
        self.assertFalse(review["ready"])
        self.assertTrue(any("missing from the pooled metadata: S3" in text for text in review["problems"]))
        review, result = self.run_files(self.split_files(metadata), options={"allow_unmatched_metadata": True})
        self.assertTrue(review["ready"], review["problems"])
        samples = result["tables"]["canonical_samples.csv"].set_index("Sample_ID")
        self.assertTrue(pd.isna(samples.loc["S3", "Collection_Date"]))

    def test_near_matches_are_listed_but_never_joined(self):
        metadata = self.metadata().assign(Sample_ID=["s1", "S2", "S3"])
        review, _ = self.run_files(self.split_files(metadata))
        self.assertEqual(review["matching_hints"], {"S1": "s1"})
        self.assertFalse(review["ready"])

    def test_conflicts_need_a_choice_and_each_choice_does_what_it_says(self):
        variants = long_table().drop(columns=["Location_Name"])
        metadata = self.metadata()
        metadata.loc[metadata["Sample_ID"].eq("S2"), "Collection_Date"] = "2026-02-02"
        files = [{"path": self.write("v.csv", variants), "role": "snv"},
                 {"path": self.write("m.csv", metadata), "role": "pooled_metadata"}]
        review, _ = self.run_files(files)
        self.assertFalse(review["ready"])
        self.assertEqual(len(review["metadata_conflicts"]), 1)
        expected = {"metadata": "2026-02-02", "variant_table": "2026-01-12", "unresolved": None}
        for choice, date in expected.items():
            review, result = self.run_files(files, options={"conflict_resolution": choice})
            self.assertTrue(review["ready"], (choice, review["problems"]))
            value = result["tables"]["canonical_samples.csv"].set_index("Sample_ID").loc["S2", "Collection_Date"]
            self.assertEqual(None if pd.isna(value) or value == "" else str(value)[:10], date, choice)

    def test_duplicate_metadata_ids(self):
        identical = pd.concat([self.metadata(), self.metadata().iloc[[0]]], ignore_index=True)
        review, _ = self.run_files(self.split_files(identical))
        self.assertTrue(review["ready"], review["problems"])
        differing = identical.copy()
        differing.loc[3, "Location_Name"] = "Elsewhere"
        review, _ = self.run_files(self.split_files(differing), options={"conflict_resolution": "metadata"})
        self.assertFalse(review["ready"])  # only "unresolved" can settle duplicate rows
        review, result = self.run_files(self.split_files(differing), options={"conflict_resolution": "unresolved"})
        self.assertTrue(review["ready"], review["problems"])
        self.assertEqual(result["tables"]["canonical_samples.csv"].set_index("Sample_ID").loc["S1", "Collection_Site"], "")

    def test_rows_without_ids_are_reported(self):
        metadata = pd.concat([self.metadata(), pd.DataFrame({"Sample_ID": [""], "Collection_Date": ["average:"]})])
        review, _ = self.run_files(self.split_files(metadata))
        self.assertIn("Metadata rows without a Sample_ID", {row["Check"] for row in review["matching"]})
        self.assertTrue(review["ready"], review["problems"])


class ReplicateReadinessTest(Workspace):
    def test_documented_and_undocumented_replicates(self):
        table = long_table().assign(Specimen_ID=["SP1", "SP1", "SP1", "SP2", "SP1", "SP3"])
        review, _ = self.run_files([{"path": self.write("unverified.csv", table), "role": "snv"}])
        readiness = {row["Analysis"]: row for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Empirical beta-binomial calibration"]["Status"], UNAVAILABLE)
        self.assertIn("Technical replicates were not verified", readiness["Empirical beta-binomial calibration"]["Reason"])
        example = ROOT / "examples" / "integrated"
        review, _ = self.run_files([{"path": example / "synthetic_replicate_snvs.csv", "role": "snv"}])
        readiness = {row["Analysis"]: row["Status"] for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Empirical beta-binomial calibration"], AVAILABLE)
        self.assertEqual(readiness["Experimental mutation-frequency grouping"], EXPERIMENTAL)

    def test_paper_format_counts_are_reconstructable_not_measured(self):
        review, _ = self.run_files([{"path": self.write("frequency_depth.csv", long_table().drop(columns=["Alt_Depth"])),
                                     "role": "snv"}])
        readiness = {row["Analysis"]: row for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Count-based VAF calculations"]["Status"], PARTIAL)
        self.assertIn("alternative read counts were not directly provided",
                      readiness["Count-based VAF calculations"]["Reason"])
        no_dates = self.write("no_dates.csv", long_table().drop(columns=["Collection_Date"]))
        review, _ = self.run_files([{"path": no_dates, "role": "snv"}])
        readiness = {row["Analysis"]: row["Status"] for row in review["analysis_readiness"]}
        self.assertEqual(readiness["Temporal mutation-frequency analysis"], UNAVAILABLE)


class MappingProfileTest(Workspace):
    def test_profile_holds_no_records_and_applies_only_to_matching_schemas(self):
        renamed = long_table(Sample_ID="SampleName", Alt_Frequency="AF", Read_Depth="DP")
        path = self.write("first.csv", renamed)
        review, _ = self.run_files([{"path": path, "role": "snv", "options": {"depth_definition": "total"}}])
        profile = build_profile(review, {"date_order": "auto"})
        text = str(profile)
        for value in ["S1", "Site North", "2026-01-05", "0.1"]:
            self.assertNotIn(value, text)  # mappings only, no sample values
        profile = read_profile(__import__("json").dumps(profile).encode())
        from genomics.input_inspection import inspect_file
        same = inspect_file(self.write("second.csv", renamed.iloc[:2]), "second.csv")
        other = inspect_file(self.write("third.csv", long_table()), "third.csv")
        roles, choices, options, _, report = apply_profile(profile, [other, same])
        self.assertEqual(roles, {1: "snv"})
        self.assertEqual(choices[1]["Alt_Frequency"], "AF")
        self.assertEqual(options[1]["depth_definition"], "total")
        roles, *_ , report = apply_profile(profile, [other])
        self.assertEqual(roles, {})
        self.assertIn("missing column(s)", report[0])
        with self.assertRaises(ValueError):
            read_profile(b'{"files": []}')


class LargeDatasetTest(Workspace):
    def test_large_wide_matrix_is_read_quickly(self):
        rng = np.random.default_rng(4)
        samples, mutations = 120, 1500
        data = {"samplename": [f"W{i}" for i in range(samples)]}
        values = np.where(rng.random((samples, mutations)) < 0.05, rng.random((samples, mutations)).round(4), np.nan)
        for j in range(mutations):
            data[f"A{j + 1}G"] = ["NA" if np.isnan(v) else str(v) for v in values[:, j]]
        table = pd.DataFrame(data)
        start = time.perf_counter()
        long, counts = wide_to_long(table, detect_wide_matrix(table.columns))
        self.assertLess(time.perf_counter() - start, 10)
        self.assertEqual(len(long), int(np.isfinite(values).sum()))
        self.assertEqual(counts["cells"], samples * mutations)

    def test_too_large_a_dataset_is_refused_with_an_explanation(self):
        os.environ["MAX_OBSERVATION_CELLS"] = "5"
        self.addCleanup(os.environ.pop, "MAX_OBSERVATION_CELLS", None)
        review, _ = self.run_files([{"path": self.write("std.csv", long_table()), "role": "snv"}])
        self.assertFalse(review["ready"])
        self.assertTrue(any("exceeds the available memory" in text and "No samples or variants are removed" in text
                            for text in review["problems"]))


if __name__ == "__main__":
    unittest.main()
