"""Tests for file formats, VCF import, column detection, role detection, and review."""

import gzip
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomics import file_formats  # noqa: E402
from genomics.column_detection import (  # noqa: E402
    apply_mapping,
    detect_columns,
    resolve_choices,
    unresolved,
    value_warnings,
)
from genomics.file_formats import (  # noqa: E402
    FASTA_KIND,
    FASTQ_KIND,
    detect_kind,
    file_stem,
    plain_copy,
    split_name,
    vcf_to_snv_table,
)
from genomics.input_inspection import (  # noqa: E402
    inspect_file,
    prepare_inputs,
    review_dataset,
    suggest_roles,
)
from validate_inputs import validate_metadata  # noqa: E402
from wastewater.processed_snv import VCF_SNV_FORMAT, load_processed_snvs  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "integrated"

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=SYNREF.1,length=120>\n"
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWW_1\tWW_2\n"
)


class TemporaryFolderTest(unittest.TestCase):
    def setUp(self):
        self.folder = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.folder, ignore_errors=True)

    def write(self, name, text):
        path = self.folder / name
        path.write_text(text, encoding="utf-8")
        return path


class FileFormatTest(TemporaryFolderTest):
    def test_compressed_names_and_stems(self):
        self.assertEqual(split_name("a.fasta.gz"), (".fasta", True))
        self.assertEqual(split_name("table.CSV"), (".csv", False))
        self.assertEqual(file_stem("S1.vcf.gz"), "S1")
        self.assertEqual(file_stem("WW 1.2.tsv"), "WW_1.2")

    def test_gzip_fasta_is_decompressed_and_recognized(self):
        source = self.folder / "seqs.fa.gz"
        with gzip.open(source, "wt") as handle:
            handle.write(">A\nACGT\n>B\nACGA\n")
        plain = plain_copy(source, self.folder / "seqs.fasta")
        self.assertEqual(detect_kind(plain, "seqs.fa.gz"), FASTA_KIND)
        info = inspect_file(plain, "seqs.fa.gz")
        self.assertEqual(info["records"], 2)
        self.assertTrue(info["compressed"])
        self.assertEqual(info["candidates"], ["fasta"])

    def test_decompression_has_a_size_limit(self):
        source = self.folder / "big.csv.gz"
        with gzip.open(source, "wb") as handle:
            handle.write(b"a,b\n" + b"1,2\n" * 5000)
        with patch.object(file_formats, "MAX_DECOMPRESSED_BYTES", 1000):
            with self.assertRaisesRegex(ValueError, "expands to more than"):
                plain_copy(source, self.folder / "big.csv")
        self.assertFalse((self.folder / "big.csv").exists())

    def test_fastq_is_refused_not_treated_as_a_table(self):
        path = self.write("reads.fastq", "@read1\nACGT\n+\nIIII\n")
        info = inspect_file(path, "reads.fastq")
        self.assertEqual(info["kind"], FASTQ_KIND)
        self.assertTrue(info["refused"])
        self.assertIn("offline workflow", info["problems"][0])

    def test_text_file_with_fasta_content_is_noted(self):
        path = self.write("sequences.txt", ">A\nACGT\n>B\nACGG\n")
        info = inspect_file(path, "sequences.txt")
        self.assertEqual(info["kind"], FASTA_KIND)
        self.assertTrue(any("content is FASTA" in note for note in info["notes"]))


class VcfTest(TemporaryFolderTest):
    def test_multi_allelic_allele_depths_are_split_per_allele(self):
        path = self.write("multi.vcf", VCF_HEADER
                          + "SYNREF.1\t15\t.\tT\tC,G\t50\tPASS\tDP=200\tAD:DP\t80,15,5:100\t90,0,10:100\n")
        table, details = vcf_to_snv_table(path, "multi")
        rows = {(row["Sample_ID"], row["Alt"]): row for row in table.to_dict("records")}
        self.assertAlmostEqual(rows[("WW_1", "C")]["Alt_Frequency"], 0.15)
        self.assertAlmostEqual(rows[("WW_1", "G")]["Alt_Frequency"], 0.05)
        self.assertEqual(rows[("WW_1", "C")]["Alt_Depth"], 15)
        self.assertEqual(rows[("WW_1", "C")]["Read_Depth"], 100)
        # WW_2 has 0 reads for C: not a variant call, and not converted to 0%.
        self.assertNotIn(("WW_2", "C"), rows)
        self.assertAlmostEqual(rows[("WW_2", "G")]["Alt_Frequency"], 0.10)
        self.assertEqual(details["fields_used"], ["FORMAT/AD"])
        self.assertEqual(sum(details["set_aside"].values()), 1)

    def test_genotype_only_and_indels_are_set_aside(self):
        path = self.write("gt.vcf", VCF_HEADER.replace("##FORMAT=<ID=AD", "##FORMAT=<ID=GT")
                          + "SYNREF.1\t20\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\n"
                          + "SYNREF.1\t30\t.\tAT\tA\t.\tPASS\t.\tGT\t1/1\t0/1\n")
        table, details = vcf_to_snv_table(path, "gt")
        self.assertTrue(table.empty)
        reasons = " ".join(details["set_aside"])
        self.assertIn("GT only", reasons)
        self.assertIn("indel", reasons)

    def test_sites_only_lofreq_style_vcf_uses_info_and_file_name(self):
        path = self.write("WW_7.vcf", "##fileformat=VCFv4.0\n"
                          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                          "SYNREF.1\t42\t.\tT\tA\t300\tPASS\tDP=500;AF=0.32;DP4=1,2,3,4\n")
        table, details = vcf_to_snv_table(path, file_stem("WW_7.vcf"))
        self.assertEqual(list(table["Sample_ID"]), ["WW_7"])
        self.assertAlmostEqual(table["Alt_Frequency"].iloc[0], 0.32)
        self.assertEqual(table["Read_Depth"].iloc[0], 500)
        self.assertIn("file name", details["sample_id_source"])

    def test_info_af_is_not_used_for_one_of_several_samples(self):
        header = VCF_HEADER.replace("AD,Number=R", "GT,Number=1")
        path = self.write("cohort.vcf", header + "SYNREF.1\t42\t.\tT\tA\t.\tPASS\tAF=0.4\tGT\t0/1\t0/0\n")
        table, details = vcf_to_snv_table(path, "cohort")
        self.assertTrue(table.empty)  # the cohort-level INFO/AF describes no single sample

    def test_vcf_loads_as_processed_snvs_with_format_and_notes(self):
        path = self.write("calls.vcf", VCF_HEADER
                          + "SYNREF.1\t15\t.\tT\tC\t50\tPASS\t.\tAD\t80,20\t70,30\n")
        data, warnings = load_processed_snvs(path)
        self.assertEqual(data.attrs["input_format"], VCF_SNV_FORMAT)
        self.assertEqual(sorted(data["Sample_ID"]), ["WW_1", "WW_2"])
        self.assertTrue(any("FORMAT/AD" in warning for warning in warnings))
        self.assertTrue(any("collection dates" in warning for warning in warnings))
        self.assertTrue(data["Collection_Date"].isna().all())  # never invented


class ColumnDetectionTest(unittest.TestCase):
    def test_aliases_ignore_case_spacing_and_abbreviations(self):
        matches = detect_columns(["sample name", "DATE", "Lat", "long", "Site", "Coverage Percent"],
                                 "pooled_metadata")
        found = {match["field"]: match["column"] for match in matches if match["status"] == "matched"}
        self.assertEqual(found["Sample_ID"], "sample name")
        self.assertEqual(found["Collection_Date"], "DATE")
        self.assertEqual(found["Latitude"], "Lat")
        self.assertEqual(found["Longitude"], "long")
        self.assertEqual(found["Location_Name"], "Site")
        self.assertEqual(found["Genome_Coverage_Percent"], "Coverage Percent")

    def test_several_candidates_are_ambiguous_until_the_user_chooses(self):
        matches = detect_columns(["Sample", "Date", "Collection_Date", "AF", "DP", "POS", "REF",
                                  "ALT", "CHROM"], "snv")
        date = next(match for match in matches if match["field"] == "Collection_Date")
        self.assertEqual(date["status"], "ambiguous")
        self.assertEqual(date["candidates"], ["Date", "Collection_Date"])
        self.assertFalse(unresolved(matches) == [])
        chosen = resolve_choices(matches, {"Collection_Date": "Collection_Date"})
        self.assertEqual(unresolved(chosen), [])
        table = pd.DataFrame(columns=["Sample", "Date", "Collection_Date", "AF", "DP", "POS",
                                      "REF", "ALT", "CHROM"])
        renamed = apply_mapping(table, chosen)
        self.assertIn("Ignored_Date", renamed.columns)
        self.assertIn("Alt_Frequency", renamed.columns)
        self.assertIn("Reference_Genome", renamed.columns)

    def test_units_are_never_rescaled_silently(self):
        table = pd.DataFrame({"Sample_ID": ["A"], "Genome_Coverage_Percent": ["0.95"]})
        problems = value_warnings(table, detect_columns(table.columns, "pooled_metadata"))
        self.assertIn("fraction", problems[0])
        table = pd.DataFrame({"Sample": ["A"], "Reference": ["R"], "Position": ["1"], "Ref": ["A"],
                              "Alt": ["G"], "Frequency": ["45"], "Depth": ["30"]})
        problems = value_warnings(table, detect_columns(table.columns, "snv"))
        self.assertIn("percentage", problems[0])


class RoleDetectionTest(TemporaryFolderTest):
    def inspect_example(self, names):
        infos = [inspect_file(EXAMPLE / name, name) for name in names]
        suggest_roles(infos)
        return {info["name"]: info for info in infos}

    def test_example_files_get_roles_with_stated_evidence(self):
        infos = self.inspect_example([
            "synthetic_individual_sequences.fasta", "synthetic_individual_metadata.csv",
            "synthetic_reference.fasta", "synthetic_pooled_snvs.csv",
            "synthetic_pooled_base_counts.csv", "synthetic_variant_signatures.csv"])
        expected = {"synthetic_individual_sequences.fasta": "fasta",
                    "synthetic_individual_metadata.csv": "metadata",
                    "synthetic_reference.fasta": "reference",
                    "synthetic_pooled_snvs.csv": "snv",
                    "synthetic_pooled_base_counts.csv": "positions",
                    "synthetic_variant_signatures.csv": "signatures"}
        for name, role in expected.items():
            self.assertEqual(infos[name]["suggested_role"], role, name)
            self.assertEqual(infos[name]["role_status"], "clear", name)
        self.assertIn("reference named in", infos["synthetic_reference.fasta"]["role_evidence"])

    def test_single_sequence_without_evidence_needs_a_choice(self):
        path = self.write("one.fasta", ">X1\n" + "ACGT" * 10 + "\n")
        info = inspect_file(path, "one.fasta")
        suggest_roles([info])
        self.assertIsNone(info["suggested_role"])
        self.assertEqual(info["role_status"], "needs choice")

    def test_metadata_role_comes_from_matching_ids(self):
        fasta = self.write("seqs.fasta", ">S1\nACGT\n>S2\nACGA\n")
        metadata = self.write("meta.csv", "Sample,Date,Site\nS1,2020-01-01,X\nS2,2020-01-02,Y\n")
        infos = [inspect_file(fasta, "seqs.fasta"), inspect_file(metadata, "meta.csv")]
        suggest_roles(infos)
        self.assertEqual(infos[1]["suggested_role"], "metadata")
        self.assertIn("2 of 2 values", infos[1]["role_evidence"])

    def test_unmatched_metadata_without_header_hints_is_not_guessed(self):
        metadata = self.write("meta.csv", "Sample,Date\nQ1,2020-01-01\n")
        info = inspect_file(metadata, "meta.csv")
        suggest_roles([info])
        self.assertIsNone(info["suggested_role"])


class ReviewTest(TemporaryFolderTest):
    def review(self, files, roles=None, choices=None, options=None):
        infos = []
        staged = {}
        for index, (name, text) in enumerate(files):
            staged[index] = self.write(name, text)
            infos.append(inspect_file(staged[index], name))
        suggest_roles(infos)
        return staged, review_dataset(staged, infos, roles or {}, choices or {}, options or {})

    def test_fasta_without_metadata_needs_a_stated_pathogen(self):
        files = [("seqs.fasta", ">S1\nACGTACGT\n>S2\nACGTACGA\n")]
        _, review = self.review(files)
        self.assertFalse(review["ready"])
        self.assertTrue(any("pathogen" in problem for problem in review["problems"]))
        staged, review = self.review(files, options={"declared_pathogen": "Test virus",
                                                     "complete_genomes": True})
        self.assertTrue(review["ready"])
        plan = {row["Analysis"]: row["Status"] for row in review["analyses"]}
        self.assertEqual(plan["Sequence comparison (MAFFT, p-distance)"], "Will run")
        self.assertEqual(plan["Temporal screening of sequence pairs"], "Skipped")
        paths, names = prepare_inputs(staged, review, self.folder,
                                      {"declared_pathogen": "Test virus"})
        metadata = pd.read_csv(paths["metadata"])
        self.assertEqual(list(metadata["Accession"]), ["S1", "S2"])
        self.assertEqual(set(metadata["Pathogen"]), {"Test virus"})
        self.assertNotIn("Collection_Date", metadata.columns)  # nothing invented

    def test_fasta_annotations_fill_absent_metadata_columns_only(self):
        files = [("seqs.fasta", ">S1 [collection_date=2020-06-01] [organism=Virus X]\nACGT\n"
                                ">S2 [collection_date=2020-06-03] [organism=Virus X]\nACGA\n"),
                 ("meta.csv", "Accession,Location\nS1,Tempe\nS2,Mesa\n")]
        staged, review = self.review(files)
        self.assertTrue(any("FASTA description annotations" in note for note in review["annotation_notes"]))
        paths, _ = prepare_inputs(staged, review, self.folder)
        metadata = pd.read_csv(paths["metadata"])
        self.assertEqual(list(metadata["Collection_Date"]), ["2020-06-01", "2020-06-03"])
        self.assertEqual(set(metadata["Pathogen_Source"]), {"FASTA description annotation"})
        self.assertEqual(list(metadata["Location_Name"]), ["Tempe", "Mesa"])

    def test_mismatched_ids_and_reference_labels_are_reported(self):
        files = [("seqs.fasta", ">S1\nACGT\n>S2\nACGA\n"),
                 ("meta.csv", "Accession,Pathogen\nS1,V\nS3,V\n"),
                 ("ref.fasta", ">REF.1 reference\nACGTACGT\n"),
                 ("calls.csv", "Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth\n"
                               "W1,OTHER.1,2,C,T,0.5,40\n")]
        _, review = self.review(files, roles={1: "metadata"})
        problems = " ".join(review["problems"])
        self.assertIn("without a metadata row", problems)
        self.assertIn("without a FASTA record", problems)
        self.assertTrue(any("never assumed" in warning for warning in review["warnings"]))
        plan = {row["Analysis"]: row["Status"] for row in review["analyses"]}
        self.assertEqual(plan["Cross-source comparison"], "Will run")  # labels checked at run time

    def test_ambiguous_columns_block_the_run_until_chosen(self):
        files = [("calls.csv", "Sample,Date,Collection_Date,Reference,Position,Ref,Alt,AF,DP\n"
                               "W1,2020-01-01,2020-01-01,R.1,5,A,G,0.2,50\n")]
        _, review = self.review(files)
        self.assertFalse(review["ready"])
        self.assertTrue(any("Choose which column is Collection_Date" in problem
                            for problem in review["problems"]))
        _, review = self.review(files, choices={0: {"Collection_Date": "Collection_Date"}})
        self.assertTrue(review["ready"])

    def test_two_files_in_a_single_role_conflict(self):
        files = [("ref1.fasta", ">R1 reference\nACGT\n"), ("ref2.fasta", ">R2 reference\nACGA\n")]
        _, review = self.review(files)
        self.assertTrue(any("only one is allowed" in problem for problem in review["problems"]))

    def test_several_vcf_files_are_combined_into_one_pooled_table(self):
        single = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                  "R.1\t5\t.\tA\tG\t.\tPASS\tDP=100;AF=0.25\n")
        staged, review = self.review([("WW_1.vcf", single), ("WW_2.vcf", single.replace("0.25", "0.5"))])
        self.assertEqual(review["pooled_ids"], ["WW_1", "WW_2"])
        paths, names = prepare_inputs(staged, review, self.folder)
        data, _ = load_processed_snvs(paths["snv"])
        self.assertEqual(sorted(data["Sample_ID"]), ["WW_1", "WW_2"])
        self.assertEqual(sorted(data["Alt_Frequency"]), [0.25, 0.5])


class OptionalMetadataTest(unittest.TestCase):
    def frame(self, **columns):
        base = {"Accession": ["A", "B"], "Pathogen": ["V", "V"]}
        base.update(columns)
        return pd.DataFrame(base, dtype="string")

    def test_missing_date_location_and_host_columns_are_warnings(self):
        warnings = []
        data, errors, _ = validate_metadata(self.frame(), warnings)
        self.assertEqual(errors, [])
        self.assertTrue(data["Collection_Date"].isna().all())
        self.assertEqual(len(warnings), 3)
        self.assertIn("No Collection_Date column", warnings[0])

    def test_blank_dates_are_warnings_but_invalid_dates_stay_errors(self):
        warnings = []
        _, errors, _ = validate_metadata(self.frame(Collection_Date=["2020-01-05", ""]), warnings)
        self.assertEqual(errors, [])
        self.assertTrue(any("Missing Collection_Date" in warning for warning in warnings))
        _, errors, _ = validate_metadata(self.frame(Collection_Date=["2020-01-05", "not a date"]), [])
        self.assertTrue(any("Invalid Collection_Date" in error for error in errors))

    def test_partial_dates_are_not_given_a_day(self):
        warnings = []
        data, errors, _ = validate_metadata(self.frame(Collection_Date=["2020", "2020-06"]), warnings)
        self.assertEqual(errors, [])
        self.assertTrue(data["Collection_Date"].isna().all())  # not 2020-01-01 / 2020-06-01
        self.assertTrue(any("only a year or month" in warning for warning in warnings))

    def test_pathogen_remains_required(self):
        data = pd.DataFrame({"Accession": ["A"]}, dtype="string")
        _, errors, _ = validate_metadata(data, [])
        self.assertTrue(any("Pathogen" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
