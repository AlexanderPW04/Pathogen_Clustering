"""Sprint 5 Phase 4: upload-to-results tests for mutation annotation.

Real Flask routes (/review, /run, result downloads). Annotation files are the
official NCBI RefSeq records in examples/annotation/ (real). Variant tables are
SYNTHETIC: small tables written here at real SARS-CoV-2 or influenza positions
(reference bases taken from the real sequences), plus the synthetic integrated
and mixture demos. MAFFT is mocked; GeoNames and Claude are off.
"""

import io
import os
import re
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import STAGED, app  # noqa: E402
from genomics.integrated_analysis import DEFAULT_SETTINGS  # noqa: E402
from genomics.mutation_annotation import MISSING_ANNOTATION, SECTION_TITLE  # noqa: E402
from genomics.references import load_reference_fasta  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ANNOTATION = ROOT / "examples" / "annotation"
EXAMPLES = ROOT / "examples" / "integrated"
REFERENCE = (EXAMPLES / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
SARS2 = load_reference_fasta(ANNOTATION / "NC_045512.2.fasta")["Sequence"]
ANNOTATION_TABLES = ["sprint5_annotated_mutations.csv", "sprint5_mutation_annotations.csv",
                     "sprint5_annotated_observations.csv", "sprint5_annotation_sources.csv"]
GROUP_TABLES = ["sprint5_mutation_groups.csv", "sprint5_mutation_group_summary.csv"]
# (position, alt) on NC_045512.2; Ref comes from the real reference sequence.
SITES = [(241, "T"), (3037, "T"), (14408, "T"), (23403, "G"), (28881, "A"), (28882, "A"), (27972, "T")]


def html_text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)).replace("&#39;", "'").replace("&quot;", '"')


def sars2_table(sites=SITES, wrong_ref_at=None):
    """A SYNTHETIC count-based SNV table: three dated samples at one site."""

    rows = []
    for number, (date, depth) in enumerate([("2026-01-05", 200), ("2026-01-12", 400), ("2026-01-19", 300)]):
        for index, (position, alt) in enumerate(sites):
            ref = SARS2[position - 1]
            if position == wrong_ref_at:
                ref = "A" if ref != "A" else "C"
            alt_depth = depth // (2 + index % 3) - number * 5
            rows.append({"Sample_ID": f"WW{number + 1}", "Reference_Genome": "NC_045512.2", "Position": position,
                         "Ref": ref, "Alt": alt, "Alt_Frequency": round(alt_depth / depth, 6), "Read_Depth": depth,
                         "Alt_Depth": alt_depth, "Collection_Date": date, "Location_Name": "Synthetic Site"})
    return pd.DataFrame(rows).to_csv(index=False).encode()


class Phase4WebTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": "", "AI_ENRICHMENT_ENABLED": "false"}),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(REFERENCE))]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def upload(self, files, mode="pooled", **form):
        uploads = [(io.BytesIO(content), name) for name, content in files]
        review = self.client.post("/review", data={"files": uploads, "intended_mode": mode, **form},
                                  content_type="multipart/form-data")
        self.assertEqual(review.status_code, 200)
        page = review.get_data(as_text=True)
        token = re.search(r'name="token" value="([^"]+)"', page)
        self.assertIsNotNone(token, html_text(page)[:2000])
        with redirect_stdout(io.StringIO()):
            response = self.client.post("/run", data={"token": token.group(1), "intended_mode": mode, **form})
        self.assertEqual(response.status_code, 200)
        result = response.get_data(as_text=True)
        download = re.search(r"/results/([^/]+)/", result)
        self.assertIsNotNone(download, html_text(result)[:3000])
        return page, result, download.group(1)

    def download(self, download_id, name, expected=200):
        response = self.client.get(f"/results/{download_id}/{name}")
        self.assertEqual(response.status_code, expected, name)
        if expected != 200:
            return None
        return (pd.read_csv(io.BytesIO(response.data), keep_default_na=False, na_values=[""])
                if name.endswith(".csv") else response.data)

    @staticmethod
    def file(folder, name):
        return name, (folder / name).read_bytes()

    def test_pooled_counts_with_genbank(self):
        review, page, download_id = self.upload([("snvs.csv", sars2_table()),
                                                 self.file(ANNOTATION, "NC_045512.2.gb")])
        self.assertIn("Genome annotation (GenBank or GFF3) for coding consequences", review)
        text = html_text(page)
        self.assertIn(SECTION_TITLE, text)
        self.assertIn("Predicted from the reference annotation", text)
        self.assertIn("Distinct mutations evaluated for annotation: 7", text)
        self.assertIn("Successfully annotated: 7", text)
        self.assertIn("does not show an effect on transmissibility", text)
        self.assertNotIn("could not be drawn", page)
        for name in ANNOTATION_TABLES:
            self.download(download_id, name)
        self.download(download_id, "mutation_annotation_overview.png")
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv").set_index("Variant_Key")
        spike = mutations.loc["NC_045512.2:23403:A>G"]
        self.assertEqual((spike["Amino_Acid_Changes"], spike["Consequence_Category"]), ("S:D614G", "Missense"))
        self.assertEqual(spike["Pooled_Samples_Observed"], 3)
        self.assertEqual(spike["Count_Basis"], "Supplied (3)")
        self.assertEqual(mutations.loc["NC_045512.2:241:C>T", "Consequence_Category"], "Noncoding")
        self.assertEqual(mutations.loc["NC_045512.2:27972:C>T", "Amino_Acid_Changes"], "ORF8:Q27*")
        self.assertEqual(mutations.loc["NC_045512.2:28881:G>A", "Shared_Codon_With"], "NC_045512.2:28882:G>A")
        self.assertEqual(mutations.loc["NC_045512.2:23403:A>G", "Series_With_Two_Or_More_Dates"], 1)
        # VAF, depth, and intervals are copied from the Phase 1-2 evidence, never recalculated.
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        observed = self.download(download_id, "sprint5_annotated_observations.csv")
        self.assertEqual(len(observed), len(evidence))
        for column in ["Sample_ID", "Variant_Key", "VAF", "Total_Depth", "Wilson_Lower", "Wilson_Upper",
                       "Alt_Count_Source", "Reference_Status", "QC_Status"]:
            pd.testing.assert_series_equal(observed[column], evidence[column], check_names=False)
        spike_vafs = evidence.loc[evidence["Variant_Key"].eq("NC_045512.2:23403:A>G"), "VAF"]
        self.assertAlmostEqual(spike["Median_VAF"], spike_vafs.median())
        self.assertTrue(observed.loc[observed["Variant_Key"].eq("NC_045512.2:23403:A>G"),
                                     "Predicted_Amino_Acid_Change"].eq("S:D614G").all())
        status = self.download(download_id, "analysis_status.csv").set_index("Analysis")
        self.assertEqual(status.loc["Mutation annotation and protein consequences (Sprint 5 Phase 4)", "Status"],
                         "Completed")
        provenance = self.download(download_id, "provenance.csv")
        self.assertIn("GenBank flat file: NC_045512.2.gb",
                      set(provenance.loc[provenance["Category"].eq("Annotation (Phase 4)"), "Value"]))
        # Phase 3 stays off by default: profiles only, no group tables.
        for name in GROUP_TABLES:
            self.download(download_id, name, expected=404)
        self.assertNotIn("Exploratory_Group", mutations.columns)

    def test_missing_annotation_explains_and_everything_else_runs(self):
        _, page, download_id = self.upload([("snvs.csv", sars2_table())])
        text = html_text(page)
        self.assertIn(SECTION_TITLE, text)
        self.assertIn(MISSING_ANNOTATION, text)
        for name in ANNOTATION_TABLES:
            self.download(download_id, name, expected=404)
        for name in ["sprint5_vaf_observations.csv", "sprint5_temporal_summary.csv", "sprint5_mutation_profiles.csv",
                     "pooled_variants_with_qc.csv"]:
            self.download(download_id, name)
        status = self.download(download_id, "analysis_status.csv").set_index("Analysis")
        self.assertEqual(status.loc["Mutation annotation and protein consequences (Sprint 5 Phase 4)", "Status"],
                         "Skipped")

    def test_gff3_uses_only_a_matching_reference_fasta(self):
        # GFF3 states no genetic code, so the standard code is confirmed as a user would (final-polish 1C).
        _, page, download_id = self.upload([("snvs.csv", sars2_table()), self.file(ANNOTATION, "NC_045512.2.gff3"),
                                            self.file(ANNOTATION, "NC_045512.2.fasta")], genetic_code="1")
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv").set_index("Variant_Key")
        self.assertEqual(mutations.loc["NC_045512.2:14408:C>T", "Amino_Acid_Changes"], "ORF1ab:P4715L")
        self.assertTrue(mutations["Annotation_Status"].eq("Annotated").all())
        self.assertIn("matching reference FASTA is used", html_text(page))
        # Without any sequence, codons cannot be checked: nothing is annotated, the run still completes.
        _, page, download_id = self.upload([("snvs.csv", sars2_table()), self.file(ANNOTATION, "NC_045512.2.gff3")],
                                           genetic_code="1")
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv")
        self.assertTrue(mutations["Annotation_Status"].eq("Not annotated").all())
        self.assertTrue(mutations["Reference_Check"].eq("Not checked (no reference sequence)").all())
        self.download(download_id, "sprint5_vaf_observations.csv")

    def test_annotation_of_another_pathogen_is_never_applied(self):
        review, page, download_id = self.upload([("snvs.csv", sars2_table()),
                                                 self.file(ANNOTATION, "influenza_A_PR8_segments.gb")])
        self.assertIn("never transferred between references", html_text(review))
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv")
        self.assertTrue(mutations["Annotation_Status"].eq("Not annotated").all())
        self.assertTrue(mutations["Genes"].isna().all())
        self.assertIn("Unverified or unannotatable: 7", html_text(page))
        self.download(download_id, "mutation_annotation_overview.png", expected=404)

    def test_wrong_reference_base_is_not_annotated(self):
        _, page, download_id = self.upload([("snvs.csv", sars2_table(wrong_ref_at=23403)),
                                            self.file(ANNOTATION, "NC_045512.2.gb")])
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv").set_index("Position")
        self.assertTrue(mutations.loc[23403, "Reference_Check"].startswith("Mismatch: reference has A"))
        self.assertEqual(mutations.loc[23403, "Annotation_Status"], "Not annotated")
        self.assertEqual(mutations.loc[14408, "Annotation_Status"], "Annotated")

    def test_segmented_influenza_reference(self):
        from genomics.genome_annotation import load_genbank
        records = load_genbank(ANNOTATION / "influenza_A_PR8_segments.gb")["records"]
        rows = []
        for reference_id, position in [("NC_002016.1", 760), ("NC_002016.1", 300), ("NC_002017.1", 500)]:
            ref = records[reference_id]["Sequence"][position - 1]
            for sample, date in [("FLU1", "2026-02-01"), ("FLU2", "2026-02-08")]:
                rows.append({"Sample_ID": sample, "Reference_Genome": reference_id, "Position": position, "Ref": ref,
                             "Alt": "C" if ref != "C" else "T", "Alt_Frequency": 0.25, "Read_Depth": 400,
                             "Alt_Depth": 100, "Collection_Date": date, "Location_Name": "Synthetic Site"})
        table = pd.DataFrame(rows).to_csv(index=False).encode()
        _, page, download_id = self.upload([("flu.csv", table), self.file(ANNOTATION, "influenza_A_PR8_segments.gb")])
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv")
        self.assertTrue(mutations["Annotation_Status"].eq("Annotated").all())
        self.assertEqual(set(mutations["Segment"].astype(str)), {"7", "4"})
        by_position = mutations.set_index("Position")
        self.assertEqual(by_position.loc[760, "Consequence_Category"] in {"Missense", "Synonymous",
                                                                          "Differs between overlapping coding regions"},
                         True)
        self.assertIn("M2", by_position.loc[760, "Genes"])
        self.assertEqual(by_position.loc[500, "Genes"], "HA")
        self.assertIn("segment 7", html_text(page))

    def test_sprint2_and_3_combined_with_a_test_annotation(self):
        # SYNTHETIC test-only GFF3 for the invented SYNREF.1: one CDS at 73..99 (an open reading frame
        # in the invented sequence), used only to exercise the combined path.
        self.assertEqual(REFERENCE[72:99], "ATGGAGGCGGGCAATGCACTCGGGTAG")
        gff = ("##gff-version 3\n##sequence-region SYNREF.1 1 120\n"
               "SYNREF.1\ttest\tCDS\t73\t99\t.\t+\t0\tID=cds-SYN1;gene=synthetic_gene_B;protein_id=SYN1\n")
        _, page, download_id = self.upload(
            [self.file(EXAMPLES, "synthetic_individual_sequences.fasta"),
             self.file(EXAMPLES, "synthetic_individual_metadata.csv"), self.file(EXAMPLES, "synthetic_reference.fasta"),
             self.file(EXAMPLES, "synthetic_pooled_snvs.csv"), self.file(EXAMPLES, "synthetic_pooled_base_counts.csv"),
             ("synthetic_test_annotation.gff3", gff.encode())], mode="combined", genetic_code="1")
        self.assertIn("Individual Sequence Relationships", page)
        self.download(download_id, "cross_source_comparison.csv")
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv").set_index("Variant_Key")
        start = mutations.loc["SYNREF.1:75:G>A"]
        self.assertEqual((start["Consequence_Category"], start["Amino_Acid_Changes"]),
                         ("Start lost", "synthetic_gene_B:M1I"))
        self.assertEqual(mutations.loc["SYNREF.1:15:T>C", "Consequence_Category"], "Noncoding")
        self.assertIn("Individual_Sequences_Carrying", mutations.columns)
        self.assertTrue(mutations["Evidence_Sources"].str.contains("Individual").any())
        self.assertTrue(mutations["Evidence_Sources"].str.contains("Pooled").any())

    def test_individual_sequences_only_with_annotation(self):
        gff = ("##gff-version 3\n##sequence-region SYNREF.1 1 120\n"
               "SYNREF.1\ttest\tCDS\t73\t99\t.\t+\t0\tID=cds-SYN1;gene=synthetic_gene_B;protein_id=SYN1\n")
        _, page, download_id = self.upload(
            [self.file(EXAMPLES, "synthetic_individual_sequences.fasta"),
             self.file(EXAMPLES, "synthetic_individual_metadata.csv"), self.file(EXAMPLES, "synthetic_reference.fasta"),
             ("synthetic_test_annotation.gff3", gff.encode())], mode="individual", genetic_code="1")
        text = html_text(page)
        self.assertIn(SECTION_TITLE, text)
        self.assertIn("Predicted from the reference annotation", text)
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv")
        self.assertTrue(len(mutations))
        self.assertTrue(mutations["Evidence_Sources"].eq("Individual").all())
        self.assertNotIn("Pooled_Samples_Observed", mutations.columns)  # no pooled evidence is invented
        self.download(download_id, "sprint5_annotated_observations.csv", expected=404)

    def test_annotation_shown_with_experimental_groups_without_validating_them(self):
        gff = "##gff-version 3\n##sequence-region SYNMIX.1 1 {}\n"
        reference = (EXAMPLES / "synthetic_mixture_reference.fasta").read_text(encoding="utf-8")
        length = len("".join(reference.split("\n", 1)[1].split()))
        files = [self.file(EXAMPLES, name) for name in ["synthetic_mixture_snvs.csv",
                                                         "synthetic_mixture_base_counts.csv",
                                                         "synthetic_mixture_reference.fasta"]]
        files.append(("mixture_test.gff3", gff.format(length).encode()))
        _, page, download_id = self.upload(files, grouping_mode="experimental")
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv")
        grouped = mutations["Exploratory_Group"].dropna()
        self.assertTrue(len(grouped))
        self.assertTrue(grouped.str.endswith("(experimental, not validated)").all())
        self.assertIn("Experimental (not validated)", html_text(page))
        self.assertIn("does not validate the exploratory mutation-frequency groups", html_text(page))
        _, _, default_id = self.upload(files)
        self.assertNotIn("Exploratory_Group", self.download(default_id, "sprint5_annotated_mutations.csv").columns)
        for name in GROUP_TABLES:
            self.download(default_id, name, expected=404)

    def test_grouping_default_is_unchanged(self):
        self.assertEqual(DEFAULT_SETTINGS["grouping_mode"], "profiles")

    @unittest.skipUnless(os.environ.get("FONTENELE_MMC2"), "set FONTENELE_MMC2 to the published mmc2.xlsx")
    def test_published_workbook_with_ncbi_annotation(self):
        path = Path(os.environ["FONTENELE_MMC2"])
        _, page, download_id = self.upload([(path.name, path.read_bytes()), self.file(ANNOTATION, "NC_045512.2.gb")])
        text = html_text(page)
        self.assertIn("Distinct mutations evaluated for annotation: 5677", text)
        self.assertIn("Successfully annotated: 5677", text)
        comparison = self.download(download_id, "sprint5_annotation_comparison.csv")
        self.assertEqual(comparison["Agreement"].value_counts().to_dict(),
                         {"Agrees": 5613, "Agrees (no codon in either)": 64})
        self.assertEqual(int(comparison["Note"].fillna("").str.contains("labels this noncoding SNV 'Syn'").sum()), 64)
        mutations = self.download(download_id, "sprint5_annotated_mutations.csv")
        self.assertEqual(mutations["Consequence_Category"].value_counts().to_dict(),
                         {"Missense": 4109, "Synonymous": 1243, "Stop gained": 243, "Noncoding": 64,
                          "Start lost": 12, "Stop lost": 3, "Stop retained": 2,
                          "Differs between overlapping coding regions": 1})
        self.assertEqual(len(self.download(download_id, "sprint5_annotated_observations.csv")), 7973)
        self.assertTrue(mutations["Count_Basis"].dropna().str.startswith("Reconstructed").all())
        for name in GROUP_TABLES:
            self.download(download_id, name, expected=404)


if __name__ == "__main__":
    unittest.main()
