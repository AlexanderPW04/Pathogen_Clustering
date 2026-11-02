"""Flask tests for the guided upload -> review -> run workflow and the new figures.

MAFFT is mocked with tests/fake_mafft.py and GeoNames is disabled, so these
are not evidence that the external tools work; see test_sequence_variants.py
for the real-MAFFT test that runs only when MAFFT is installed.
"""

import gzip
import io
import os
import re
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import app  # noqa: E402
from genomics import figures  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "integrated"
REFERENCE = (EXAMPLE / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
ALL_FILES = ["synthetic_individual_sequences.fasta", "synthetic_individual_metadata.csv",
             "synthetic_reference.fasta", "synthetic_pooled_snvs.csv",
             "synthetic_pooled_base_counts.csv", "synthetic_variant_signatures.csv"]
PNG = b"\x89PNG"


def example(name, content=None, upload_name=None):
    data = (EXAMPLE / name).read_bytes() if content is None else content
    return (io.BytesIO(data), upload_name or name)


class GuidedWorkflowTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": ""}),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(REFERENCE))]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.client = app.test_client()

    def review(self, files, **form):
        response = self.client.post("/review", data={"files": files, **form},
                                    content_type="multipart/form-data")
        page = response.get_data(as_text=True)
        token = re.search(r'name="token" value="([^"]+)"', page)
        return page, token.group(1) if token else None

    def run_reviewed(self, token, **form):
        response = self.client.post("/run", data={"token": token, **form})
        return response.get_data(as_text=True)

    def test_review_detects_every_role_and_run_builds_the_integrated_dashboard(self):
        page, token = self.review([example(name) for name in ALL_FILES])
        self.assertIn("Step 2 · Review your dataset", page)
        self.assertEqual(page.count("badge completed\">Detected"), 6)
        for text in ["Reference compatibility", "Analyses available with these files",
                     "SYNREF_B.1 versus SYNREF.1", "Will run"]:
            self.assertIn(text, page)
        self.assertNotIn("disabled>Run analysis", page)
        result = self.run_reviewed(token)
        self.assertIn("<h2>2. Cross-source evidence</h2>", result)
        self.assertIn("Variants observed in both sources:</strong> 2", result)
        self.assertIn("How the two sources were connected", result)
        self.assertIn("Variant-overlap matrix across both sources", result)
        self.assertIn("Potential clusters:</strong>", result)  # Sprint 2 kept
        self.assertIn("Analyses predicted in Step 2", result)

    def test_compressed_uploads_are_recognized(self):
        fasta = gzip.compress((EXAMPLE / "synthetic_individual_sequences.fasta").read_bytes())
        table = gzip.compress((EXAMPLE / "synthetic_pooled_snvs.csv").read_bytes())
        page, _ = self.review([example("x", fasta, "individual.fasta.gz"),
                               example("y", table, "pooled.csv.gz")])
        self.assertIn("gzip-compressed", page)
        self.assertIn("FASTA (5 records)", page)
        self.assertEqual(page.count("badge completed\">Detected"), 2)

    def test_raw_reads_are_refused_and_nothing_can_run(self):
        page, token = self.review([example("r", b"@r1\nACGT\n+\nIIII\n", "reads.fastq")])
        self.assertIn("Not accepted", page)
        self.assertIn("offline workflow", page)
        self.assertIn("disabled>Run analysis", page)
        result = self.run_reviewed(token)
        self.assertIn("Resolve the problems listed in Step 2", result)

    def test_expired_or_unknown_token_asks_for_a_new_upload(self):
        self.assertIn("expired", self.run_reviewed("not-a-token"))

    def test_sequences_without_metadata_still_get_genetic_comparison(self):
        page, token = self.review([example("synthetic_individual_sequences.fasta")])
        self.assertIn("Pathogen of all individual sequences", page)
        form = {"declared_pathogen": "Synthetic respiratory virus",
                "genomic_region_mode": "complete_genomes", "confirm_complete_genomes": "on"}
        page, token = self.review([], token=token, **form)
        self.assertNotIn("disabled>Run analysis", page)
        result = self.run_reviewed(token, **form)
        self.assertIn("Pairwise comparisons:</strong> 10", result)
        self.assertIn("No Collection_Date column", result)
        self.assertIn("Built from FASTA record IDs", result)
        self.assertRegex(result, r"Temporal screening of sequence pairs</td><td><span class=\"badge "
                                 r"skipped\">Skipped")
        self.assertIn("Potential clusters:</strong> 0", result)

    def test_pooled_vcf_alone_runs_pooled_analyses_only(self):
        vcf = (b"##fileformat=VCFv4.2\n##contig=<ID=SYNREF.1,length=120>\n"
               b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWW_1\tWW_2\tWW_3\n"
               b"SYNREF.1\t15\t.\tT\tC\t.\tPASS\t.\tAD:DP\t70,30:100\t60,40:100\t50,50:100\n"
               b"SYNREF.1\t42\t.\tT\tA,G\t.\tPASS\t.\tAD:DP\t80,10,10:100\t90,10,0:100\t75,25,0:100\n")
        page, token = self.review([example("v", vcf, "pooled.vcf")])
        self.assertIn("VCF files carry no collection dates", page)
        result = self.run_reviewed(token)
        self.assertIn("Pooled samples:</strong> 3", result)
        self.assertIn("VCF (single-nucleotide entries converted", result)
        self.assertIn("Cross-source comparison needs individual sequences", result)
        self.assertIn("Upload individual FASTA sequences and metadata", result)

    def test_different_reference_label_is_never_forced_together(self):
        snvs = (EXAMPLE / "synthetic_pooled_snvs.csv").read_text(encoding="utf-8")
        snvs = snvs.replace("SYNREF.1", "OTHER.1").replace("SYNREF_B.1", "OTHER.1")
        files = [example(name) for name in ALL_FILES[:3]] + [example("s", snvs.encode(), "pooled.csv")]
        page, token = self.review(files)
        self.assertIn("OTHER.1 versus SYNREF.1", page)
        result = self.run_reviewed(token)
        self.assertIn("Variants observed in both sources:</strong> 0", result)
        self.assertIn("Incomparable: reference not verified as the same", result)

    def test_ambiguous_columns_must_be_chosen_before_running(self):
        table = ("Sample,Date,Collection_Date,Reference,Position,Ref,Alt,AF,DP\n"
                 "W1,2026-03-01,2026-03-01,SYNREF.1,15,T,C,0.2,50\n"
                 "W2,2026-03-08,2026-03-08,SYNREF.1,15,T,C,0.4,60\n")
        page, token = self.review([example("t", table.encode(), "calls.csv")])
        self.assertIn('name="map_0_Collection_Date"', page)
        self.assertIn("Choose which column is Collection_Date", page)
        self.assertIn("Resolve the problems", self.run_reviewed(token))
        result = self.run_reviewed(token, map_0_Collection_Date="Collection_Date")
        self.assertIn("Pooled samples:</strong> 2", result)


class NewFigureTest(unittest.TestCase):
    def test_completeness_marks_unavailable_fields_instead_of_zero(self):
        png = figures.completeness_figure([
            {"Source": "Individual sequences", "Samples": 3, "With_Collection_Date": 0,
             "With_Location": 2, "With_Coordinates": 1, "With_Genome_Coverage": "needs a reference FASTA"},
            {"Source": "Pooled samples", "Samples": 4, "With_Collection_Date": 4, "With_Location": 4,
             "With_Coordinates": "needs pooled metadata", "With_Genome_Coverage": "not supplied"}])
        self.assertTrue(png.startswith(PNG))

    def test_geography_without_coordinates_draws_only_the_precision_panel(self):
        samples = pd.DataFrame({
            "Sample_Type": ["Individual sequence", "Pooled sample"],
            "Precision_Category": ["Broad administrative area", "Named place or site, no coordinates"],
            "Latitude": [33.4, None], "Longitude": [-112.0, None]})
        self.assertTrue(figures.geography_figure(samples).startswith(PNG))

    def test_relationship_figures(self):
        pairs = pd.DataFrame({"Sample_ID_1": ["A", "A", "B"], "Sample_ID_2": ["B", "C", "C"],
                              "P_Distance": [0.01, 0.2, 0.03], "Days_Apart": [3, None, 30],
                              "Km_Apart": [5.0, None, None], "Same_Location": [None, None, "False"]})
        self.assertTrue(figures.distance_figure(pairs, [0.05, 14, 50]).startswith(PNG))
        network = figures.relationship_network_figure(["A", "B", "C"], [("A", "B")], [("B", "C")])
        self.assertTrue(network.startswith(PNG))

    def test_depth_and_qc_figures_count_unverifiable_checks_separately(self):
        variants = pd.DataFrame({"Read_Depth": [10, 50, None], "Analysis_Included": [False, True, True],
                                 "QC_Variant_Depth": ["Fail (depth 10 < 20)", "Pass",
                                                      "Not evaluable (depth not supplied)"]})
        self.assertTrue(figures.depth_distribution_figure(variants, 20).startswith(PNG))
        self.assertTrue(figures.qc_outcome_figure(variants).startswith(PNG))


if __name__ == "__main__":
    unittest.main()
