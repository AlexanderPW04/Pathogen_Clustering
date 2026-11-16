"""Sprint 5 Phase 3: upload-to-results tests for exploratory mutation-frequency groups.

Real Flask routes (/review, /run, result downloads). Data are the SYNTHETIC
files in examples/integrated/ (the mixture demo has planted groups and
documented technical-replicate libraries; see make_mixture_example.py).
MAFFT is mocked; GeoNames and Claude are off.
"""

import io
import os
import re
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import STAGED, app  # noqa: E402
from genomics import mutation_grouping as mg  # noqa: E402
from genomics.mutation_group_analysis import (  # noqa: E402
    DETECTION_ONLY_REASON,
    EXPLANATION,
    NOT_SHOWN_REASON,
)

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "integrated"
REFERENCE = (EXAMPLES / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
MIXTURE = ["synthetic_mixture_snvs.csv", "synthetic_mixture_base_counts.csv", "synthetic_mixture_reference.fasta"]
PHASE3_TABLES = ["sprint5_mutation_profiles.csv", "sprint5_mutation_eligibility.csv",
                 "sprint5_grouping_parameters.csv"]
GROUP_TABLES = ["sprint5_mutation_groups.csv", "sprint5_mutation_group_summary.csv",
                "sprint5_mutation_group_trajectories.csv", "sprint5_grouping_comparison.csv"]


def html_text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)).replace("&#39;", "'").replace("&quot;", '"')


class Phase3WebTest(unittest.TestCase):
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
        if token is None:
            return page, None, None
        with redirect_stdout(io.StringIO()):
            response = self.client.post("/run", data={"token": token.group(1), "intended_mode": mode, **form})
        self.assertEqual(response.status_code, 200)
        result = response.get_data(as_text=True)
        download = re.search(r"/results/([^/]+)/", result)
        return page, result, download.group(1) if download else None

    def download(self, download_id, name, expected=200):
        response = self.client.get(f"/results/{download_id}/{name}")
        self.assertEqual(response.status_code, expected, name)
        if expected != 200:
            return None
        return (pd.read_csv(io.BytesIO(response.data), keep_default_na=False, na_values=[""])
                if name.endswith(".csv") else response.data)

    @staticmethod
    def example(name):
        return name, (EXAMPLES / name).read_bytes()

    def mixture(self, **form):
        return self.upload([self.example(name) for name in MIXTURE], **form)

    # ------------------------------------------------------------------ default
    def test_default_mode_shows_profiles_but_no_groups(self):
        _, page, download_id = self.mixture()
        text = html_text(page)
        self.assertIn("Exploratory Mutation-Frequency Groups", text)
        self.assertIn(EXPLANATION, text)
        self.assertIn("Profiles only: no groups shown", text)
        self.assertIn(NOT_SHOWN_REASON, text)
        self.assertIn("Mutation-frequency profiles (not grouped)", text)
        self.assertNotIn("could not be drawn", page)
        for name in PHASE3_TABLES:
            self.download(download_id, name)
        for name in GROUP_TABLES:
            self.download(download_id, name, expected=404)
        # Phase 2 calibration is unchanged by Phase 3: all 12 documented specimens calibrated.
        calibration = self.download(download_id, "sprint5_overdispersion_calibration.csv")
        self.assertEqual(int(calibration["Status"].str.startswith("Estimated").sum()), 12)

    def test_profiles_preserve_coverage_and_provenance(self):
        _, _, download_id = self.mixture()
        profiles = self.download(download_id, "sprint5_mutation_profiles.csv")
        eligibility = self.download(download_id, "sprint5_mutation_eligibility.csv")
        self.assertEqual(profiles["Variant_Key"].nunique(), 29)
        self.assertTrue(profiles["Unit_Type"].str.startswith("Specimen: 2 verified").all())
        self.assertEqual(profiles["Unit_ID"].nunique(), 12)  # 24 libraries -> 12 specimens
        self.assertTrue(profiles["Count_Basis"].eq("Measured").all())
        # A combined unit's counts are the sum of its two libraries.
        snvs = pd.read_csv(EXAMPLES / "synthetic_mixture_snvs.csv")
        row = profiles[profiles["State"].eq("Detected")].iloc[0]
        libraries = row["Sample_IDs"].split("; ")
        position = int(row["Position"])
        base = pd.read_csv(EXAMPLES / "synthetic_mixture_base_counts.csv")
        cells = base[base["Sample_ID"].isin(libraries) & base["Position"].eq(position)]
        depth = int(cells[["A_Count", "C_Count", "G_Count", "T_Count"]].sum(axis=1).sum())
        self.assertEqual(int(row["Total_Depth"]), depth)
        self.assertIn(row["Sample_IDs"].split("; ")[0], set(snvs["Sample_ID"]) | set(base["Sample_ID"]))
        # Every sample-level state is counted; nothing is unknown here because base counts cover everything.
        states = eligibility[["Detected", "Supported_Non_Detection", "Counted_Not_Called", "Unknown_Coverage",
                              "Low_Or_No_Coverage", "Covered_Value_Unknown", "Rejected_By_QC"]].sum(axis=1)
        self.assertTrue((states == eligibility["Samples_On_Reference"]).all())

    # -------------------------------------------------------------- experimental
    def test_experimental_groups_on_planted_mixture(self):
        _, page, download_id = self.mixture(grouping_mode="experimental")
        text = html_text(page)
        self.assertIn("Experimental (not validated)", text)
        self.assertIn("Empirically calibrated", text)
        self.assertIn("not lineages", text)
        for name in GROUP_TABLES:
            self.download(download_id, name)
        groups = self.download(download_id, "sprint5_mutation_groups.csv")
        truth = pd.read_csv(EXAMPLES / "synthetic_mixture_truth.csv").set_index("Variant_Key")["Planted_Group"]
        depth = groups[groups["Method"].eq(mg.DEPTH_AWARE_METHOD)].fillna({"Group_ID": ""})
        grouped = depth[depth["Group_ID"].ne("")]
        self.assertGreaterEqual(depth["Group_ID"].replace("", np.nan).nunique(), 2)
        for _, members in grouped.groupby("Group_ID"):
            self.assertEqual(truth[members["Variant_Key"]].nunique(), 1)  # every demo group is pure
            self.assertTrue(truth[members["Variant_Key"]].iloc[0].startswith("Planted"))
        self.assertTrue(depth["Model"].eq(mg.CALIBRATED_MODE).all())
        comparison = self.download(download_id, "sprint5_grouping_comparison.csv")
        self.assertTrue(comparison["Comparison"].str.startswith("Hypothetical").any())
        self.assertTrue(comparison.loc[comparison["Comparison"].str.startswith("Hypothetical"), "Note"]
                        .str.startswith("Hypothetical; not calibrated").all())
        for figure in ["mutation_group_heatmap.png", "mutation_group_trajectories.png", "grouping_comparison.png"]:
            self.download(download_id, figure)

    def test_no_group_structure_gives_no_depth_aware_groups(self):
        files = [self.example(name) for name in ["synthetic_replicate_snvs.csv", "synthetic_replicate_base_counts.csv",
                                                 "synthetic_replicate_reference.fasta"]]
        _, page, download_id = self.upload(files, grouping_mode="experimental")
        groups = self.download(download_id, "sprint5_mutation_groups.csv")
        depth = groups[groups["Method"].eq(mg.DEPTH_AWARE_METHOD)]
        self.assertTrue(depth["Group_ID"].isna().all())
        self.assertIn("valid result, not an error", html_text(page))

    def test_detection_only_input_abstains_with_a_reason(self):
        _, page, download_id = self.upload([self.example("synthetic_mixture_snvs.csv")],
                                           grouping_mode="experimental")
        self.assertIn(DETECTION_ONLY_REASON, html_text(page))
        for name in GROUP_TABLES:
            self.download(download_id, name, expected=404)
        self.download(download_id, "sprint5_mutation_profiles.csv")
        self.download(download_id, "sprint5_vaf_observations.csv")  # Phase 1-2 still run

    def test_frequency_only_input_explains_instead_of_failing(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_mixture_snvs.csv").assign(Read_Depth="", Alt_Depth="")
        files = [("frequency_only.csv", snvs.to_csv(index=False).encode()),
                 self.example("synthetic_mixture_base_counts.csv")]
        _, page, download_id = self.upload(files, grouping_mode="experimental")
        self.assertNotIn("could not be drawn", page)
        self.assertIn("Exploratory Mutation-Frequency Groups", page)
        eligibility = self.download(download_id, "sprint5_mutation_eligibility.csv")
        self.assertIn("Frequency_Only_Units", eligibility)

    def test_mutations_on_different_references_are_never_grouped_together(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_mixture_snvs.csv")
        counts = pd.read_csv(EXAMPLES / "synthetic_mixture_base_counts.csv")
        moved = snvs["Position"] > 500
        snvs.loc[moved, "Reference_Genome"] = "SYNMIX.2"
        counts.loc[counts["Position"] > 500, "Reference_ID"] = "SYNMIX.2"
        files = [("two_references.csv", snvs.to_csv(index=False).encode()),
                 ("two_reference_counts.csv", counts.to_csv(index=False).encode())]
        _, page, download_id = self.upload(files, grouping_mode="experimental")
        groups = self.download(download_id, "sprint5_mutation_groups.csv").fillna({"Group_ID": ""})
        for _, members in groups[groups["Group_ID"].ne("")].groupby(["Method", "Group_ID"]):
            self.assertEqual(members["Reference_ID"].nunique(), 1)

    def test_single_time_point_runs(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_mixture_snvs.csv")
        first = snvs[snvs["Collection_Date"].eq(snvs["Collection_Date"].min())]
        counts = pd.read_csv(EXAMPLES / "synthetic_mixture_base_counts.csv")
        counts = counts[counts["Sample_ID"].isin(first["Sample_ID"])]
        files = [("one_date.csv", first.to_csv(index=False).encode()),
                 ("one_date_counts.csv", counts.to_csv(index=False).encode())]
        _, page, download_id = self.upload(files, grouping_mode="experimental")
        self.assertNotIn("could not be drawn", page)
        self.assertIn("Exploratory Mutation-Frequency Groups", page)
        self.download(download_id, "sprint5_mutation_eligibility.csv")

    def test_off_mode_and_invalid_setting(self):
        _, page, download_id = self.mixture(grouping_mode="off")
        self.assertIn("switched off", html_text(page))
        self.download(download_id, "sprint5_mutation_profiles.csv", expected=404)
        self.download(download_id, "sprint5_vaf_observations.csv")
        _, page, _ = self.mixture(grouping_minimum_shared="1")
        self.assertIn("minimum number of shared samples", html_text(page))

    def test_sprint2_and_3_combined(self):
        _, page, download_id = self.upload(
            [self.example("synthetic_individual_sequences.fasta"), self.example("synthetic_individual_metadata.csv"),
             self.example("synthetic_reference.fasta"), self.example("synthetic_pooled_snvs.csv"),
             self.example("synthetic_pooled_base_counts.csv")], mode="combined", grouping_mode="experimental")
        self.assertIn("Exploratory Mutation-Frequency Groups", page)
        self.assertIn("Individual sequence relationships", page)
        self.download(download_id, "cross_source_comparison.csv")
        self.download(download_id, "sprint5_mutation_profiles.csv")
        overlap = self.client.get(f"/results/{download_id}/sprint5_group_individual_overlap.csv")
        if overlap.status_code == 200:
            table = pd.read_csv(io.BytesIO(overlap.data))
            self.assertTrue(table["Interpretation"].str.startswith("Descriptive only").all())

    @unittest.skipUnless(os.environ.get("FONTENELE_MMC2"), "set FONTENELE_MMC2 to the published mmc2.xlsx")
    def test_published_workbook_gets_no_groups_or_zeros(self):
        path = Path(os.environ["FONTENELE_MMC2"])
        _, page, download_id = self.upload([(path.name, path.read_bytes())], grouping_mode="experimental")
        text = html_text(page)
        self.assertIn(DETECTION_ONLY_REASON, text)
        self.assertNotIn("could not be drawn", page)
        profiles = self.download(download_id, "sprint5_mutation_profiles.csv")
        self.assertEqual(len(profiles), 7958)
        self.assertTrue(profiles["State"].eq("Detected").all())  # no invented zeros
        self.assertTrue(profiles["Count_Basis"].eq("Reconstructed (approximate)").all())
        self.assertTrue(profiles["Rho"].isna().all())  # no invented calibration
        eligibility = self.download(download_id, "sprint5_mutation_eligibility.csv")
        self.assertEqual(len(eligibility), 5669)
        for name in GROUP_TABLES:
            self.download(download_id, name, expected=404)
        self.download(download_id, "sprint5_vaf_observations.csv")


if __name__ == "__main__":
    unittest.main()
