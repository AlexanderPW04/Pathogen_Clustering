"""Sprint 5 Phase 2: upload-to-results tests for overdispersion calibration and sensitivity.

Real Flask routes (/review, /run, result downloads). The technical-replicate
data are the SYNTHETIC files from examples/integrated/make_replicate_example.py
(true rho 0.01 for SPEC_A, 0.002 for SPEC_B). MAFFT is mocked; GeoNames and
Claude are off.
"""

import io
import os
import re
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402

from app import STAGED, app  # noqa: E402
from genomics import overdispersion as od  # noqa: E402
from genomics.frequency_uncertainty import wilson_interval  # noqa: E402
from genomics.overdispersion_analysis import NO_CALIBRATION_TEXT  # noqa: E402
from wastewater.processed_snv import load_processed_snvs  # noqa: E402

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "integrated"
REFERENCE = (EXAMPLES / "synthetic_reference.fasta").read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", "")
REPLICATE_FILES = ["synthetic_replicate_snvs.csv", "synthetic_replicate_base_counts.csv",
                   "synthetic_replicate_reference.fasta"]


def html_text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)).replace("&#39;", "'").replace("&quot;", '"')


class Phase2WebTest(unittest.TestCase):
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

    def test_technical_replicates_are_calibrated(self):
        _, page, download_id = self.upload([self.example(name) for name in REPLICATE_FILES])
        text = html_text(page)
        self.assertIn("Sequencing Uncertainty and Variant Frequency", text)
        self.assertIn("Empirically calibrated", text)
        self.assertNotIn("could not be drawn", page)
        calibration = self.download(download_id, "sprint5_overdispersion_calibration.csv").set_index("Specimen_ID")
        self.assertEqual(calibration.loc["SPEC_A", "Status"], od.ESTIMATED)
        self.assertGreaterEqual(calibration.loc["SPEC_A", "Sites_Used"], od.MINIMUM_SITES)
        self.assertTrue(0.003 < calibration.loc["SPEC_A", "Rho_Estimate"] < 0.03)  # true value 0.01
        self.assertGreater(calibration.loc["SPEC_A", "Rho_Upper_95"], calibration.loc["SPEC_A", "Rho_Estimate"])
        self.assertIn(calibration.loc["SPEC_B", "Status"], od.USABLE_STATUSES)
        self.assertLess(calibration.loc["SPEC_B", "Rho_Estimate"], calibration.loc["SPEC_A", "Rho_Estimate"])
        self.assertEqual(calibration.loc["SPEC_C", "Status"], od.TOO_FEW_SITES)
        self.assertEqual(calibration.loc["SPEC_D", "Status"], od.NOT_REPLICATE_PAIR)
        self.assertIn("different collection dates", calibration.loc["SPEC_D", "Status_Detail"])

        sites = self.download(download_id, "sprint5_calibration_sites.csv")
        multi = sites[sites["Position"] == 597]
        self.assertTrue(len(multi) == 2 and not multi["Used"].astype(bool).any())

        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        status = evidence.set_index("Sample_ID")["Calibration_Status"]
        self.assertTrue(status.loc[["REP_A1", "REP_A2", "REP_B1", "REP_B2"]].eq("Calibrated (technical replicates)").all())
        self.assertTrue(status.loc[["REP_C1", "REP_D1"]].str.contains("replicate estimate unusable").all())
        self.assertTrue(status.loc[["WW_SINGLE"]].str.contains("sample not replicated").all())
        # The Phase 1 Wilson baseline is unchanged: recompute it from the downloaded counts.
        lower, upper = wilson_interval(evidence["Alt_Read_Count"].astype(float), evidence["Total_Depth"].astype(float))
        np.testing.assert_allclose(evidence["Wilson_Lower"], lower, rtol=1e-9)
        np.testing.assert_allclose(evidence["Wilson_Upper"], upper, rtol=1e-9)
        calibrated = evidence[evidence["Rho"].notna()]
        self.assertTrue((calibrated["BetaBinomial_Lower"] <= calibrated["Wilson_Lower"] + 1e-12).all())
        self.assertTrue((calibrated["BetaBinomial_Upper"] >= calibrated["Wilson_Upper"] - 1e-12).all())
        self.assertTrue(evidence.loc[evidence["Sample_ID"].eq("WW_SINGLE"), "Rho"].isna().all())
        self.download(download_id, "replicate_concordance.png")
        self.download(download_id, "uncertainty_models.png")
        self.assertIn("Overdispersion calibration (Sprint 5 Phase 2)", text)

    def test_no_replicates_shows_the_required_statement(self):
        _, page, download_id = self.upload([self.example("synthetic_pooled_snvs.csv"),
                                            self.example("synthetic_pooled_base_counts.csv"),
                                            self.example("synthetic_reference.fasta")])
        text = html_text(page)
        self.assertIn(NO_CALIBRATION_TEXT, text)
        self.assertIn("Not calibrated", text)
        sensitivity = self.download(download_id, "sprint5_overdispersion_sensitivity.csv")
        self.assertEqual(set(sensitivity["Rho_Source"]), {"Hypothetical; not calibrated"})
        self.assertEqual(sorted(sensitivity["Assumed_Rho"].unique()), [0.002, 0.02, 0.05])
        self.download(download_id, "sprint5_overdispersion_calibration.csv", expected=404)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertTrue(evidence["Rho"].isna().all())
        self.assertTrue(evidence["BetaBinomial_Lower"].isna().all())
        models = self.download(download_id, "sprint5_uncertainty_models.csv")
        self.assertEqual(models.loc[0, "Model"], "Binomial (Wilson)")

    def test_user_values_and_empty_setting(self):
        files = [self.example("synthetic_pooled_snvs.csv"), self.example("synthetic_pooled_base_counts.csv")]
        _, _, download_id = self.upload(files, overdispersion_rho_values="0.01")
        sensitivity = self.download(download_id, "sprint5_overdispersion_sensitivity.csv")
        self.assertEqual(sorted(sensitivity["Assumed_Rho"].unique()), [0.01])
        _, page, download_id = self.upload(files, overdispersion_rho_values="")
        self.download(download_id, "sprint5_overdispersion_sensitivity.csv", expected=404)
        self.assertIn("No hypothetical overdispersion values were entered", html_text(page))

    def test_invalid_setting_gives_a_message(self):
        files = [self.example("synthetic_pooled_snvs.csv")]
        review, page, _ = self.upload(files, overdispersion_rho_values="2")
        combined = html_text(review + (page or ""))
        self.assertIn("Hypothetical overdispersion values must be numbers", combined)

    def test_missing_depth_gives_explanations_not_errors(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_pooled_snvs.csv", dtype=str, keep_default_na=False)
        snvs = snvs.assign(Read_Depth="", Alt_Depth="")
        _, page, download_id = self.upload([("no_depth.csv", snvs.to_csv(index=False).encode())])
        text = html_text(page)
        self.assertIn("Read depth is unavailable", text)
        self.assertIn(NO_CALIBRATION_TEXT, text)
        self.assertNotIn("could not be drawn", page)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertTrue(evidence["Calibration_Status"].str.startswith("No interval").all())

    def test_replicates_without_dates_still_calibrate(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_replicate_snvs.csv", dtype=str, keep_default_na=False)
        snvs = snvs[~snvs["Specimen_ID"].eq("SPEC_D")].assign(Collection_Date="")
        files = [("replicates_no_dates.csv", snvs.to_csv(index=False).encode()),
                 self.example("synthetic_replicate_base_counts.csv")]
        _, page, download_id = self.upload(files)
        self.assertIn("Temporal trends unavailable because collection dates were not provided", page)
        calibration = self.download(download_id, "sprint5_overdispersion_calibration.csv").set_index("Specimen_ID")
        self.assertEqual(calibration.loc["SPEC_A", "Status"], od.ESTIMATED)

    def test_sprint2_and_3_combined_keeps_the_panel(self):
        _, page, download_id = self.upload(
            [self.example("synthetic_individual_sequences.fasta"), self.example("synthetic_individual_metadata.csv"),
             self.example("synthetic_reference.fasta"), self.example("synthetic_pooled_snvs.csv"),
             self.example("synthetic_pooled_base_counts.csv")], mode="combined")
        self.assertIn("Sequencing Uncertainty and Variant Frequency", page)
        self.assertIn("Individual Sequence Relationships", page)
        self.download(download_id, "cross_source_comparison.csv")
        self.download(download_id, "sprint5_uncertainty_models.csv")

    @unittest.skipUnless(os.environ.get("FONTENELE_MMC2"), "set FONTENELE_MMC2 to the published mmc2.xlsx")
    def test_published_workbook_gets_no_invented_rho(self):
        path = Path(os.environ["FONTENELE_MMC2"])
        _, page, download_id = self.upload([(path.name, path.read_bytes())])
        self.assertIn(NO_CALIBRATION_TEXT, html_text(page))
        self.assertNotIn("could not be drawn", page)
        evidence = self.download(download_id, "sprint5_vaf_observations.csv")
        self.assertEqual(len(evidence), 7973)
        self.assertTrue(evidence["Rho"].isna().all())
        self.assertEqual(int(evidence["Wilson_Lower"].notna().sum()), 7970)
        sensitivity = self.download(download_id, "sprint5_overdispersion_sensitivity.csv")
        self.assertEqual(len(sensitivity), 7970 * 3)
        self.assertTrue(sensitivity["Rho_Source"].eq("Hypothetical; not calibrated").all())
        self.download(download_id, "sprint5_overdispersion_calibration.csv", expected=404)


class SpecimenLoaderTest(unittest.TestCase):
    def write(self, folder, name, table):
        path = Path(folder) / name
        table.to_csv(path, index=False)
        return path

    def test_specimen_id_is_kept_and_checked(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_replicate_snvs.csv", dtype=str, keep_default_na=False)
        with tempfile.TemporaryDirectory() as folder:
            data, _ = load_processed_snvs(self.write(folder, "ok.csv", snvs))
            self.assertEqual(data.groupby("Sample_ID")["Specimen_ID"].first()["REP_A1"], "SPEC_A")
            bad = snvs.copy()
            bad.loc[bad.index[0], "Specimen_ID"] = "OTHER"
            with self.assertRaisesRegex(ValueError, "one Specimen_ID"):
                load_processed_snvs(self.write(folder, "bad.csv", bad))

    def test_metadata_can_supply_specimen_id(self):
        snvs = pd.read_csv(EXAMPLES / "synthetic_replicate_snvs.csv", dtype=str, keep_default_na=False)
        metadata = snvs[["Sample_ID", "Specimen_ID"]].drop_duplicates()
        with tempfile.TemporaryDirectory() as folder:
            data, _ = load_processed_snvs(self.write(folder, "snvs.csv", snvs.drop(columns=["Specimen_ID"])),
                                          metadata_path=self.write(folder, "meta.csv", metadata))
            self.assertEqual(data.groupby("Sample_ID")["Specimen_ID"].first()["REP_B2"], "SPEC_B")


if __name__ == "__main__":
    unittest.main()
