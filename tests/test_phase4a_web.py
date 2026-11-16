"""Sprint 5 Phase 4A website flow: review controls, mapping profiles, readiness, and downloads.

SYNTHETIC data only, except the optional real-workbook test (FONTENELE_MMC2).
"""

import io
import json
import os
import re
import unittest
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd

from app import STAGED, app
from tests.test_phase4a_import import DATES, REFERENCE_SEQUENCE, SAMPLES, long_table

MMC2 = os.environ.get("FONTENELE_MMC2")


def wide_bytes():
    wide = {"samplename": SAMPLES, "dates": [DATES[s] for s in SAMPLES], "location_name": ["N", "N", "S"]}
    for index, (sample, mutation, value) in enumerate([("S1", "C2T", "0.1"), ("S2", "C2T", "0.5"),
                                                       ("S1", "G7A", "0.25"), ("S3", "G7A", "0.3"),
                                                       ("S2", "G17T", "0.1")]):
        wide.setdefault(f"GENE:X{index}Y::{mutation}", ["NA"] * 3)[SAMPLES.index(sample)] = value
    buffer = io.BytesIO()
    pd.DataFrame(wide).to_excel(buffer, index=False)
    return buffer.getvalue()


def text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page))


class Phase4AWebTest(unittest.TestCase):
    def setUp(self):
        os.environ["AI_ENRICHMENT_ENABLED"] = "false"
        self.client = app.test_client()
        self.addCleanup(STAGED.clear)

    def upload(self, files):
        page = self.client.post("/review", data={"files": files}, content_type="multipart/form-data")
        page = page.get_data(as_text=True)
        return page, re.search(r'name="token" value="([^"]+)"', page).group(1)

    def wide_files(self):
        return [(io.BytesIO(wide_bytes()), "wide.xlsx"),
                (io.BytesIO((">SYN.1 SYNTHETIC\n" + REFERENCE_SEQUENCE + "\n").encode()), "ref.fasta")]

    def confirmed_form(self, token):
        return {"token": token, "role_0": "snv", "role_1": "reference", "opt_0_layout_confirmed": "on",
                "opt_0_declared_reference": "SYN.1"}

    def test_review_asks_only_what_headers_cannot_settle_then_shows_step_3(self):
        page, token = self.upload(self.wide_files())
        words = text(page)
        self.assertIn("How to read each table", words)
        self.assertIn("Confirm the wide layout", words)
        self.assertIn("no reference-genome column", words)
        self.assertIn("disabled>Run analysis", page)
        page = self.client.post("/review", data=self.confirmed_form(token)).get_data(as_text=True)
        words = text(page)
        self.assertNotIn("disabled>Run analysis", page)
        self.assertIn("Step 3 · Analyses these data support", words)
        for status in ["Available", "Unavailable", "Experimental"]:
            self.assertIn(f'">{status}</span>', page)
        self.assertIn("Mutation grouping is experimental and disabled by default", words)
        self.assertIn("No read depth was provided", words)

    def test_mapping_profile_round_trip(self):
        _, token = self.upload(self.wide_files())
        response = self.client.post("/review/profile", data=self.confirmed_form(token))
        self.assertEqual(response.status_code, 200)
        profile = json.loads(response.data)
        self.assertEqual(profile["files"][0]["options"]["declared_reference"], "SYN.1")
        for value in ["S1", "0.25", "2026-01-05"]:
            self.assertNotIn(value, response.get_data(as_text=True))
        _, fresh = self.upload(self.wide_files())  # a new upload of the same format
        page = self.client.post("/review", data={"token": fresh,
                                                 "mapping_profile": (io.BytesIO(response.data), "profile.json")},
                                content_type="multipart/form-data").get_data(as_text=True)
        self.assertIn("Profile entry 1 (snv) applied to wide.xlsx", text(page))
        self.assertNotIn("disabled>Run analysis", page)
        page = self.client.post("/review", data={"token": fresh,
                                                 "mapping_profile": (io.BytesIO(b'{"a": 1}'), "x.json")},
                                content_type="multipart/form-data").get_data(as_text=True)
        self.assertIn("not a mapping profile saved by this application", text(page))

    def test_results_offer_xlsx_and_pdf_with_the_same_values(self):
        _, token = self.upload(self.wide_files())
        with redirect_stdout(io.StringIO()):
            page = self.client.post("/run", data=self.confirmed_form(token)).get_data(as_text=True)
        self.assertIn("Results workbook (XLSX)", page)
        self.assertIn("Analysis report (PDF)", page)
        links = dict((name, link) for link, name in re.findall(r'href="(/results/[^"]+/([^"/]+))"', page))
        workbook = pd.read_excel(io.BytesIO(self.client.get(links["analysis_results.xlsx"]).data), sheet_name=None)
        csv = pd.read_csv(io.BytesIO(self.client.get(links["canonical_variant_observations.csv"]).data))
        sheet = workbook["Observations"]
        self.assertEqual(sheet.shape, csv.shape)
        pd.testing.assert_series_equal(sheet["Reported_VAF"], csv["Reported_VAF"], check_dtype=False)
        self.assertIn("Summary", workbook)
        report = self.client.get(links["analysis_report.pdf"]).data
        self.assertTrue(report.startswith(b"%PDF"))
        self.assertIn("How each uploaded file was read (import provenance)", text(page))
        self.assertIn("Analyses the data supported", text(page))

    def test_two_uploads_keep_their_own_choices(self):
        csv = long_table().to_csv(index=False)
        noted = "Notes line\n" + csv
        _, first = self.upload([(io.BytesIO(noted.encode()), "noted.csv")])
        _, second = self.upload([(io.BytesIO(csv.encode()), "plain.csv")])
        page = self.client.post("/review", data={"token": first, "role_0": "snv", "opt_0_header_row": "2"})
        self.assertIn("Column headers are read from row 2, as you chose", text(page.get_data(as_text=True)))
        page = self.client.post("/review", data={"token": second, "role_0": "snv"}).get_data(as_text=True)
        self.assertNotIn("as you chose", text(page))
        self.assertNotIn("disabled>Run analysis", page)

    def test_conflict_choice_on_the_review_page(self):
        variants = long_table().drop(columns=["Location_Name"]).to_csv(index=False)
        metadata = pd.DataFrame({"Sample_ID": SAMPLES, "Collection_Date": ["2026-01-05", "2026-02-02", "2026-01-19"]})
        page, token = self.upload([(io.BytesIO(variants.encode()), "variants.csv"),
                                   (io.BytesIO(metadata.to_csv(index=False).encode()), "meta.csv")])
        form = {"token": token, "role_0": "snv", "role_1": "pooled_metadata"}
        page = self.client.post("/review", data=form).get_data(as_text=True)
        self.assertIn("Sample metadata matching", text(page))
        self.assertIn("metadata value(s) conflict", text(page))
        self.assertIn("disabled>Run analysis", page)
        page = self.client.post("/review", data={**form, "conflict_resolution": "unresolved"}).get_data(as_text=True)
        self.assertNotIn("disabled>Run analysis", page)

    @unittest.skipUnless(MMC2 and Path(MMC2).exists(), "Set FONTENELE_MMC2 to the published mmc2.xlsx")
    def test_real_workbook_exports(self):
        page, token = self.upload([(io.BytesIO(Path(MMC2).read_bytes()), "mmc2.xlsx")])
        self.assertIn("Step 3 · Analyses these data support", text(page))
        with redirect_stdout(io.StringIO()):
            page = self.client.post("/run", data={"token": token, "role_0": "snv"}).get_data(as_text=True)
        links = dict((name, link) for link, name in re.findall(r'href="(/results/[^"]+/([^"/]+))"', page))
        workbook = pd.read_excel(io.BytesIO(self.client.get(links["analysis_results.xlsx"]).data),
                                 sheet_name=["Observations", "Samples"])
        self.assertEqual(len(workbook["Observations"]), 7973)
        self.assertEqual(len(workbook["Samples"]), 52)
        self.assertTrue(self.client.get(links["analysis_report.pdf"]).data.startswith(b"%PDF"))


if __name__ == "__main__":
    unittest.main()
