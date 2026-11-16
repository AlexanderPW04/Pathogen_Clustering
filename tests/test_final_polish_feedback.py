"""Final polish, Milestone 2: feedback during long steps and on-demand report downloads.

The browser behaviour (banner, elapsed time, ignored repeat clicks, report status
messages) was checked in a real browser; these tests pin the server side and the
page markup it depends on. SYNTHETIC data only.
"""

import os
import threading
import time
import unittest
import zlib
from time import monotonic
from unittest.mock import patch

import app as app_module
from app import RESULT_DOWNLOADS, RESULT_LOCK, app
from genomics.report_exports import SUMMARY_NAME


def store(files):
    """Put a fake stored analysis into the download cache; returns its id."""

    download_id = "test-" + str(time.perf_counter_ns())
    packed = {name: (zlib.compress(content), mimetype) for name, (content, mimetype) in files.items()}
    with RESULT_LOCK:
        RESULT_DOWNLOADS[download_id] = (monotonic(), packed)
    return download_id


class ReportDownloadTest(unittest.TestCase):
    def setUp(self):
        os.environ["AI_ENRICHMENT_ENABLED"] = "false"
        self.client = app.test_client()

    def tearDown(self):
        with RESULT_LOCK:
            for key in [key for key in RESULT_DOWNLOADS if str(key).startswith("test-")]:
                del RESULT_DOWNLOADS[key]

    def test_repeated_clicks_build_the_workbook_once(self):
        download_id = store({SUMMARY_NAME: (b"{}", "application/json")})
        calls = []

        def slow_build(name, read):
            calls.append(name)
            time.sleep(0.6)
            return b"workbook-bytes", []

        responses = {}

        def click(index):
            responses[index] = app.test_client().get(f"/results/{download_id}/analysis_results.xlsx")

        with patch.object(app_module, "build_from_stored", side_effect=slow_build):
            threads = [threading.Thread(target=click, args=(index,)) for index in range(3)]
            for thread in threads:
                thread.start()
                time.sleep(0.05)
            for thread in threads:
                thread.join()
        self.assertEqual(calls, ["analysis_results.xlsx"])  # built once, then reused
        self.assertEqual({response.status_code for response in responses.values()}, {200})
        self.assertEqual({response.data for response in responses.values()}, {b"workbook-bytes"})
        self.assertEqual(app_module.REPORT_BUILDS, {})  # no lock is kept after the builds

    def test_a_failed_build_is_reported_and_can_be_retried(self):
        download_id = store({SUMMARY_NAME: (b"{}", "application/json")})
        with patch.object(app_module, "build_from_stored", side_effect=RuntimeError("synthetic failure")):
            failed = self.client.get(f"/results/{download_id}/analysis_report.pdf")
        self.assertEqual(failed.status_code, 500)
        self.assertIn(b"CSV downloads hold the same values", failed.data)
        with patch.object(app_module, "build_from_stored", return_value=(b"%PDF-synthetic", [])):
            retried = self.client.get(f"/results/{download_id}/analysis_report.pdf")
        self.assertEqual((retried.status_code, retried.data), (200, b"%PDF-synthetic"))

    def test_expired_results_explain_themselves(self):
        for name in ["canonical_samples.csv", "analysis_results.xlsx"]:
            response = self.client.get(f"/results/no-such-results/{name}")
            self.assertEqual(response.status_code, 404)
            self.assertIn(b"kept in memory for 10 minutes", response.data)

    def test_page_markup_for_feedback(self):
        page = self.client.get("/").get_data(as_text=True)
        self.assertIn('id="working"', page)
        self.assertIn("Processing genomic data. Large sequencing datasets may take several minutes.", page)
        self.assertNotIn("up to two minutes", page)
        self.assertIn("No progress percentage is shown", page)
        template = (app_module.Path(app_module.__file__).parent / "templates" / "_review_import.html").read_text(
            encoding="utf-8")
        self.assertIn("data-no-busy", template)  # the profile download leaves the page usable
        index = (app_module.Path(app_module.__file__).parent / "templates" / "index.html").read_text(encoding="utf-8")
        self.assertIn('class="button-link report-link"', index)
        self.assertIn("Preparing ", index)


if __name__ == "__main__":
    unittest.main()
