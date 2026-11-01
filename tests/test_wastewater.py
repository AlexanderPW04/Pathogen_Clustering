"""Focused tests for the separate processed-wastewater mode."""

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from wastewater.processed_snv import load_processed_snvs, summarize_snvs


STANDARD_ROWS = (
    "Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth,Collection_Date,Location_Name\n"
    "W1,NC_045512.2,241,C,T,0.2,100,2020-05-01,Site A\n"
    "W2,NC_045512.2,241,C,T,0.4,200,2020-05-02,Site B\n"
)


class WastewaterLoaderTest(unittest.TestCase):
    def test_published_column_mapping_and_reference(self):
        source = pd.DataFrame({
            "SAMPLE": ["253"], "COLLECTION_DATE": ["30_May_2020"],
            "STATE": ["AZ"], "LOCATION": ["Guadalupe"],
            "POSITION": ["13469"], "REF": ["g"], "ALT": ["T"],
            "ALT_FREQ": ["0.0284464"], "TOTAL_DP": ["1371"],
        })
        source.attrs["fontenele_2021_table2"] = True
        with patch("wastewater.processed_snv.read_table", return_value=source):
            data, warnings = load_processed_snvs("mmc2.xlsx")
        self.assertEqual(warnings, [])
        self.assertEqual(data.loc[0, "Sample_ID"], "253")
        self.assertEqual(data.loc[0, "Reference_Genome"], "NC_045512.2")
        self.assertEqual(data.loc[0, "Collection_Date"].date().isoformat(), "2020-05-30")
        self.assertEqual(data.loc[0, "Location_Name"], "AZ: Guadalupe")
        self.assertEqual(data.loc[0, "Ref"], "G")
        self.assertEqual(data.loc[0, "Read_Depth"], 1371)

    def test_paper_like_columns_without_source_note_need_reference(self):
        source = pd.DataFrame({
            "SAMPLE": ["W1"], "COLLECTION_DATE": ["30_May_2020"],
            "STATE": ["AZ"], "LOCATION": ["Site A"],
            "POSITION": ["241"], "REF": ["C"], "ALT": ["T"],
            "ALT_FREQ": ["0.2"], "TOTAL_DP": ["100"],
        })
        with patch("wastewater.processed_snv.read_table", return_value=source):
            with self.assertRaisesRegex(ValueError, "Reference_Genome"):
                load_processed_snvs("unknown.xlsx")

    def test_summary_counts_distinct_samples_not_duplicate_rows(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            path.write_text(STANDARD_ROWS + "W1,NC_045512.2,241,C,T,0.2,100,2020-05-01,Site A\n", encoding="utf-8")
            data, warnings = load_processed_snvs(path)
        summary = summarize_snvs(data)
        self.assertEqual(summary["samples"], 2)
        self.assertEqual(summary["records"], 3)
        self.assertEqual(summary["unique_snvs"], 1)
        self.assertEqual(summary["top_snvs"][0]["Samples_Observed"], 2)
        self.assertIn("repeated", warnings[0])

    def test_invalid_numeric_and_nucleotide_values(self):
        bad = STANDARD_ROWS.replace(",241,C,T,0.2,100,", ",0,C,C,1.2,bad,")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.csv"
            path.write_text(bad, encoding="utf-8")
            with self.assertRaises(ValueError) as caught:
                load_processed_snvs(path)
        message = str(caught.exception)
        self.assertIn("position", message)
        self.assertIn("frequency", message)
        self.assertIn("read depth", message)
        self.assertIn("identical", message)

    def test_metadata_join_and_missing_sample_detection(self):
        snvs = STANDARD_ROWS.replace(",2020-05-01,Site A", ",,").replace(",2020-05-02,Site B", ",,")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            metadata = Path(folder) / "metadata.csv"
            path.write_text(snvs, encoding="utf-8")
            metadata.write_text("Sample_ID,Collection_Date,Location_Name\nW1,2020-05-01,Site A\nW2,2020-05-02,Site B\n", encoding="utf-8")
            data, _ = load_processed_snvs(path, metadata)
            self.assertEqual(data.loc[0, "Location_Name"], "Site A")
            self.assertEqual(data.loc[1, "Collection_Date"].date().isoformat(), "2020-05-02")
            metadata.write_text("Sample_ID,Collection_Date,Location_Name\nW1,2020-05-01,Site A\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "W2"):
                load_processed_snvs(path, metadata)

    def test_missing_depth_date_and_location_stay_missing(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            path.write_text(
                "Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth,Collection_Date,Location_Name\n"
                "W1,NC_045512.2,241,C,T,0.2,,,\n",
                encoding="utf-8",
            )
            data, warnings = load_processed_snvs(path)
        summary = summarize_snvs(data)
        self.assertTrue(pd.isna(data.loc[0, "Read_Depth"]))
        self.assertIsNone(summary["date_start"])
        self.assertEqual(summary["locations"], [])
        self.assertEqual(summary["depth_count"], 0)
        self.assertEqual(len(warnings), 3)

    def test_conflicting_metadata_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "snvs.csv"
            metadata = Path(folder) / "metadata.csv"
            path.write_text(STANDARD_ROWS, encoding="utf-8")
            metadata.write_text("Sample_ID,Collection_Date,Location_Name\nW1,2020-05-01,Other Site\nW2,2020-05-02,Site B\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "conflict"):
                load_processed_snvs(path, metadata)

class WastewaterFlaskTest(unittest.TestCase):
    def test_wastewater_page_and_download_do_not_run_sprint_2(self):
        from app import app

        with patch("app.run_pipeline") as sprint2:
            with app.test_client() as client:
                home = client.get("/")
                page = client.get("/wastewater")
                response = client.post(
                    "/wastewater",
                    data={"snv_file": (io.BytesIO(STANDARD_ROWS.encode()), "calls.csv")},
                    content_type="multipart/form-data",
                )
                html = response.get_data(as_text=True)
                self.assertIn("Wastewater Surveillance", page.get_data(as_text=True))
                self.assertIn("Wastewater Surveillance", home.get_data(as_text=True))
                self.assertIn("Unique samples:</strong> 2", html)
                self.assertIn("SNV records:</strong> 2", html)
                self.assertIn("not</strong> confirmed absent", html)
                self.assertIn("calls.csv", html)
                link = html.split('href="/wastewater/download/')[1].split('"')[0]
                download = client.get("/wastewater/download/" + link)
        sprint2.assert_not_called()
        self.assertEqual(download.status_code, 200)
        self.assertIn("text/csv", download.content_type)
        self.assertIn(b"W1,NC_045512.2,241,C,T", download.data)

    def test_wastewater_page_shows_validation_error(self):
        from app import app

        with app.test_client() as client:
            response = client.post(
                "/wastewater",
                data={"snv_file": (io.BytesIO(b"Sample_ID,Position\nW1,241\n"), "bad.csv")},
                content_type="multipart/form-data",
            )
        self.assertIn("Missing processed SNV column", response.get_data(as_text=True))


if __name__ == "__main__":
    unittest.main()
