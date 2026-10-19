"""Tests for the separate processed-measurement analysis mode."""

import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from environmental_analysis import load_measurements, make_charts, summarize_measurements


CSV_ROWS = (
    "Sample ID,Date,Site,Feature Identifier,Frequency,Quality Flag\n"
    "S1,2025-01-01,North,feature_a,0.1,Pass\n"
    "S1,2025-01-01,North,feature_b,0.4,Pass\n"
    "S2,2025-01-08,North,feature_a,0.3,Review\n"
    "S2,2025-01-08,North,feature_b,,\n"
    "S3,,South,feature_a,0.2,Pass\n"
)


class EnvironmentalLoaderTest(unittest.TestCase):
    def write_csv(self, folder, contents=CSV_ROWS):
        path = Path(folder) / "measurements.csv"
        path.write_text(contents, encoding="utf-8")
        return path

    def test_aliases_missing_values_and_summaries(self):
        with tempfile.TemporaryDirectory() as folder:
            data, metric = load_measurements(self.write_csv(folder))

        self.assertEqual(metric, "Frequency")
        self.assertEqual(data.loc[0, "Sample_ID"], "S1")
        self.assertEqual(data.loc[0, "Location_Name"], "North")
        self.assertEqual(data["Sample_ID"].nunique(), 3)
        self.assertEqual(int(data["Collection_Date"].isna().sum()), 1)
        self.assertEqual(int(data["Frequency"].isna().sum()), 1)

        summary = summarize_measurements(data, metric)
        self.assertEqual(summary["counts"]["samples"], 3)
        self.assertEqual(summary["counts"]["features"], 2)
        self.assertEqual(summary["counts"]["rows"], 5)
        self.assertEqual(len(summary["tables"]["temporal_changes"]), 4)
        changes = summary["tables"]["temporal_changes"]
        missing_feature = changes.loc[
            changes["Feature_ID"].eq("feature_b")
            & changes["Collection_Date"].eq(pd.Timestamp("2025-01-08"))
        ]
        self.assertTrue(missing_feature["Mean_Value"].isna().all())
        alpha_change = changes.loc[
            changes["Feature_ID"].eq("feature_a")
            & changes["Collection_Date"].eq(pd.Timestamp("2025-01-08")),
            "Change_From_Previous",
        ].iloc[0]
        self.assertAlmostEqual(alpha_change, 0.2)
        per_sample = summary["tables"]["samples"].set_index("Sample_ID")
        self.assertEqual(per_sample.loc["S2", "Features_Listed"], 2)
        self.assertEqual(per_sample.loc["S2", "Missing_Values"], 1)
        locations = summary["tables"]["location_summary"]
        north_alpha = locations.loc[
            locations["Feature_ID"].eq("feature_a")
            & locations["Location_Name"].eq("North"), "Mean_Value"
        ].iloc[0]
        self.assertAlmostEqual(north_alpha, 0.2)
        missing = summary["tables"]["missing_data"].set_index("Field")
        self.assertEqual(missing.loc["Collection_Date", "Missing_Rows"], 1)
        self.assertEqual(missing.loc["Frequency", "Missing_Rows"], 1)
        self.assertEqual(summary["tables"]["quality_summary"]["Rows"].sum(), 5)

    def test_measurement_value_column_is_supported(self):
        rows = CSV_ROWS.replace("Frequency", "Measurement Value")
        with tempfile.TemporaryDirectory() as folder:
            data, metric = load_measurements(self.write_csv(folder, rows))
        self.assertEqual(metric, "Measurement_Value")
        self.assertTrue(data[metric].dropna().between(0, 1).all())

    def test_excel_input_is_supported(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "measurements.xlsx"
            pd.DataFrame({
                "Sample_ID": ["S1"], "Collection_Date": ["2025-01-01"],
                "Location_Name": ["North"], "Feature_ID": ["feature_a"],
                "Measurement_Value": [12.5],
            }).to_excel(path, index=False)
            data, metric = load_measurements(path)
        self.assertEqual(metric, "Measurement_Value")
        self.assertEqual(data.loc[0, metric], 12.5)

    def test_invalid_measurement_is_reported(self):
        rows = CSV_ROWS.replace("feature_a,0.1,Pass", "feature_a,unknown,Pass")
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "non-numeric"):
                load_measurements(self.write_csv(folder, rows))

    def test_frequency_must_be_between_zero_and_one(self):
        rows = CSV_ROWS.replace("feature_a,0.1,Pass", "feature_a,25,Pass")
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "fraction from 0 to 1"):
                load_measurements(self.write_csv(folder, rows))

    def test_chart_is_a_png_image(self):
        with tempfile.TemporaryDirectory() as folder:
            data, metric = load_measurements(self.write_csv(folder))
        summary = summarize_measurements(data, metric)
        image = make_charts(data, metric, summary["tables"])
        self.assertTrue(image.startswith(b"\x89PNG\r\n\x1a\n"))

    def test_features_only_input_keeps_missing_values_and_counts_all_features(self):
        rows = "Sample_ID,Feature_ID\n001,feature_a\n001,feature_b\n002,feature_c\n"
        with tempfile.TemporaryDirectory() as folder:
            data, metric = load_measurements(self.write_csv(folder, rows))
        summary = summarize_measurements(data, metric)
        self.assertIsNone(metric)
        self.assertEqual(summary["counts"]["samples"], 2)
        self.assertEqual(summary["counts"]["features"], 3)
        self.assertEqual(summary["counts"]["frequency_values"], 0)
        self.assertEqual(summary["counts"]["coverage_values"], 0)
        self.assertEqual(data.loc[0, "Sample_ID"], "001")
        self.assertTrue(data["Collection_Date"].isna().all())
        self.assertTrue(summary["tables"]["date_counts"].empty)
        self.assertTrue(summary["tables"]["location_counts"].empty)
        self.assertTrue(make_charts(data, metric, summary["tables"]).startswith(b"\x89PNG"))

    def test_missing_frequency_is_distinct_from_reported_zero_and_coverage(self):
        rows = (
            "Sample_ID,Feature_ID,Frequency,Coverage,Read_Depth\n"
            "S1,feature_a,0,25.5,100\n"
            "S2,feature_b,,,\n"
        )
        with tempfile.TemporaryDirectory() as folder:
            data, metric = load_measurements(self.write_csv(folder, rows))
        summary = summarize_measurements(data, metric)
        self.assertEqual(data.loc[0, "Frequency"], 0)
        self.assertTrue(pd.isna(data.loc[1, "Frequency"]))
        self.assertTrue(pd.isna(data.loc[1, "Coverage"]))
        self.assertEqual(summary["counts"]["features"], 2)
        self.assertEqual(summary["counts"]["frequency_values"], 1)
        self.assertEqual(summary["counts"]["coverage_values"], 1)
        self.assertEqual(summary["counts"]["read_depth_values"], 1)


class EnvironmentalFlaskTest(unittest.TestCase):
    def test_synthetic_example_renders_and_downloads_without_sprint2(self):
        from app import app

        with patch("app.run_pipeline") as sprint2:
            with app.test_client() as client:
                page = client.get("/environmental")
                response = client.post("/environmental", data={"run_example": "yes"})
                html = response.get_data(as_text=True)
                self.assertEqual(page.status_code, 200)
                self.assertIn("Environmental Data Analysis", page.get_data(as_text=True))
                self.assertIn("SYNTHETIC EXAMPLE", html)
                self.assertIn("Unique genomic features:</strong> 2", html)
                self.assertIn("data:image/png;base64,", html)
                link = html.split('href="/environmental/download/')[1].split('"')[0]
                download = client.get("/environmental/download/" + link)
        sprint2.assert_not_called()
        self.assertEqual(download.status_code, 200)
        self.assertIn(b"Sample_ID,Collection_Date", download.data)

    def test_uploaded_csv_renders_sample_summary(self):
        from app import app

        with patch("app.run_pipeline") as sprint2:
            with app.test_client() as client:
                response = client.post(
                    "/environmental",
                    data={"measurement_file": (io.BytesIO(CSV_ROWS.encode()), "measurements.csv")},
                    content_type="multipart/form-data",
                )
        sprint2.assert_not_called()
        html = response.get_data(as_text=True)
        self.assertEqual(response.status_code, 200)
        self.assertIn("measurements.csv", html)
        self.assertIn("Total samples:</strong> 3", html)
        self.assertIn("Measurements by sample", html)

    def test_uploaded_csv_displays_validation_error(self):
        from app import app

        with app.test_client() as client:
            response = client.post(
                "/environmental",
                data={"measurement_file": (io.BytesIO(b"Sample_ID,Frequency\nS1,0.1\n"), "bad.csv")},
                content_type="multipart/form-data",
            )
        self.assertIn("Missing column(s)", response.get_data(as_text=True))

    def test_wastewater_measurement_page_upload_and_download(self):
        from app import app

        rows = (
            "Sample_ID,Feature_ID,Frequency\n"
            "S1,feature_a,0\nS2,feature_b,\n"
        )
        with patch("app.run_pipeline") as sprint2:
            with app.test_client() as client:
                page = client.get("/wastewater-data")
                response = client.post(
                    "/wastewater-data",
                    data={"measurement_file": (io.BytesIO(rows.encode()), "measurements.csv")},
                    content_type="multipart/form-data",
                )
                html = response.get_data(as_text=True)
                self.assertIn("Wastewater Data Analysis", page.get_data(as_text=True))
                self.assertIn("Total samples:</strong> 2", html)
                self.assertIn("Unique genomic features:</strong> 2", html)
                self.assertIn("Coverage information unavailable", html)
                link = html.split('href="/wastewater-data/download/')[1].split('"')[0]
                download = client.get("/wastewater-data/download/" + link)
        sprint2.assert_not_called()
        self.assertEqual(response.status_code, 200)
        self.assertEqual(download.status_code, 200)
        organized = pd.read_csv(io.BytesIO(download.data))
        self.assertEqual(organized.loc[0, "Frequency"], 0)
        self.assertTrue(pd.isna(organized.loc[1, "Frequency"]))

    def test_wastewater_excel_without_frequencies_remains_usable(self):
        from app import app

        workbook = io.BytesIO()
        pd.DataFrame({"Sample_ID": ["001"], "Feature_ID": ["feature_a"]}).to_excel(
            workbook, index=False
        )
        workbook.seek(0)
        with app.test_client() as client:
            response = client.post(
                "/wastewater-data",
                data={"measurement_file": (workbook, "features.xlsx")},
                content_type="multipart/form-data",
            )
        self.assertEqual(response.status_code, 200)
        html = response.get_data(as_text=True)
        self.assertIn("Frequency information unavailable", html)
        self.assertIn("Coverage information unavailable", html)
        self.assertIn("data:image/png;base64,", html)

    def test_missing_or_invalid_upload_and_expired_download_are_readable(self):
        from app import ENVIRONMENTAL_DOWNLOADS, app

        with app.test_client() as client:
            missing = client.post("/wastewater-data")
            invalid = client.post(
                "/wastewater-data",
                data={"measurement_file": (io.BytesIO(b"bad"), "bad.txt")},
                content_type="multipart/form-data",
            )
            with patch.dict(ENVIRONMENTAL_DOWNLOADS, {"old": (0, b"data")}):
                expired = client.get("/wastewater-data/download/old")
        self.assertIn("Choose a measurement", missing.get_data(as_text=True))
        self.assertIn("must use the .csv or .xlsx", invalid.get_data(as_text=True))
        self.assertEqual(expired.status_code, 404)


if __name__ == "__main__":
    unittest.main()
