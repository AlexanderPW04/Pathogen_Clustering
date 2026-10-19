import unittest

import pandas as pd

from add_temporal_geographic import enrich_pairs
from detect_potential_clusters import find_related_pairs


class TemporalGeographicTests(unittest.TestCase):
    def make_pair(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "Sample_ID_1": ["S1"],
                "Sample_ID_2": ["S2"],
                "Pathogen": ["Pathogen A"],
                "Gene_Segment": ["Segment 1"],
                "Differences": [0],
                "Comparable_Positions": [4],
                "P_Distance": [0.0],
                "Percent_Identity": [100.0],
            }
        )

    def make_metadata(
        self,
        date_1: str = "2026-01-01",
        date_2: str = "2026-01-01",
        latitude_1: float | None = 0.0,
        longitude_1: float | None = 0.0,
        latitude_2: float | None = 0.0,
        longitude_2: float | None = 0.0,
    ) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "Sample_ID": ["S1", "S2"],
                "Collection_Date": pd.to_datetime([date_1, date_2]),
                "Latitude": [latitude_1, latitude_2],
                "Longitude": [longitude_1, longitude_2],
                "Location_Name": ["Phoenix", "  phoenix  "],
            }
        )

    def test_two_samples_on_same_date(self) -> None:
        result = enrich_pairs(self.make_metadata(), self.make_pair())

        self.assertEqual(result.loc[0, "Days_Apart"], 0)
        self.assertTrue(result.loc[0, "Same_Location"])

    def test_two_samples_with_different_dates(self) -> None:
        metadata = self.make_metadata(date_2="2026-01-11")
        result = enrich_pairs(metadata, self.make_pair())

        self.assertEqual(result.loc[0, "Days_Apart"], 10)

    def test_two_samples_at_same_coordinates(self) -> None:
        result = enrich_pairs(self.make_metadata(), self.make_pair())

        self.assertEqual(result.loc[0, "Km_Apart"], 0.0)

    def test_two_samples_with_different_coordinates(self) -> None:
        metadata = self.make_metadata(longitude_2=1.0)
        result = enrich_pairs(metadata, self.make_pair())

        self.assertAlmostEqual(result.loc[0, "Km_Apart"], 111.195, places=3)

    def test_missing_coordinates(self) -> None:
        metadata = self.make_metadata(latitude_2=None)
        result = enrich_pairs(metadata, self.make_pair())

        self.assertTrue(pd.isna(result.loc[0, "Km_Apart"]))

    def test_coordinate_sources_are_kept_in_pairwise_output(self) -> None:
        metadata = self.make_metadata()
        metadata.loc[1, "Latitude"] = 33.45
        metadata.loc[1, "Longitude"] = -112.07
        metadata["Coordinate_Source"] = ["Provided", "GeoNames"]

        result = enrich_pairs(metadata, self.make_pair())

        self.assertEqual(result.loc[0, "Coordinate_Source_1"], "Provided")
        self.assertEqual(result.loc[0, "Coordinate_Source_2"], "GeoNames")
        self.assertEqual(result.loc[0, "Latitude_2"], 33.45)
        self.assertEqual(result.loc[0, "Longitude_2"], -112.07)

    def test_state_names_cannot_pass_as_fifty_kilometer_sites(self) -> None:
        metadata = self.make_metadata()
        metadata["Location_Name"] = ["USA: ND", " usa: nd "]
        result = enrich_pairs(metadata, self.make_pair())

        self.assertTrue(pd.isna(result.loc[0, "Km_Apart"]))
        self.assertTrue(pd.isna(result.loc[0, "Same_Location"]))
        self.assertIn("too coarse", result.loc[0, "Geographic_Precision_1"])
        self.assertEqual(result.loc[0, "Latitude_1"], 0.0)
        related = find_related_pairs(result, 0.02, 14, 50)
        self.assertTrue(related.empty)

    def test_administrative_geonames_centroid_cannot_pass(self) -> None:
        metadata = self.make_metadata()
        metadata["Geographic_Precision"] = [
            "Administrative area (too coarse for distance cutoff)",
            "GeoNames place (approximate)",
        ]
        result = enrich_pairs(metadata, self.make_pair())

        self.assertTrue(pd.isna(result.loc[0, "Km_Apart"]))
        self.assertTrue(pd.isna(result.loc[0, "Same_Location"]))


if __name__ == "__main__":
    unittest.main()
