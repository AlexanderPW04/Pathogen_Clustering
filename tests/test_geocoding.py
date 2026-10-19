import os
import unittest
from unittest.mock import Mock, patch

import pandas as pd
import requests

from geocode_locations import fill_missing_coordinates


class GeocodingTests(unittest.TestCase):
    def test_provided_coordinates_are_kept_without_a_request(self):
        data = pd.DataFrame(
            {"Location_Name": ["Phoenix"], "Latitude": [33.4], "Longitude": [-112.0]}
        )

        with patch.dict(os.environ, {"GEONAMES_USERNAME": "student"}):
            with patch("geocode_locations.requests.get") as get:
                result = fill_missing_coordinates(data)

        get.assert_not_called()
        self.assertEqual(result.loc[0, "Latitude"], 33.4)
        self.assertEqual(result.loc[0, "Longitude"], -112.0)
        self.assertEqual(result.loc[0, "Coordinate_Source"], "Provided")

    def test_missing_username_uses_location_fallback(self):
        data = pd.DataFrame(
            {"Location_Name": ["Phoenix"], "Latitude": [None], "Longitude": [None]}
        )

        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
            with patch("geocode_locations.requests.get") as get:
                result = fill_missing_coordinates(data)

        get.assert_not_called()
        self.assertTrue(pd.isna(result.loc[0, "Latitude"]))
        self.assertTrue(pd.isna(result.loc[0, "Coordinate_Source"]))

    def test_repeated_locations_share_one_lookup(self):
        data = pd.DataFrame(
            {
                "Location_Name": ["Phoenix", " phoenix  ", "Phoenix"],
                "Latitude": [None, None, 33.0],
                "Longitude": [None, None, None],
            }
        )
        response = Mock()
        response.json.return_value = {"geonames": [{"lat": "33.45", "lng": "-112.07"}]}

        with patch.dict(os.environ, {"GEONAMES_USERNAME": "student"}):
            with patch("geocode_locations.requests.get", return_value=response) as get:
                result = fill_missing_coordinates(data)

        self.assertEqual(get.call_count, 1)
        self.assertEqual(result["Latitude"].tolist(), [33.45] * 3)
        self.assertEqual(result["Longitude"].tolist(), [-112.07] * 3)
        self.assertEqual(result["Coordinate_Source"].tolist(), ["GeoNames"] * 3)
        self.assertEqual(
            get.call_args.args[0], "https://secure.geonames.org/searchJSON"
        )
        self.assertIs(get.call_args.kwargs["verify"], True)

    def test_state_only_location_is_not_sent_for_a_centroid(self):
        data = pd.DataFrame(
            {"Location_Name": ["USA: ND"], "Latitude": [None], "Longitude": [None]}
        )

        with patch.dict(os.environ, {"GEONAMES_USERNAME": "student"}):
            with patch("geocode_locations.requests.get") as get:
                result = fill_missing_coordinates(data)

        get.assert_not_called()
        self.assertTrue(pd.isna(result.loc[0, "Latitude"]))
        self.assertIn("too coarse", result.loc[0, "Geographic_Precision"])

    def test_administrative_geonames_result_is_marked_coarse(self):
        data = pd.DataFrame(
            {"Location_Name": ["North Dakota"], "Latitude": [None], "Longitude": [None]}
        )
        response = Mock()
        response.json.return_value = {
            "geonames": [{"lat": "47.5", "lng": "-100.5", "fcode": "ADM1"}]
        }

        with patch.dict(os.environ, {"GEONAMES_USERNAME": "student"}):
            with patch("geocode_locations.requests.get", return_value=response):
                result = fill_missing_coordinates(data)

        self.assertEqual(result.loc[0, "Coordinate_Source"], "GeoNames")
        self.assertIn("too coarse", result.loc[0, "Geographic_Precision"])

    def test_service_failure_leaves_coordinates_missing(self):
        data = pd.DataFrame(
            {"Location_Name": ["Phoenix"], "Latitude": [None], "Longitude": [None]}
        )

        with patch.dict(os.environ, {"GEONAMES_USERNAME": "student"}):
            with patch(
                "geocode_locations.requests.get", side_effect=requests.Timeout()
            ):
                result = fill_missing_coordinates(data)

        self.assertTrue(pd.isna(result.loc[0, "Latitude"]))
        self.assertTrue(pd.isna(result.loc[0, "Longitude"]))


if __name__ == "__main__":
    unittest.main()
