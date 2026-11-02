"""NCBI Virus geographic columns: Geo_Location, Country, and the state column "USA".

The fixture tests/fixtures/ncbi_virus_structure_* copies the 17-column
structure of a real NCBI Virus download with SYNTHETIC values (no real
records). MAFFT is mocked and GeoNames disabled unless a test says otherwise.
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
from genomics.column_detection import detect_columns  # noqa: E402
from genomics.location_hierarchy import apply_location_hierarchy, location_hierarchy  # noqa: E402
from geocode_locations import fill_missing_coordinates  # noqa: E402
from validate_inputs import load_metadata, validate_metadata  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"
METADATA = FIXTURES / "ncbi_virus_structure_metadata.csv"
SEQUENCES = FIXTURES / "ncbi_virus_structure_sequences.fasta"
FIRST = SEQUENCES.read_text(encoding="utf-8").split("\n")[1]
COLUMNS = ("Accession,Organism_Name,Species,Genotype,Isolate,Segment,Length,Nuc_Completeness,Geo_Location,"
           "Country,USA,Host,Tissue_Specimen_Source,Submitters,Collection_Date,Release_Date,Molecule_type").split(",")


def matches(columns):
    return {match["field"]: match for match in detect_columns(columns, "metadata")}


class ColumnDetectionTest(unittest.TestCase):
    def test_the_17_ncbi_columns_map_without_ambiguity(self):
        found = matches(COLUMNS)
        self.assertEqual(found["Location_Name"]["column"], "Geo_Location")
        self.assertEqual(found["Country"]["column"], "Country")
        self.assertEqual(found["State_Province"]["column"], "USA")
        self.assertEqual(found["Pathogen"]["column"], "Organism_Name")
        self.assertEqual(found["Host_Source"]["column"], "Host")
        self.assertEqual(found["Accession"]["column"], "Accession")
        self.assertEqual(found["Collection_Date"]["column"], "Collection_Date")
        self.assertFalse([name for name, match in found.items() if match["status"] == "ambiguous"])

    def test_header_spelling_variants(self):
        found = matches(["accession", "GEO LOCATION", "country-name", "State", "organism name"])
        self.assertEqual(found["Location_Name"]["column"], "GEO LOCATION")
        self.assertEqual(found["Country"]["column"], "country-name")
        self.assertEqual(found["State_Province"]["column"], "State")
        found = matches(["Accession", "geo_loc_name", "Country"])
        self.assertEqual(found["Location_Name"]["column"], "geo_loc_name")

    def test_two_specific_location_columns_are_still_a_real_choice(self):
        found = matches(["Accession", "Geo_Location", "Location"])
        self.assertEqual(found["Location_Name"]["status"], "ambiguous")


class HierarchyTest(unittest.TestCase):
    def rows(self, **columns):
        result, summary = location_hierarchy(pd.DataFrame(columns), location_header="Geo_Location")
        return result, summary

    def test_state_level_ncbi_values(self):
        result, summary = self.rows(Location_Name=["USA: SD"], Country=["USA"], USA=["SD"])
        self.assertEqual(result.loc[0, "Location_Level"], "state")
        self.assertEqual(summary["state_level"], 1)

    def test_missing_geo_location_uses_the_broader_columns(self):
        result, _ = self.rows(Location_Name=["", "not provided"], Country=["USA", "Brazil"], USA=["SD", ""])
        self.assertEqual(list(result["Location_Name"]), ["USA: SD", "Brazil"])
        self.assertEqual(list(result["Location_Level"]), ["state", "country"])
        self.assertTrue(result["Location_Name_Source"].str.startswith("Built").all())

    def test_country_only_dataset(self):
        result, summary = self.rows(Country=["Brazil", "Brazil"])
        self.assertEqual(list(result["Location_Name"]), ["Brazil", "Brazil"])
        self.assertEqual(summary["country_level"], 2)

    def test_more_specific_location_is_kept_and_not_coarse(self):
        result, _ = self.rows(Location_Name=["USA: SD, Sioux Falls"], Country=["USA"], USA=["SD"])
        self.assertEqual(result.loc[0, "Location_Name"], "USA: SD, Sioux Falls")
        self.assertEqual(result.loc[0, "Location_Level"], "")

    def test_contradictory_columns_are_flagged(self):
        result, summary = self.rows(Location_Name=["USA: SD", "USA: SD"], Country=["Canada", "USA"], USA=["SD", "ND"])
        self.assertEqual(list(result["Location_Level"]), ["conflict", "conflict"])
        self.assertIn("different country", result.loc[0, "Location_Conflict"])
        self.assertIn("different state", result.loc[1, "Location_Conflict"])
        self.assertEqual(summary["conflicts"], 2)

    def test_explicit_location_is_never_overwritten(self):
        table = apply_location_hierarchy(pd.DataFrame({"Location_Name": ["Tempe, Arizona"], "Country": ["USA"]}))
        self.assertEqual(table.loc[0, "Location_Name"], "Tempe, Arizona")


class ValidationAndPrecisionTest(unittest.TestCase):
    def test_validated_metadata_keeps_source_fields_and_levels(self):
        validated, errors, _ = validate_metadata(load_metadata(METADATA), [])
        self.assertEqual(errors, [])
        self.assertTrue(validated["Location_Name"].eq("USA: SD").all())
        self.assertTrue(validated["Location_Level"].eq("state").all())
        for column in ["Country", "USA", "Isolate", "Genotype", "Species", "Submitters"]:
            self.assertIn(column, validated.columns)
        self.assertEqual(validated["Isolate"].nunique(), 2)
        self.assertEqual(sorted(validated["Segment"].astype(int).unique()), list(range(1, 9)))

    def test_country_and_conflict_levels_are_too_coarse_and_not_geocoded(self):
        samples = pd.DataFrame({"Location_Name": ["Brazil", "USA: SD", "Peru: Lima"], "Latitude": [pd.NA] * 3,
                                "Longitude": [pd.NA] * 3, "Location_Level": ["country", "conflict", pd.NA]})
        with patch.dict(os.environ, {"GEONAMES_USERNAME": "test-user"}), \
                patch("geocode_locations.requests.get") as get, redirect_stdout(io.StringIO()):
            get.return_value.json.return_value = {"geonames": []}
            result = fill_missing_coordinates(samples)
        self.assertEqual(result.loc[0, "Geographic_Precision"], "Country-level location (too coarse for distance cutoff)")
        self.assertIn("Conflicting location columns", result.loc[1, "Geographic_Precision"])
        queried = [call.kwargs["params"]["q"] for call in get.call_args_list]
        self.assertEqual(queried, ["Peru: Lima"])  # only the specific place is looked up


class ReviewAndRunTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": ""}),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(FIRST))]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def review(self, metadata_bytes=None, metadata_name="sequences (3).csv"):
        files = [(io.BytesIO(SEQUENCES.read_bytes()), "sequences.fasta"),
                 (io.BytesIO(metadata_bytes or METADATA.read_bytes()), metadata_name)]
        page = self.client.post("/review", data={"files": files}, content_type="multipart/form-data").get_data(as_text=True)
        return page, re.search(r'name="token" value="([^"]+)"', page).group(1)

    def test_review_recognizes_ncbi_geography_and_runs_two_within_segment_pairs(self):
        page, token = self.review()
        self.assertNotIn("Choose which column is Location_Name", page)
        self.assertIn("Geographic columns recognized: &#39;Geo_Location&#39; was selected for Location_Name; "
                      "&#39;Country&#39; and &#39;USA&#39; were kept separately", page)
        self.assertIn("10 of 10 record(s) have only country- or state-level locations", page)
        self.assertNotIn("disabled>Run analysis", page)
        with redirect_stdout(io.StringIO()):
            result = self.client.post("/run", data={"token": token}).get_data(as_text=True)
        self.assertIn("Unique samples analyzed:</strong> 10", result)
        self.assertIn("Pairwise comparisons:</strong> 2", result)  # segments 1 and 2 only
        self.assertIn("Strict related pairs:</strong> 0", result)
        self.assertIn("Genomic-temporal candidate pairs with unresolved geography:</strong> 2", result)
        link = re.search(r'href="(/results/[^"]+temporal_geographic\.csv)"', result).group(1)
        pairs = pd.read_csv(io.StringIO(self.client.get(link).get_data(as_text=True)))
        self.assertEqual(sorted(pairs["Gene_Segment"]), ["Segment 1", "Segment 2"])
        self.assertTrue(pairs["Km_Apart"].isna().all())
        self.assertTrue(pairs["Same_Location"].isna().all())

    def test_contradictions_are_shown_for_review(self):
        text = METADATA.read_text(encoding="utf-8").replace("USA: SD,USA,SD", "USA: SD,Canada,SD", 1)
        page, _ = self.review(text.encode())
        self.assertIn("Location values to review", page)
        self.assertIn("TST000001: Geo_Location &#39;USA: SD&#39; names a different country than &#39;Canada&#39;", page)

    def test_missing_geo_location_values_are_filled_with_their_level(self):
        text = METADATA.read_text(encoding="utf-8").replace("partial,USA: SD,USA,SD", "partial,,USA,SD", 1)
        page, token = self.review(text.encode())
        self.assertIn("1 record(s) without a &#39;Geo_Location&#39; value use the broader country/state", page)
        with redirect_stdout(io.StringIO()):
            result = self.client.post("/run", data={"token": token}).get_data(as_text=True)
        self.assertIn("Pairwise comparisons:</strong> 2", result)

    def test_tsv_and_xlsx_versions_are_recognized(self):
        table = pd.read_csv(METADATA, dtype=str)
        tsv = table.to_csv(sep="\t", index=False).encode()
        page, _ = self.review(tsv, "sequences.tsv")
        self.assertNotIn("Choose which column is Location_Name", page)
        workbook = io.BytesIO()
        table.to_excel(workbook, index=False)
        page, _ = self.review(workbook.getvalue(), "sequences.xlsx")
        self.assertNotIn("Choose which column is Location_Name", page)
        self.assertIn("Geographic columns recognized", page)


if __name__ == "__main__":
    unittest.main()
