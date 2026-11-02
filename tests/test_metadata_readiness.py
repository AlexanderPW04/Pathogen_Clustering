"""Analysis-specific metadata readiness on the dataset review step.

Uses the synthetic NCBI-style fixture (tests/fixtures/ncbi_virus_structure_*)
with columns removed or blanked to make partially and highly incomplete
files. MAFFT is mocked and GeoNames disabled; no external request is sent.
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
from literature.preanalysis import gap_resolution  # noqa: E402
from literature.workflow import INDIVIDUAL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
INTEGRATED = ROOT / "examples" / "integrated"
METADATA = pd.read_csv(FIXTURES / "ncbi_virus_structure_metadata.csv", dtype=str, keep_default_na=False)
SEQUENCES = (FIXTURES / "ncbi_virus_structure_sequences.fasta").read_text(encoding="utf-8")
FIRST = SEQUENCES.split("\n")[1]


def text(page):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page)).replace("&#39;", "'").replace("&gt;", ">")


def section(page, start, end):
    return text(page[page.index(start):page.index(end)])


def partial_metadata():
    """Dates for 7 of 10 records; two records have a city-level location, two have none."""

    table = METADATA.drop(columns=["USA"]).copy()
    table.loc[7:9, "Collection_Date"] = ""
    table.loc[0:1, "Geo_Location"] = "USA: SD, Sioux Falls"
    table.loc[8:9, ["Geo_Location", "Country"]] = ""
    return table


def sparse_metadata():
    return METADATA[["Accession", "Organism_Name", "Segment"]]


class ReadinessTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {"GEONAMES_USERNAME": ""}),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(FIRST))]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def review(self, metadata=None, sequences=SEQUENCES, files=None, **form):
        files = files or [(io.BytesIO(sequences.encode()), "sequences.fasta"),
                          (io.BytesIO(metadata.to_csv(index=False).encode()), "metadata.csv")]
        page = self.client.post("/review", data={"files": files, **form},
                                content_type="multipart/form-data").get_data(as_text=True)
        return page, re.search(r'name="token" value="([^"]+)"', page).group(1)

    def readiness(self, page):
        return section(page, 'id="readiness"', "Analyses available with these files")

    def plan(self, page):
        return section(page, "Analyses available with these files", "Sprint 4 · Metadata enrichment (optional)")

    def run_analysis(self, token, **form):
        with redirect_stdout(io.StringIO()):
            return self.client.post("/run", data={"token": token, **form}).get_data(as_text=True)

    def pairs(self, result):
        link = re.search(r'href="(/results/[^"]+temporal_geographic\.csv)"', result).group(1)
        return pd.read_csv(io.StringIO(self.client.get(link).get_data(as_text=True)))

    # Complete, partial, and highly incomplete NCBI-style files -------------------------
    def test_complete_ncbi_file(self):
        page, token = self.review(METADATA, intended_mode="individual")
        readiness = self.readiness(page)
        self.assertIn("Metadata readiness for: Individual sequences", readiness)
        self.assertIn("All required fields were found", readiness)
        self.assertIn("Required Sample identifier Detected Column 'Accession': 10 of 10 rows", readiness)
        self.assertIn("Required Pathogen or organism Detected Column 'Organism_Name'", readiness)
        self.assertIn("Required Genomic region (segment or gene) Detected Column 'Segment'", readiness)
        self.assertIn("Recommended Collection date Detected", readiness)
        self.assertIn("Recommended Collection location Imprecise Column 'Geo_Location': 10 of 10 sequences "
                      "have a location; 10 only at state level (for example 'USA: SD')", readiness)
        self.assertIn("Collection locations are only available at the state level. Genetic comparisons can "
                      "proceed, but precise geographic-distance screening cannot be performed", readiness)
        self.assertNotIn("may be inaccurate", page)
        # Source columns no field uses are listed as kept, never dropped; Segment is used.
        kept = readiness[readiness.index("Other columns kept unchanged"):]
        for column in ["Species", "Genotype", "Isolate", "Submitters"]:
            self.assertIn(column, kept)
        self.assertNotIn("Segment", kept.split("What the missing")[0])
        plan = self.plan(page)
        self.assertIn("Temporal screening of sequence pairs Will run", plan)
        self.assertIn("Geographic-distance screening Unavailable", plan)
        self.assertNotIn("Pooled variant QC", plan)  # another mode's analysis, no files here
        result = self.run_analysis(token)
        self.assertIn("Pairwise comparisons:</strong> 2", result)

    def test_partially_complete_ncbi_file(self):
        page, token = self.review(partial_metadata(), intended_mode="individual")
        readiness = self.readiness(page)
        self.assertIn("All required fields were found", readiness)
        self.assertIn("Recommended Collection date Partly filled Column 'Collection_Date': 7 of 10 sequences "
                      "have a date", readiness)
        self.assertIn("3 of 10 sequences have no collection date; pairs involving them get no days-apart "
                      "value and cannot meet the temporal criterion", readiness)
        self.assertIn("Collection location Imprecise", readiness)
        self.assertIn("6 of 10 sequences have only state-level locations", readiness)
        self.assertIn("2 of 10 sequences have no collection location", readiness)
        plan = self.plan(page)
        self.assertIn("Temporal screening of sequence pairs Limited", plan)
        self.assertIn("Geographic-distance screening Limited", plan)
        self.assertIn("2 specific place name(s), located with GeoNames during the run; 6 only at state or "
                      "country level", plan)
        result = self.run_analysis(token)
        self.assertIn("Pairwise comparisons:</strong> 2", result)

    def test_highly_incomplete_ncbi_file_still_runs_the_genetic_comparison(self):
        page, token = self.review(sparse_metadata(), intended_mode="individual")
        readiness = self.readiness(page)
        self.assertIn("All required fields were found", readiness)
        self.assertIn("Recommended Collection date Missing", readiness)
        self.assertIn("Recommended Collection location Missing", readiness)
        self.assertIn("Optional Latitude and longitude Missing", readiness)
        self.assertIn("No collection dates were found for the individual sequences. Genetic comparisons can "
                      "proceed, but temporal screening (days apart) cannot be performed", readiness)
        self.assertIn("No collection locations were found for the individual sequences. Genetic and date "
                      "comparisons can proceed, but geographic-distance screening cannot be performed",
                      readiness)
        self.assertNotIn("disabled>Run analysis", page)
        plan = self.plan(page)
        self.assertIn("Sequence comparison (MAFFT, p-distance) Will run", plan)
        self.assertIn("Temporal screening of sequence pairs Unavailable", plan)
        result = self.run_analysis(token)
        self.assertIn("Pairwise comparisons:</strong> 2", result)
        pairs = self.pairs(result)
        # Missing values stay missing: no invented dates, zero distances, or coordinates.
        self.assertTrue(pairs["Days_Apart"].isna().all())
        self.assertTrue(pairs["Km_Apart"].isna().all())

    # Required information ---------------------------------------------------------------
    def test_missing_pathogen_is_the_only_blocker_and_can_be_stated(self):
        metadata = METADATA[["Accession", "Segment", "Collection_Date"]]
        page, token = self.review(metadata, intended_mode="individual")
        readiness = self.readiness(page)
        self.assertIn("Required information missing: Pathogen or organism (individual sequences)", readiness)
        self.assertIn("Required Pathogen or organism Missing Enter it under 'Pathogen of all individual "
                      "sequences' below", readiness)
        self.assertIn("disabled>Run analysis", page)
        page, _ = self.review(files=None, metadata=metadata, declared_pathogen="Influenza A virus",
                              intended_mode="individual")
        readiness = self.readiness(page)
        self.assertIn("Pathogen or organism Stated by you 'Influenza A virus'", readiness)
        self.assertIn("All required fields were found", readiness)

    def test_unknown_genomic_region_names_only_the_sequence_comparison(self):
        bare = re.sub(r"^(>\S+).*$", r"\1", SEQUENCES, flags=re.M)
        metadata = METADATA[["Accession", "Organism_Name", "Collection_Date"]]
        page, _ = self.review(metadata, sequences=bare, intended_mode="individual")
        readiness = self.readiness(page)
        self.assertIn("Required for one analysis: Genomic region (segment or gene): without it, sequence "
                      "comparison cannot run", readiness)
        self.assertIn("confirm complete genomes under Advanced settings", readiness)
        page, _ = self.review(metadata, sequences=bare, intended_mode="individual",
                              genomic_region_mode="complete_genomes")
        self.assertIn("Genomic region (segment or gene) Complete genomes confirmed", self.readiness(page))

    # Requirements depend on the analysis mode -------------------------------------------
    def test_pooled_mode_lists_only_pooled_requirements(self):
        files = [(io.BytesIO((INTEGRATED / "synthetic_pooled_snvs.csv").read_bytes()), "pooled.csv")]
        page, _ = self.review(files=files, intended_mode="pooled")
        readiness = self.readiness(page)
        self.assertIn("Metadata readiness for: Pooled sequencing results", readiness)
        self.assertIn("Required Variant frequency Detected", readiness)
        self.assertIn("Required Read depth Detected", readiness)
        self.assertIn("Optional Per-position A/C/G/T counts (separate table) Missing", readiness)
        self.assertIn("the study's PCoA cannot be reproduced", readiness)
        self.assertNotIn("Sample identifier", readiness.replace("Pooled sample identifier", ""))
        self.assertNotIn("Pathogen or organism", readiness)
        self.assertNotIn("Sequence comparison (MAFFT", self.plan(page))

    def test_combined_mode_needs_a_shared_reference_only_for_the_cross_source_comparison(self):
        names = ["synthetic_individual_sequences.fasta", "synthetic_individual_metadata.csv",
                 "synthetic_pooled_snvs.csv"]
        files = [(io.BytesIO((INTEGRATED / name).read_bytes()), name) for name in names]
        page, _ = self.review(files=files, intended_mode="combined")
        readiness = self.readiness(page)
        self.assertIn("Metadata readiness for: Combined individual and pooled analysis", readiness)
        self.assertIn("Required Shared reference genome (FASTA) Both sources Missing", readiness)
        self.assertIn("Shared reference genome (FASTA): without it, cross-source comparison cannot run. "
                      "Everything else can still run", readiness)
        self.assertNotIn("disabled>Run analysis", page)
        files = [(io.BytesIO((INTEGRATED / name).read_bytes()), name)
                 for name in names + ["synthetic_reference.fasta"]]
        page, _ = self.review(files=files, intended_mode="combined")
        self.assertIn("Shared reference genome (FASTA) Both sources Detected SYNREF.1", self.readiness(page))

    # Sprint 4 -----------------------------------------------------------------------------
    def test_sprint4_column_names_fields_publications_might_improve(self):
        page, _ = self.review(sparse_metadata(), intended_mode="individual", enrich="on")
        readiness = self.readiness(page)
        self.assertIn("Could Sprint 4 help?", readiness)
        self.assertIn("Collection date Missing No column was recognized Temporal screening (days apart) "
                      "Possibly. Publications may state collection dates or date ranges.", readiness)
        self.assertIn("Possibly. Publications may name a more specific collection site.", readiness)
        page, _ = self.review(sparse_metadata(), intended_mode="individual")
        readiness = self.readiness(page)
        self.assertNotIn("Could Sprint 4 help?", readiness)
        self.assertIn("could look in linked publications for: collection date, collection location, latitude "
                      "and longitude, host or sample source", readiness)

    def test_gap_resolution_counts_only_gaps_that_are_gone(self):
        session = {"gaps": [{"Sample_ID": "A", "Gaps": ["No collection date"]},
                            {"Sample_ID": "B", "Gaps": ["No collection date"]},
                            {"Sample_ID": "C", "Gaps": ["State-only location"]}]}
        derived = pd.DataFrame({
            "Sample_ID": ["A", "B", "C"], "Sample_Source": [INDIVIDUAL] * 3,
            "Location_Name": ["Tempe, Arizona", "Tempe, Arizona", "USA: SD"],
            "Latitude": [33.4, 33.4, pd.NA], "Longitude": [-111.9, -111.9, pd.NA],
            "Geographic_Precision": ["Provided coordinates"] * 2 + ["State-level location (too coarse)"],
            "Collection_Date": ["2026-03-05", pd.NA, pd.NA], "Host_Source": ["Human"] * 3,
            "Enrichment_Sources": ["Collection_Date from A1 (P1)", "Collection_Date_Range from A1 (P2)", pd.NA]})
        rows = {row["Field"]: row for row in gap_resolution(session, derived)}
        self.assertEqual((rows["Collection date"]["Resolved"], rows["Collection date"]["Improved, still limited"],
                          rows["Collection date"]["Unresolved"]), (1, 1, 0))
        self.assertIn("B (still: no collection date)", rows["Collection date"]["Improved samples"])
        self.assertEqual(rows["Collection location"]["Unresolved samples"], "C")


if __name__ == "__main__":
    unittest.main()
