"""Final polish, Milestone 3: upload controls first, a review summary first, and result navigation.

Positions on screen were measured in a real browser (see docs/FINAL_POLISH_PROGRESS.md); these
tests pin the order of the page parts that produce them. SYNTHETIC example files only.
"""

import io
import os
import re
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from app import STAGED, app

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "integrated"


class NavigationTest(unittest.TestCase):
    def setUp(self):
        os.environ["AI_ENRICHMENT_ENABLED"] = "false"
        self.client = app.test_client()
        self.addCleanup(STAGED.clear)

    def test_home_puts_the_file_guide_and_upload_before_the_extras(self):
        page = self.client.get("/").get_data(as_text=True)
        order = [page.index(text) for text in ['class="quick-guide"', 'id="files"', 'fieldset class="enrich-toggle"',
                                               "Review files</button>", 'id="file-details"']]
        self.assertEqual(order, sorted(order))
        guide = page[page.index('class="quick-guide"'):page.index('id="files"')]
        for text in ["Which files do I need?", "you never need every type", "Individual sequences",
                     "Pooled SNV data", "complete A/C/G/T nucleotide counts", "reference genome FASTA",
                     "gene annotation (GenBank/GFF3)", "supplementary sample metadata"]:
            self.assertIn(text, guide)
        how = page[page.index('id="how-it-works"'):]
        self.assertTrue(how.split(">", 2)[1].strip().startswith("<details") and "<details>" in how[:200])

    def upload_example(self):
        names = ["synthetic_individual_sequences.fasta", "synthetic_individual_metadata.csv",
                 "synthetic_reference.fasta", "synthetic_pooled_snvs.csv", "synthetic_pooled_base_counts.csv"]
        files = [(io.BytesIO((EXAMPLES / name).read_bytes()), name) for name in names]
        page = self.client.post("/review", data={"files": files}, content_type="multipart/form-data")
        page = page.get_data(as_text=True)
        return page, re.search(r'name="token" value="([^"]+)"', page).group(1)

    def test_review_summary_comes_first_and_links_to_run(self):
        page, token = self.upload_example()
        self.assertLess(page.index('id="review-status"'), page.index("1. Files recognized"))
        self.assertIn("Ready to run.", page)
        self.assertIn('href="#run-analysis"', page)
        self.assertIn('id="run-analysis"', page)
        # A problem is listed at the top, once, before the long details.
        broken = self.client.post("/review", data={"token": token, "role_0": "ignore", "role_1": "ignore",
                                                   "role_2": "ignore", "role_3": "ignore", "role_4": "ignore"})
        text = broken.get_data(as_text=True)
        top = text[text.index('id="review-status"'):text.index("1. Files recognized")]
        self.assertIn("Required information to fix before running", top)
        self.assertEqual(text.count("Required information to fix before running"), 1)

    def test_results_navigation_reaches_every_part(self):
        _, token = self.upload_example()
        with redirect_stdout(io.StringIO()):
            page = self.client.post("/run", data={"token": token}).get_data(as_text=True)
        nav = page[page.index('<nav class="sections"'):page.index("</nav>", page.index('<nav class="sections"'))]
        for target in ["#downloads", "#overview", "#cross-source", "#depth-aware", "#sequencing-uncertainty",
                       "#mutation-frequency-groups", "#provenance"]:
            self.assertIn(f'href="{target}"', nav)
        self.assertIn('id="downloads"', page)
        # Jump targets inside Sprint 5 are panels, not sections: they need the same offset below the
        # sticky bar, or their heading lands under it (found on staging, 2026-11-16).
        self.assertIn(".uncertainty-panel[id] { scroll-margin-top: calc(var(--jump-bar-height", page)
        self.assertIn('setProperty("--jump-bar-height"', page)
        self.assertLess(page.index("Major findings"), page.index("3. Genomic variation</h2>"))


if __name__ == "__main__":
    unittest.main()
