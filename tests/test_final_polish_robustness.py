"""Final polish, Milestone 4: mapping-profile validation and the shared result cache.

4A  An uploaded profile is untrusted: unknown keys, roles, fields, or option values the review form
    would refuse make the whole profile refused with a specific message; legitimate profiles work.
4B  The result cache is bounded by memory and count, evicts the least recently used result set first,
    keeps result sets apart, and expires them after ten minutes. SYNTHETIC data only.
"""

import copy
import json
import os
import threading
import unittest
from unittest.mock import patch

import app as app_module
from app import RESULT_DOWNLOADS, RESULT_LAST_USED, RESULT_LOCK, app, store_result_files
from genomics.mapping_profiles import PROFILE_FORMAT, PROFILE_VERSION, read_profile

VALID = {
    "profile_format": PROFILE_FORMAT, "profile_version": PROFILE_VERSION, "created": "2026-11-16",
    "description": "How to read these files", "options": {"date_order": "auto"},
    "files": [{"role": "snv", "kind": "Table", "format": None, "layout": "Long table: one row per sample and variant",
               "sheet": None, "fields": {"Sample_ID": "sample", "Position": "POS", "Alt_Frequency": "AF"},
               "options": {"frequency_scale": "percent", "position_base": "0", "header_row": 2,
                           "layout_confirmed": True, "base_count_minimum": 2, "declared_reference": "NC_045512.2"}},
              {"role": "reference", "kind": "FASTA sequences", "fields": {}, "options": {}}],
}


def encoded(profile):
    return json.dumps(profile).encode("utf-8")


class ProfileValidationTest(unittest.TestCase):
    def changed(self, edit):
        profile = copy.deepcopy(VALID)
        edit(profile)
        return encoded(profile)

    def refused(self, edit, words):
        with self.assertRaises(ValueError) as raised:
            read_profile(self.changed(edit))
        self.assertIn(words, str(raised.exception))

    def test_a_legitimate_profile_is_accepted_unchanged(self):
        self.assertEqual(read_profile(encoded(VALID)), VALID)

    def test_invalid_option_values_are_refused_with_the_rule(self):
        def option(name, value):
            return lambda profile: profile["files"][0]["options"].__setitem__(name, value)

        self.refused(option("header_row", -5), "'header_row' is -5; it must be a whole number from 1 to 50")
        self.refused(option("header_row", 0), "from 1 to 50")
        self.refused(option("header_row", "abc"), "from 1 to 50")
        self.refused(option("header_row", 2.5), "invalid value")
        self.refused(option("position_base", "2"), "must be 0 or 1")
        self.refused(option("frequency_scale", "ratio"), "'fraction' or 'percent'")
        self.refused(option("base_count_minimum", 5000), "from 1 to 1000")
        self.refused(option("depth_definition", "made-up"), "a value the review page offers")
        self.refused(option("run_shell", "rm"), "unknown option(s): run_shell")

    def test_structure_roles_fields_and_text_are_checked(self):
        self.refused(lambda p: p.__setitem__("extra", 1), "unknown field(s): extra")
        self.refused(lambda p: p["files"][0].__setitem__("role", "../../evil"), "unknown role")
        self.refused(lambda p: p["files"][0].__setitem__("role", "ignore"), "unknown role")
        self.refused(lambda p: p["files"][0]["fields"].__setitem__("Shell", "x"), "unknown field 'Shell'")
        self.refused(lambda p: p["files"][0]["fields"].__setitem__("Position", ""), "empty column name")
        self.refused(lambda p: p["files"][0].__setitem__("sheet", "x" * 300), "plain text")
        self.refused(lambda p: p["files"].__setitem__(0, "a string"), "must be a JSON object")
        self.refused(lambda p: p["files"][0].__setitem__("fields", ["Position"]), "JSON object")
        self.refused(lambda p: p.__setitem__("files", []), "1 to 20 files")
        self.refused(lambda p: p.__setitem__("files", [VALID["files"][1]] * 21), "1 to 20 files")
        self.refused(lambda p: p["options"].__setitem__("date_order", "guess"), "date order")
        self.refused(lambda p: p.__setitem__("created", "yesterday"), "creation date")

    def test_web_upload_of_an_invalid_profile_explains_it(self):
        os.environ["AI_ENRICHMENT_ENABLED"] = "false"
        client = app.test_client()
        import io
        csv = ("Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth\n"
               "S1,SYN.1,2,C,T,0.1,100\nS2,SYN.1,7,G,A,0.2,100\n")
        page = client.post("/review", data={"files": [(io.BytesIO(csv.encode()), "snv.csv")]},
                           content_type="multipart/form-data").get_data(as_text=True)
        import re
        token = re.search(r'name="token" value="([^"]+)"', page).group(1)
        bad = self.changed(lambda p: p["files"][0]["options"].__setitem__("header_row", -5))
        page = client.post("/review", data={"token": token, "mapping_profile": (io.BytesIO(bad), "profile.json")},
                           content_type="multipart/form-data").get_data(as_text=True)
        self.assertIn("it must be a whole number from 1 to 50", page)
        self.assertIn(token, page)  # the review stays open; nothing was applied
        app_module.STAGED.clear()


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


class ResultCacheTest(unittest.TestCase):
    def setUp(self):
        with RESULT_LOCK:
            self.saved = (dict(RESULT_DOWNLOADS), dict(RESULT_LAST_USED))
            RESULT_DOWNLOADS.clear()
            RESULT_LAST_USED.clear()
        self.clock = FakeClock()
        patcher = patch.object(app_module, "monotonic", self.clock)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.client = app.test_client()

    def tearDown(self):
        with RESULT_LOCK:
            RESULT_DOWNLOADS.clear()
            RESULT_LAST_USED.clear()
            RESULT_DOWNLOADS.update(self.saved[0])
            RESULT_LAST_USED.update(self.saved[1])

    def store(self, text, size=10):
        return store_result_files({"table.csv": ((text * size).encode(), "text/csv")})

    def test_an_active_visitor_keeps_results_while_another_runs_many_analyses(self):
        with patch.object(app_module, "RESULT_DOWNLOAD_LIMIT", 3):
            mine = self.store("mine")
            others = []
            for number in range(5):                    # another visitor's analyses, about one a minute and a half
                self.clock.now += 100
                self.assertEqual(self.client.get(f"/results/{mine}/table.csv").status_code, 200)  # I keep using mine
                self.clock.now += 1
                others.append(self.store(f"other{number}"))
            self.assertIn(mine, RESULT_DOWNLOADS)      # still downloadable while it is used (within 10 minutes)
            self.assertEqual(self.client.get(f"/results/{mine}/table.csv").data, b"mine" * 10)
            self.assertEqual(len(RESULT_DOWNLOADS), 3)  # the count limit holds
            self.assertNotIn(others[0], RESULT_DOWNLOADS)  # their own oldest, unused set went first
            self.assertIn(others[-1], RESULT_DOWNLOADS)    # their newest set (just created) is kept

    def test_a_just_finished_analysis_is_never_the_one_evicted(self):
        with patch.object(app_module, "RESULT_DOWNLOAD_LIMIT", 3):
            first = self.store("first")
            self.clock.now += 10
            self.client.get(f"/results/{first}/table.csv")   # downloaded a moment ago
            self.clock.now += 1
            second = self.store("second")
            self.clock.now += 1
            third = self.store("third")
            self.clock.now += 1
            newest = self.store("newest")
            self.assertIn(newest, RESULT_DOWNLOADS)
            self.assertIn(third, RESULT_DOWNLOADS)
            self.assertNotIn(first, RESULT_DOWNLOADS)       # least recently used of the three older sets
            self.assertIn(second, RESULT_DOWNLOADS)

    def test_memory_budget_is_never_exceeded(self):
        import os as operating_system
        with patch.object(app_module, "RESULT_CACHE_BYTES", 3000):
            for number in range(10):
                self.clock.now += 1
                store_result_files({"x.bin": (operating_system.urandom(900), "application/octet-stream")})
                total = sum(sum(len(c) for c, _ in item[1].values()) for item in RESULT_DOWNLOADS.values())
                self.assertLessEqual(total, 3000)

    def test_results_expire_and_stay_separate(self):
        first, second = self.store("alpha"), self.store("beta")
        self.assertNotEqual(first, second)
        self.assertGreaterEqual(len(first), 20)  # 128-bit random identifiers
        self.assertEqual(self.client.get(f"/results/{second}/table.csv").data, b"beta" * 10)
        self.clock.now += 601
        gone = self.client.get(f"/results/{first}/table.csv")
        self.assertEqual(gone.status_code, 404)
        self.assertIn(b"kept in memory for 10 minutes", gone.data)
        self.store("gamma")  # storing clears expired sets
        self.assertNotIn(first, RESULT_DOWNLOADS)
        self.assertNotIn(first, RESULT_LAST_USED)

    def test_concurrent_stores_stay_bounded(self):
        identifiers = []

        def worker(number):
            identifiers.append(self.store(f"t{number}", size=50))

        threads = [threading.Thread(target=worker, args=(number,)) for number in range(30)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(len(set(identifiers)), 30)
        self.assertLessEqual(len(RESULT_DOWNLOADS), app_module.RESULT_DOWNLOAD_LIMIT)
        self.assertEqual(set(RESULT_DOWNLOADS), set(RESULT_LAST_USED))


if __name__ == "__main__":
    unittest.main()
