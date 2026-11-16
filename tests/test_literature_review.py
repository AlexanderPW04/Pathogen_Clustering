"""Sprint 4 literature review on NCBI-shaped data: batched requests, outcomes, and batches.

Every external service is replaced: a fake HTTP session answers like NCBI
E-utilities and PubTator3 did for the real records QB228092.1-QB228101.1
(BioProject 980729 linked to PMID 35302933, two BioSamples), an in-memory
source supplies invented articles, and the extraction "model" is the offline
mock provider with fixed answers. No NCBI, PubTator, or Claude request is
sent, and nothing is billed.
"""

import io
import json
import os
import re
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fake_mafft import make_fake_mafft  # noqa: E402
from test_literature_eligibility import OFFLINE_AI, FakeSource, passage, ten_records  # noqa: E402

import literature_view  # noqa: E402
from app import STAGED, app  # noqa: E402
from literature.http_client import ApiClient, RequestLimitReached, ServiceUnavailable  # noqa: E402
from literature.providers import MockProvider  # noqa: E402
from literature.review import FAILED, PREDATES, REVIEWED_NONE, review_summary  # noqa: E402
from literature.sources import LiveSource  # noqa: E402
from literature.workflow import find_literature, new_session, run_analysis  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REVIEWED_NO_CHANGE = ("One publication was found and reviewed. It predates the uploaded samples and "
                      "does not provide verifiable collection locations for them. No metadata was changed.")


class Response:
    def __init__(self, payload, status=200, headers=None):
        self.status_code, self.payload, self.headers = status, payload, headers or {}

    def json(self):
        return self.payload


class NcbiShapedSession:
    """Answers like E-utilities and PubTator3 did for the ten real records; records every call."""

    def __init__(self, pubtator_status=200, ncbi_status=200):
        self.calls = []
        self.pubtator_status, self.ncbi_status = pubtator_status, ncbi_status

    def get(self, url, params=None, timeout=None):
        params = dict(params or {})
        self.calls.append((url.rsplit("/", 1)[-1] or "search", params.get("dbfrom"), params.get("db"), params.get("id")))
        if "pubtator3-api" in url:
            if self.pubtator_status != 200:
                return Response({}, self.pubtator_status)
            return Response({"PubTator3": [{"pmid": 35302933, "passages": [
                {"infons": {"type": "title"}, "text": "Invented title about wild birds in 2021.", "offset": 0},
                {"infons": {"type": "abstract"}, "offset": 40,
                 "text": "Wild bird samples collected in December 2021 were screened for influenza."}]}]})
        if self.ncbi_status != 200:
            return Response({}, self.ncbi_status)
        if url.endswith("esummary.fcgi") and params["db"] == "nuccore":
            accessions = params["id"].split(",")
            uids = [str(3417148500 + int(item.split(".")[0][-3:])) for item in accessions]
            return Response({"result": {"uids": uids, **{uid: {"caption": item.split(".")[0],
                                                               "accessionversion": item.split(".")[0] + ".1"}
                                                         for uid, item in zip(uids, accessions)}}})
        if url.endswith("elink.fcgi"):
            sets = []
            for uid in params["id"]:
                links = []
                if params["dbfrom"] == "nuccore" and params["db"] == "biosample":
                    links = ["63447206" if int(uid) % 1000 >= 600 else "63447205"]
                elif params["dbfrom"] == "nuccore" and params["db"] == "bioproject":
                    links = ["980729"]
                elif params["dbfrom"] == "bioproject" and params["db"] == "pubmed":
                    links = ["35302933"]
                sets.append({"ids": [uid], "linksetdbs": [{"linkname": f"{params['dbfrom']}_{params['db']}",
                                                            "links": links}] if links else []})
            return Response({"linksets": sets})
        if url.endswith("esummary.fcgi") and params["db"] == "biosample":
            return Response({"result": {uid: {"accession": f"SAMN{uid}"} for uid in params["id"].split(",")}})
        if url.endswith("esummary.fcgi") and params["db"] == "pubmed":
            return Response({"result": {"35302933": {"title": "Invented study of wild birds, 2021",
                                                     "pubdate": "2022 May", "epubdate": "2022 Mar 18",
                                                     "articleids": [{"idtype": "pmc", "value": "PMC0000001"}]}}})
        return Response({}, 404)


def real_shaped_records():
    """Accession without a version (as in the NCBI CSV) and the FASTA ID with one."""

    records = ten_records()
    records["Accession"] = records["Accession"].str.split(".").str[0]
    return records


def live_review(session_http):
    session = new_session(individual=real_shaped_records())
    with patch("literature.http_client.time.sleep"):
        find_literature(session, LiveSource(session=session_http), max_articles=5, max_retrieved=20)
    return session


class HttpClientTest(unittest.TestCase):
    prefix = "https://example.org/"

    def client(self, responses, **options):
        class Session:
            def __init__(self):
                self.count = 0

            def get(self, url, params=None, timeout=None):
                self.count += 1
                return responses.pop(0)
        session = Session()
        return ApiClient("Test service", [self.prefix], 0, session=session, **options), session

    def test_http_429_is_retried_a_bounded_number_of_times(self):
        client, session = self.client([Response({}, 429) for _ in range(5)], retries=2)
        with patch("literature.http_client.time.sleep") as sleep, self.assertRaises(ServiceUnavailable) as raised:
            client.get(self.prefix + "a")
        self.assertEqual(session.count, 3)  # first try + 2 retries
        self.assertEqual(client.requests_sent, 3)  # retries count toward the request budget
        self.assertIn("rate limit reached (HTTP 429) after 3 attempt(s)", str(raised.exception))
        self.assertLessEqual(max(call.args[0] for call in sleep.call_args_list), 4)  # 1 s, 2 s backoff

    def test_retry_after_is_honored_and_long_waits_stop_the_search(self):
        client, session = self.client([Response({}, 429, {"Retry-After": "3"}), Response({"ok": 1})])
        with patch("literature.http_client.time.sleep") as sleep:
            self.assertEqual(client.get(self.prefix + "b"), {"ok": 1})
        self.assertGreaterEqual(max(call.args[0] for call in sleep.call_args_list), 2.9)
        client, session = self.client([Response({}, 503, {"Retry-After": "Wed, 21 Oct 2099 07:28:00 GMT"})])
        with patch("literature.http_client.time.sleep"), self.assertRaises(ServiceUnavailable) as raised:
            client.get(self.prefix + "c")
        self.assertEqual(session.count, 1)
        self.assertIn("asked to wait", str(raised.exception))

    def test_request_budget_is_finite(self):
        client, _ = self.client([Response({}, 503), Response({"ok": 1}), Response({"ok": 2})], max_requests=2)
        with patch("literature.http_client.time.sleep"):
            self.assertEqual(client.get(self.prefix + "d"), {"ok": 1})  # 503, then success: 2 requests
            with self.assertRaises(RequestLimitReached):
                client.get(self.prefix + "e")

    def test_successful_responses_are_cached_and_shared_between_searches(self):
        url = self.prefix + f"shared-{os.getpid()}-{id(self)}"
        with patch("literature.http_client.requests.Session.get", return_value=Response({"ok": 1})) as get, \
                patch("literature.http_client.time.sleep"):
            first = ApiClient("Shared test service", [self.prefix], 0)
            second = ApiClient("Shared test service", [self.prefix], 0)
            self.assertEqual(first.get(url, {"id": ["1", "2"]}), {"ok": 1})
            self.assertEqual(second.get(url, {"id": ["1", "2"]}), {"ok": 1})
        self.assertEqual(get.call_count, 1)
        self.assertEqual(second.cache_hits, 1)
        self.assertIs(first.limiter, second.limiter)  # one rate limiter per service


class NcbiShapedReviewTest(unittest.TestCase):
    def test_ten_accessions_sharing_one_pmid_need_few_requests_and_one_retrieval(self):
        http = NcbiShapedSession()
        session = live_review(http)
        table = session["publications"]
        self.assertEqual(len(table), 10)  # one association row per accession, kept separate
        self.assertEqual(table["Article_ID"].unique().tolist(), ["PMID:35302933"])
        self.assertIn("QB228092.1 belongs to BioProject uid 980729", table["Association_Evidence"].iloc[0])
        names = [call[0] for call in http.calls]
        self.assertLessEqual(len(http.calls), 9)
        self.assertEqual(names.count("biocjson"), 1)  # article text fetched once
        nuccore = [call for call in http.calls if call[2] == "nuccore"]
        self.assertEqual(len(nuccore), 1)
        self.assertEqual(len(nuccore[0][3].split(",")), 10)  # accession and accession.version deduplicated
        project = [call for call in http.calls if call[1] == "bioproject"]
        self.assertEqual([call[3] for call in project], [["980729"]])  # one shared BioProject, looked up once
        usage = {item["Service"]: item for item in session["request_usage"]}
        self.assertEqual(usage["NCBI E-utilities"]["Requests sent"], 8)
        self.assertEqual(usage["PubTator3"]["Requests sent"], 1)
        self.assertIn("SAMN63447205", session["aliases"]["QB228092.1"])
        self.assertIn("SAMN63447206", session["aliases"]["QB228100.1"])

        summary = review_summary(session)
        card = summary["cards"][0]
        self.assertEqual(card["outcome"], PREDATES)
        self.assertEqual((card["published"], card["epub"]), ("2022 May", "2022 Mar 18"))
        self.assertEqual(card["retrieval"], "Title and abstract only (via PubTator3)")  # no full text
        self.assertEqual(summary["headline"], REVIEWED_NO_CHANGE)
        stats = dict(summary["stats"])
        self.assertEqual((stats["Accession records inspected"], stats["Distinct isolates"],
                          stats["Unique publications found"], stats["Publications retrieved"],
                          stats["Publications reviewed"], stats["Study-level links"],
                          stats["Samples with proposed improvements"]), (10, 2, 1, 1, 1, 10, 0))
        self.assertEqual({row["Status"] for row in summary["samples"]},
                         {"Linked publication predates this sample's collection"})
        self.assertEqual(session["plan"], [])  # nothing to send to Claude

    def test_text_retrieval_failure_is_an_incomplete_review(self):
        session = live_review(NcbiShapedSession(pubtator_status=503))
        summary = review_summary(session)
        self.assertEqual(summary["cards"][0]["outcome"], FAILED)
        self.assertIn("The evidence review is incomplete", summary["headline"])
        self.assertIn("This does not mean the literature has no evidence", summary["headline"])
        self.assertNotIn("No verifiable", summary["headline"])
        self.assertTrue(summary["failures"])
        self.assertEqual(summary["cards"][0]["published"], "2022 May")  # details still retrieved separately

    def test_ncbi_outage_is_reported_as_an_incomplete_search(self):
        session = live_review(NcbiShapedSession(ncbi_status=503))
        summary = review_summary(session)
        self.assertIn("The publication search could not be completed", summary["headline"])
        self.assertIn("does not mean that no publication exists", summary["headline"])
        self.assertEqual({row["Status"] for row in summary["samples"]},
                         {"Search incomplete (a service was unavailable)"})

    def test_exact_accession_in_a_study_article_makes_it_usable_for_that_sample(self):
        article = {"Title": "Invented study", "Publication_Date": "2026 Dec", "Passages": [
            passage("PMID:9200001", 1, "RESULTS", "Sequence QB228093 came from turkeys on a farm in Kingsbury "
                                                  "County, South Dakota.")]}
        source = FakeSource({"PMID:9200001": article}, {accession: [("PMID:9200001", "study")]
                                                       for accession in ["QB228092", "QB228093"]})
        session = new_session(individual=ten_records(["QB228092", "QB228093"]))
        find_literature(session, source)
        table = session["publications"].set_index("Sample_ID")
        self.assertTrue(table.loc["QB228093.1", "Eligible_For_Enrichment"])
        self.assertFalse(table.loc["QB228092.1", "Eligible_For_Enrichment"])

    def test_reviewed_article_without_supported_improvements(self):
        article = {"Title": "Invented report", "Publication_Date": "2026 Dec", "Passages": [
            passage("PMID:9200002", 1, "METHODS", "Sample QB228092 was sequenced on one run; where it was "
                                                  "collected is not reported.")]}
        source = FakeSource({"PMID:9200002": article}, {"QB228092": [("PMID:9200002", "direct")]})
        session = new_session(individual=ten_records(["QB228092"]))
        find_literature(session, source)
        with tempfile.TemporaryDirectory() as folder:
            mock = Path(folder) / "mock.json"
            mock.write_text(json.dumps({"extraction": {"PMID:9200002": {
                "facts": [], "insufficient_evidence_samples": ["QB228092.1"], "notes": ""}}}), encoding="utf-8")
            settings = {**literature_view.ai_settings(), "max_articles": 5, "max_requests": 10,
                        "max_input_tokens": 60000}
            run_analysis(session, MockProvider(mock), settings)
        summary = review_summary(session)
        self.assertEqual(summary["cards"][0]["outcome"], REVIEWED_NONE)
        self.assertEqual(summary["samples"][0]["Status"], "Reviewed; no supported improvement")
        self.assertIn("It does not provide verifiable collection locations for the uploaded samples",
                      summary["headline"])
        self.assertEqual(summary["cards"][0]["claude"], "Yes")


# --- Ten samples, ten publications, mixed outcomes -------------------------------------------

CITIES = ["Tempe", "Mesa", "Chandler", "Gilbert", "Scottsdale"]
SAMPLES = [f"MXD{number:02d}" for number in range(1, 11)]
BASES = "ACGT"


def mixed_articles():
    """Five usable, three without collection details, one ambiguous, one contradictory (all invented)."""

    articles, responses = {}, {}
    for index, sample in enumerate(SAMPLES):
        article_id = f"PMID:{9300001 + index}"
        if index < 5:
            sentence = (f"Sample {sample} was collected at the wastewater sampling site in {CITIES[index]}, "
                        f"Arizona, on 2026-03-{index + 1:02d}.")
            facts = [{"Sample_ID": sample, "Field": "Location", "Value": f"{CITIES[index]}, Arizona, USA",
                      "Location_Role": "Sample collection site", "Geographic_Level": "City or town",
                      "Passage_ID": f"{article_id}#p1", "Evidence_Quote": sentence,
                      "Evidence_Type": "Directly stated", "Ambiguity": "", "Reasoning_Summary": "Named site."}]
            insufficient = []
        elif index < 8:
            sentence = f"Sample {sample} was sequenced with a tiled amplicon protocol; collection details are not given."
            facts, insufficient = [], [sample]
        elif index == 8:
            sentence = f"Samples including {sample} were collected at one of several plants in the Phoenix area."
            facts = [{"Sample_ID": sample, "Field": "Location", "Value": "Phoenix, Arizona, USA",
                      "Location_Role": "Sample collection site", "Geographic_Level": "City or town",
                      "Passage_ID": f"{article_id}#p1", "Evidence_Quote": sentence, "Evidence_Type": "Ambiguous",
                      "Ambiguity": "Several plants; this sample's plant is not stated.",
                      "Reasoning_Summary": "Area only."}]
            insufficient = []
        else:
            sentence = f"Sample {sample} was collected in Denver, Colorado."
            facts = [{"Sample_ID": sample, "Field": "Location", "Value": "Denver, Colorado, USA",
                      "Location_Role": "Sample collection site", "Geographic_Level": "City or town",
                      "Passage_ID": f"{article_id}#p1", "Evidence_Quote": sentence,
                      "Evidence_Type": "Directly stated", "Ambiguity": "", "Reasoning_Summary": "Named city."}]
            insufficient = []
        articles[article_id] = {"Title": f"Invented article {index + 1}", "Publication_Date": "2026 Dec",
                                "Passages": [passage(article_id, 1, "METHODS", sentence)]}
        responses[article_id] = {"facts": facts, "insufficient_evidence_samples": insufficient, "notes": ""}
    return articles, {"extraction": responses, "association": {}}


class CountingSource(FakeSource):
    def __init__(self, *arguments, **options):
        super().__init__(*arguments, **options)
        self.fetched = []

    def fetch(self, article_ids):
        self.fetched += list(article_ids)
        return super().fetch(article_ids)


def mixed_files():
    sequences, rows = [], []
    reference = "".join(BASES[(index * 7) % 4] for index in range(48))
    for index, sample in enumerate(SAMPLES):
        sequence = reference[:index] + BASES[(BASES.index(reference[index]) + 1) % 4] + reference[index + 1:]
        sequences.append(f">{sample} |Synthetic test virus segment 4\n{sequence}\n")
        rows.append({"Accession": sample, "Pathogen": "Synthetic test virus", "Gene_Segment": "Segment 4",
                     "Collection_Date": f"2026-03-{index + 1:02d}", "Location_Name": "USA: AZ",
                     "Host_Source": "Synthetic source"})
    return "".join(sequences), pd.DataFrame(rows).to_csv(index=False), reference


class MixedOutcomeTest(unittest.TestCase):
    def setUp(self):
        patches = [patch.dict(os.environ, {**OFFLINE_AI, "AI_MAX_ARTICLES_PER_RUN": "5"}),
                   patch("literature.config.ENV_FILE", ROOT / "missing.env"),
                   patch("align_sequences.shutil.which", return_value="mafft"),
                   patch.dict(app.config, {"LITERATURE_INLINE_JOBS": True})]
        for item in patches:
            item.start()
            self.addCleanup(item.stop)
        self.addCleanup(STAGED.clear)
        self.client = app.test_client()

    def post(self, url, data):
        with redirect_stdout(io.StringIO()):
            response = self.client.post(url, data=data)
        if response.status_code == 302:
            response = self.client.get(response.headers["Location"])
        return response.get_data(as_text=True)

    def test_ten_publications_with_mixed_outcomes(self):
        articles, responses = mixed_articles()
        source = CountingSource(articles, {sample: [(f"PMID:{9300001 + index}", "direct")]
                                           for index, sample in enumerate(SAMPLES)})
        sequences, metadata, reference = mixed_files()
        files = [(io.BytesIO(sequences.encode()), "mixed.fasta"), (io.BytesIO(metadata.encode()), "mixed.csv")]
        page = self.client.post("/review", data={"files": files, "enrich": "on", "intended_mode": "individual"},
                                content_type="multipart/form-data").get_data(as_text=True)
        staged = re.search(r'name="token" value="([^"]+)"', page).group(1)
        with tempfile.TemporaryDirectory() as folder:
            mock = Path(folder) / "mixed_mock.json"
            mock.write_text(json.dumps(responses), encoding="utf-8")
            with patch("literature_view.literature_source", return_value=source), \
                    patch("literature_view.choose_provider", return_value=(MockProvider(mock), "Mock (test)")), \
                    patch("align_sequences.subprocess.run", side_effect=make_fake_mafft(reference)):
                page = self.post("/enrich/start", {"token": staged, "enrich": "on", "enrich_shown": "1"})
                token = re.search(r'name="token" value="([^"]+)"', page[page.index('id="literature"'):]).group(1)
                found = page  # the free search ran when Sprint 4 started
                estimate = self.post("/literature/analyze", {"token": token})
                first = self.post("/literature/analyze", {"token": token, "confirm": "1", "batch": "1"})
                second = self.post("/literature/analyze", {"token": token, "confirm": "1", "batch": "2"})
                session = literature_view.get_session(token)
                settings = literature_view.run_settings(session)
                facts = {fact["Sample_ID"]: fact for fact in session["facts"]}
                # Try to approve every fact, including the ambiguous and contradictory ones.
                choices = {f"review_{fact['Fact_ID']}": "approve" for fact in session["facts"]}
                result = self.post("/enrich/run", {"token": token, **choices})

        # Each unique publication is retrieved once, and the ten associations stay separate.
        self.assertEqual(sorted(source.fetched), sorted(articles))
        self.assertIn("Unique publications found:</strong> 10", found)
        self.assertIn("Publications with a direct sample link:</strong> 10", found)
        self.assertIn("analyzed in 2 batches of at most 5 (the per-run limit)", found)
        # Batches: each estimated and run separately; the paid limit is never raised.
        self.assertIn("Batch 1:</strong> PMID:9300001, PMID:9300002, PMID:9300003, PMID:9300004, PMID:9300005",
                      estimate)
        self.assertEqual(settings["max_articles"], 5)
        self.assertEqual(len(re.findall(r'name="review_(F\d+)" value="approve"', first)), 5)
        self.assertEqual(len(re.findall(r'name="review_(F\d+)" value="approve"', second)), 5)  # kept from batch 1
        self.assertLessEqual(session["run"]["usage"]["Requests sent"], settings["max_requests"])
        # Valid improvements only for the matching samples; ambiguous and contradictory need review.
        self.assertEqual(sorted(sample for sample, fact in facts.items() if fact["Approvable"]), SAMPLES[:5])
        self.assertEqual(facts["MXD09"]["Evidence_Status"], "Ambiguous")
        self.assertEqual(facts["MXD10"]["Evidence_Status"], "Contradictory")
        self.assertNotIn("MXD06", facts)
        self.assertEqual(session["decisions"][facts["MXD09"]["Fact_ID"]], "reject")  # approval was refused
        # The analysis ran with five enriched and five original locations.
        self.assertIn("Pairwise comparisons:</strong> 45", result)
        report = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", result[result.index('id="enrichment-report"'):]))
        for index, sample in enumerate(SAMPLES[:5]):
            self.assertIn(f"{sample} Location USA: AZ {CITIES[index]}, Arizona, USA Approved "
                          f"{CITIES[index]}, Arizona, USA PMID:{9300001 + index}", report)
        for sample in SAMPLES[5:8]:
            self.assertIn(f"{sample} — unchanged — No proposal Original kept —", report)
        self.assertIn("MXD09 Location USA: AZ Phoenix, Arizona, USA Not approvable (Ambiguous)", report)
        self.assertIn("MXD10 Location USA: AZ Denver, Colorado, USA Not approvable (Contradictory)", report)
        enriched = pd.read_csv(io.StringIO(self.client.get(re.search(
            r'href="(/results/[^"]+/literature_enriched_metadata\.csv)"', result).group(1)).get_data(as_text=True)))
        enriched = enriched.set_index("Sample_ID")
        self.assertEqual(enriched.loc["MXD01", "Location_Name"], "Tempe, Arizona, USA")
        self.assertEqual(enriched.loc["MXD10", "Location_Name"], "USA: AZ")
        self.assertEqual(enriched.loc["MXD06", "Location_Name"], "USA: AZ")


if __name__ == "__main__":
    unittest.main()
