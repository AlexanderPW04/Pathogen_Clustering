"""Literature sources with one interface: live NCBI + PubTator3, or a local corpus.

The local corpus holds the clearly labeled synthetic articles used by the
offline demonstration and the automated tests. Both sources answer the same
four questions, so discovery and retrieval code does not care which is used.
"""

import json
from pathlib import Path

from literature.config import ncbi_settings
from literature.identifiers import PATTERNS, is_standard, mentions
from literature.ncbi import NcbiClient
from literature.pubtator import PubTatorClient


class LiveSource:
    label = "NCBI E-utilities and PubTator3 (live public services)"
    synthetic = False

    def __init__(self, session=None):
        settings = ncbi_settings()
        self.ncbi = NcbiClient(settings["api_key"], settings["email"], settings["max_requests"], session)
        self.pubtator = PubTatorClient(max(10, settings["max_requests"] // 2), session)

    def record_links(self, kind, value):
        """[(article_id, relationship, evidence, method)] from one identifier."""

        if not is_standard(kind, value):
            return []  # project-specific identifiers cannot be resolved by NCBI
        if kind == "PMID":
            return [(f"PMID:{value}", "provided", f"PMID {value} is listed in the sample metadata.",
                     "Metadata publication field")]
        if kind in {"PMCID", "DOI"}:
            pmid = self.ncbi.pmid_for(kind, value)
            if pmid is None:
                return []
            return [(f"PMID:{pmid}", "provided", f"{kind} {value} in the sample metadata resolves to PMID {pmid}.",
                     "Metadata publication field + PubMed search")]
        if kind in {"Nucleotide", "BioSample", "BioProject", "SRA"}:
            return [(f"PMID:{pmid}", relationship, evidence, "NCBI E-utilities record links")
                    for pmid, relationship, evidence in self.ncbi.pubmed_for_record(kind, value)]
        return []

    def searchable(self, alias):
        return any(PATTERNS[kind].match(alias) for kind in ["Nucleotide", "SRA", "BioSample"])

    def text_search(self, alias):
        return [f"PMID:{pmid}" for pmid in self.pubtator.search(alias, limit=3)]

    def describe(self, article_ids):
        pmids = [article_id.split(":", 1)[1] for article_id in article_ids]
        return {f"PMID:{pmid}": summary for pmid, summary in self.ncbi.summaries(pmids).items()}

    def fetch(self, article_ids):
        pmids = [article_id.split(":", 1)[1] for article_id in article_ids]
        return {f"PMID:{pmid}": article for pmid, article in self.pubtator.articles(pmids).items()}

    def failures(self):
        return self.ncbi.client.failures + self.pubtator.client.failures

    def requests_sent(self):
        return self.ncbi.client.requests_sent + self.pubtator.client.requests_sent


class LocalCorpusSource:
    """Articles from a JSON file (the synthetic demonstration corpus)."""

    synthetic = True

    def __init__(self, path):
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        self.label = data.get("label", "Local article corpus")
        self.articles = {article["Article_ID"]: article for article in data["articles"]}

    def record_links(self, kind, value):
        links = []
        for article_id, article in self.articles.items():
            record_links = article.get("Record_Links", {})
            if kind == "Publication" and value == article_id:
                links.append((article_id, "provided", f"{value} is listed in the sample metadata.",
                              "Metadata publication field"))
            elif kind == "BioProject" and value in record_links.get("study", []):
                links.append((article_id, "study", f"Study record {value} links {article_id}.",
                              "Study record link (synthetic)"))
            elif value in record_links.get("direct", []):
                links.append((article_id, "direct", f"Sample record {value} links {article_id}.",
                              "Sample record link (synthetic)"))
        return links

    def searchable(self, alias):
        return True

    def text_search(self, alias):
        # Loose matching on purpose, like a search engine: callers must verify.
        needle = alias.casefold()
        return [article_id for article_id, article in self.articles.items()
                if any(needle in passage["Text"].casefold() for passage in article.get("Passages", []))]

    def describe(self, article_ids):
        return {article_id: {key: self.articles[article_id].get(key, "") for key in
                             ["Title", "Publication_Date", "PMID", "PMCID", "DOI", "Source_URL"]}
                for article_id in article_ids if article_id in self.articles}

    def fetch(self, article_ids):
        fetched = {}
        for article_id in article_ids:
            article = self.articles.get(article_id)
            if article is None:
                continue
            passages = [{"Article_ID": article_id, "Passage_ID": f"{article_id}#{item['Passage_ID']}",
                         "Section": item["Section"], "Type": item.get("Type", "paragraph"),
                         "Offset": 0, "Text": item["Text"], "Annotations": []}
                        for item in article.get("Passages", [])]
            fetched[article_id] = {"Passages": passages,
                                   "Text_Status": article.get("Text_Status", "Synthetic text")
                                   if passages else "No text available"}
        return fetched

    def failures(self):
        return []

    def requests_sent(self):
        return 0


def text_mentions(article, alias):
    return any(mentions(passage["Text"], alias) for passage in article.get("Passages", []))
