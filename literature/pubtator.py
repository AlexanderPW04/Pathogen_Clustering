"""PubTator3: biomedical literature search and annotated article text.

PubTator3 returns titles, abstracts, and (for PMC Open Access articles,
which permit text mining) full text, split into passages with section names
and biomedical concept annotations (species, diseases, genes, variants).
It does not mark sample-collection locations. API:
https://www.ncbi.nlm.nih.gov/research/pubtator3/api
"""

from literature.http_client import ApiClient

PUBTATOR = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/"
ABSTRACT_SECTIONS = {"TITLE", "ABSTRACT"}


class PubTatorClient:
    def __init__(self, max_requests=40, session=None):
        # PubTator asks for no more than three requests per second.
        self.client = ApiClient("PubTator3", [PUBTATOR], 0.34, max_requests, timeout=30, session=session)

    def search(self, text, limit=5):
        """PMIDs of the top search results (candidates only, never accepted unchecked)."""

        result = self.client.get(PUBTATOR + "search/", {"text": text, "page": 1})
        return [str(item["pmid"]) for item in result.get("results", [])[:limit] if item.get("pmid")]

    def articles(self, pmids):
        """Annotated passages for PubMed IDs: {pmid: article} (missing PMIDs are absent)."""

        if not pmids:
            return {}
        result = self.client.get(PUBTATOR + "publications/export/biocjson",
                                 {"pmids": ",".join(pmids), "full": "true"})
        documents = result.get("PubTator3", []) if isinstance(result, dict) else []
        return {document_pmid(document): parse_document(document) for document in documents
                if document_pmid(document)}


def document_pmid(document):
    """The PMID of a BioC document (full-text documents use the PMC number as 'id')."""

    return str(document.get("pmid") or document.get("id") or "")


def section_name(infons):
    section = str(infons.get("section_type") or "").upper()
    if section:
        return section
    kind = str(infons.get("type", "")).casefold()
    return "TITLE" if kind in {"title", "front"} else "ABSTRACT" if kind == "abstract" else "OTHER"


def parse_document(document):
    """Turn a BioC JSON document into the passage list used by this package."""

    article_id = f"PMID:{document_pmid(document)}"
    passages = []
    for index, passage in enumerate(document.get("passages", [])):
        text = str(passage.get("text") or "")
        if not text.strip():
            continue
        infons = passage.get("infons", {})
        annotations = [{"type": item.get("infons", {}).get("type", ""), "text": item.get("text", "")}
                       for item in passage.get("annotations", [])[:30]]
        passages.append({
            "Article_ID": article_id, "Passage_ID": f"{article_id}#p{index}",
            "Section": section_name(infons), "Type": str(infons.get("type", "")),
            "Offset": passage.get("offset", 0), "Text": text, "Annotations": annotations,
        })
    full_text = any(passage["Section"] not in ABSTRACT_SECTIONS for passage in passages)
    status = ("Full text (PMC Open Access, via PubTator3)" if full_text
              else "Title and abstract only (via PubTator3)" if passages else "No text available")
    return {"Passages": passages, "Text_Status": status}


def check_pubtator():
    """One small free search, run only when the diagnostics button is pressed."""

    client = PubTatorClient(max_requests=1)
    client.client.retries = 0
    try:
        found = client.search("SARS-CoV-2 wastewater sequencing", limit=1)
    except Exception:
        return {"status": "Unavailable", "message": "Could not get a search result from PubTator3."}
    if not found:
        return {"status": "Check failed", "message": "PubTator3 answered but returned no results."}
    return {"status": "Operational", "message": "A PubTator3 literature search returned results."}
