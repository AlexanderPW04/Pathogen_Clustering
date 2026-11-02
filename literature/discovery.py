"""Find publications associated with samples and classify each association.

Only traceable links count as associations that can support sample-specific
metadata: a record link from the sample, a publication listed in its
metadata, or the sample's exact identifier inside the article text. A search
hit on its own is never accepted, and the top result is never chosen
automatically.
"""

import pandas as pd

from literature.http_client import RequestLimitReached, ServiceUnavailable
from literature.identifiers import sample_identifiers
from literature.sources import text_mentions

DIRECT = "Direct sample-publication association"
STUDY = "Verified study-publication association"
IN_TEXT = "Exact accession mentioned in article text"
POSSIBLE = "Possible association (context only)"
UNRELATED = "Unrelated or unsupported association"
CATEGORY_RANK = {DIRECT: 1, STUDY: 2, IN_TEXT: 3, POSSIBLE: 4, UNRELATED: 5}
ELIGIBLE = {DIRECT, STUDY, IN_TEXT}

PUBLICATION_COLUMNS = [
    "Sample_ID", "Sample_Source", "Article_ID", "PMID", "PMCID", "DOI", "Title",
    "Publication_Date", "Source_URL", "Relationship_Type", "Association_Evidence",
    "Verification_Status", "Retrieval_Method", "Identifier_In_Text", "Text_Status",
    "Ambiguity_Notes", "Eligible_For_Enrichment",
]


def discover_publications(samples, source, max_samples=25, max_articles=5):
    """Return (publication table, fetched articles, notes) for the selected samples."""

    notes = []
    selected = samples.head(max_samples)
    if len(samples) > max_samples:
        notes.append(f"Only the first {max_samples} selected samples were searched (request limit).")
    links = []  # one dictionary per (sample, article) candidate
    aliases_by_sample = {}
    stopped = False
    for _, row in selected.iterrows():
        sample_id = str(row["Sample_ID"])
        identifiers, aliases = sample_identifiers(row)
        aliases_by_sample[sample_id] = aliases
        found = []
        for kind, value, column in identifiers + [("Sample", sample_id, "Sample_ID")]:
            if stopped:
                break
            try:
                for article_id, relationship, evidence, method in source.record_links(kind, value):
                    found.append({"Article_ID": article_id, "Relationship": relationship,
                                  "Evidence": evidence, "Method": f"{method} ({column})"})
            except (ServiceUnavailable, RequestLimitReached) as problem:
                notes.append(str(problem))
                stopped = isinstance(problem, RequestLimitReached)
        if not found and not stopped:
            for alias in aliases:
                if not source.searchable(alias):
                    continue
                try:
                    for article_id in source.text_search(alias):
                        found.append({"Article_ID": article_id, "Relationship": "search",
                                      "Evidence": f"Text search for '{alias}' returned this article.",
                                      "Method": "Identifier text search (candidate only)"})
                except (ServiceUnavailable, RequestLimitReached) as problem:
                    notes.append(str(problem))
                    stopped = isinstance(problem, RequestLimitReached)
                    break
        seen = set()
        for item in found:
            if item["Article_ID"] in seen:
                continue
            seen.add(item["Article_ID"])
            links.append({"Sample_ID": sample_id, "Sample_Source": row.get("Sample_Source", ""), **item})

    # Retrieve only a few articles, strongest associations first.
    strength = {"direct": 0, "provided": 0, "study": 1, "search": 2}
    ranking = {}
    for link in links:
        best, count = ranking.get(link["Article_ID"], (9, 0))
        ranking[link["Article_ID"]] = (min(best, strength[link["Relationship"]]), count - 1)
    chosen = sorted(ranking, key=lambda article_id: (ranking[article_id], article_id))[:max_articles]
    if len(ranking) > max_articles:
        notes.append(f"{len(ranking) - max_articles} candidate article(s) were not retrieved "
                     f"(limit {max_articles} articles per run).")
    articles, summaries = {}, {}
    try:
        summaries = source.describe(chosen)
        articles = source.fetch(chosen)
    except (ServiceUnavailable, RequestLimitReached) as problem:
        notes.append(f"Article retrieval stopped: {problem}")

    rows = [classify_link(link, aliases_by_sample[link["Sample_ID"]], summaries.get(link["Article_ID"], {}),
                          articles.get(link["Article_ID"]), link["Article_ID"] in chosen,
                          sum(other["Sample_ID"] == link["Sample_ID"] for other in links))
            for link in links]
    for article_id in chosen:
        article = articles.get(article_id)
        if article is not None:
            article.update(summaries.get(article_id, {}))
            article["Article_ID"] = article_id
    table = pd.DataFrame(rows, columns=PUBLICATION_COLUMNS)
    unmatched = [sample_id for sample_id in aliases_by_sample if sample_id not in set(table["Sample_ID"])]
    return table, articles, notes, unmatched


def classify_link(link, aliases, summary, article, retrieved, links_for_sample):
    """One row of the publication-candidate table."""

    in_text = article is not None and any(text_mentions(article, alias) for alias in aliases)
    ambiguity = []
    if link["Relationship"] in {"direct", "provided"}:
        category = DIRECT
        status = ("Verified (NCBI record link)" if link["Relationship"] == "direct"
                  else "Provided in sample metadata (not independently verified)")
    elif link["Relationship"] == "study":
        category = STUDY
        status = "Verified (study record link)"
        ambiguity.append("Study-level link: sample-specific facts need the sample named in the text.")
    elif not retrieved or article is None:
        category = POSSIBLE
        status = "Unverified (text not retrieved)"
    elif in_text:
        category = IN_TEXT
        status = "Verified (exact identifier found in the retrieved text)"
    else:
        category = UNRELATED
        status = "Rejected (identifier not in the retrieved text)"
    if links_for_sample > 1:
        ambiguity.append(f"{links_for_sample} candidate publications for this sample.")
    if not retrieved:
        ambiguity.append("Not retrieved (article limit).")
    elif article is None or not article.get("Passages"):
        ambiguity.append("No article text was available.")
    elif not in_text and category in {DIRECT, STUDY}:
        ambiguity.append("The sample identifier does not appear in the available text.")
    text_status = article.get("Text_Status", "No text available") if article else "Not retrieved"
    return {
        "Sample_ID": link["Sample_ID"], "Sample_Source": link["Sample_Source"],
        "Article_ID": link["Article_ID"], "PMID": summary.get("PMID", ""),
        "PMCID": summary.get("PMCID", ""), "DOI": summary.get("DOI", ""),
        "Title": summary.get("Title", ""), "Publication_Date": summary.get("Publication_Date", ""),
        "Source_URL": summary.get("Source_URL", ""), "Relationship_Type": category,
        "Association_Evidence": link["Evidence"], "Verification_Status": status,
        "Retrieval_Method": link["Method"], "Identifier_In_Text": in_text,
        "Text_Status": text_status, "Ambiguity_Notes": " ".join(ambiguity),
        "Eligible_For_Enrichment": category in ELIGIBLE and bool(article and article.get("Passages")),
    }
