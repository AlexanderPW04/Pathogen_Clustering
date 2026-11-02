"""Find publications associated with samples and classify each association.

Only traceable links count as associations that can support sample-specific
metadata: a record link from the sample, a publication listed in its
metadata, or the sample's exact identifier inside the article text. A search
hit on its own is never accepted, and the top result is never chosen
automatically.

A link is not yet evidence. An article is sent for analysis about a sample
only when it could document that sample's collection:
- it was not published before the sample was collected (literature/chronology.py);
- a study-level link (shared BioProject or dataset identifier) also needs the
  sample named in the article text, for example in a sample-to-site table.
A general study location is never assigned to every sample of a study.
"""

import pandas as pd

from literature.chronology import chronology
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
    "Publication_Date", "Epub_Date", "Source_URL", "Relationship_Type", "Association_Evidence",
    "Verification_Status", "Retrieval_Method", "Identifier_In_Text", "Text_Status",
    "Ambiguity_Notes", "Retrieval_Status", "Chronology", "Predates_Collection", "Eligibility_Reason",
    "Eligible_For_Enrichment",
]


def dataset_links(source, identifiers, notes):
    """Study-level links from identifiers entered for the whole dataset."""

    links = []
    for kind, value in identifiers:
        try:
            for article_id, _, evidence, method in source.record_links(kind, value):
                links.append({"Article_ID": article_id, "Relationship": "study",
                              "Evidence": f"{kind} {value} was entered for the whole dataset. {evidence}",
                              "Method": "Dataset identifier entered by the user"})
        except (ServiceUnavailable, RequestLimitReached) as problem:
            notes.append(str(problem))
    if identifiers and not links:
        notes.append("The dataset identifier(s) entered did not resolve to a publication.")
    return links


def identifier_key(kind, value):
    """Identifiers that name the same record share a key (an accession with or without its version)."""

    text = str(value).strip()
    if kind == "Nucleotide":
        return (kind, text.split(".")[0].upper())
    return (kind, text)


def batched(source):
    """True for sources that answer many identifiers at once (LiveSource)."""

    return getattr(source, "supports_batches", False) is True


def lookup_links(source, identifiers, notes):
    """{(kind, value): [(article_id, relationship, evidence, method)]} for unique identifiers.

    Sources with a batch method (LiveSource) answer all identifiers in a few
    requests; others are asked once per unique identifier.
    """

    if batched(source):
        try:
            return source.record_links_batch(identifiers)
        except (ServiceUnavailable, RequestLimitReached) as problem:
            notes.append(f"Record links could not be retrieved: {problem}")
            return {}
    found = {}
    for kind, value in identifiers:
        try:
            found[(kind, value)] = source.record_links(kind, value)
        except (ServiceUnavailable, RequestLimitReached) as problem:
            notes.append(str(problem))
            if isinstance(problem, RequestLimitReached):
                break
    return found


def lookup_aliases(source, identifiers, notes):
    """{(kind, value): [other identifiers of the same record]} (BioSample accessions)."""

    try:
        if batched(source):
            return source.record_aliases_batch(identifiers)
        return {(kind, value): source.record_aliases(kind, value) for kind, value in identifiers}
    except (ServiceUnavailable, RequestLimitReached) as problem:
        notes.append(f"Additional sample identifiers could not be looked up: {problem}")
        return {}


def discover_publications(samples, source, max_samples=25, max_articles=5, dataset_identifiers=(),
                          alias_store=None):
    """Return (publication table, retrieved articles, notes, unmatched samples).

    Stage A (discovery): every unique identifier is looked up once, so ten
    accessions of one study cost a few requests, and each sample keeps its own
    association rows. Stage B (retrieval): each unique publication is
    described and fetched once (at most `max_articles`, strongest links
    first); the same text is then assessed separately for every linked sample.

    `dataset_identifiers` are (type, value) pairs entered for the whole
    dataset; they link the study's publication to every selected sample at
    the study level, never as sample-specific evidence. `alias_store`, when
    given, receives each sample's names (including BioSample accessions
    found in NCBI records) for the later passage selection.
    """

    notes = []
    shared = dataset_links(source, dataset_identifiers, notes)
    selected = samples.head(max_samples)
    if len(samples) > max_samples:
        notes.append(f"Only the first {max_samples} selected samples were searched (request limit).")
    sample_rows, unique = [], {}
    for _, row in selected.iterrows():
        sample_id = str(row["Sample_ID"])
        identifiers, aliases = sample_identifiers(row)
        identifiers = identifiers + [("Sample", sample_id, "Sample_ID")]
        keys = []
        for kind, value, column in identifiers:
            key = identifier_key(kind, value)
            keys.append((key, column))
            if key not in unique or ("." in value and "." not in unique[key][1]):
                unique[key] = (kind, value)  # the versioned accession is the one looked up
        sample_rows.append((sample_id, row, keys, aliases))
    found_by_key = lookup_links(source, list(unique.values()), notes)
    stopped = any("limit of" in note for note in notes)

    links, aliases_by_sample, collected, searched = [], {}, {}, {}
    linked_nucleotides = {}
    for sample_id, row, keys, aliases in sample_rows:
        aliases_by_sample[sample_id] = aliases
        collected[sample_id] = row.get("Collection_Date")
        found = []
        for key, column in keys:
            for article_id, relationship, evidence, method in found_by_key.get(unique[key], []):
                found.append({"Article_ID": article_id, "Relationship": relationship,
                              "Evidence": evidence, "Method": f"{method} ({column})"})
        if found:
            linked_nucleotides[sample_id] = [unique[key] for key, _ in keys if key[0] == "Nucleotide"]
        found += [dict(link) for link in shared]
        if not found and not stopped:
            for alias in aliases:
                if not source.searchable(alias):
                    continue
                if alias not in searched:
                    try:
                        searched[alias] = source.text_search(alias)
                    except (ServiceUnavailable, RequestLimitReached) as problem:
                        notes.append(str(problem))
                        stopped = isinstance(problem, RequestLimitReached)
                        break
                found += [{"Article_ID": article_id, "Relationship": "search",
                           "Evidence": f"Text search for '{alias}' returned this article.",
                           "Method": "Identifier text search (candidate only)"} for article_id in searched[alias]]
        seen = set()
        for item in found:
            if item["Article_ID"] not in seen:
                seen.add(item["Article_ID"])
                links.append({"Sample_ID": sample_id, "Sample_Source": row.get("Sample_Source", ""), **item})

    # Other names of linked samples (BioSample accessions), looked up together.
    wanted = list(dict.fromkeys(item for items in linked_nucleotides.values() for item in items))
    extra = lookup_aliases(source, wanted, notes) if wanted else {}
    for sample_id, items in linked_nucleotides.items():
        names = list(aliases_by_sample[sample_id])
        names += [name for item in items for name in extra.get(item, []) if name not in names]
        aliases_by_sample[sample_id] = sorted(names)

    # Retrieve each unique publication once, strongest associations first.
    strength = {"direct": 0, "provided": 0, "study": 1, "search": 2}
    ranking = {}
    for link in links:
        best, count = ranking.get(link["Article_ID"], (9, 0))
        ranking[link["Article_ID"]] = (min(best, strength[link["Relationship"]]), count - 1)
    chosen = sorted(ranking, key=lambda article_id: (ranking[article_id], article_id))[:max_articles]
    if len(ranking) > max_articles:
        notes.append(f"{len(ranking) - max_articles} candidate article(s) were not retrieved "
                     f"(limit {max_articles} articles per run).")
    articles, retrieval = retrieve_articles(source, chosen, notes)
    for article_id in ranking:
        retrieval.setdefault(article_id, f"Not retrieved (limit of {max_articles} articles per search)")

    rows = [classify_link(link, aliases_by_sample[link["Sample_ID"]], articles.get(link["Article_ID"], {}),
                          articles.get(link["Article_ID"]) if articles.get(link["Article_ID"], {}).get("Text_Retrieved")
                          else None, link["Article_ID"] in chosen,
                          sum(other["Sample_ID"] == link["Sample_ID"] for other in links),
                          collected.get(link["Sample_ID"]), retrieval[link["Article_ID"]])
            for link in links]
    if alias_store is not None:
        alias_store.update(aliases_by_sample)
    table = pd.DataFrame(rows, columns=PUBLICATION_COLUMNS)
    unmatched = [sample_id for sample_id in aliases_by_sample if sample_id not in set(table["Sample_ID"])]
    return table, articles, notes, unmatched


def retrieve_articles(source, article_ids, notes):
    """({article_id: article}, {article_id: retrieval status}) with one request per batch.

    Publication details (title, dates) and text are requested separately, so
    a failure of one does not hide the other. Every requested article gets an
    entry, even without text, so its dates and status can be reported.
    """

    summaries, texts, retrieval = {}, {}, {}
    summary_problem = text_problem = None
    if article_ids:
        try:
            summaries = source.describe(article_ids)
        except (ServiceUnavailable, RequestLimitReached) as problem:
            summary_problem = str(problem)
            notes.append(f"Publication details could not be retrieved: {problem}")
        try:
            texts = source.fetch(article_ids)
        except (ServiceUnavailable, RequestLimitReached) as problem:
            text_problem = str(problem)
            notes.append(f"Article text could not be retrieved: {problem}")
    articles = {}
    for article_id in article_ids:
        text = texts.get(article_id)
        article = {"Article_ID": article_id, "Passages": [], **summaries.get(article_id, {})}
        if text is not None:
            article.update(text)
        article["Text_Retrieved"] = text is not None
        if text_problem:
            status = f"Retrieval incomplete: {text_problem}"
        elif text is None:
            status = "Not available from the text service"
        else:
            status = text.get("Text_Status", "Retrieved")
        if summary_problem and article_id not in summaries:
            status += f" Publication details unavailable: {summary_problem}"
        article["Retrieval_Status"] = status
        article["Retrieval_Failed"] = bool(text_problem or (summary_problem and article_id not in summaries))
        retrieval[article_id] = status
        articles[article_id] = article
    return articles, retrieval


def eligibility(category, in_text, has_text, predates, chronology_note, retrieval_status=""):
    """(eligible, reason): can this article document this sample's collection?"""

    reasons = []
    if category not in ELIGIBLE:
        reasons.append("No verified link between this sample and the article.")
    if not has_text:
        reasons.append(f"No article text is available ({retrieval_status})." if retrieval_status
                       else "No article text is available.")
    if predates:
        reasons.append(chronology_note)
    if category == STUDY and not in_text:
        reasons.append("Linked only through the broader study, and no identifier of this sample "
                       "(accession, isolate, BioSample, or other name) appears in the available text.")
    if reasons:
        return False, " ".join(reasons)
    if category == DIRECT:
        return True, "Usable: the sample's own record links this article."
    return True, "Usable: the article text names this sample."


def classify_link(link, aliases, summary, article, retrieved, links_for_sample, collection_date=None,
                  retrieval_status=""):
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
    predates, chronology_note = chronology(
        [summary.get("Publication_Date", ""), summary.get("Epub_Date", "")], collection_date)
    eligible, reason = eligibility(category, in_text, bool(article and article.get("Passages")), predates,
                                   chronology_note, retrieval_status)
    return {
        "Sample_ID": link["Sample_ID"], "Sample_Source": link["Sample_Source"],
        "Article_ID": link["Article_ID"], "PMID": summary.get("PMID", ""),
        "PMCID": summary.get("PMCID", ""), "DOI": summary.get("DOI", ""),
        "Title": summary.get("Title", ""), "Publication_Date": summary.get("Publication_Date", ""),
        "Epub_Date": summary.get("Epub_Date", ""),
        "Source_URL": summary.get("Source_URL", ""), "Relationship_Type": category,
        "Association_Evidence": link["Evidence"], "Verification_Status": status,
        "Retrieval_Method": link["Method"], "Identifier_In_Text": in_text,
        "Text_Status": text_status, "Ambiguity_Notes": " ".join(ambiguity), "Retrieval_Status": retrieval_status,
        "Chronology": chronology_note, "Predates_Collection": bool(predates), "Eligibility_Reason": reason,
        "Eligible_For_Enrichment": eligible,
    }


def article_overview(table, inspected):
    """Counts that keep links, study associations, and usable evidence apart."""

    if table is None or table.empty:
        return {"articles": 0, "direct": 0, "study_only": 0, "named": 0, "usable": 0, "unsuitable": 0,
                "samples_without": inspected}
    groups = list(table.groupby("Article_ID"))
    usable = {article for article, group in groups if group["Eligible_For_Enrichment"].any()}
    return {
        "articles": len(groups),
        "direct": int(table.loc[table["Relationship_Type"] == DIRECT, "Article_ID"].nunique()),
        "study_only": sum(1 for _, group in groups if set(group["Relationship_Type"]) == {STUDY}),
        "named": int(table.loc[table["Identifier_In_Text"].astype(bool), "Article_ID"].nunique()),
        "usable": len(usable),
        "unsuitable": len(groups) - len(usable),
        "samples_without": inspected - table.loc[table["Eligible_For_Enrichment"], "Sample_ID"].nunique(),
    }


def no_evidence_summary(table):
    """What to tell the user when no article can document any selected sample's collection."""

    groups = list(table.groupby("Article_ID", sort=False))
    count = len(groups)
    study_only = all(set(group["Relationship_Type"]) <= {STUDY} for _, group in groups)
    kind = "study-associated" if study_only else "linked"
    if count == 1:
        opening = f"One {kind} publication was found, but it does not"
    else:
        opening = f"{count} {kind} publications were found, but they do not"
    message = (f"{opening} provide verifiable collection information for these samples. No metadata changes "
               "are proposed, nothing is sent to Claude, and nothing is billed. You can continue with the "
               "original metadata.")
    reasons = []
    for article_id, group in groups:
        title = group.iloc[0]["Title"]
        details = (["Linked to these samples, but no passage of the available text describes sample "
                    "collection."] if group["Eligible_For_Enrichment"].any()
                   else list(dict.fromkeys(group["Eligibility_Reason"])))
        if len(group) > 1:  # the article-level summary speaks about all linked samples
            details = [detail.replace("this sample was collected", "these samples were collected")
                       .replace("identifier of this sample", "identifier of these samples") for detail in details]
        links = f"{len(group)} sample link" + ("s" if len(group) != 1 else "")
        reasons.append({"article": article_id, "title": title, "links": links, "reasons": details})
    return {"message": message, "reasons": reasons}
