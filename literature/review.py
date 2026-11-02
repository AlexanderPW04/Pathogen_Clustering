"""What the literature review found, per publication and per sample.

The workflow has five stages:
A. discovery: which publications are linked to which records, and how;
B. retrieval: each unique publication's details and text, fetched once;
C. evidence assessment (free): can the article document each linked
   sample's collection (dates, sample named in the text, link strength)?
D. Claude extraction (paid, only after cost approval, in batches);
E. verification and review: checked proposals, a person's decisions, and
   the values applied.

This module only summarizes the session. Each publication gets one outcome
category, each sample one status, so "nothing could be improved" is a
completed result, and a failed retrieval is reported as an incomplete review,
never as an absence of evidence.
"""

import pandas as pd

from literature.discovery import DIRECT, STUDY

PROPOSED = "Metadata improvements proposed"
REVIEWED_NONE = "Article reviewed; no supported improvements"
STUDY_ONLY = "Article associated with study only; no sample-specific evidence"
PREDATES = "Article predates the sample collection"
UNAVAILABLE = "Article retrieved, but relevant information was unavailable"
FAILED = "Article retrieval incomplete due to API failure"
NOT_RETRIEVED = "Not retrieved (article limit for one search)"
PENDING = "Eligible evidence; Claude extraction not run yet"
NO_PUBLICATION = "No associated publication found"
REVIEWED = {PROPOSED, REVIEWED_NONE, STUDY_ONLY, PREDATES, UNAVAILABLE}
NUMBER_WORDS = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven", 8: "Eight",
                9: "Nine", 10: "Ten"}


def counted(number, singular, plural=None):
    word = NUMBER_WORDS.get(number, str(number))
    return f"{word} {singular if number == 1 else (plural or singular + 's')}"


def approvable_facts(session, article_id=None, sample_id=None):
    return [fact for fact in session.get("facts", []) if fact["Approvable"]
            and (article_id is None or fact["Article_ID"] == article_id)
            and (sample_id is None or str(fact["Sample_ID"]) == str(sample_id))]


def article_outcome(session, article_id, group):
    """The outcome category of one publication, from its rows in the candidate table."""

    article = session["articles"].get(article_id)
    if article is None:
        return NOT_RETRIEVED
    if article.get("Retrieval_Failed"):
        return FAILED
    analyzed = session.get("analyzed", {})
    if group["Eligible_For_Enrichment"].any():
        planned = {item["Article_ID"] for item in session.get("plan", [])}
        if article_id not in planned:
            return UNAVAILABLE  # no passage describes sample collection
        if article_id in analyzed:
            return PROPOSED if approvable_facts(session, article_id) else REVIEWED_NONE
        return PENDING
    if group["Predates_Collection"].all():
        return PREDATES
    if not article.get("Passages"):
        return UNAVAILABLE
    if set(group["Relationship_Type"]) <= {STUDY}:
        return STUDY_ONLY
    return REVIEWED_NONE


def article_cards(session):
    """One entry per unique publication: details, retrieval, outcome, and its linked samples."""

    table = session.get("publications")
    if table is None or table.empty:
        return []
    cards = []
    for article_id, group in table.groupby("Article_ID", sort=False):
        first = group.iloc[0]
        article = session["articles"].get(article_id) or {}
        reasons = list(dict.fromkeys(group["Eligibility_Reason"]))
        if len(group) > 1:  # one card speaks about all linked samples
            reasons = [reason.replace("this sample was collected", "these samples were collected")
                       .replace("identifier of this sample", "identifier of these samples") for reason in reasons]
        cards.append({
            "article": article_id, "title": first["Title"],
            "url": first["Source_URL"] if str(first["Source_URL"]).startswith("https://") else "",
            "published": first["Publication_Date"], "epub": first.get("Epub_Date", ""),
            "retrieval": article.get("Retrieval_Status", first.get("Retrieval_Status", "")),
            "text": first["Text_Status"], "outcome": article_outcome(session, article_id, group),
            "direct": int((group["Relationship_Type"] == DIRECT).sum()),
            "study": int((group["Relationship_Type"] == STUDY).sum()),
            "samples": list(group["Sample_ID"]),
            "usable_samples": list(group.loc[group["Eligible_For_Enrichment"], "Sample_ID"]),
            "chronology": list(dict.fromkeys(group["Chronology"])),
            "reasons": reasons, "evidence": list(dict.fromkeys(group["Association_Evidence"]))[:3],
        })
    return cards


def distinct_isolates(samples):
    """Number of distinct isolate (or strain) names, or None when no such column exists."""

    for column in samples.columns:
        if str(column).casefold() in {"isolate", "strain", "isolate_name"}:
            values = samples[column].dropna().astype(str).str.strip()
            values = values[values != ""]
            return int(values.nunique()) if len(values) else None
    return None


def sample_statuses(session, cards, search_failed=False):
    """One compact row per inspected sample."""

    table = session.get("publications")
    outcome_of = {card["article"]: card["outcome"] for card in cards}
    rows = []
    for sample_id in session.get("selected") or list(session["samples"]["Sample_ID"]):
        group = table[table["Sample_ID"] == sample_id] if table is not None and not table.empty else pd.DataFrame()
        facts = approvable_facts(session, sample_id=sample_id)
        decisions = session.get("decisions", {})
        if group.empty:
            status = "Search incomplete (a service was unavailable)" if search_failed else NO_PUBLICATION
        elif facts:
            approved = sum(decisions.get(fact["Fact_ID"]) == "approve" for fact in facts)
            status = (f"Improvements proposed ({len(facts)}; {approved} approved)" if decisions
                      else f"Improvements proposed ({len(facts)}), awaiting your review")
        else:
            outcomes = [outcome_of.get(article) for article in group["Article_ID"]]
            usable = group["Eligible_For_Enrichment"].any()
            if FAILED in outcomes or NOT_RETRIEVED in outcomes:
                status = "Review incomplete: a linked publication was not retrieved"
            elif usable and PENDING in outcomes:
                status = "Usable evidence found; Claude extraction not run yet"
            elif usable:
                status = "Reviewed; no supported improvement"
            elif all(outcome == PREDATES for outcome in outcomes):
                status = "Linked publication predates this sample's collection"
            elif all(outcome in {STUDY_ONLY, PREDATES} for outcome in outcomes):
                status = "Study-level link only; no sample-specific evidence"
            else:
                status = "Reviewed; no supported improvement"
        rows.append({"Sample_ID": sample_id, "Publications": ", ".join(group["Article_ID"]) if not group.empty else "—",
                     "Link": ", ".join(sorted(set(group["Relationship_Type"]))) if not group.empty else "—",
                     "Status": status})
    return rows


def review_summary(session):
    """Counts, a headline sentence, publication cards, and per-sample statuses for the page."""

    table = session.get("publications")
    cards = article_cards(session)
    search_failures = [note for note in session.get("discovery_notes", [])
                       if "could not" in note or "limit of" in note or "HTTP" in note]
    samples = sample_statuses(session, cards, bool(search_failures)) if table is not None else []
    inspected = len(session.get("selected") or [])
    outcomes = [card["outcome"] for card in cards]
    failures = [f"{card['article']}: {card['retrieval']}" for card in cards if card["outcome"] == FAILED]
    failures += [note for note in search_failures if not any(note in item for item in failures)]
    applied = session.get("applied_count")
    if applied is None and session.get("applied"):
        applied = sum(item["Outcome"].startswith("Applied") for item in session["applied"]["log"])
    stats = [
        ("Accession records inspected", inspected),
        ("Distinct isolates", distinct_isolates(session["samples"]) if inspected else None),
        ("Unique publications found", len(cards)),
        ("Publications with a direct sample link", sum(card["direct"] > 0 for card in cards)),
        ("Publications with study-level links only", sum(card["direct"] == 0 and card["study"] > 0 for card in cards)),
        ("Direct sample links", sum(card["direct"] for card in cards)),
        ("Study-level links", sum(card["study"] for card in cards)),
        ("Publications retrieved", sum(card["outcome"] not in {FAILED, NOT_RETRIEVED} for card in cards)),
        ("Publications reviewed", sum(outcome in REVIEWED for outcome in outcomes)),
        ("Publications with usable sample-specific evidence", sum(bool(card["usable_samples"]) for card in cards)),
        ("Samples with proposed improvements", sum(row["Status"].startswith("Improvements") for row in samples)),
        ("Samples without a supported improvement", sum(not row["Status"].startswith("Improvements")
                                                        for row in samples)),
        ("Retrieval failures", len(failures)),
        ("Changes approved", sum(choice == "approve" for choice in session.get("decisions", {}).values())),
        ("Changes applied", applied if applied is not None else 0),
    ]
    stats = [(label, "Not determinable" if value is None else value) for label, value in stats]
    return {"headline": headline(inspected, cards, samples, search_failures), "stats": stats, "cards": cards,
            "samples": samples, "failures": failures}


def headline(inspected, cards, samples, search_failures=()):
    """One plain sentence about the review, true to what was actually retrieved and assessed."""

    if not cards and search_failures:
        return (f"The publication search could not be completed ({search_failures[0]}). The evidence review is "
                "incomplete; this does not mean that no publication exists. The original metadata remains "
                "unchanged, and the genomic analysis can still run.")
    if not cards:
        return (f"{counted(inspected, 'accession record')} {'was' if inspected == 1 else 'were'} searched; "
                "no associated publication was found. The original metadata remains unchanged.")
    kinds = "study " if all(card["direct"] == 0 for card in cards) else ""
    linked = (f"{counted(inspected, 'accession record')} {'was' if inspected == 1 else 'were'} associated "
              f"with {counted(len(cards), kinds + 'publication').lower()}.")
    incomplete = [card for card in cards if card["outcome"] in {FAILED, NOT_RETRIEVED}]
    if incomplete:
        return (f"{linked} The evidence review is incomplete: {counted(len(incomplete), 'publication').lower()} "
                f"could not be retrieved ({incomplete[0]['retrieval']}). This does not mean the literature has "
                "no evidence. The original metadata remains unchanged.")
    proposed = sum(row["Status"].startswith("Improvements") for row in samples)
    pending = [card for card in cards if card["outcome"] == PENDING]
    once = ("The publication was retrieved and reviewed once." if len(cards) == 1
            else "Each publication was retrieved and reviewed once.")
    if proposed:
        return (f"{linked} {once} Metadata improvements were proposed for {counted(proposed, 'sample').lower()}; "
                "nothing changes until you approve them.")
    if pending:
        return (f"{linked} {once} {counted(len(pending), 'publication')} {'contains' if len(pending) == 1 else 'contain'} "
                "sample-specific evidence that Claude can analyze after you approve the cost. "
                "The original metadata is unchanged so far.")
    return (f"{linked} {once} No verifiable collection-site improvements were found for these records. "
            "The original metadata remains unchanged.")


def enrichment_summary(session, log=None):
    """Per-sample rows: original, proposed, approval, applied value, and supporting publication.

    Samples without any approvable proposal get one row saying their original
    values were kept, so the summary lists every inspected sample.
    """

    log = log if log is not None else ((session.get("applied") or {}).get("log") or [])
    applied = {item["Fact_ID"]: item for item in log}
    decisions = session.get("decisions", {})
    rows = []
    for sample_id in session.get("selected") or list(session["samples"]["Sample_ID"]):
        facts = [fact for fact in session.get("facts", []) if str(fact["Sample_ID"]) == str(sample_id)]
        if not facts:
            rows.append({"Sample": sample_id, "Field": "—", "Original": "unchanged", "Proposed": "—",
                         "Approval": "No proposal", "Applied": "Original kept", "Publication": "—"})
            continue
        for fact in facts:
            decision = decisions.get(fact["Fact_ID"])
            if not fact["Approvable"]:
                approval = f"Not approvable ({fact['Evidence_Status']}); needs no action"
            else:
                approval = {"approve": "Approved", "reject": "Rejected"}.get(decision, "Awaiting review")
            item = applied.get(fact["Fact_ID"])
            applied_value = (item["Value"] if item and item["Outcome"].startswith("Applied")
                             else item["Outcome"] if item else "Original kept")
            rows.append({"Sample": sample_id, "Field": fact["Field"].replace("_", " "),
                         "Original": fact.get("Original_Value", ""), "Proposed": fact["Value"],
                         "Approval": approval, "Applied": applied_value, "Publication": fact["Article_ID"]})
    return rows
