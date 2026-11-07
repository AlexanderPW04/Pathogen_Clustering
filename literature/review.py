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

import re

import pandas as pd

from literature.discovery import DIRECT, STUDY
from literature.geography import STATE_ABBREVIATIONS
from literature.identifiers import mentions
from literature.passages import COORDINATES, EXCLUDED_SECTIONS

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


COUNTRIES = ["United States", "USA", "Canada", "Mexico", "Brazil", "Argentina", "Chile", "Peru", "Colombia",
             "United Kingdom", "England", "Scotland", "Ireland", "France", "Germany", "Spain", "Italy",
             "Netherlands", "Belgium", "Poland", "Russia", "China", "Japan", "South Korea", "India", "Pakistan",
             "Bangladesh", "Thailand", "Vietnam", "Indonesia", "Australia", "New Zealand", "South Africa",
             "Nigeria", "Kenya", "Egypt", "Israel", "Iran", "Turkey"]
STATE_NAMES = sorted({name.title() for name in STATE_ABBREVIATIONS.values()}, key=len, reverse=True)
PLACE_NAMES = re.compile(r"\b(" + "|".join(STATE_NAMES + COUNTRIES) + r")\b|\b((?:[A-Z][a-z]+ ){1,2}County)\b")


def counted(number, singular, plural=None):
    word = NUMBER_WORDS.get(number, str(number))
    return f"{word} {singular if number == 1 else (plural or singular + 's')}"


def article_places(article):
    """Geographic names found in an article's text (deterministic, not interpreted).

    Finding a name only means a place is mentioned somewhere; it says nothing
    about where any particular sample was collected. Affiliations, references,
    and acknowledgements are not read, as they never describe collection.
    """

    places = []
    for passage in (article or {}).get("Passages", []):
        if passage["Section"] in EXCLUDED_SECTIONS or "affil" in passage.get("Type", "").casefold():
            continue
        places += [first or second for first, second in PLACE_NAMES.findall(passage["Text"])]
        if COORDINATES.search(passage["Text"]):
            places.append("coordinates")
    return list(dict.fromkeys(places))


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


def plain_explanation(outcome, usable, facts, article_count=1):
    """One sentence a reader needs no background for."""

    if outcome == PREDATES:
        return ("This article was published before these samples were collected and does not identify "
                "them. Its collection locations cannot be applied to your samples.")
    if outcome == STUDY_ONLY:
        return ("This article belongs to the same study, but it does not name your samples, so its "
                "locations cannot be applied to them.")
    if outcome == FAILED:
        return "The article could not be retrieved, so its location information could not be checked."
    if outcome == NOT_RETRIEVED:
        return "This article was not retrieved (search limit), so it could not be checked."
    if outcome == UNAVAILABLE:
        return "The article's text is not available, or no passage describes sample collection."
    if outcome == PENDING:
        return (f"The article is linked to {len(usable)} of your samples and may say where they were "
                "collected. Claude can check this after you approve the cost.")
    if outcome == PROPOSED:
        return f"Verified collection details were found ({facts} proposal(s) to review below)."
    return "The article was checked; it gives no verifiable collection details for your samples."


def article_answers(session, article, outcome, usable):
    """The direct questions: is there location text, does it apply, can metadata improve?"""

    facts = len(approvable_facts(session, article["Article_ID"])) if article else 0
    if outcome in {FAILED, NOT_RETRIEVED}:
        location, places = "Could not check", []
    elif not (article or {}).get("Passages"):
        location, places = "Could not check", []
    else:
        places = article_places(article)
        location = "Yes" if places else "No"
    if outcome == PROPOSED:
        applies, improvement = "Yes", "Yes"
    elif outcome == PENDING:
        applies, improvement = "Uncertain", "Needs Claude review"
    elif outcome in {FAILED, NOT_RETRIEVED}:
        applies, improvement = "Uncertain", "Could not check"
    else:
        applies, improvement = "No", "No"
    return {"location": location, "places": places[:6], "applies": applies, "improvement": improvement,
            "plain": plain_explanation(outcome, usable, facts)}


def article_cards(session):
    """One entry per unique publication: details, retrieval, outcome, and its linked samples."""

    table = session.get("publications")
    if table is None or table.empty:
        return []
    cards = []
    for article_id, group in table.groupby("Article_ID", sort=False):
        first = group.iloc[0]
        article = session["articles"].get(article_id) or {}
        outcome = article_outcome(session, article_id, group)
        usable = list(group.loc[group["Eligible_For_Enrichment"], "Sample_ID"])
        year = re.search(r"\b(19|20)\d{2}\b", str(first["Publication_Date"]))
        reasons = list(dict.fromkeys(group["Eligibility_Reason"]))
        if len(group) > 1:  # one card speaks about all linked samples
            reasons = [reason.replace("this sample was collected", "these samples were collected")
                       .replace("identifier of this sample", "identifier of these samples") for reason in reasons]
        cards.append({
            "article": article_id, "title": first["Title"],
            "url": first["Source_URL"] if str(first["Source_URL"]).startswith("https://") else "",
            "published": first["Publication_Date"], "epub": first.get("Epub_Date", ""),
            "retrieval": article.get("Retrieval_Status", first.get("Retrieval_Status", "")),
            "text": first["Text_Status"], "outcome": outcome, "year": year.group() if year else "Unknown",
            "claude": "Yes" if article_id in session.get("analyzed", {}) else "No",
            **article_answers(session, session["articles"].get(article_id), outcome, usable),
            "direct": int((group["Relationship_Type"] == DIRECT).sum()),
            "study": int((group["Relationship_Type"] == STUDY).sum()),
            "samples": list(group["Sample_ID"]),
            "usable_samples": usable,
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
    claude = claude_state(session, cards, bool(search_failures))
    return {"headline": headline(inspected, cards, samples, search_failures, applied or 0),
            "why": why_lines(session, inspected, cards, claude, applied or 0), "stats": stats, "cards": cards,
            "samples": samples, "failures": failures, "claude": claude,
            "stages": stage_statuses(cards, samples, claude, bool(search_failures)),
            "same_status": len({row["Status"] for row in samples}) == 1}


def claude_state(session, cards, search_failed):
    """(label, explanation, tone): exactly what happened with Claude, never implied."""

    run = session.get("run") or {}
    state = run.get("state")
    pending = [card for card in cards if card["outcome"] == PENDING]
    if state == "running":
        return {"label": "Claude analysis running", "why": run.get("progress", ""), "tone": "running"}
    if state == "completed":
        if session.get("synthetic"):
            return {"label": "Claude analysis completed",
                    "why": "Offline demonstration: fixed answers, no AI model was called.", "tone": "done"}
        return {"label": "Claude analysis completed", "why": "", "tone": "done"}
    if state in {"failed", "stopped early", "cancelled"}:
        return {"label": "Claude analysis failed" if state == "failed" else f"Claude analysis {state}",
                "why": " ".join(run.get("messages", [])[:2]), "tone": "failed"}
    if pending:
        why = "Claude has not been run. It reads only the usable passages, after you approve the cost."
    elif search_failed or any(card["outcome"] in {FAILED, NOT_RETRIEVED} for card in cards):
        why = "Claude was not used because the evidence check could not be completed."
    elif not cards:
        why = "Claude was not used because no publication was found."
    elif len(cards) == 1:
        why = "Claude was not used because this publication cannot provide sample-specific collection information."
    else:
        why = "Claude was not used because these publications cannot provide sample-specific collection information."
    return {"label": "Claude not run", "why": why, "tone": "skipped"}


def stage_statuses(cards, samples, claude, search_failed):
    """What actually ran, in order: (label, 'done' | 'failed' | 'skipped')."""

    retrieved = [card for card in cards if card["outcome"] not in {FAILED, NOT_RETRIEVED}]
    proposed = any(row["Status"].startswith("Improvements") for row in samples)
    undecided = any(card["outcome"] == PENDING for card in cards)
    incomplete = search_failed or len(retrieved) < len(cards)
    stages = [("Literature search completed" if not search_failed else "Literature search failed",
               "failed" if search_failed else "done")]
    if cards:
        stages.append(("Article retrieved" if len(retrieved) == len(cards) else "Article retrieval failed",
                       "done" if len(retrieved) == len(cards) else "failed"))
    else:
        stages.append(("No article to retrieve", "skipped"))
    stages.append(("Location evidence checked" if not incomplete and cards else "Location evidence not checked",
                   "done" if not incomplete and cards else "skipped"))
    stages.append((claude["label"], claude["tone"]))
    if proposed:
        stages.append(("Metadata changes proposed", "done"))
    elif undecided or incomplete:
        stages.append(("No changes proposed yet", "skipped"))
    else:
        stages.append(("No changes supported", "done"))
    return stages


def why_lines(session, inspected, cards, claude, applied):
    """A short explanation under the conclusion; details stay in the expandable section."""

    lines = []
    isolates = distinct_isolates(session["samples"]) if inspected else None
    if isolates:
        lines.append(f"{counted(inspected, 'accession').rstrip('.')} {'represents' if inspected == 1 else 'represent'} "
                     f"{counted(isolates, 'distinct isolate').lower()}.")
    if cards:
        retrieved = sum(card["outcome"] not in {FAILED, NOT_RETRIEVED} for card in cards)
        lines.append(f"{counted(retrieved, 'associated publication')} {'was' if retrieved == 1 else 'were'} retrieved.")
    if len(cards) == 1 and cards[0]["year"] != "Unknown":
        dates = pd.to_datetime(session["samples"].get("Collection_Date"), errors="coerce").dropna()
        collected = (f"; the samples were collected in {dates.min():%B %Y}"
                     if len(dates) and dates.min().to_period("M") == dates.max().to_period("M") else "")
        lines.append(f"The publication is from {cards[0]['year']}{collected}.")
    if cards and not any(card["usable_samples"] for card in cards):
        lines.append("No usable direct sample-specific evidence was established.")
    for card in cards[:1] if len(cards) == 1 else []:
        lines.append(card["plain"])
    if claude["why"]:
        lines.append(claude["why"])
    elif claude["label"] == "Claude analysis completed":
        lines.append("Claude read the usable passages; every proposal below was checked against the article.")
    lines.append(f"{counted(applied, 'change')} applied." if applied else "No metadata was changed.")
    return lines


def headline(inspected, cards, samples, search_failures=(), applied=0):
    """The concise conclusion, true to what was actually retrieved and assessed."""

    searched = f"{inspected} accession record{'' if inspected == 1 else 's'} searched."
    unchanged = f"{counted(applied, 'change')} applied." if applied else "No metadata was changed."
    if not cards and search_failures:
        return (f"The publication search could not be completed, so locations could not be checked "
                f"({search_failures[0]}) This does not mean that no publication exists. {unchanged}")
    if not cards:
        return f"{searched} No publication found. {unchanged}"
    found = f"{counted(len(cards), 'publication')} found"
    incomplete = [card for card in cards if card["outcome"] in {FAILED, NOT_RETRIEVED}]
    if incomplete:
        subject = ("it" if len(cards) == 1 else "all of them" if len(incomplete) == len(cards)
                   else counted(len(incomplete), "publication").lower())
        return (f"{searched} {found}, but {subject} could not be retrieved, so location information could not "
                "be checked. The evidence review is incomplete. This does not mean the literature has no "
                f"evidence. {unchanged}")
    proposed = sum(row["Status"].startswith("Improvements") for row in samples)
    pending = [card for card in cards if card["outcome"] == PENDING]
    if proposed:
        return (f"{searched} {found}. Verified collection details were proposed for "
                f"{counted(proposed, 'sample').lower()}. Nothing changes until you approve them."
                + (f" {unchanged}" if applied else ""))
    if pending:
        return (f"{searched} {found}. {counted(len(pending), 'publication')} may hold collection details for "
                f"your samples; a Claude review is needed to confirm. {unchanged}")
    return f"{searched} {found}. No verified collection locations could be assigned to these samples. {unchanged}"


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
