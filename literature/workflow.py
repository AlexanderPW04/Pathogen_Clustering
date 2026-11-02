"""The enrichment workflow: find literature, analyze it, review, and apply.

A session is a plain dictionary kept in server memory for a short time (no
database). Each step reads what the previous step stored. Nothing here
changes the original metadata or the genomic results.
"""

import json
import time

import pandas as pd

from literature.apply import apply_approved, precision_counts, rescreen
from literature.budget import RunBudget
from literature.discovery import CATEGORY_RANK, discover_publications, no_evidence_summary
from literature.identifiers import column_key, linkable_columns, sample_identifiers
from literature.passages import select_passages
from literature.prompts import (ASSOCIATION_SYSTEM, EXTRACTION_SYSTEM, association_message,
                                extraction_message)
from literature.providers import ProviderError
from literature.schema import ASSOCIATION_SCHEMA, EXTRACTION_SCHEMA, PROMPT_VERSION
from literature.verification import AMBIGUOUS, CONTRADICTORY, verify_fact

INDIVIDUAL = "Individual sequence"
POOLED = "Pooled sample"
DROP_COLUMNS = {"Sequence", "Aligned_Sequence", "Normalized_Accession", "Description", "Record_Number"}
MAX_SAMPLES_PER_REQUEST = 20


def sample_table(individual=None, pooled=None):
    """One metadata row per sample from either source, without sequences."""

    frames = []
    for data, source in [(individual, INDIVIDUAL), (pooled, POOLED)]:
        if data is None or data.empty:
            continue
        frame = data.drop(columns=[column for column in data.columns if column in DROP_COLUMNS]).copy()
        frame["Sample_ID"] = frame["Sample_ID"].astype(str)
        frame["Sample_Source"] = source
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["Sample_ID", "Sample_Source"])
    table = pd.concat(frames, ignore_index=True)
    return table.drop_duplicates("Sample_ID", keep="first").reset_index(drop=True)


def attach_identifiers(samples, raw):
    """Add identifier columns (BioProject, PMID, ...) that validation did not keep.

    raw is the uploaded metadata table as text; rows are matched on the
    sample's accession (individual sequences) or Sample_ID (pooled samples).
    """

    if raw is None or raw.empty or samples.empty:
        return samples
    keys = {column_key(column): column for column in raw.columns}
    key = keys.get("accession") or keys.get("sampleid") or keys.get("sample")
    columns = [column for column in linkable_columns(raw)
               if column != key and column not in samples.columns]
    if key is None or not columns:
        return samples
    extra = raw[[key] + columns].copy()
    extra["_match"] = extra[key].astype(str).str.strip()
    extra = extra.drop(columns=[key]).drop_duplicates("_match")
    match = samples["Accession"] if "Accession" in samples else samples["Sample_ID"]
    merged = samples.assign(_match=match.astype(str).str.strip()).merge(extra, on="_match", how="left")
    return merged.drop(columns=["_match"])


def new_session(individual=None, pooled=None, pairs=None, thresholds=None, synthetic=False, label="",
                individual_metadata=None, pooled_metadata=None):
    individual_rows = sample_table(individual, None)
    pooled_rows = sample_table(None, pooled)
    samples = pd.concat([attach_identifiers(individual_rows, individual_metadata),
                         attach_identifiers(pooled_rows, pooled_metadata)], ignore_index=True)
    samples = samples.drop_duplicates("Sample_ID", keep="first").reset_index(drop=True)
    return {"samples": samples, "pairs": pairs, "thresholds": thresholds,
            "synthetic": synthetic, "label": label, "created": time.time(), "source_label": "",
            "publications": None, "articles": {}, "aliases": {}, "discovery_notes": [], "unmatched": [],
            "plan": [], "plan_notes": [], "no_evidence": None, "facts": [], "insufficient": [],
            "run": {"state": "not started"}, "decisions": {}, "applied": None}


def find_literature(session, source, selected_ids=None, max_samples=25, max_articles=5):
    """Step A: publication candidates, retrieved text, and the analysis plan."""

    samples = session["samples"]
    if selected_ids:
        samples = samples[samples["Sample_ID"].isin(selected_ids)]
    session["selected"] = list(samples["Sample_ID"])
    aliases = {str(row["Sample_ID"]): sample_identifiers(row)[1] for _, row in samples.iterrows()}
    table, articles, notes, unmatched = discover_publications(
        samples, source, max_samples, max_articles, session.get("dataset_identifiers", ()), alias_store=aliases)
    session.update(publications=table, articles=articles, discovery_notes=notes, unmatched=unmatched,
                   source_label=source.label, run={"state": "not started"}, facts=[],
                   insufficient=[], decisions={}, applied=None, estimate=None)
    session["aliases"] = aliases
    build_plan(session, max_articles)
    return session


def original_text(row, column):
    value = row.get(column)
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return ""
    if column == "Collection_Date":
        value = pd.to_datetime(value, errors="coerce")
        return "" if pd.isna(value) else value.date().isoformat()
    return str(value)


def sample_row(session, sample_id):
    rows = session["samples"].loc[session["samples"]["Sample_ID"] == str(sample_id)]
    return rows.iloc[0] if not rows.empty else None


def build_plan(session, max_articles):
    """One extraction request per eligible article, with its selected passages.

    Articles that cannot document any selected sample's collection are never
    sent: no request, no cost. session["no_evidence"] then says why.
    """

    plan, notes = [], []
    session["no_evidence"] = None
    table = session["publications"]
    eligible = table[table["Eligible_For_Enrichment"]] if table is not None and not table.empty else table
    if eligible is None or eligible.empty:
        session["plan"], session["plan_notes"] = [], []
        if table is not None and not table.empty:
            session["no_evidence"] = no_evidence_summary(table)
        return
    for article_id, group in table[~table["Article_ID"].isin(eligible["Article_ID"])].groupby("Article_ID"):
        notes.append(f"{article_id}: not sent for analysis. {group['Eligibility_Reason'].iloc[0]}")
    for article_id, group in eligible.groupby("Article_ID", sort=False):
        article = session["articles"].get(article_id)
        sample_ids = list(dict.fromkeys(group["Sample_ID"]))[:MAX_SAMPLES_PER_REQUEST]
        aliases = sorted({alias for sample_id in sample_ids for alias in session["aliases"].get(sample_id, [])})
        passages = select_passages(article, aliases)
        if not passages:
            notes.append(f"{article_id}: no passage describes sample collection; it is not sent for analysis.")
            continue
        excluded = sorted(set(table.loc[table["Article_ID"] == article_id, "Sample_ID"]) - set(group["Sample_ID"]))
        if excluded:
            notes.append(f"{article_id}: sent only for {len(sample_ids)} sample(s) it can document; "
                         f"{len(excluded)} other linked sample(s) are left out ({', '.join(excluded[:5])}"
                         f"{', …' if len(excluded) > 5 else ''}).")
        samples = []
        for sample_id in sample_ids:
            row = sample_row(session, sample_id)
            samples.append({"Sample_ID": sample_id, "Aliases": session["aliases"].get(sample_id, [sample_id]),
                            "Original_Location": original_text(row, "Location_Name"),
                            "Original_Collection_Date": original_text(row, "Collection_Date")})
        linkage = {sample_id: min(group.loc[group["Sample_ID"] == sample_id, "Relationship_Type"],
                                  key=CATEGORY_RANK.get) for sample_id in sample_ids}
        best = min(linkage.values(), key=CATEGORY_RANK.get)
        plan.append({"Article_ID": article_id, "Title": article.get("Title", ""),
                     "Source_URL": article.get("Source_URL", ""), "Samples": samples, "Linkage": linkage,
                     "Passages": passages, "User": extraction_message(samples, article, passages, best)})
    session["plan"], session["plan_notes"] = plan[:max_articles], notes
    if not plan:
        session["no_evidence"] = no_evidence_summary(table)


def estimate_plan(session, provider, settings):
    """Token counts and the worst-case cost of the planned extraction requests."""

    budget = RunBudget(settings, provider.is_live)
    total_tokens, method = 0, "none"
    for item in session["plan"]:
        tokens, method = provider.count_tokens(EXTRACTION_SYSTEM, item["User"], EXTRACTION_SCHEMA)
        item["Estimated_Input_Tokens"] = tokens
        total_tokens += tokens
    worst = sum(budget.worst_case_cost(item["Estimated_Input_Tokens"]) for item in session["plan"])
    # Association checks (at most one per article) run afterwards only while
    # the same limits still allow them, so they can never exceed the caps.
    return {"requests": len(session["plan"]), "input_tokens": total_tokens, "token_method": method,
            "max_output_tokens": settings["max_output_tokens"], "worst_case_cost": round(worst, 4),
            "cost_limit": settings["max_cost_usd"],
            "fits": worst <= settings["max_cost_usd"] and total_tokens <= settings["max_input_tokens"]
            and len(session["plan"]) <= settings["max_requests"]}


def call_provider(provider, budget, cache, run, task, system, user, schema, article_id, tokens=None):
    """One request within the limits; None when it was not sent or failed."""

    key = cache.key(provider.model, task, system, user, provider.name) if cache else None
    stored = cache.get(key) if cache else None
    if stored is not None:
        budget.record(stored)
        return stored
    if tokens is None:
        tokens, _ = provider.count_tokens(system, user, schema)
    reason = budget.refusal(tokens)
    if reason:
        budget.stopped = reason
        run["messages"].append(f"Stopped before {article_id}: {reason}.")
        return None
    try:
        outcome = provider.run(task, system, user, schema, budget.max_output_tokens, article_id)
    except ProviderError as problem:
        run["messages"].append(f"{article_id}: {problem}")
        run["fatal"] = run.get("fatal") or problem.fatal
        return None
    budget.record(outcome)
    if outcome["error"]:
        run["messages"].append(f"{article_id}: {outcome['error']}")
        return None
    if cache:
        cache.put(key, outcome)
    return outcome


def run_analysis(session, provider, settings, cache=None, cancel=None):
    """Step B: extraction (stage A), verification, and selective association checks (stage B)."""

    budget = RunBudget(settings, provider.is_live)
    run = {"state": "running", "progress": "Starting", "provider": provider.name, "model": provider.model,
           "prompt_version": PROMPT_VERSION, "started": time.strftime("%Y-%m-%d %H:%M:%S"),
           "messages": [], "budget": budget, "fatal": False}
    session["run"] = run
    facts, insufficient = [], []
    sent = {}
    for number, item in enumerate(session["plan"], start=1):
        if cancel is not None and cancel.is_set():
            run["state"] = "cancelled"
            break
        if run["fatal"] or budget.stopped:
            break
        run["progress"] = f"Extraction: article {number} of {len(session['plan'])} ({item['Article_ID']})"
        outcome = call_provider(provider, budget, cache, run, "extraction", EXTRACTION_SYSTEM, item["User"],
                                EXTRACTION_SCHEMA, item["Article_ID"], item.get("Estimated_Input_Tokens"))
        if outcome is None:
            continue
        sent[item["Article_ID"]] = item
        passages = {passage["Passage_ID"]: passage for passage in item["Passages"]}
        request_samples = {sample["Sample_ID"]: {"Aliases": sample["Aliases"],
                                                 "Row": sample_row(session, sample["Sample_ID"])}
                           for sample in item["Samples"]}
        for fact in outcome["data"]["facts"]:
            fact = {**fact, "Fact_ID": f"F{len(facts) + 1}", "Article_ID": item["Article_ID"],
                    "Article_Title": item["Title"], "Article_URL": item["Source_URL"],
                    "Sample_Article_Linkage": item["Linkage"].get(str(fact["Sample_ID"]), "Study-level fact"),
                    "Model_Name": outcome["model"], "Prompt_Version": PROMPT_VERSION, "Provider": provider.name,
                    "Second_Stage": "Not run"}
            checked = verify_fact(fact, passages, request_samples, session["aliases"])
            row = request_samples.get(str(fact["Sample_ID"]), {}).get("Row")
            checked["Original_Value"] = original_value(row, fact["Field"]) if row is not None else ""
            facts.append(checked)
        for sample_id in outcome["data"]["insufficient_evidence_samples"]:
            insufficient.append({"Sample_ID": sample_id, "Article_ID": item["Article_ID"]})

    # Stage B: a second, separate check of the facts that could change metadata.
    if run["state"] == "running" and not run["fatal"]:
        for article_id, item in sent.items():
            if (cancel is not None and cancel.is_set()) or budget.stopped:
                break
            candidates = [fact for fact in facts if fact["Article_ID"] == article_id and fact["Approvable"]]
            if not candidates:
                continue
            run["progress"] = f"Association check: {article_id}"
            outcome = call_provider(provider, budget, cache, run, "association", ASSOCIATION_SYSTEM,
                                    association_message(item["Samples"], candidates, item, item["Passages"]),
                                    ASSOCIATION_SCHEMA, article_id)
            if outcome is not None:
                apply_association(candidates, outcome["data"]["decisions"])
    for fact in facts:
        if fact["Approvable"] and fact["Second_Stage"] == "Not run":
            fact["Second_Stage"] = "Not run (stopped before this check)"
        fact["Review_Status"] = "Pending review" if fact["Approvable"] else "Not approvable"
    session.update(facts=facts, insufficient=insufficient)
    if run["state"] == "running":
        run["state"] = "stopped early" if (run["fatal"] or budget.stopped) else "completed"
    run["progress"] = ""
    run["usage"] = budget.summary()
    return session


def apply_association(candidates, decisions):
    """Stage B can only lower a fact's status; agreement is not independent verification."""

    by_id = {fact["Fact_ID"]: fact for fact in candidates}
    for decision in decisions:
        fact = by_id.get(decision["Fact_ID"])
        if fact is None:
            continue
        fact["Second_Stage"] = decision["Linkage"]
        if decision["Linkage"] == "Specific to this sample":
            continue
        fact["Evidence_Status"] = CONTRADICTORY if decision["Linkage"].startswith("Contradicts") else AMBIGUOUS
        fact["Approvable"] = False
        fact["Check_Notes"] = (f"{fact['Check_Notes']} Second-stage check: {decision['Linkage']} — "
                               f"{decision['Explanation']}").strip()


def original_value(row, field):
    if field in {"Location", "Site_Identifier"}:
        return original_text(row, "Location_Name") or "(missing)"
    if field == "Coordinates":
        latitude, longitude = original_text(row, "Latitude"), original_text(row, "Longitude")
        return f"{latitude}, {longitude}" if latitude and longitude else "(missing)"
    if field in {"Collection_Date", "Collection_Date_Range"}:
        return original_text(row, "Collection_Date") or "(missing)"
    return original_text(row, field) or "(missing)"


def record_decisions(session, choices):
    """Store review choices; only approvable facts can be approved (enforced here)."""

    decisions = {}
    for fact in session["facts"]:
        choice = choices.get(fact["Fact_ID"], "unresolved")
        if choice not in {"approve", "reject", "unresolved"} or (choice == "approve" and not fact["Approvable"]):
            choice = "unresolved" if fact["Approvable"] else "reject"
        decisions[fact["Fact_ID"]] = choice
        fact["Review_Status"] = {"approve": "Approved", "reject": "Rejected",
                                 "unresolved": "Left unresolved"}[choice]
    session["decisions"] = decisions
    return decisions


def apply_session(session, geocoder=None):
    """Step D: derived metadata and before/after screening for approved facts."""

    original = session["samples"]
    derived, log = apply_approved(original, session["facts"], session["decisions"], geocoder)
    applied = {"metadata": derived, "log": log, "comparison": None, "summary": {},
               "precision_before": precision_counts(original), "precision_after": precision_counts(derived)}
    individual = derived["Sample_Source"] == INDIVIDUAL
    if session["pairs"] is not None and session["thresholds"] is not None and individual.any():
        table, summary = rescreen(original[individual], derived[individual], session["pairs"],
                                  session["thresholds"])
        applied["comparison"], applied["summary"] = table, summary
    session["applied"] = applied
    return applied


def downloads(session):
    """CSV and JSON evidence files for the result store: {name: (bytes, mimetype)}."""

    def csv(frame):
        return frame.to_csv(index=False).encode("utf-8"), "text/csv"

    files = {}
    if session["publications"] is not None:
        files["literature_publication_candidates.csv"] = csv(session["publications"])
    passages = [{key: passage[key] for key in ["Article_ID", "Passage_ID", "Section", "Hash", "Text"]}
                for item in session["plan"] for passage in item["Passages"]]
    if passages:
        files["literature_evidence_passages.csv"] = csv(pd.DataFrame(passages))
    if session["facts"]:
        files["literature_ai_proposals.csv"] = csv(pd.DataFrame(session["facts"]))
    applied = session.get("applied")
    if applied:
        files["literature_enriched_metadata.csv"] = csv(applied["metadata"])
        files["literature_application_log.csv"] = csv(pd.DataFrame(applied["log"], columns=[
            "Fact_ID", "Sample_ID", "Field", "Value", "Outcome"]))
        if applied["comparison"] is not None:
            files["literature_before_after_pairs.csv"] = csv(applied["comparison"])
    run = session["run"]
    report = {"prompt_version": PROMPT_VERSION, "literature_source": session["source_label"],
              "synthetic_demonstration": session["synthetic"], "provider": run.get("provider"),
              "model": run.get("model"), "started": run.get("started"), "state": run.get("state"),
              "usage": run.get("usage"), "messages": run.get("messages", []),
              "discovery_notes": session["discovery_notes"], "decisions": session["decisions"],
              "note": "Model judgments are not calibrated probabilities. Approved values came from quoted, "
                      "verified passages and a person's review; they are not confirmed collection sites."}
    files["literature_enrichment_report.json"] = (json.dumps(report, indent=2, default=str).encode("utf-8"),
                                                  "application/json")
    return files
