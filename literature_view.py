"""Web-facing helpers for the optional literature-enrichment section.

Keeps enrichment sessions in memory for 30 minutes (no database), chooses the
extraction provider, enforces who may start paid requests, runs the analysis
in a background thread so it can be cancelled between requests, and turns a
session into the plain values the template shows.
"""

import hmac
import os
import threading
from pathlib import Path
from secrets import token_urlsafe
from time import monotonic

import pandas as pd

from geocode_locations import lookup_place
from literature.config import ai_settings, ai_status, ncbi_settings
from literature.discovery import DIRECT, IN_TEXT, STUDY
from literature.providers import ClaudeProvider, ExtractionCache, LocalModelProvider, MockProvider
from literature.sources import LiveSource, LocalCorpusSource
from literature.workflow import (apply_session, downloads, estimate_plan, find_literature,
                                 record_decisions, run_analysis)

LITERATURE_FOLDER = Path(__file__).resolve().parent / "examples" / "literature"
LITERATURE_EXAMPLE_FILES = {"fasta": "synthetic_literature_sequences.fasta",
                            "metadata": "synthetic_literature_metadata.csv"}
SYNTHETIC_CORPUS = LITERATURE_FOLDER / "synthetic_articles.json"
MOCK_RESPONSES = LITERATURE_FOLDER / "mock_responses.json"
SESSIONS = {}
SESSION_LOCK = threading.Lock()
SESSION_SECONDS = 1800
SESSION_LIMIT = 4
LOCAL_ADDRESSES = {"127.0.0.1", "::1"}
SHARED_DATA = ("sample IDs and their alternative names, each sample's original location label and "
               "collection date, and short passages from the selected articles. Sequences, variant "
               "tables, coordinates, and other metadata columns are not sent.")


def store_session(session):
    with SESSION_LOCK:
        now = monotonic()
        for key in [key for key, item in SESSIONS.items() if now - item["stored"] > SESSION_SECONDS]:
            del SESSIONS[key]
        while len(SESSIONS) >= SESSION_LIMIT:
            del SESSIONS[min(SESSIONS, key=lambda key: SESSIONS[key]["stored"])]
        token = token_urlsafe(16)
        session["stored"] = now
        SESSIONS[token] = session
    return token


def get_session(token):
    with SESSION_LOCK:
        session = SESSIONS.get(token or "")
        if session is None or monotonic() - session["stored"] > SESSION_SECONDS:
            return None
        return session


def run_settings(session):
    """Limits for this session; the free synthetic demo may analyze every invented article."""

    settings = ai_settings()
    if session["synthetic"]:
        settings.update(max_articles=12, max_requests=30, max_input_tokens=200000)
    return settings


def literature_source(session):
    if session["synthetic"]:
        return LocalCorpusSource(session.get("corpus") or SYNTHETIC_CORPUS)
    return LiveSource()


def choose_provider(session, settings):
    """(provider, explanation) — provider is None when analysis is not possible."""

    if session["synthetic"]:
        return MockProvider(session.get("mock") or MOCK_RESPONSES), ("Offline demonstration: fixture responses written for the "
                                              "synthetic articles. No AI model is called and nothing is billed.")
    if not settings["enabled"]:
        return None, ("AI analysis is disabled on this server. An administrator can enable it by setting "
                      "AI_ENRICHMENT_ENABLED=true and ANTHROPIC_API_KEY in the server environment "
                      "(never in the browser or in Git).")
    if settings["local_model_url"]:
        return LocalModelProvider(settings["local_model_url"]), "A locally hosted model (LOCAL_MODEL_URL)."
    if not settings["key_configured"]:
        return None, "AI analysis is enabled but ANTHROPIC_API_KEY is not set on this server."
    if settings["prices"] is None:
        return None, (f"No price is known for {settings['model']}; set AI_PRICE_INPUT_PER_MTOK and "
                      "AI_PRICE_OUTPUT_PER_MTOK so the spending limit can be enforced.")
    return ClaudeProvider(settings), f"Claude API, model {settings['model']} (billable)."


def live_request_allowed(remote_address, code):
    """Paid requests: from this computer, or with the server's access code."""

    expected = os.getenv("AI_ACCESS_CODE", "").strip()
    if expected:
        return hmac.compare_digest(str(code or "").encode(), expected.encode())
    return remote_address in LOCAL_ADDRESSES


def find(session, selected_ids):
    settings = run_settings(session)
    limits = ncbi_settings()
    find_literature(session, literature_source(session), selected_ids, limits["max_samples"],
                    settings["max_articles"])


def prepare(session):
    """Provider, explanation, and cost estimate shown before any paid request."""

    settings = run_settings(session)
    provider, explanation = choose_provider(session, settings)
    estimate = estimate_plan(session, provider, settings) if provider and session["plan"] else None
    session["estimate"] = estimate
    return provider, explanation, estimate


def start(session, provider, inline=False):
    """Run stage A/B in a background thread (or inline for tests); cancellable between requests."""

    settings = run_settings(session)
    cancel = threading.Event()
    session["cancel"] = cancel
    session["run"] = {"state": "running", "progress": "Starting", "messages": []}
    cache = ExtractionCache(settings["cache_dir"] or None) if provider.is_live else None

    def work():
        try:
            run_analysis(session, provider, settings, cache, cancel)
        except Exception as problem:  # report, never crash the web process
            session["run"] = {"state": "failed", "progress": "", "messages": [
                f"The analysis stopped unexpectedly ({type(problem).__name__})."]}

    if inline:
        work()
    else:
        threading.Thread(target=work, daemon=True).start()


def review_and_apply(session, choices):
    record_decisions(session, choices)
    return apply_session(session, lookup_place)


def check_claude(remote_address, code):
    """Explicit Claude API check with the free token-counting endpoint (no text is generated)."""

    status = ai_status()
    if status["status"] != "Configured":
        return status
    if not live_request_allowed(remote_address, code):
        return {"status": "Not checked", "message": "The Claude check can only be run on the server "
                                                    "computer or with the access code."}
    import anthropic
    settings = ai_settings()
    try:
        anthropic.Anthropic(timeout=20.0, max_retries=0).messages.count_tokens(
            model=settings["model"], messages=[{"role": "user", "content": "Connection check."}])
    except (anthropic.AuthenticationError, anthropic.PermissionDeniedError):
        return {"status": "Check failed", "message": "The Claude API rejected the configured key."}
    except anthropic.NotFoundError:
        return {"status": "Check failed", "message": f"Model {settings['model']} is not available to this key."}
    except anthropic.APIError:
        return {"status": "Unavailable", "message": "Could not complete the Claude API check."}
    return {"status": "Operational", "message": (f"The Claude API accepted the key and model "
                                                 f"{settings['model']} (token-counting request; no text "
                                                 "generated and no charge).")}


def safe_url(url):
    """Only https links are shown; article metadata could contain anything else."""

    text = str(url or "").strip()
    return text if text.startswith("https://") else ""


def short(text, limit=220):
    text = "" if text is None or (not isinstance(text, str) and pd.isna(text)) else str(text)
    return text if len(text) <= limit else text[:limit - 1] + "…"


def view(token, session, provider=None, explanation="", confirm=False, download_id=None, files=()):
    """Plain values for templates/_literature.html."""

    publications = session["publications"]
    stage = "start"
    stats, publication_rows = [], []
    if publications is not None:
        stage = "found"
        inspected = len(session.get("selected", []))
        linked = set(publications.loc[publications["Eligible_For_Enrichment"], "Sample_ID"])
        stats = [
            ("Samples inspected", inspected),
            ("Publications found", publications["Article_ID"].nunique()),
            ("Direct sample links", int((publications["Relationship_Type"] == DIRECT).sum())),
            ("Study-level links", int((publications["Relationship_Type"] == STUDY).sum())),
            ("Identifier found in text", int((publications["Relationship_Type"] == IN_TEXT).sum())),
            ("Samples without a usable publication", inspected - len(linked)),
        ]
        publication_rows = [{
            "sample": row["Sample_ID"], "article": row["Article_ID"], "title": short(row["Title"], 120),
            "url": safe_url(row["Source_URL"]), "date": row["Publication_Date"], "relationship": row["Relationship_Type"],
            "status": row["Verification_Status"], "text": row["Text_Status"], "notes": row["Ambiguity_Notes"],
            "eligible": bool(row["Eligible_For_Enrichment"]),
        } for _, row in publications.iterrows()]
    run = session["run"]
    if run.get("state") not in {None, "not started"}:
        stage = "running" if run["state"] == "running" else "analyzed"
    if session.get("applied"):
        stage = "applied"
    settings = run_settings(session)
    usage = run.get("usage") or {}
    facts = [{
        "id": fact["Fact_ID"], "sample": fact["Sample_ID"], "field": fact["Field"].replace("_", " "),
        "original": fact.get("Original_Value", ""), "proposed": fact["Value"], "quote": fact["Evidence_Quote"],
        "section": fact.get("Evidence_Section", ""), "passage": fact["Passage_ID"], "article": fact["Article_ID"],
        "url": safe_url(fact.get("Article_URL")), "linkage": fact.get("Sample_Article_Linkage", ""),
        "status": fact["Evidence_Status"], "comparison": fact.get("Metadata_Comparison", ""),
        "approvable": fact["Approvable"], "notes": fact.get("Check_Notes", ""),
        "second_stage": fact.get("Second_Stage", ""), "reasoning": fact.get("Reasoning_Summary", ""),
        "decision": session["decisions"].get(fact["Fact_ID"], "unresolved" if fact["Approvable"] else "reject"),
    } for fact in session["facts"]]
    applied = session.get("applied")
    applied_view = None
    if applied:
        comparison = applied["comparison"]
        changed = []
        if comparison is not None:
            changed = comparison[comparison["Change"] != "Unchanged"].fillna("").astype(str).to_dict("records")
        applied_view = {"summary": applied["summary"], "log": applied["log"], "changed": changed[:40],
                        "changed_total": len(changed),
                        "precision": sorted(set(applied["precision_before"]) | set(applied["precision_after"])),
                        "before": applied["precision_before"], "after": applied["precision_after"]}
    return {
        "token": token, "stage": stage, "synthetic": session["synthetic"], "label": session["label"],
        "source_label": session["source_label"], "samples": list(session["samples"]["Sample_ID"]),
        "sample_sources": dict(zip(session["samples"]["Sample_ID"], session["samples"]["Sample_Source"])),
        "selected": set(session.get("selected") or (session.get("eligible") if session.get("stage") == "pre-analysis"
                                                     else session["samples"]["Sample_ID"])),
        "pre_analysis": session.get("stage") == "pre-analysis",
        "gaps": [row for row in session.get("gaps", []) if row["Gaps"]],
        "gap_total": len(session.get("gaps", [])),
        "dataset_identifiers": session.get("dataset_identifiers", []),
        "stats": stats, "publications": publication_rows, "notes": session["discovery_notes"] + session["plan_notes"],
        "plan": [{"article": item["Article_ID"], "title": short(item["Title"], 120),
                  "samples": ", ".join(sample["Sample_ID"] for sample in item["Samples"]),
                  "passages": len(item["Passages"]),
                  "characters": sum(len(passage["Text"]) for passage in item["Passages"])}
                 for item in session["plan"]],
        "provider_name": provider.name if provider else "", "explanation": explanation,
        "estimate": session.get("estimate"), "confirm": confirm, "settings": settings,
        "shared_data": SHARED_DATA, "live": bool(provider and provider.is_live),
        "run": run, "usage": usage, "facts": facts,
        "approvable_count": sum(fact["approvable"] for fact in facts),
        "insufficient": session["insufficient"], "applied": applied_view,
        "download_id": download_id, "files": sorted(files),
    }


def enrichment_downloads(session):
    return downloads(session)
