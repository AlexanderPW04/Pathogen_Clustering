"""Evaluate literature enrichment on the versioned synthetic labeled set.

Default (free, offline): the mock provider, which measures the pipeline's
deterministic checks on planted model mistakes, not any model's accuracy.

    python validation/evaluate_literature_enrichment.py

Live Claude evaluation on the same invented articles (billable, opt-in):

    set AI_ENRICHMENT_ENABLED=true and ANTHROPIC_API_KEY in the environment, then
    python validation/evaluate_literature_enrichment.py --provider claude --allow-paid

A future locally hosted model can be compared the same way with
--provider local (LOCAL_MODEL_URL). Results on synthetic articles are not
evidence of real-world accuracy.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomics.file_formats import read_table  # noqa: E402
from literature.config import ai_settings, ai_status  # noqa: E402
from literature.evaluation import evaluate  # noqa: E402
from literature.providers import ClaudeProvider, ExtractionCache, LocalModelProvider, MockProvider  # noqa: E402
from literature.sources import LocalCorpusSource  # noqa: E402
from literature.workflow import find_literature, new_session, run_analysis  # noqa: E402

EXAMPLE = ROOT / "examples" / "literature"


def synthetic_samples():
    """The synthetic metadata in the shape of the validated Sprint 2 sample table."""

    raw = read_table(EXAMPLE / "synthetic_literature_metadata.csv")
    samples = pd.DataFrame({
        "Sample_ID": raw["Accession"], "Accession": raw["Accession"],
        "Collection_Date": pd.to_datetime(raw["Collection_Date"].replace("", None), errors="coerce"),
        "Location_Name": raw["Location_Name"], "Latitude": pd.NA, "Longitude": pd.NA,
        "Coordinate_Source": pd.NA,
        "Geographic_Precision": "State-level location (too coarse for distance cutoff)",
    })
    return samples, raw


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--provider", choices=["mock", "claude", "local"], default="mock")
    parser.add_argument("--allow-paid", action="store_true", help="Required for --provider claude.")
    parser.add_argument("--output", type=Path, help="Folder for per-sample results and proposals (CSV).")
    args = parser.parse_args()

    settings = {**ai_settings(), "max_articles": 12, "max_requests": 24}
    cache = None
    if args.provider == "mock":
        provider = MockProvider(EXAMPLE / "mock_responses.json")
    elif args.provider == "claude":
        status = ai_status()
        if not args.allow_paid or status["status"] != "Configured":
            print("Live Claude evaluation needs --allow-paid and a configured server environment:")
            print(f"  {status['status']}: {status['message']}")
            return 2
        provider = ClaudeProvider(settings)
        cache = ExtractionCache(settings["cache_dir"] or None)
    else:
        if not settings["local_model_url"]:
            print("Set LOCAL_MODEL_URL to evaluate a locally hosted model.")
            return 2
        provider = LocalModelProvider(settings["local_model_url"])

    samples, raw = synthetic_samples()
    session = new_session(individual=samples, synthetic=True, individual_metadata=raw)
    find_literature(session, LocalCorpusSource(EXAMPLE / "synthetic_articles.json"), max_articles=12)
    run_analysis(session, provider, settings, cache)
    metrics, table = evaluate(session, EXAMPLE / "expected_labels.json")

    print(f"Provider: {provider.name} ({provider.model})")
    print(f"Run: {session['run']['state']}; usage: {session['run']['usage']}")
    for message in session["run"]["messages"]:
        print(f"  note: {message}")
    for label, value in metrics.items():
        print(f"{label}: {value}")
    if args.provider == "mock":
        print("Mock responses were written by hand with planted errors; these numbers test the checks only.")
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.output / f"evaluation_{args.provider}.csv", index=False)
        pd.DataFrame(session["facts"]).to_csv(args.output / f"proposals_{args.provider}.csv", index=False)
        print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
