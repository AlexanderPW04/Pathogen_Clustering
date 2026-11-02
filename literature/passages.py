"""Choose the few article passages worth sending to the extraction model.

Methods, sample tables, and passages naming a sample or describing sample
collection are preferred. References, acknowledgements, author
contributions, and affiliations are never sent: they cannot describe where a
sample was collected, and leaving them out keeps requests small.
"""

import hashlib
import re

from literature.identifiers import mentions

SECTION_WEIGHT = {"METHODS": 3.0, "TABLE": 3.0, "SUPPL": 2.0, "RESULTS": 2.0, "ABSTRACT": 2.0,
                  "DISCUSS": 1.0, "CONCL": 1.0, "FIG": 1.0, "INTRO": 0.5, "TITLE": 0.5, "OTHER": 0.5}
EXCLUDED_SECTIONS = {"REF", "ACK_FUND", "AUTH_CONT", "COMP_INT", "ABBR", "REVIEW_INFO",
                     "AFFILIATION", "KEYWORD"}
COLLECTION_WORDS = re.compile(
    r"\b(collected|collection|sampled|sampling|specimens?|samples?|wastewater|sewershed|"
    r"treatment plant|influent|site|county|city|town|province|state|located|obtained)\b",
    re.IGNORECASE)
DATE_WORDS = re.compile(
    r"\b(\d{4}-\d{2}-\d{2}|(January|February|March|April|May|June|July|August|September|October|"
    r"November|December)\b|\d{1,2}/\d{1,2}/\d{2,4})", re.IGNORECASE)
COORDINATES = re.compile(r"-?\d{1,3}\.\d{2,}\s*°?\s*[NS]?\s*,\s*-?\d{1,3}\.\d{2,}")
SURVEILLANCE_WORDS = re.compile(
    r"\b(increase[ds]?|decrease[ds]?|trend|cluster|prevalence|detections?|emerg\w*|limitation)",
    re.IGNORECASE)


def score_passage(passage, aliases):
    """Higher scores mean more likely to describe sample collection."""

    section = passage["Section"]
    if section in EXCLUDED_SECTIONS or "affil" in passage.get("Type", "").casefold():
        return -1.0
    text = passage["Text"]
    score = SECTION_WEIGHT.get(section, 0.5)
    if any(mentions(text, alias) for alias in aliases):
        score += 5
    score += min(len(COLLECTION_WORDS.findall(text)), 3)
    score += 1 if DATE_WORDS.search(text) else 0
    score += 2 if COORDINATES.search(text) else 0
    score += 0.5 if SURVEILLANCE_WORDS.search(text) else 0
    return score


def select_passages(article, aliases, max_passages=6, max_characters=6000, min_score=3.0):
    """Relevant passages of one article, in reading order, without duplicates."""

    candidates = []
    seen = set()
    for index, passage in enumerate(article.get("Passages", [])):
        digest = hashlib.sha256(" ".join(passage["Text"].split()).casefold().encode()).hexdigest()
        if digest in seen:
            continue
        seen.add(digest)
        score = score_passage(passage, aliases)
        if score >= min_score:
            candidates.append((score, index, passage))
    chosen = sorted(candidates, key=lambda item: (-item[0], item[1]))[:max_passages]
    selected, used = [], 0
    for _, index, passage in sorted(chosen, key=lambda item: item[1]):
        text = passage["Text"]
        if used + len(text) > max_characters:
            remaining = max_characters - used
            if remaining < 300:
                break
            text = text[:remaining]  # keep the start; quotations are checked against this text
        selected.append({**passage, "Text": text, "Hash": hashlib.sha256(text.encode()).hexdigest()[:16]})
        used += len(text)
    return selected
