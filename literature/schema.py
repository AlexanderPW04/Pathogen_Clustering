"""Structured-output contracts for the two extraction stages, and their checker.

The same JSON Schemas are sent to Claude (structured outputs) and checked
again here in Python, so downstream code never depends on free-form text.
A future locally hosted model must return the same shapes (see
SPRINT4_AI_NLP_INTEGRATION_REPORT.md, "Extraction contract").
"""

PROMPT_VERSION = "literature-extraction-v1"

FIELDS = ["Location", "Coordinates", "Site_Identifier", "Collection_Date", "Collection_Date_Range",
          "Sample_Type", "Collection_Method", "Study_Identifier", "Study_Context"]
LOCATION_ROLES = ["Sample collection site", "General study area", "Sequencing or processing laboratory",
                  "Author affiliation", "Background or cited location", "Not a location"]
GEOGRAPHIC_LEVELS = ["Coordinates", "Facility or site", "City or town", "County",
                     "State or province", "Country", "Not geographic"]
EVIDENCE_TYPES = ["Directly stated", "Sample-to-site table", "Inferred from broader context", "Ambiguous"]
CONTEXT_CATEGORIES = ["Increased detections", "Genetic changes over time", "Multiple collection locations",
                      "Change in sampling prevalence", "Reported cluster", "Population-level trend",
                      "Study limitation"]
CONTEXT_LABELS = ["Explicitly reported in the publication", "Interpretation or hypothesis",
                  "Unsupported or unclear"]
LINKAGE_DECISIONS = ["Specific to this sample", "Applies to the study, not shown to be this sample",
                     "Contradicts the original metadata", "Not supported"]


def strict_object(properties):
    return {"type": "object", "properties": properties, "required": list(properties),
            "additionalProperties": False}


STRING = {"type": "string"}

EXTRACTION_SCHEMA = strict_object({
    "facts": {"type": "array", "items": strict_object({
        "Sample_ID": STRING,
        "Field": {"type": "string", "enum": FIELDS},
        "Value": STRING,
        "Location_Role": {"type": "string", "enum": LOCATION_ROLES},
        "Geographic_Level": {"type": "string", "enum": GEOGRAPHIC_LEVELS},
        "Passage_ID": STRING,
        "Evidence_Quote": STRING,
        "Evidence_Type": {"type": "string", "enum": EVIDENCE_TYPES},
        "Ambiguity": STRING,
        "Reasoning_Summary": STRING,
    })},
    "context_statements": {"type": "array", "items": strict_object({
        "Category": {"type": "string", "enum": CONTEXT_CATEGORIES},
        "Statement": STRING,
        "Label": {"type": "string", "enum": CONTEXT_LABELS},
        "Passage_ID": STRING,
        "Evidence_Quote": STRING,
    })},
    "insufficient_evidence_samples": {"type": "array", "items": STRING},
    "notes": STRING,
})

ASSOCIATION_SCHEMA = strict_object({
    "decisions": {"type": "array", "items": strict_object({
        "Fact_ID": STRING,
        "Linkage": {"type": "string", "enum": LINKAGE_DECISIONS},
        "Explanation": STRING,
    })},
})

SCHEMAS = {"extraction": EXTRACTION_SCHEMA, "association": ASSOCIATION_SCHEMA}

PYTHON_TYPES = {"object": dict, "array": list, "string": str, "boolean": bool,
                "integer": int, "number": (int, float)}


def schema_problems(value, schema, path="$"):
    """List of problems with value against the subset of JSON Schema used above."""

    expected = PYTHON_TYPES[schema["type"]]
    if not isinstance(value, expected) or (schema["type"] in {"integer", "number"} and isinstance(value, bool)):
        return [f"{path}: expected {schema['type']}"]
    problems = []
    if "enum" in schema and value not in schema["enum"]:
        problems.append(f"{path}: '{value}' is not an allowed value")
    if schema["type"] == "object":
        for key in schema.get("required", []):
            if key not in value:
                problems.append(f"{path}: missing '{key}'")
        for key, item in value.items():
            if key not in schema["properties"]:
                problems.append(f"{path}: unexpected '{key}'")
            else:
                problems.extend(schema_problems(item, schema["properties"][key], f"{path}.{key}"))
    if schema["type"] == "array":
        for index, item in enumerate(value):
            problems.extend(schema_problems(item, schema["items"], f"{path}[{index}]"))
    return problems
