"""Read explicit genomic-region labels without guessing from sequence length."""

GENE_STOP_WORDS = {"complete", "partial", "full", "whole", "coding",
                   "sequence", "genome", "segment"}


def metadata_region(value, column="Gene_Segment"):
    """Turn one metadata region value into a comparable label and clues."""

    text = str(value).strip()
    words = text.split()
    lower_words = [word.casefold() for word in words]

    if len(words) == 1 and words[0].isdigit():
        segment = words[0]
        return {"label": f"Segment {segment}", "segment": segment,
                "genes": set(), "whole": False}

    if len(words) >= 2 and lower_words[0] == "segment":
        segment = words[1].strip(".,;:()[]{}").upper()
        return {"label": f"Segment {segment}", "segment": segment,
                "genes": set(), "whole": False}

    if lower_words in (["complete", "genome"], ["whole", "genome"],
                       ["full", "genome"]):
        return {"label": "Complete genome", "segment": None,
                "genes": set(), "whole": True}

    if column in {"Gene", "Gene_Name"} or lower_words[:1] == ["gene"]:
        gene = words[-1].strip(".,;:()[]{}").upper()
        return {"label": f"Gene {gene}", "segment": None,
                "genes": {gene.casefold()}, "whole": False}

    # A custom explicit region name is kept as the user's own grouping label.
    return {"label": text, "segment": None, "genes": set(), "whole": False}


def fasta_region(description):
    """Read explicit segment, gene, or whole-genome words from a FASTA header."""

    # Biopython's description starts with the accession; it is not annotation.
    parts = str(description).split(maxsplit=1)
    annotation = parts[1] if len(parts) == 2 else ""
    words = annotation.replace("|", " ").replace(",", " ").replace(";", " ").split()
    clean_words = [word.strip(".,;:()[]{}").casefold() for word in words]

    segment = None
    gene = None
    gene_aliases = set()
    whole = False

    for index, word in enumerate(clean_words):
        if word == "segment" and index + 1 < len(words):
            candidate = clean_words[index + 1]
            if candidate.isalnum() and len(candidate) <= 5:
                segment = candidate.upper()

        if word in {"complete", "whole", "full"} and index + 1 < len(words):
            if clean_words[index + 1] == "genome":
                whole = True

        if word == "gene":
            if (index > 0 and clean_words[index - 1].isalpha()
                    and clean_words[index - 1] not in GENE_STOP_WORDS):
                gene = clean_words[index - 1].upper()
                gene_aliases.add(clean_words[index - 1])
            elif (index + 1 < len(words) and clean_words[index + 1].isalpha()
                    and clean_words[index + 1] not in GENE_STOP_WORDS):
                gene = clean_words[index + 1].upper()
                gene_aliases.add(clean_words[index + 1])

        # An explicit name plus acronym, such as "hemagglutinin (HA)".
        if index > 0 and words[index].startswith("(") and words[index].endswith(")"):
            acronym = clean_words[index]
            name = clean_words[index - 1]
            if (2 <= len(acronym) <= 5 and acronym.isalpha()
                    and words[index].strip("()").isupper()
                    and name.isalpha() and len(name) >= 5
                    and name[0] == acronym[0]):
                gene = acronym.upper()
                gene_aliases.update({acronym, name})

    if whole and (segment or gene):
        return {"label": None, "segment": segment, "genes": gene_aliases,
                "whole": True}
    if segment:
        label = f"Segment {segment}"
    elif whole:
        label = "Complete genome"
    elif gene:
        label = f"Gene {gene}"
    else:
        label = None

    return {"label": label, "segment": segment, "genes": gene_aliases,
            "whole": whole}


def regions_conflict(first, second):
    """Find contradictions in annotations that describe the same sample."""

    if first["whole"] and (second["segment"] or second["genes"]):
        return True
    if second["whole"] and (first["segment"] or first["genes"]):
        return True
    if first["segment"] and second["segment"]:
        if first["segment"] != second["segment"]:
            return True
    if first["genes"] and second["genes"]:
        if first["genes"].isdisjoint(second["genes"]):
            return True
    return False
