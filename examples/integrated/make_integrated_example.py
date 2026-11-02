"""Write the SYNTHETIC paired individual + pooled demonstration files.

Everything here is invented for a classroom demonstration: the 120-base
reference, the "Synthetic respiratory virus", the sample names, dates, sites,
read counts, and frequencies. None of it is real surveillance data.

The data are designed so the correct answer is known in advance (see
EXPECTED_CATEGORIES). Run from the project folder:

    python examples/integrated/make_integrated_example.py
"""

from pathlib import Path

import pandas as pd


FOLDER = Path(__file__).resolve().parent

REFERENCE_ID = "SYNREF.1"
OTHER_REFERENCE_ID = "SYNREF_B.1"  # a different invented reference
REFERENCE = (
    "GTCATCGTTGTGCTTCTGTCGGAAGATGAGGGTGCAGTGCATCGGTAGATCATACGTACTCGTACATTGG"
    "AGATGGAGGCGGGCAATGCACTCGGGTAGTGATTACGTCTGTGCGTCTTC"
)

# Planted substitutions (1-based reference positions).
V1 = (15, "T", "C")   # shared: individuals A and B, pooled North and South
V2 = (42, "T", "A")   # shared: individual C, pooled South
V3 = (60, "T", "G")   # individual only: pooled samples had adequate depth here
V4 = (75, "G", "A")   # pooled only: every individual sequence called G here
V5 = (100, "T", "C")  # unresolved: pooled depth is only 8 reads at 98-103
V6 = (30, "G", "T")   # unresolved: no individual sequence is callable at 30

EXPECTED_CATEGORIES = {
    f"{REFERENCE_ID}:15:T>C": "Shared",
    f"{REFERENCE_ID}:42:T>A": "Shared",
    f"{REFERENCE_ID}:60:T>G": "Individual only",
    f"{REFERENCE_ID}:75:G>A": "Pooled only",
    f"{REFERENCE_ID}:100:T>C": "Unresolved: pooled",
    f"{REFERENCE_ID}:30:G>T": "Unresolved: individual",
    f"{OTHER_REFERENCE_ID}:15:T>C": "Incomparable",
    f"{REFERENCE_ID}:110:G>A": "Incomparable",  # wrong Ref base (reference has T)
}

GENES = [("synthetic_gene_A", 1, 40), ("synthetic_gene_B", 41, 90), ("synthetic_gene_C", 91, 120)]


def substitute(sequence, *changes):
    """Apply 1-based (position, ref, alt) substitutions after checking ref."""

    bases = list(sequence)
    for position, ref, alt in changes:
        assert bases[position - 1] == ref, (position, bases[position - 1], ref)
        bases[position - 1] = alt
    return "".join(bases)


def mask(sequence, start, end, code="N"):
    """Replace 1-based positions start..end with an uncallable code."""

    return sequence[: start - 1] + code * (end - start + 1) + sequence[end:]


def individual_sequences():
    """Five invented individual sequences, each testing one situation."""

    sample_a = mask(substitute(REFERENCE, V1, V3), 26, 34)            # N run over 30
    sample_b = mask(substitute(REFERENCE, V1, V5), 30, 30, "R")       # IUPAC code at 30
    sample_c = substitute(REFERENCE, V2)[35:]                         # partial: 36-120
    with_n = mask(REFERENCE, 30, 30)
    insertion = "CC" if "C" not in REFERENCE[49:51] else "TT"
    sample_d = with_n[:50] + insertion + with_n[50:87] + with_n[90:]  # +2 after 50, -3 at 88-90
    unrelated = ("ATTA" * 30)                                         # a different "pathogen"
    return {
        "IND_A": sample_a,
        "IND_B": sample_b,
        "IND_C": sample_c,
        "IND_D": sample_d,
        "IND_E": unrelated,
    }


def individual_metadata():
    rows = [
        ("IND_A", "Synthetic respiratory virus", "2026-03-02", "Tempe, Arizona", 33.4255, -111.9400),
        ("IND_B", "Synthetic respiratory virus", "2026-03-09", "Tempe, Arizona", 33.4300, -111.9300),
        ("IND_C", "Synthetic respiratory virus", "2026-03-20", "Phoenix, Arizona", 33.4484, -112.0740),
        ("IND_D", "Synthetic respiratory virus", "2026-03-23", "USA: AZ", None, None),
        ("IND_E", "Synthetic unrelated virus", "2026-03-25", "Tucson, Arizona", 32.2226, -110.9747),
    ]
    return pd.DataFrame(
        [{"Accession": sample, "Pathogen": pathogen, "Gene_Segment": "Complete genome",
          "Collection_Date": date, "Location_Name": place, "Latitude": latitude,
          "Longitude": longitude, "Host_Source": "Synthetic individual specimen"}
         for sample, pathogen, date, place, latitude, longitude in rows]
    )


# Pooled samples: site, date, sequencing depth, and planted frequencies.
POOLED = [
    ("WW_N1", "North Synthetic Sewershed", "2026-03-01", 140, {V1: 0.10, V4: 0.20, V6: 0.05}),
    ("WW_N2", "North Synthetic Sewershed", "2026-03-08", 180, {V1: 0.25, V4: 0.18, V6: 0.08}),
    ("WW_N3", "North Synthetic Sewershed", "2026-03-15", 220, {V1: 0.45, V4: 0.22}),
    ("WW_N4", "North Synthetic Sewershed", "2026-03-22", 260, {V1: 0.70, V4: 0.15, V6: 0.12}),
    ("WW_S1", "South Synthetic Sewershed", "2026-03-03", 120, {V2: 0.30, V4: 0.25}),
    ("WW_S2", "South Synthetic Sewershed", "2026-03-10", 160, {V1: 0.05, V2: 0.35}),
    ("WW_S3", "South Synthetic Sewershed", "2026-03-17", 200, {V1: 0.10, V2: 0.40, V4: 0.30}),
    ("WW_S4", "South Synthetic Sewershed", "2026-03-24", 240, {V2: 0.50, V4: 0.28}),
]
DROPOUT = (98, 103)       # every pooled sample has only 8 reads here
DROPOUT_DEPTH = 8
S4_NO_COVERAGE = (1, 20)  # WW_S4 has zero reads here: 78% genome coverage


def depth_for(sample_id, position, depth):
    if DROPOUT[0] <= position <= DROPOUT[1]:
        return DROPOUT_DEPTH
    if sample_id == "WW_S4" and S4_NO_COVERAGE[0] <= position <= S4_NO_COVERAGE[1]:
        return 0
    return depth


def gene_at(position):
    for name, start, end in GENES:
        if start <= position <= end:
            return name
    return ""


def pooled_tables():
    snv_rows = []
    count_rows = []
    for sample_id, site, date, depth, frequencies in POOLED:
        for position in range(1, len(REFERENCE) + 1):
            reads = depth_for(sample_id, position, depth)
            counts = {"A": 0, "C": 0, "G": 0, "T": 0}
            counts[REFERENCE[position - 1]] = reads
            for (variant_position, ref, alt), frequency in frequencies.items():
                if variant_position == position and reads:
                    alt_reads = round(frequency * reads)
                    counts[alt] = alt_reads
                    counts[ref] = reads - alt_reads
            count_rows.append({
                "Sample_ID": sample_id, "Reference_ID": REFERENCE_ID, "Position": position,
                "A_Count": counts["A"], "C_Count": counts["C"],
                "G_Count": counts["G"], "T_Count": counts["T"],
            })
        for (position, ref, alt), frequency in sorted(frequencies.items()):
            reads = depth_for(sample_id, position, depth)
            # WW_N1's call at position 30 is given a non-significant p-value, so
            # the paper preset (Fisher p < 0.05) excludes it.
            p_value = 0.2 if (sample_id == "WW_N1" and position == 30) else 0.001
            snv_rows.append({
                "Sample_ID": sample_id, "Reference_Genome": REFERENCE_ID, "Position": position,
                "Ref": ref, "Alt": alt, "Alt_Frequency": frequency, "Read_Depth": reads,
                "Alt_Depth": round(frequency * reads), "Collection_Date": date,
                "Location_Name": site, "Source_Fisher_P": p_value,
                "Source_Gene_Feature": gene_at(position),
            })

    # A record on a different reference, and one whose Ref base disagrees
    # with SYNREF.1 (the reference has T at 110). Neither may be matched.
    snv_rows.append({
        "Sample_ID": "WW_X1", "Reference_Genome": OTHER_REFERENCE_ID, "Position": 15,
        "Ref": "T", "Alt": "C", "Alt_Frequency": 0.40, "Read_Depth": 150, "Alt_Depth": 60,
        "Collection_Date": "2026-03-12", "Location_Name": "North Synthetic Sewershed",
        "Source_Fisher_P": 0.001, "Source_Gene_Feature": "",
    })
    snv_rows.append({
        "Sample_ID": "WW_N3", "Reference_Genome": REFERENCE_ID, "Position": 110,
        "Ref": "G", "Alt": "A", "Alt_Frequency": 0.12, "Read_Depth": 220, "Alt_Depth": 26,
        "Collection_Date": "2026-03-15", "Location_Name": "North Synthetic Sewershed",
        "Source_Fisher_P": 0.001, "Source_Gene_Feature": "synthetic_gene_C",
    })
    return pd.DataFrame(snv_rows), pd.DataFrame(count_rows)


def signatures():
    source = "Invented profiles for the BMI 461 synthetic demonstration"
    rows = []
    for name, variants in [("Synthetic profile 1", [V1, V3, V5]), ("Synthetic profile 2", [V2, V4])]:
        for position, ref, alt in variants:
            rows.append({"Signature_Name": name, "Reference_ID": REFERENCE_ID,
                         "Position": position, "Ref": ref, "Alt": alt, "Source": source})
    return pd.DataFrame(rows)


def write_all(folder=FOLDER):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "synthetic_reference.fasta").write_text(
        f">{REFERENCE_ID} SYNTHETIC invented reference genome (not a real organism)\n{REFERENCE}\n",
        encoding="utf-8",
    )
    fasta_lines = []
    for sample_id, sequence in individual_sequences().items():
        fasta_lines.append(f">{sample_id} SYNTHETIC complete genome\n{sequence}")
    (folder / "synthetic_individual_sequences.fasta").write_text(
        "\n".join(fasta_lines) + "\n", encoding="utf-8"
    )
    individual_metadata().to_csv(folder / "synthetic_individual_metadata.csv", index=False)
    snvs, counts = pooled_tables()
    snvs.to_csv(folder / "synthetic_pooled_snvs.csv", index=False)
    counts.to_csv(folder / "synthetic_pooled_base_counts.csv", index=False)
    signatures().to_csv(folder / "synthetic_variant_signatures.csv", index=False)


if __name__ == "__main__":
    write_all()
    print(f"Synthetic integrated example written to {FOLDER}")
