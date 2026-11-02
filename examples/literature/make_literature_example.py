"""Write the SYNTHETIC sequences and metadata for the literature-enrichment demo.

Every sample, sequence, date, place, identifier, and article in this folder is
invented. The articles (synthetic_articles.json), the expected answers
(expected_labels.json), and the mock model responses (mock_responses.json)
are written by hand; this script only produces the FASTA and metadata files.

Sequence design: SYNLIT_01, _09 and _10 differ by one or two bases (a
genetically close group); SYNLIT_03 and _04 form a second close pair; the
other samples differ from both groups at 15 or more of 60 positions.
"""

import random
from pathlib import Path

FOLDER = Path(__file__).resolve().parent
BASE = "ATGGCTAGCTTACGGATCCGTAGCTAGGCTTAACGGATCGATCGTACGGCTAGCTAAGCT"
PATHOGEN = "Synthetic example pathogen"

SAMPLES = [
    # id, date, location, publication, bioproject
    ("SYNLIT_01", "2026-03-03", "USA: Arizona", "SYNPUB-001", ""),
    ("SYNLIT_02", "2026-03-04", "USA: Arizona", "SYNPUB-002", ""),
    ("SYNLIT_03", "2026-03-06", "USA: Arizona", "SYNPUB-003", ""),
    ("SYNLIT_04", "2026-03-07", "USA: Arizona", "SYNPUB-003", ""),
    ("SYNLIT_05", "2026-03-02", "USA: Arizona", "SYNPUB-004", ""),
    ("SYNLIT_06", "2026-03-09", "USA: Arizona", "SYNPUB-005", ""),
    ("SYNLIT_07", "2026-03-10", "USA: Arizona", "", ""),
    ("SYNLIT_08", "2026-03-05", "USA: Arizona", "", "SYNPRJ-08"),
    ("SYNLIT_09", "", "USA: Arizona", "SYNPUB-001", ""),
    ("SYNLIT_10", "", "USA: Arizona", "SYNPUB-009", ""),
    ("SYNLIT_11", "2026-03-02", "USA: Arizona", "SYNPUB-008", ""),
    ("SYNLIT_12", "2026-03-12", "USA: Arizona", "SYNPUB-010", ""),
    ("SYNLIT_13", "2026-03-11", "USA", "", ""),
]


def mutate(sequence, positions):
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    bases = list(sequence)
    for position in positions:
        bases[position] = swap[bases[position]]
    return "".join(bases)


def sequences():
    generator = random.Random(461)
    close_pair = mutate(BASE, [5, 17, 29, 41, 53, 7, 19, 31, 43, 55, 11, 23, 35, 47, 59, 2])
    designed = {
        "SYNLIT_01": BASE,
        "SYNLIT_09": mutate(BASE, [20]),
        "SYNLIT_10": mutate(BASE, [20, 40]),
        "SYNLIT_03": close_pair,
        "SYNLIT_04": mutate(close_pair, [25]),
    }
    for sample_id, *_ in SAMPLES:
        if sample_id not in designed:
            designed[sample_id] = mutate(BASE, sorted(generator.sample(range(60), 18)))
    return designed


def main():
    designed = sequences()
    fasta = "".join(f">{sample_id} synthetic literature-demo sequence\n{designed[sample_id]}\n"
                    for sample_id, *_ in SAMPLES)
    (FOLDER / "synthetic_literature_sequences.fasta").write_text(fasta, encoding="utf-8")
    lines = ["Accession,Pathogen,Gene_Segment,Collection_Date,Location_Name,Latitude,Longitude,"
             "Publication_ID,BioProject,Host_Source"]
    for sample_id, collected, location, publication, project in SAMPLES:
        lines.append(f'{sample_id},{PATHOGEN},Segment 1,{collected},"{location}",,,{publication},{project},'
                     "Synthetic source")
    (FOLDER / "synthetic_literature_metadata.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
