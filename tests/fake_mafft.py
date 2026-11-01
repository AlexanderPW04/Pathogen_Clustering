"""A small, deterministic stand-in for MAFFT used only in automated tests.

Real MAFFT is not installed in every test environment. This helper aligns
short sequences with a simple Needleman-Wunsch algorithm (free end gaps), so
tests can exercise partial sequences, insertions, and deletions. Tests that
use it are mocked-MAFFT tests; they are not evidence that real MAFFT works.
"""

import subprocess
from io import StringIO
from pathlib import Path

from Bio import SeqIO


def needleman_wunsch(reference, sample, match=2, mismatch=-1, gap=-3):
    """Return (reference_row, sample_row) with gaps; end gaps cost nothing."""

    rows, columns = len(reference), len(sample)
    score = [[0] * (columns + 1) for _ in range(rows + 1)]
    move = [[""] * (columns + 1) for _ in range(rows + 1)]
    for i in range(1, rows + 1):
        move[i][0] = "up"
    for j in range(1, columns + 1):
        move[0][j] = "left"

    for i in range(1, rows + 1):
        for j in range(1, columns + 1):
            a, b = reference[i - 1], sample[j - 1]
            pair = match if a == b else (0 if b not in "ACGT" else mismatch)
            options = [
                (score[i - 1][j - 1] + pair, "diagonal"),
                (score[i - 1][j] + gap, "up"),
                (score[i][j - 1] + gap, "left"),
            ]
            score[i][j], move[i][j] = max(options, key=lambda item: item[0])

    # Free trailing gaps: start from the best cell in the last row or column.
    candidates = [(score[rows][j], rows, j) for j in range(columns + 1)]
    candidates += [(score[i][columns], i, columns) for i in range(rows + 1)]
    _, i, j = max(candidates)
    reference_row = reference[i:]
    sample_row = "-" * (rows - i)
    if j < columns:
        reference_row = reference_row + "-" * (columns - j)
        sample_row = sample_row + sample[j:]

    while i > 0 or j > 0:
        step = move[i][j]
        if step == "diagonal":
            reference_row = reference[i - 1] + reference_row
            sample_row = sample[j - 1] + sample_row
            i, j = i - 1, j - 1
        elif step == "up":
            reference_row = reference[i - 1] + reference_row
            sample_row = "-" + sample_row
            i -= 1
        else:
            reference_row = "-" + reference_row
            sample_row = sample[j - 1] + sample_row
            j -= 1
    return reference_row, sample_row


def star_alignment(reference, sequences):
    """Align several sequences through the reference (a multiple alignment)."""

    pairwise = {name: needleman_wunsch(reference, sequence) for name, sequence in sequences.items()}
    # Longest insertion after each reference position (0 = before position 1).
    longest = [0] * (len(reference) + 1)
    pieces = {}
    for name, (reference_row, sample_row) in pairwise.items():
        bases = [""] * len(reference)
        inserted = [""] * (len(reference) + 1)
        position = 0
        for reference_base, sample_base in zip(reference_row, sample_row):
            if reference_base == "-":
                inserted[position] += sample_base
            else:
                bases[position] = sample_base
                position += 1
        for index, text in enumerate(inserted):
            longest[index] = max(longest[index], len(text))
        pieces[name] = (bases, inserted)

    aligned = {}
    for name, (bases, inserted) in pieces.items():
        text = inserted[0].ljust(longest[0], "-")
        for index, base in enumerate(bases):
            text += base + inserted[index + 1].ljust(longest[index + 1], "-")
        aligned[name] = text
    return aligned


def make_fake_mafft(reference):
    """Build a replacement for subprocess.run that behaves like `mafft --auto`."""

    def fake_run(command, **kwargs):
        text = Path(command[-1]).read_text(encoding="utf-8")
        records = {str(record.id): str(record.seq).upper()
                   for record in SeqIO.parse(StringIO(text), "fasta")}
        if "REFERENCE" in records:
            reference_row, sample_row = needleman_wunsch(records["REFERENCE"], records["SAMPLE"])
            output = f">REFERENCE\n{reference_row}\n>SAMPLE\n{sample_row}\n"
        else:
            aligned = star_alignment(reference, records)
            output = "".join(f">{name}\n{row.lower()}\n" for name, row in aligned.items())
        return subprocess.CompletedProcess(command, 0, stdout=output)

    return fake_run
