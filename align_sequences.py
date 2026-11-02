"""Align each pathogen and gene/segment group with the MAFFT program."""

import shutil
import subprocess
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO

from validate_inputs import ALLOWED_SEQUENCE_CHARACTERS, clean_sequence


# The hosted app's worker is stopped after 180 s, so MAFFT must stop first.
MAFFT_TIMEOUT_SECONDS = 150


def find_mafft():
    """Find the external MAFFT executable before analysis starts."""

    mafft_path = shutil.which("mafft")
    if mafft_path is None:
        raise RuntimeError(
            "MAFFT was not found. Install MAFFT and add it to PATH, "
            "then run the analysis again."
        )
    return mafft_path


def run_mafft(mafft_path, input_path):
    """Run MAFFT on one FASTA file and return the aligned FASTA text."""

    input_path = Path(input_path)
    command = [mafft_path, "--auto", "--nuc", str(input_path)]
    folder = None

    # The official Windows package may expose MAFFT as a batch file. cmd.exe
    # drops the outer quotes when a command line has more than two, which
    # breaks paths with spaces, so the input is passed by its plain file name
    # from its own folder and only the MAFFT path can need quotes.
    if Path(mafft_path).suffix.lower() in {".bat", ".cmd"}:
        command = ["cmd.exe", "/c", mafft_path, "--auto", "--nuc", input_path.name]
        folder = str(input_path.parent)

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
            cwd=folder,
            timeout=MAFFT_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"MAFFT did not finish within {MAFFT_TIMEOUT_SECONDS} seconds. Use fewer or "
            "shorter sequences on the website, or run_pipeline.py on your own computer."
        ) from error
    return result.stdout


def align_samples(data, mafft_path=None):
    """Return a copy with MAFFT-aligned sequences in each comparable group."""

    if mafft_path is None:
        mafft_path = find_mafft()

    aligned = data.copy()
    groups = data.groupby(["Pathogen", "Gene_Segment"], sort=True)

    for group_names, group in groups:
        pathogen, gene_segment = group_names

        # A single sample has no other sequence to align against.
        if len(group) < 2:
            continue

        with TemporaryDirectory() as folder:
            input_path = Path(folder) / "group.fasta"

            with input_path.open("w", encoding="utf-8") as fasta_file:
                for _, sample in group.iterrows():
                    fasta_file.write(
                        f">{sample['Sample_ID']}\n{sample['Sequence']}\n"
                    )

            try:
                aligned_text = run_mafft(mafft_path, input_path)
            except (OSError, subprocess.CalledProcessError) as error:
                details = getattr(error, "stderr", "") or str(error)
                raise RuntimeError(
                    f"MAFFT failed for {pathogen} / {gene_segment}: "
                    f"{details.strip()}"
                ) from error

        records = list(SeqIO.parse(StringIO(aligned_text), "fasta"))
        aligned_by_id = {}
        for record in records:
            sample_id = str(record.id)
            if sample_id in aligned_by_id:
                raise RuntimeError(
                    f"MAFFT returned a duplicate sample ID: {sample_id}"
                )
            sequence = clean_sequence(record.seq)
            if set(sequence) - ALLOWED_SEQUENCE_CHARACTERS:
                raise RuntimeError(
                    f"MAFFT returned invalid sequence characters for {sample_id}."
                )
            aligned_by_id[sample_id] = sequence

        expected_ids = set(group["Sample_ID"].astype(str))
        if set(aligned_by_id) != expected_ids:
            raise RuntimeError(
                f"MAFFT changed or omitted sample IDs for "
                f"{pathogen} / {gene_segment}."
            )

        aligned_lengths = {len(sequence) for sequence in aligned_by_id.values()}
        if len(aligned_lengths) != 1:
            raise RuntimeError(
                f"MAFFT returned unequal aligned lengths for "
                f"{pathogen} / {gene_segment}."
            )

        for index, sample in group.iterrows():
            sample_id = str(sample["Sample_ID"])
            aligned.at[index, "Sequence"] = aligned_by_id[sample_id]

    return aligned
