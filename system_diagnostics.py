"""Small, on-demand checks for the server's external dependencies."""

import subprocess
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO

from align_sequences import find_mafft


def check_mafft():
    """Align two invented sequences to confirm MAFFT actually runs."""

    try:
        mafft_path = find_mafft()
    except RuntimeError:
        return {"status": "Unavailable", "message": "MAFFT was not found on the server."}

    with TemporaryDirectory() as folder:
        input_path = Path(folder) / "diagnostic.fasta"
        input_path.write_text(
            ">DIAG_A\nACGTACGTACGTACGTACGT\n"
            ">DIAG_B\nACGTACGTACGACGTACGT\n",
            encoding="utf-8",
        )
        command = [mafft_path, "--auto", "--nuc", str(input_path)]
        if Path(mafft_path).suffix.lower() in {".bat", ".cmd"}:
            command = ["cmd.exe", "/c"] + command

        try:
            result = subprocess.run(
                command, capture_output=True, text=True, timeout=5
            )
        except subprocess.TimeoutExpired:
            return {"status": "Check failed", "message": "The MAFFT test timed out."}
        except OSError:
            return {"status": "Unavailable", "message": "MAFFT could not be started."}

    if result.returncode != 0:
        return {"status": "Check failed", "message": "MAFFT returned an execution error."}

    try:
        records = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
    except ValueError:
        records = []
    sequences = {record.id: str(record.seq) for record in records}
    if (
        len(records) != 2
        or set(sequences) != {"DIAG_A", "DIAG_B"}
        or not sequences["DIAG_A"]
        or len(sequences["DIAG_A"]) != len(sequences["DIAG_B"])
    ):
        return {
            "status": "Check failed",
            "message": "MAFFT did not return both sequences at equal aligned lengths.",
        }

    version = None
    version_command = [mafft_path, "--version"]
    if Path(mafft_path).suffix.lower() in {".bat", ".cmd"}:
        version_command = ["cmd.exe", "/c"] + version_command
    try:
        version_result = subprocess.run(
            version_command, capture_output=True, text=True, timeout=1
        )
        if version_result.returncode == 0:
            for word in (version_result.stdout + " " + version_result.stderr).split():
                candidate = word.strip(" ,;()")
                if (
                    candidate.startswith("v")
                    and any(character.isdigit() for character in candidate)
                    and all(character in "v0123456789." for character in candidate)
                ):
                    version = candidate
                    break
    except (OSError, subprocess.TimeoutExpired):
        pass  # Alignment succeeded; an unavailable version is not a failure.

    return {
        "status": "Operational",
        "message": "Two synthetic sequences aligned successfully to equal lengths.",
        "version": version,
    }
