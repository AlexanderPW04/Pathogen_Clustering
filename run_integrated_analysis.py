"""Run the integrated individual + pooled analysis from the command line.

Use this for processed datasets that are too large for the hosted website
(for example a full per-position A/C/G/T table for the paper's PCoA). It runs
the same code as the website and writes every table (CSV) and figure (PNG)
into one output folder.

Examples:

    python run_integrated_analysis.py --example --output-folder results/example

    python run_integrated_analysis.py --snv mmc2.xlsx --output-folder results/paper

    python run_integrated_analysis.py --fasta seqs.fasta --metadata meta.csv \
        --reference NC_045512.2.fasta --snv pooled.csv --positions counts.csv \
        --output-folder results/combined
"""

import argparse
import math
import shutil
import sys
from pathlib import Path

from genomics.integrated_analysis import run_integrated_analysis
from run_pipeline import run_pipeline


PROJECT = Path(__file__).resolve().parent
EXAMPLE_FOLDER = PROJECT / "examples" / "integrated"
EXAMPLE_FILES = {
    "fasta": "synthetic_individual_sequences.fasta",
    "metadata": "synthetic_individual_metadata.csv",
    "reference": "synthetic_reference.fasta",
    "snv": "synthetic_pooled_snvs.csv",
    "positions": "synthetic_pooled_base_counts.csv",
    "signatures": "synthetic_variant_signatures.csv",
}


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(description="Run the integrated genomic surveillance analysis.")
    parser.add_argument("--example", action="store_true",
                        help="Use the SYNTHETIC paired example in examples/integrated/")
    parser.add_argument("--fasta", type=Path, help="Individual sequences (.fasta/.fa/.fna)")
    parser.add_argument("--metadata", type=Path, help="Individual metadata (.csv/.xlsx)")
    parser.add_argument("--reference", type=Path, help="Reference FASTA with one sequence")
    parser.add_argument("--snv", type=Path, help="Pooled processed SNV table")
    parser.add_argument("--pooled-metadata", type=Path, help="Pooled-sample metadata")
    parser.add_argument("--positions", type=Path, help="Per-position depth or A/C/G/T counts")
    parser.add_argument("--masked-sites", type=Path, help="Problematic-site list")
    parser.add_argument("--signatures", type=Path, help="Documented variant profiles")
    parser.add_argument("--max-p-distance", type=float, default=0.05)
    parser.add_argument("--max-days-apart", type=float, default=14)
    parser.add_argument("--max-km-apart", type=float, default=50)
    parser.add_argument("--complete-genomes", action="store_true",
                        help="Confirm all individual sequences are complete genomes of one pathogen")
    parser.add_argument("--qc-mode", choices=["paper", "exploratory"], default="paper")
    parser.add_argument("--date-window-days", type=float, default=14)
    parser.add_argument("--output-folder", type=Path, required=True)
    return parser.parse_args(arguments)


def main(arguments=None):
    options = parse_arguments(arguments)
    inputs = {
        "fasta": options.fasta, "metadata": options.metadata, "reference": options.reference,
        "snv": options.snv, "pooled_metadata": options.pooled_metadata,
        "positions": options.positions, "masked_sites": options.masked_sites,
        "signatures": options.signatures,
    }
    if options.example:
        inputs.update({key: EXAMPLE_FOLDER / name for key, name in EXAMPLE_FILES.items()})
        print("Using the SYNTHETIC integrated example (invented data).")

    supplied = {key: path for key, path in inputs.items() if path is not None}
    for key, path in supplied.items():
        if not Path(path).is_file():
            print(f"ERROR: {key} file not found: {path}", file=sys.stderr)
            return 2
    if bool(inputs["fasta"]) != bool(inputs["metadata"]):
        print("ERROR: --fasta and --metadata must be given together.", file=sys.stderr)
        return 2
    if not inputs["fasta"] and not inputs["snv"]:
        print("ERROR: give individual sequences, a pooled SNV table, or both.", file=sys.stderr)
        return 2
    for value in [options.max_p_distance, options.max_days_apart, options.max_km_apart,
                  options.date_window_days]:
        if value < 0 or not math.isfinite(value):
            print("ERROR: thresholds and the date window must be finite and non-negative.",
                  file=sys.stderr)
            return 2

    output = options.output_folder
    output.mkdir(parents=True, exist_ok=True)

    individual_data = None
    mafft_path = None
    if inputs["fasta"]:
        # Work on copies so pipeline CSVs are written to the output folder.
        fasta_copy = output / Path(inputs["fasta"]).name
        metadata_copy = output / Path(inputs["metadata"]).name
        shutil.copyfile(inputs["fasta"], fasta_copy)
        shutil.copyfile(inputs["metadata"], metadata_copy)
        shared = {}
        status = run_pipeline(fasta_copy, metadata_copy, options.max_p_distance,
                              options.max_days_apart, options.max_km_apart,
                              complete_genomes=options.complete_genomes, analysis_data=shared)
        if status != 0:
            return status
        if inputs["reference"]:
            individual_data = shared["samples"]
            mafft_path = shared["mafft_path"]

    if not inputs["reference"] and not inputs["snv"]:
        print("Sequence-relationship results written; no reference or pooled table was given.")
        return 0

    try:
        result = run_integrated_analysis(
            individual_data=individual_data, mafft_path=mafft_path,
            reference_path=inputs["reference"], snv_path=inputs["snv"],
            pooled_metadata_path=inputs["pooled_metadata"], position_path=inputs["positions"],
            masked_sites_path=inputs["masked_sites"], signature_path=inputs["signatures"],
            settings={"qc_mode": options.qc_mode, "date_window_days": options.date_window_days},
            input_labels={key: Path(path).name for key, path in supplied.items()},
        )
    except (ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    for name, table in result["tables"].items():
        table.to_csv(output / name, index=False)
    for name, figure in result["figures"].items():
        (output / f"{name}.png").write_bytes(figure["png"])
    for section, notes in result["notes"].items():
        for text in notes:
            print(f"[{section}] {text}")
    print("Summary:")
    for key, value in result["summary"].items():
        print(f"  {key}: {value}")
    print(f"Wrote {len(result['tables'])} tables and {len(result['figures'])} figures to {output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
