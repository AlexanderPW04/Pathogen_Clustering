"""Small Flask interface for the existing analysis pipeline."""

from contextlib import redirect_stderr
from io import StringIO
from math import isfinite
from pathlib import Path
from shutil import copyfile
from tempfile import TemporaryDirectory

import pandas as pd
from flask import Flask, render_template, request

from run_pipeline import output_paths, region_assignment_path, run_pipeline
from validate_inputs import FASTA_EXTENSIONS, METADATA_EXTENSIONS


app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 10 * 1024 * 1024  # 10 MB per request.
EXAMPLE_FOLDER = Path(__file__).resolve().parent / "examples"


@app.route("/health")
def health():
    """Let the hosting service check that the website is responding."""

    return "OK", 200


@app.errorhandler(413)
def upload_too_large(error):
    """Show a clear size-limit message instead of a server error page."""

    return "Upload too large. The FASTA and metadata files must total less than 10 MB.", 413


@app.errorhandler(500)
def internal_error(error):
    """Keep unexpected internal details out of the public response."""

    return "The analysis could not be completed. Please check your files and try again.", 500


def read_thresholds(values):
    """Convert the three threshold fields to valid non-negative numbers."""

    try:
        max_p_distance = float(values["max_p_distance"])
        max_days_apart = float(values["max_days_apart"])
        max_km_apart = float(values["max_km_apart"])
    except ValueError:
        return None, "All thresholds must be non-negative numbers."

    thresholds = [max_p_distance, max_days_apart, max_km_apart]

    for threshold in thresholds:
        if threshold < 0 or not isfinite(threshold):
            return None, "All thresholds must be non-negative numbers."

    return thresholds, None


def make_cluster_list(membership):
    """Convert the membership DataFrame into simple dictionaries for HTML."""

    clusters = []
    group_columns = ["Cluster_ID", "Pathogen", "Gene_Segment"]
    grouped_membership = membership.groupby(group_columns, sort=True)

    for group_names, group in grouped_membership:
        cluster_id, pathogen, gene_segment = group_names
        samples = group["Sample_ID"].astype(str).tolist()

        cluster = {
            "cluster_id": cluster_id,
            "pathogen": pathogen,
            "gene_segment": gene_segment,
            "samples": samples,
        }
        clusters.append(cluster)

    return clusters


def find_geography_unresolved_candidates(pairs, thresholds):
    """List genetic-and-date matches that lack usable geographic evidence."""

    candidates = []
    for _, pair in pairs.iterrows():
        if pd.isna(pair["P_Distance"]) or pd.isna(pair["Days_Apart"]):
            continue
        if pair["P_Distance"] > thresholds[0] or pair["Days_Apart"] > thresholds[1]:
            continue
        if pd.notna(pair["Km_Apart"]):
            # Known distances, even those above the cutoff, are not unresolved.
            continue
        if str(pair["Same_Location"]).strip().casefold() == "true":
            # The existing strict stage already accepts this named-place fallback.
            continue

        candidates.append(
            {
                "sample_1": str(pair["Sample_ID_1"]),
                "sample_2": str(pair["Sample_ID_2"]),
                "pathogen": pair["Pathogen"],
                "gene_segment": pair["Gene_Segment"],
                "p_distance": pair["P_Distance"],
                "days_apart": pair["Days_Apart"],
            }
        )
    return candidates


@app.route("/", methods=["GET", "POST"])
def index():
    """Show the upload form and process it when the user submits it."""

    error = None
    notice = None
    analysis_complete = False
    clusters = []
    related_pair_count = 0
    candidates = []
    summary = {}
    input_names = None
    is_synthetic_example = False

    values = {
        "max_p_distance": request.form.get("max_p_distance", "0.05"),
        "max_days_apart": request.form.get("max_days_apart", "14"),
        "max_km_apart": request.form.get("max_km_apart", "50"),
        "genomic_region_mode": request.form.get("genomic_region_mode", "automatic"),
        "confirm_complete_genomes": request.form.get("confirm_complete_genomes") == "on",
    }

    if request.method == "POST":
        use_example = "run_example" in request.form
        if use_example:
            is_synthetic_example = True
            values["genomic_region_mode"] = "automatic"
            fasta_upload = None
            metadata_upload = None
        else:
            fasta_upload = request.files.get("fasta_file")
            metadata_upload = request.files.get("metadata_file")

        # Check both uploaded filenames before saving anything.
        if values["genomic_region_mode"] not in {"automatic", "complete_genomes"}:
            error = "Choose Automatic detection or Complete genomes."
        elif (
            values["genomic_region_mode"] == "complete_genomes"
            and not values["confirm_complete_genomes"]
        ):
            error = "Please confirm that all sequences are complete genomes of one pathogen."
        elif use_example and not all(
            (EXAMPLE_FOLDER / name).is_file()
            for name in ("synthetic_sequences.fasta", "synthetic_metadata.csv")
        ):
            error = "The bundled synthetic example files are missing."
        elif not use_example and (fasta_upload is None or fasta_upload.filename == ""):
            error = "Please choose a FASTA sequence file."
        elif not use_example and Path(fasta_upload.filename).suffix.lower() not in FASTA_EXTENSIONS:
            error = "The sequence file must use .fasta, .fa, or .fna."
        elif not use_example and (metadata_upload is None or metadata_upload.filename == ""):
            error = "Please choose a metadata file."
        elif not use_example and Path(metadata_upload.filename).suffix.lower() not in METADATA_EXTENSIONS:
            error = "The metadata file must use .xlsx or .csv."

        thresholds = None
        if error is None:
            thresholds, error = read_thresholds(values)

        if error is None:
            # Every request gets its own folder, which is deleted automatically.
            with TemporaryDirectory() as temporary_folder:
                if use_example:
                    fasta_name = "synthetic_sequences.fasta"
                    metadata_name = "synthetic_metadata.csv"
                    fasta_path = Path(temporary_folder) / fasta_name
                    metadata_path = Path(temporary_folder) / metadata_name
                    copyfile(EXAMPLE_FOLDER / fasta_name, fasta_path)
                    copyfile(EXAMPLE_FOLDER / metadata_name, metadata_path)
                else:
                    fasta_name = Path(fasta_upload.filename).name
                    metadata_name = Path(metadata_upload.filename).name
                    fasta_path = Path(temporary_folder) / f"sequences{Path(fasta_name).suffix.lower()}"
                    metadata_path = Path(temporary_folder) / f"metadata{Path(metadata_name).suffix.lower()}"
                    fasta_upload.save(fasta_path)
                    metadata_upload.save(metadata_path)

                input_names = {"fasta": fasta_name, "metadata": metadata_name}

                # Capture backend errors so they can be shown on the webpage.
                pipeline_errors = StringIO()
                with redirect_stderr(pipeline_errors):
                    status = run_pipeline(
                        fasta_path,
                        metadata_path,
                        thresholds[0],
                        thresholds[1],
                        thresholds[2],
                        complete_genomes=(
                            values["genomic_region_mode"] == "complete_genomes"
                        ),
                        run_summary=summary,
                    )

                if status != 0:
                    error = pipeline_errors.getvalue().strip()
                    error = error.replace(str(temporary_folder), "[temporary upload]")
                    if "Traceback" in error:
                        error = "The analysis could not be completed. Please check your files."
                    if not error:
                        error = "The analysis pipeline could not complete."
                else:
                    notice = pipeline_errors.getvalue().strip() or None
                    output_files = output_paths(metadata_path)
                    related_path = output_files[2]
                    membership_path = output_files[3]
                    enriched_path = output_files[1]

                    try:
                        related_pairs = pd.read_csv(related_path)
                        enriched_pairs = pd.read_csv(
                            enriched_path,
                            dtype={"Sample_ID_1": "string", "Sample_ID_2": "string"},
                        )
                        membership = pd.read_csv(
                            membership_path, dtype={"Sample_ID": "string"}
                        )
                        assignments = pd.read_csv(
                            region_assignment_path(metadata_path),
                            dtype={"Sample_ID": "string"},
                        )
                    except Exception:
                        app.logger.exception("Could not read pipeline results")
                        error = "Could not read the analysis results. Please try again."
                    else:
                        related_pair_count = len(related_pairs)
                        clusters = make_cluster_list(membership)
                        candidates = find_geography_unresolved_candidates(
                            enriched_pairs, thresholds
                        )
                        candidate_sample_ids = set()
                        for pair in candidates:
                            candidate_sample_ids.add(pair["sample_1"])
                            candidate_sample_ids.add(pair["sample_2"])

                        summary["unique_samples"] = int(
                            assignments["Sample_ID"].nunique()
                        )
                        summary["pairwise_comparisons"] = len(enriched_pairs)
                        summary["candidate_unique_samples"] = len(
                            candidate_sample_ids
                        )
                        summary["valid_geographic_comparisons"] = int(
                            enriched_pairs["Km_Apart"].notna().sum()
                        )
                        region_by_sample = assignments.set_index("Sample_ID")
                        for cluster in clusters:
                            cluster["sample_regions"] = []
                            for sample_id in cluster["samples"]:
                                region = region_by_sample.loc[sample_id]
                                cluster["sample_regions"].append(
                                    f"{sample_id}: {region['Gene_Segment']} "
                                    f"({region['Gene_Segment_Source']})"
                                )
                        analysis_complete = True

    return render_template(
        "index.html",
        error=error,
        notice=notice,
        analysis_complete=analysis_complete,
        clusters=clusters,
        related_pair_count=related_pair_count,
        candidates=candidates,
        summary=summary,
        input_names=input_names,
        is_synthetic_example=is_synthetic_example,
        values=values,
    )


if __name__ == "__main__":
    app.run(debug=False)
