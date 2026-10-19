"""Small Flask interface for the existing analysis pipeline."""

from base64 import b64encode
from contextlib import redirect_stderr
from io import BytesIO, StringIO
from math import isfinite
from pathlib import Path
from secrets import token_urlsafe
from shutil import copyfile
from tempfile import TemporaryDirectory
from threading import Lock
from time import monotonic

import pandas as pd
from flask import Flask, abort, render_template, request, send_file

from geocode_locations import check_geonames
from run_pipeline import output_paths, region_assignment_path, run_pipeline
from system_diagnostics import check_mafft
from validate_inputs import FASTA_EXTENSIONS, METADATA_EXTENSIONS
from environmental_analysis import load_measurements, make_charts, summarize_measurements, table_records
from wastewater.processed_snv import load_processed_snvs, summarize_snvs


app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 10 * 1024 * 1024  # 10 MB per request.
EXAMPLE_FOLDER = Path(__file__).resolve().parent / "examples"
STATUS_CACHE_SECONDS = 300
STATUS_CACHE = {"checked_at": 0.0, "results": None}
STATUS_LOCK = Lock()
WASTEWATER_DOWNLOADS = {}  # Small, short-lived CSV downloads; no database.
WASTEWATER_LOCK = Lock()
WASTEWATER_DOWNLOAD_SECONDS = 600
ENVIRONMENTAL_DOWNLOADS = {}  # Short-lived CSV files; no database or permanent upload storage.
ENVIRONMENTAL_LOCK = Lock()
ENVIRONMENTAL_DOWNLOAD_SECONDS = 600
SYNTHETIC_ENVIRONMENTAL_FILE = EXAMPLE_FOLDER / "synthetic_environmental_measurements.csv"


@app.route("/health")
def health():
    """Let the hosting service check that the website is responding."""

    return "OK", 200


@app.errorhandler(413)
def upload_too_large(error):
    """Show a clear size-limit message instead of a server error page."""

    return "Upload too large. Files in one request must total less than 10 MB.", 413


@app.errorhandler(500)
def internal_error(error):
    """Keep unexpected internal details out of the public response."""

    return "The analysis could not be completed. Please check your files and try again.", 500


def get_system_status():
    """Run both checks at most once every five minutes per server process."""

    with STATUS_LOCK:
        now = monotonic()
        if (
            STATUS_CACHE["results"] is not None
            and now - STATUS_CACHE["checked_at"] < STATUS_CACHE_SECONDS
        ):
            return STATUS_CACHE["results"], True

        try:
            mafft_status = check_mafft()
        except Exception:
            app.logger.warning("MAFFT diagnostic failed unexpectedly")
            mafft_status = {
                "status": "Check failed",
                "message": "The MAFFT check could not finish.",
            }

        try:
            geonames_status = check_geonames()
        except Exception:
            app.logger.warning("GeoNames diagnostic failed unexpectedly")
            geonames_status = {
                "status": "Check failed",
                "message": "The GeoNames check could not finish.",
            }

        results = {"mafft": mafft_status, "geonames": geonames_status}
        STATUS_CACHE["checked_at"] = now
        STATUS_CACHE["results"] = results
        return results, False


@app.route("/system-status", methods=["POST"])
def system_status():
    """Check external dependencies only when the visitor presses the button."""

    results, cached = get_system_status()
    values = {
        "max_p_distance": "0.05",
        "max_days_apart": "14",
        "max_km_apart": "50",
        "genomic_region_mode": "automatic",
        "confirm_complete_genomes": False,
    }
    return render_template(
        "index.html", values=values, system_status=results, status_cached=cached
    )


@app.route("/wastewater-data", methods=["GET", "POST"], endpoint="wastewater_data_page")
@app.route("/environmental", methods=["GET", "POST"])
def environmental_page():
    """Analyze processed measurements without running the pathogen pipeline."""

    error = None
    result = None

    if request.method == "POST":
        use_example = request.form.get("run_example") == "yes"
        upload = request.files.get("measurement_file")

        if use_example:
            if not SYNTHETIC_ENVIRONMENTAL_FILE.is_file():
                error = "The bundled synthetic environmental example is missing."
        elif upload is None or not upload.filename:
            error = "Choose a measurement CSV or XLSX file."
        elif Path(upload.filename).suffix.lower() not in METADATA_EXTENSIONS:
            error = "Measurement files must use the .csv or .xlsx format."

        if error is None:
            try:
                with TemporaryDirectory() as folder:
                    if use_example:
                        input_path = SYNTHETIC_ENVIRONMENTAL_FILE
                        filename = "Synthetic example (invented data)"
                    else:
                        filename = Path(upload.filename).name
                        input_path = Path(folder) / f"measurements{Path(filename).suffix.lower()}"
                        upload.save(input_path)
                        if input_path.stat().st_size > 3 * 1024 * 1024:
                            raise ValueError("Measurement files must be under 3 MB for this class prototype.")

                    data, metric = load_measurements(input_path)
                    summary = summarize_measurements(data, metric)
                    chart = b64encode(
                        make_charts(data, metric, summary["tables"])
                    ).decode("ascii")

                    organized = data.copy()
                    organized["Collection_Date"] = (
                        organized["Collection_Date"].dt.strftime("%Y-%m-%d").fillna("")
                    )
                    csv_bytes = organized.to_csv(index=False).encode("utf-8")
                    if len(csv_bytes) > 4 * 1024 * 1024:
                        raise ValueError("The organized CSV is too large for this prototype (4 MB limit).")

                    with ENVIRONMENTAL_LOCK:
                        now = monotonic()
                        expired = [
                            key for key, item in ENVIRONMENTAL_DOWNLOADS.items()
                            if now - item[0] > ENVIRONMENTAL_DOWNLOAD_SECONDS
                        ]
                        for key in expired:
                            del ENVIRONMENTAL_DOWNLOADS[key]
                        if len(ENVIRONMENTAL_DOWNLOADS) >= 4:
                            oldest = min(
                                ENVIRONMENTAL_DOWNLOADS,
                                key=lambda key: ENVIRONMENTAL_DOWNLOADS[key][0],
                            )
                            del ENVIRONMENTAL_DOWNLOADS[oldest]
                        download_id = token_urlsafe(16)
                        ENVIRONMENTAL_DOWNLOADS[download_id] = (now, csv_bytes)

                    tables = summary["tables"]
                    preview = table_records(tables["measurements"].head(100))
                    result = {
                        "filename": filename,
                        "metric": metric,
                        "summary": summary["counts"],
                        "tables": {
                            "samples": table_records(tables["samples"]),
                            "temporal_changes": table_records(tables["temporal_changes"]),
                            "location_summary": table_records(tables["location_summary"]),
                            "location_counts": table_records(tables["location_counts"]),
                            "date_counts": table_records(tables["date_counts"]),
                            "missing_data": table_records(tables["missing_data"]),
                            "quality_summary": table_records(tables["quality_summary"]),
                        },
                        "measurement_columns": list(organized.columns),
                        "preview": preview,
                        "chart": chart,
                        "download_id": download_id,
                    }
            except ValueError as exc:
                error = str(exc)
            except Exception:
                app.logger.exception("Environmental analysis failed")
                error = "Could not analyze this measurement file. Check that it is a valid CSV or XLSX workbook."

    wastewater_view = request.endpoint == "wastewater_data_page"
    return render_template(
        "environmental.html", error=error, result=result,
        page_title="Wastewater Data Analysis" if wastewater_view else "Environmental Data Analysis",
        download_endpoint="wastewater_data_download" if wastewater_view else "environmental_download",
    )


@app.route("/wastewater-data/download/<download_id>", endpoint="wastewater_data_download")
@app.route("/environmental/download/<download_id>")
def environmental_download(download_id):
    """Return the recent organized CSV while its short-lived download is available."""

    with ENVIRONMENTAL_LOCK:
        item = ENVIRONMENTAL_DOWNLOADS.get(download_id)
        if item is None or monotonic() - item[0] > ENVIRONMENTAL_DOWNLOAD_SECONDS:
            abort(404)
        csv_bytes = item[1]
    return send_file(
        BytesIO(csv_bytes),
        mimetype="text/csv",
        as_attachment=True,
        download_name=(
            "wastewater_measurements_organized.csv"
            if request.endpoint == "wastewater_data_download"
            else "environmental_measurements_organized.csv"
        ),
    )


@app.route("/wastewater", methods=["GET", "POST"])
def wastewater_page():
    """Show processed wastewater SNVs without using the Sprint 2 pipeline."""

    error = None
    result = None
    if request.method == "POST":
        snv_upload = request.files.get("snv_file")
        metadata_upload = request.files.get("sample_metadata")
        if not snv_upload or not snv_upload.filename:
            error = "Choose a processed SNV .xlsx or .csv file."
        elif Path(snv_upload.filename).suffix.lower() not in {".xlsx", ".csv"}:
            error = "Processed SNV data must be an .xlsx or .csv file."
        elif metadata_upload and metadata_upload.filename and Path(metadata_upload.filename).suffix.lower() not in {".xlsx", ".csv"}:
            error = "Optional sample metadata must be an .xlsx or .csv file."

        if error is None:
            try:
                with TemporaryDirectory() as folder:
                    snv_path = Path(folder) / f"snvs{Path(snv_upload.filename).suffix.lower()}"
                    snv_upload.save(snv_path)
                    metadata_path = None
                    if metadata_upload and metadata_upload.filename:
                        metadata_path = Path(folder) / f"sample_metadata{Path(metadata_upload.filename).suffix.lower()}"
                        metadata_upload.save(metadata_path)
                    if snv_path.stat().st_size > 3 * 1024 * 1024:
                        raise ValueError("Processed SNV files must be under 3 MB for this class prototype.")
                    if metadata_path and metadata_path.stat().st_size > 1024 * 1024:
                        raise ValueError("Optional wastewater metadata files must be under 1 MB.")

                    data, warnings = load_processed_snvs(snv_path, metadata_path)
                    summary = summarize_snvs(data)
                    display = data.copy()
                    display["Collection_Date"] = display["Collection_Date"].dt.strftime("%Y-%m-%d").fillna("")
                    csv_bytes = display.to_csv(index=False).encode("utf-8")
                    if len(csv_bytes) > 4 * 1024 * 1024:
                        raise ValueError("The standardized CSV is too large for this prototype (4 MB limit).")
                    with WASTEWATER_LOCK:
                        now = monotonic()
                        expired = [key for key, item in WASTEWATER_DOWNLOADS.items() if now - item[0] > WASTEWATER_DOWNLOAD_SECONDS]
                        for key in expired:
                            del WASTEWATER_DOWNLOADS[key]
                        if len(WASTEWATER_DOWNLOADS) >= 4:
                            oldest = min(WASTEWATER_DOWNLOADS, key=lambda key: WASTEWATER_DOWNLOADS[key][0])
                            del WASTEWATER_DOWNLOADS[oldest]
                        download_id = token_urlsafe(16)
                        WASTEWATER_DOWNLOADS[download_id] = (now, csv_bytes)

                    preview = display.head(100).copy()
                    preview["Read_Depth"] = preview["Read_Depth"].astype(object).where(preview["Read_Depth"].notna(), "")
                    result = {
                        "filename": Path(snv_upload.filename).name,
                        "metadata_filename": Path(metadata_upload.filename).name if metadata_path else None,
                        "summary": summary,
                        "warnings": warnings,
                        "preview": preview.to_dict("records"),
                        "download_id": download_id,
                    }
            except ValueError as exc:
                error = str(exc)
            except Exception:
                app.logger.exception("Wastewater analysis failed")
                error = "Could not read the processed SNV file. Check that it is a valid CSV or XLSX workbook."

    return render_template("wastewater.html", error=error, result=result)


@app.route("/wastewater/download/<download_id>")
def wastewater_download(download_id):
    """Download a recent standardized CSV; the uploaded files are already gone."""

    with WASTEWATER_LOCK:
        item = WASTEWATER_DOWNLOADS.get(download_id)
        if item is None or monotonic() - item[0] > WASTEWATER_DOWNLOAD_SECONDS:
            abort(404)
        csv_bytes = item[1]
    return send_file(
        BytesIO(csv_bytes),
        mimetype="text/csv",
        as_attachment=True,
        download_name="wastewater_standardized_snvs.csv",
    )


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
