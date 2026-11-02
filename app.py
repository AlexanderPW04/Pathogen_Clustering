"""Flask interface: one dashboard for individual and pooled genomic evidence.

The page guides the visitor through three steps: upload any mix of files,
review what each file was recognized as (and correct it), then run every
analysis the files allow. The original Sprint 2 sequence-relationship
pipeline runs for individual sequences, and the reference-based integrated
analysis (genomics.integrated_analysis) connects them with pooled samples.
The earlier Sprint 3 single-source pages remain available.
"""

from base64 import b64encode
from contextlib import redirect_stderr
from io import BytesIO, StringIO
from math import isfinite
from pathlib import Path
from secrets import token_urlsafe
from shutil import copyfile, rmtree
from tempfile import TemporaryDirectory, gettempdir, mkdtemp
from threading import Lock
from time import monotonic

import pandas as pd
from flask import Flask, abort, render_template, request, send_file

from geocode_locations import check_geonames
from run_pipeline import output_paths, region_assignment_path, run_pipeline
from system_diagnostics import check_mafft
from validate_inputs import FASTA_EXTENSIONS, METADATA_EXTENSIONS
from environmental_analysis import load_measurements, make_charts, summarize_measurements, table_records
from genomics.file_formats import (
    EXCEL_KIND, FASTA_KIND, VCF_KIND, detect_kind, plain_copy, split_name,
)
from genomics.input_inspection import (
    ROLE_LABELS, empty_info, inspect_file, prepare_inputs, review_dataset, suggest_roles,
)
from genomics.integrated_analysis import run_integrated_analysis
from integrated_view import build_dashboard, download_files
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

# The paired, explicitly synthetic individual + pooled demonstration.
INTEGRATED_FOLDER = EXAMPLE_FOLDER / "integrated"
INTEGRATED_EXAMPLE_FILES = {
    "fasta": "synthetic_individual_sequences.fasta",
    "metadata": "synthetic_individual_metadata.csv",
    "reference": "synthetic_reference.fasta",
    "snv": "synthetic_pooled_snvs.csv",
    "positions": "synthetic_pooled_base_counts.csv",
    "signatures": "synthetic_variant_signatures.csv",
}
TABLE_EXTENSIONS = {".csv", ".tsv", ".xlsx"}
# (form field, input key, allowed extensions, readable name)
OPTIONAL_UPLOADS = [
    ("reference_file", "reference", FASTA_EXTENSIONS, "reference FASTA"),
    ("snv_file", "snv", TABLE_EXTENSIONS | {".vcf"}, "pooled SNV table"),
    ("pooled_metadata_file", "pooled_metadata", TABLE_EXTENSIONS, "pooled-sample metadata"),
    ("position_file", "positions", TABLE_EXTENSIONS, "per-position depth table"),
    ("masked_sites_file", "masked_sites", TABLE_EXTENSIONS, "masked-site list"),
    ("signature_file", "signatures", TABLE_EXTENSIONS, "variant-profile table"),
]
INPUT_LABELS = {
    "fasta": "Individual sequences (FASTA)",
    "metadata": "Individual metadata",
    "reference": "Reference FASTA",
    "snv": "Pooled SNV table",
    "pooled_metadata": "Pooled-sample metadata",
    "positions": "Per-position depth/base counts",
    "masked_sites": "Masked-site list",
    "signatures": "Variant profiles",
}
RESULT_DOWNLOADS = {}  # Short-lived integrated-result files; no database.
RESULT_LOCK = Lock()
RESULT_DOWNLOAD_SECONDS = 600
RESULT_DOWNLOAD_LIMIT = 4

# Guided uploads wait on disk (not in Git, not permanent) between the review
# and run steps; each set is deleted after 30 minutes or when newer sets need room.
STAGING_ROOT = Path(gettempdir()) / "genomic_surveillance_uploads"
STAGED = {}
STAGED_LOCK = Lock()
STAGED_SECONDS = 1800
STAGED_LIMIT = 6
MAX_GUIDED_FILES = 20


def default_values(form):
    """Form values with the prototype defaults."""

    return {
        "max_p_distance": form.get("max_p_distance", "0.05"),
        "max_days_apart": form.get("max_days_apart", "14"),
        "max_km_apart": form.get("max_km_apart", "50"),
        "genomic_region_mode": form.get("genomic_region_mode", "automatic"),
        "confirm_complete_genomes": form.get("confirm_complete_genomes") == "on",
        "qc_mode": form.get("qc_mode", "paper"),
        "date_window_days": form.get("date_window_days", "14"),
    }


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
    values = default_values({})
    return render_template(
        "index.html", values=values, system_status=results, status_cached=cached, step=1
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


def read_date_window(values):
    """Read the cross-source date window as a non-negative number of days."""

    try:
        days = float(values["date_window_days"])
    except ValueError:
        days = -1
    if days < 0 or not isfinite(days):
        return None, "The cross-source date window must be a non-negative number of days."
    return days, None


def store_result_files(files):
    """Keep a few recent result files in memory for ten minutes; no database."""

    with RESULT_LOCK:
        now = monotonic()
        expired = [
            key for key, item in RESULT_DOWNLOADS.items()
            if now - item[0] > RESULT_DOWNLOAD_SECONDS
        ]
        for key in expired:
            del RESULT_DOWNLOADS[key]
        while len(RESULT_DOWNLOADS) >= RESULT_DOWNLOAD_LIMIT:
            oldest = min(RESULT_DOWNLOADS, key=lambda key: RESULT_DOWNLOADS[key][0])
            del RESULT_DOWNLOADS[oldest]
        download_id = token_urlsafe(16)
        RESULT_DOWNLOADS[download_id] = (now, files)
    return download_id


@app.route("/results/<download_id>/<filename>")
def result_download(download_id, filename):
    """Download one table or figure from a recent integrated analysis."""

    with RESULT_LOCK:
        item = RESULT_DOWNLOADS.get(download_id)
        if (item is None or monotonic() - item[0] > RESULT_DOWNLOAD_SECONDS
                or filename not in item[1]):
            abort(404)
        content, mimetype = item[1][filename]
    return send_file(BytesIO(content), mimetype=mimetype, as_attachment=True,
                     download_name=filename)


def run_relationship_analysis(fasta_path, metadata_path, thresholds, complete_genomes,
                              temporary_folder):
    """Run the original Sprint 2 pipeline and read its results for the page."""

    outcome = {
        "error": None, "notice": None, "complete": False, "clusters": [],
        "related_pair_count": 0, "candidates": [], "summary": {}, "analysis_data": {},
        "files": {},
    }
    summary = outcome["summary"]

    # Capture backend errors so they can be shown on the webpage.
    pipeline_errors = StringIO()
    with redirect_stderr(pipeline_errors):
        status = run_pipeline(
            fasta_path,
            metadata_path,
            thresholds[0],
            thresholds[1],
            thresholds[2],
            complete_genomes=complete_genomes,
            run_summary=summary,
            analysis_data=outcome["analysis_data"],
        )

    if status != 0:
        error = pipeline_errors.getvalue().strip()
        error = error.replace(str(temporary_folder), "[temporary upload]")
        if "Traceback" in error:
            error = "The analysis could not be completed. Please check your files."
        if not error:
            error = "The analysis pipeline could not complete."
        outcome["error"] = error
        return outcome

    outcome["notice"] = pipeline_errors.getvalue().strip() or None
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
        membership = pd.read_csv(membership_path, dtype={"Sample_ID": "string"})
        assignments = pd.read_csv(
            region_assignment_path(metadata_path),
            dtype={"Sample_ID": "string"},
        )
    except Exception:
        app.logger.exception("Could not read pipeline results")
        outcome["error"] = "Could not read the analysis results. Please try again."
        return outcome

    outcome["related_pair_count"] = len(related_pairs)
    clusters = make_cluster_list(membership)
    candidates = find_geography_unresolved_candidates(enriched_pairs, thresholds)
    candidate_sample_ids = set()
    for pair in candidates:
        candidate_sample_ids.add(pair["sample_1"])
        candidate_sample_ids.add(pair["sample_2"])

    summary["unique_samples"] = int(assignments["Sample_ID"].nunique())
    summary["pairwise_comparisons"] = len(enriched_pairs)
    summary["candidate_unique_samples"] = len(candidate_sample_ids)
    summary["valid_geographic_comparisons"] = int(enriched_pairs["Km_Apart"].notna().sum())
    region_by_sample = assignments.set_index("Sample_ID")
    for cluster in clusters:
        cluster["sample_regions"] = []
        for sample_id in cluster["samples"]:
            region = region_by_sample.loc[sample_id]
            cluster["sample_regions"].append(
                f"{sample_id}: {region['Gene_Segment']} "
                f"({region['Gene_Segment_Source']})"
            )

    outcome["clusters"] = clusters
    outcome["candidates"] = candidates
    outcome["pairs"] = enriched_pairs
    for path in list(output_files) + [region_assignment_path(metadata_path)]:
        outcome["files"][f"sequence_relationships_{path.name}"] = path.read_bytes()
    outcome["complete"] = True
    return outcome


def accepted_name(filename, extensions):
    """True when the base extension (ignoring a final .gz) is allowed."""

    extension, _ = split_name(filename)
    return extension in extensions


def save_upload(upload, folder, key):
    """Save one named-field upload, decompressing .gz into a plain file."""

    name = Path(upload.filename).name
    extension, _ = split_name(name)
    raw = folder / f"{key}_upload"
    upload.save(raw)
    path = folder / f"{key}{extension}"
    try:
        plain_copy(raw, path)
    finally:
        raw.unlink(missing_ok=True)
    return path, name


def save_inputs(folder, use_example, use_integrated, fasta_upload, metadata_upload, uploads):
    """Copy example files or save uploads; return (paths, display names)."""

    paths = {}
    names = {}
    if use_example:
        for key, name in [("fasta", "synthetic_sequences.fasta"),
                          ("metadata", "synthetic_metadata.csv")]:
            paths[key] = folder / name
            copyfile(EXAMPLE_FOLDER / name, paths[key])
            names[key] = name
        return paths, names
    if use_integrated:
        for key, name in INTEGRATED_EXAMPLE_FILES.items():
            paths[key] = folder / name
            copyfile(INTEGRATED_FOLDER / name, paths[key])
            names[key] = name
        return paths, names

    if fasta_upload is not None and fasta_upload.filename:
        paths["fasta"], names["fasta"] = save_upload(fasta_upload, folder, "sequences")
        paths["metadata"], names["metadata"] = save_upload(metadata_upload, folder, "metadata")
    for key, upload in uploads.items():
        paths[key], names[key] = save_upload(upload, folder, key)
    return paths, names


def new_page():
    return {
        "error": None, "notice": None, "analysis_complete": False, "clusters": [],
        "related_pair_count": 0, "candidates": [], "summary": {}, "input_names": None,
        "is_synthetic_example": False, "is_integrated_example": False, "dashboard": None,
        "integrated_error": None, "download_id": None, "downloads": [],
        "skipped_downloads": [], "review": None, "token": None, "step": 1,
        "analysis_plan": None, "declared_pathogen": "",
    }


def check_settings(values):
    """Validate thresholds and options; return (thresholds, date_window, error)."""

    if values["genomic_region_mode"] not in {"automatic", "complete_genomes"}:
        return None, None, "Choose Automatic detection or Complete genomes."
    if values["genomic_region_mode"] == "complete_genomes" and not values["confirm_complete_genomes"]:
        return None, None, "Please confirm that all sequences are complete genomes of one pathogen."
    thresholds, error = read_thresholds(values)
    if error:
        return None, None, error
    date_window, error = read_date_window(values)
    if error:
        return None, None, error
    if values["qc_mode"] not in {"paper", "exploratory"}:
        return None, None, "Choose the paper-reproduction or exploratory quality-control mode."
    return thresholds, date_window, None


def relationship_tables(relationship, thresholds):
    """The Sprint 2 pairwise results that the dashboard turns into figures."""

    return {"pairs": relationship.get("pairs"), "thresholds": thresholds,
            "candidates": relationship["candidates"], "clusters": relationship["clusters"]}


def run_analyses(paths, names, values, thresholds, date_window, temporary_folder, page):
    """Run every analysis the saved inputs allow and fill the results page."""

    page["input_names"] = names
    page["step"] = 3
    relationship = None
    if "fasta" in paths:
        relationship = run_relationship_analysis(
            paths["fasta"], paths["metadata"], thresholds,
            values["genomic_region_mode"] == "complete_genomes", temporary_folder,
        )
        page["error"] = relationship["error"]
        page["notice"] = relationship["notice"]
        page["analysis_complete"] = relationship["complete"]
        page["clusters"] = relationship["clusters"]
        page["related_pair_count"] = relationship["related_pair_count"]
        page["candidates"] = relationship["candidates"]
        page["summary"] = relationship["summary"]

    relationship_files = relationship["files"] if relationship else {}
    completed = relationship is not None and relationship["complete"]
    if not ("reference" in paths or "snv" in paths or completed):
        return page
    individual_data = None
    mafft_path = None
    if completed:
        individual_data = relationship["analysis_data"].get("samples")
        mafft_path = relationship["analysis_data"].get("mafft_path")
    try:
        result = run_integrated_analysis(
            individual_data=individual_data,
            mafft_path=mafft_path,
            reference_path=paths.get("reference"),
            snv_path=paths.get("snv"),
            pooled_metadata_path=paths.get("pooled_metadata"),
            position_path=paths.get("positions"),
            masked_sites_path=paths.get("masked_sites"),
            signature_path=paths.get("signatures"),
            settings={"qc_mode": values["qc_mode"], "date_window_days": date_window},
            input_labels=names,
            relationship=relationship_tables(relationship, thresholds) if completed else None,
        )
    except (ValueError, RuntimeError) as problem:
        page["integrated_error"] = str(problem).replace(str(temporary_folder), "[temporary upload]")
    except Exception:
        app.logger.exception("Integrated analysis failed")
        page["integrated_error"] = (
            "The integrated analysis could not be completed. Check the reference and "
            "pooled-sample files."
        )
    else:
        files, skipped = download_files(result, relationship_files)
        page["download_id"] = store_result_files(files)
        page["dashboard"] = build_dashboard(result, page["download_id"])
        page["downloads"] = sorted(files)
        page["skipped_downloads"] = skipped
        return page
    if relationship_files:
        files = {name: (content, "text/csv") for name, content in relationship_files.items()}
        page["download_id"] = store_result_files(files)
        page["downloads"] = sorted(files)
    return page


def render_page(values, page):
    return render_template("index.html", values=values, input_labels=INPUT_LABELS,
                           role_labels=ROLE_LABELS, **page)


@app.route("/", methods=["GET", "POST"])
def index():
    """Show the dashboard; POST runs the examples or a direct named-field upload.

    The named fields (fasta_file, metadata_file, reference_file, snv_file, ...)
    are kept for scripts and tests; the page itself uses the guided workflow.
    """

    values = default_values(request.form)
    page = new_page()

    if request.method == "POST":
        use_example = "run_example" in request.form
        use_integrated = "run_integrated_example" in request.form
        page["is_synthetic_example"] = use_example
        page["is_integrated_example"] = use_integrated
        uploads = {}
        fasta_upload = None
        metadata_upload = None
        if use_example or use_integrated:
            values["genomic_region_mode"] = "automatic"
        else:
            fasta_upload = request.files.get("fasta_file")
            metadata_upload = request.files.get("metadata_file")
            for field, key, _, _ in OPTIONAL_UPLOADS:
                upload = request.files.get(field)
                if upload is not None and upload.filename:
                    uploads[key] = upload
        has_fasta = fasta_upload is not None and fasta_upload.filename != ""
        has_metadata = metadata_upload is not None and metadata_upload.filename != ""

        # Check every filename before saving anything.
        error = None
        if values["genomic_region_mode"] not in {"automatic", "complete_genomes"}:
            error = "Choose Automatic detection or Complete genomes."
        elif (values["genomic_region_mode"] == "complete_genomes"
              and not values["confirm_complete_genomes"]):
            error = "Please confirm that all sequences are complete genomes of one pathogen."
        elif use_example and not all(
            (EXAMPLE_FOLDER / name).is_file()
            for name in ("synthetic_sequences.fasta", "synthetic_metadata.csv")
        ):
            error = "The bundled synthetic example files are missing."
        elif use_integrated and not all(
            (INTEGRATED_FOLDER / name).is_file() for name in INTEGRATED_EXAMPLE_FILES.values()
        ):
            error = "The bundled synthetic integrated example files are missing."
        elif not (use_example or use_integrated) and not (has_fasta or has_metadata or "snv" in uploads):
            error = ("Please choose a FASTA sequence file with its metadata, a pooled SNV "
                     "table, or both.")
        elif has_metadata and not has_fasta:
            error = "Please choose a FASTA sequence file."
        elif has_fasta and not accepted_name(fasta_upload.filename, FASTA_EXTENSIONS):
            error = "The sequence file must use .fasta, .fa, or .fna (optionally .gz)."
        elif has_fasta and not has_metadata:
            error = "Please choose a metadata file."
        elif has_metadata and not accepted_name(metadata_upload.filename, METADATA_EXTENSIONS):
            error = "The metadata file must use .xlsx or .csv."

        if error is None:
            for _, key, extensions, label in OPTIONAL_UPLOADS:
                if key in uploads and not accepted_name(uploads[key].filename, extensions):
                    error = (f"The {label} must use " + ", ".join(sorted(extensions))
                             + " (text files may be .gz).")
                    break
        if error is None and "snv" not in uploads and any(
            key in uploads for key in ("pooled_metadata", "positions", "masked_sites")
        ):
            error = "Pooled-sample metadata, position, and masked-site files need a pooled SNV table."

        thresholds = date_window = None
        if error is None:
            thresholds, date_window, error = check_settings(values)

        if error is None:
            # Every request gets its own folder, which is deleted automatically.
            with TemporaryDirectory() as temporary_folder:
                folder = Path(temporary_folder)
                try:
                    paths, names = save_inputs(folder, use_example, use_integrated,
                                               fasta_upload, metadata_upload, uploads)
                except ValueError as problem:
                    error = str(problem)
                else:
                    run_analyses(paths, names, values, thresholds, date_window,
                                 temporary_folder, page)

        if error is not None:
            page["error"] = error

    return render_page(values, page)


# ----------------------------------------------------------------------------
# Guided workflow: Step 1 upload -> Step 2 review -> Step 3 results
# ----------------------------------------------------------------------------

def remove_staged(token):
    item = STAGED.pop(token, None)
    if item is not None:
        rmtree(item["folder"], ignore_errors=True)


def clean_staged():
    """Delete expired staged uploads and keep at most STAGED_LIMIT sets."""

    now = monotonic()
    for token in [key for key, item in STAGED.items() if now - item["created"] > STAGED_SECONDS]:
        remove_staged(token)
    while len(STAGED) >= STAGED_LIMIT:
        remove_staged(min(STAGED, key=lambda key: STAGED[key]["created"]))


PLAIN_EXTENSIONS = {FASTA_KIND: ".fasta", VCF_KIND: ".vcf", EXCEL_KIND: ".xlsx"}


def stage_uploads(uploads):
    """Save, decompress, and inspect the uploaded files; return a staging token."""

    if len(uploads) > MAX_GUIDED_FILES:
        raise ValueError(f"Choose at most {MAX_GUIDED_FILES} files at a time.")
    STAGING_ROOT.mkdir(parents=True, exist_ok=True)
    folder = Path(mkdtemp(prefix="upload_", dir=STAGING_ROOT))
    staged = {}
    inspections = []
    try:
        for index, upload in enumerate(uploads):
            name = Path(upload.filename).name
            raw = folder / f"raw_{index}"
            upload.save(raw)
            plain = folder / f"file_{index}"
            try:
                plain_copy(raw, plain)
            except ValueError as problem:
                info = empty_info(name)
                info.update(kind="Unreadable", refused=True, problems=[str(problem)])
                inspections.append(info)
                continue
            finally:
                raw.unlink(missing_ok=True)
            kind = detect_kind(plain, name)
            extension = PLAIN_EXTENSIONS.get(kind)
            if extension is None:
                base, _ = split_name(name)
                extension = base if base in {".csv", ".tsv", ".txt"} else ".txt"
            path = plain.rename(folder / f"file_{index}{extension}")
            staged[index] = path
            inspections.append(inspect_file(path, name))
    except Exception:
        rmtree(folder, ignore_errors=True)
        raise
    suggest_roles(inspections)
    with STAGED_LOCK:
        clean_staged()
        token = token_urlsafe(16)
        STAGED[token] = {"folder": folder, "created": monotonic(), "files": staged,
                         "inspections": inspections}
    return token


def review_choices(form, inspections):
    """Read the role, column, and pathogen choices from the review form."""

    roles = {}
    choices = {}
    for index in range(len(inspections)):
        role = form.get(f"role_{index}")
        if role in ROLE_LABELS:
            roles[index] = role
        prefix = f"map_{index}_"
        picks = {key[len(prefix):]: value for key, value in form.items()
                 if key.startswith(prefix) and value}
        if picks:
            choices[index] = picks
    options = {"declared_pathogen": form.get("declared_pathogen", "").strip()[:120],
               "complete_genomes": form.get("genomic_region_mode") == "complete_genomes"}
    return roles, choices, options


def staged_item(token):
    with STAGED_LOCK:
        item = STAGED.get(token or "")
        if item is None or monotonic() - item["created"] > STAGED_SECONDS:
            return None
        return item


@app.route("/review", methods=["POST"])
def review_upload():
    """Step 2: show what each file is and which analyses the dataset allows."""

    values = default_values(request.form)
    page = new_page()
    page["step"] = 2
    token = request.form.get("token")
    uploads = [upload for upload in request.files.getlist("files") if upload and upload.filename]
    try:
        if uploads:
            token = stage_uploads(uploads)
        item = staged_item(token)
        if item is None:
            raise ValueError("Choose one or more files to upload. Uploaded files are kept for "
                             "30 minutes while you review them.")
        roles, choices, options = review_choices(request.form, item["inspections"])
        page["review"] = review_dataset(item["files"], item["inspections"], roles, choices, options)
        page["token"] = token
        page["declared_pathogen"] = options["declared_pathogen"]
    except ValueError as problem:
        page["error"] = str(problem)
        page["step"] = 1
    return render_page(values, page)


@app.route("/run", methods=["POST"])
def run_reviewed():
    """Step 3: run the applicable analyses on the reviewed files."""

    values = default_values(request.form)
    page = new_page()
    token = request.form.get("token")
    item = staged_item(token)
    if item is None:
        page["error"] = "The uploaded files have expired. Please upload them again."
        return render_page(values, page)
    roles, choices, options = review_choices(request.form, item["inspections"])
    review = review_dataset(item["files"], item["inspections"], roles, choices, options)
    thresholds, date_window, error = check_settings(values)
    if error is None and not review["ready"]:
        error = "Resolve the problems listed in Step 2 before running the analysis."
    if error is not None:
        page.update(error=error, review=review, token=token, step=2,
                    declared_pathogen=options["declared_pathogen"])
        return render_page(values, page)
    with TemporaryDirectory() as temporary_folder:
        try:
            paths, names = prepare_inputs(item["files"], review, temporary_folder, options)
        except ValueError as problem:
            page["error"] = str(problem)
        else:
            run_analyses(paths, names, values, thresholds, date_window, temporary_folder, page)
    page["analysis_plan"] = review["analyses"]
    page["token"] = token
    return render_page(values, page)


if __name__ == "__main__":
    app.run(debug=False)
