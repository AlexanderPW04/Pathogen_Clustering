"""Flask interface: one dashboard for individual and pooled genomic evidence.

The page guides the visitor through three steps: upload any mix of files,
review what each file was recognized as (and correct it), then run every
analysis the files allow. The original Sprint 2 sequence-relationship
pipeline runs for individual sequences, and the reference-based integrated
analysis (genomics.integrated_analysis) connects them with pooled samples.
The earlier Sprint 3 single-source pages remain available.

An optional literature-enrichment section (Sprint 4, literature/ and
literature_view.py) can propose metadata from linked publications; nothing
from it changes an analysis until a person approves each proposal.
"""

import ctypes
import gc
import zlib
from base64 import b64encode
from contextlib import redirect_stderr
from io import BytesIO, StringIO
from math import isfinite
from pathlib import Path
from secrets import token_urlsafe
from shutil import copyfile, rmtree
from tempfile import TemporaryDirectory, gettempdir, mkdtemp
from threading import BoundedSemaphore, Lock
from time import monotonic, time

import pandas as pd
from flask import Flask, abort, redirect, render_template, request, send_file, url_for
from werkzeug.datastructures import FileStorage

from compare_genetic_distances import TEXT_LABEL_COLUMNS, read_pair_table
from detect_potential_clusters import enough_overlap
from geocode_locations import check_geonames
from run_pipeline import output_paths, region_assignment_path, run_pipeline
from system_diagnostics import check_mafft
from validate_inputs import FASTA_EXTENSIONS, METADATA_EXTENSIONS
from environmental_analysis import load_measurements, make_charts, summarize_measurements, table_records
from genomics.file_formats import (
    ANNOTATION_EXTENSIONS, EXCEL_KIND, FASTA_KIND, MAX_DECOMPRESSED_BYTES, VCF_KIND, detect_kind, plain_copy,
    read_table, split_name,
)
from genomics.collection_dates import DATE_ORDERS
from genomics.overdispersion_analysis import DEFAULT_SENSITIVITY_RHOS, parse_rho_values
from genomics.mutation_group_analysis import PROFILES_MODE, parse_grouping_mode, parse_minimum_shared
from genomics.mutation_grouping import DEFAULT_MINIMUM_SHARED
from genomics.input_inspection import (
    ROLE_LABELS, ServerBusy, apply_declared_role, empty_info, inspect_file, prepare_inputs, review_dataset,
    suggest_roles,
)
from genomics.integrated_analysis import run_integrated_analysis
from genomics.import_options import DEPTH_DEFINITIONS, FREQUENCY_SCALES, POSITION_BASES, read_file_options
from genomics.genome_annotation import genetic_code_choices, parse_genetic_code
from genomics.mapping_profiles import apply_profile, profile_json, read_profile
from genomics.report_exports import MIMETYPES as REPORT_MIMETYPES, SUMMARY_NAME, build_from_stored
from genomics.metadata_matching import RESOLUTIONS
from genomics.spreadsheet_safety import escape_csv_bytes, safe_csv
from genomics.table_layouts import BASE_COUNT_LAYOUT, WIDE_LAYOUT
from dashboard_guide import GLOSSARY
from integrated_view import build_dashboard, download_files
from literature.config import ai_status
from literature.providers import ProviderError
from literature.review import enrichment_summary, review_summary
from literature.pubtator import check_pubtator
from literature.identifiers import parse_identifiers
from literature.preanalysis import (
    approved_changes, before_after, contribution, gap_resolution, pooled_overrides, preview_session,
    sections_using, write_individual_metadata,
)
from literature.workflow import new_session, record_decisions
from literature_view import (
    LITERATURE_EXAMPLE_FILES, LITERATURE_FOLDER, check_claude, choose_provider, chosen_batch,
    enrichment_downloads, find, get_session, live_request_allowed, prepare, review_and_apply, run_settings,
    start, store_session, view,
)
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
    ("annotation_file", "annotation", ANNOTATION_EXTENSIONS, "genome annotation (GenBank or GFF3)"),
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
    "annotation": "Genome annotation",
    "snv": "Pooled SNV table",
    "pooled_metadata": "Pooled-sample metadata",
    "positions": "Per-position depth/base counts",
    "masked_sites": "Masked-site list",
    "signatures": "Variant profiles",
}
# Optional "what do you want to analyze?" choice on the upload page. File roles
# are still detected automatically; the choice only lets the review step say
# which files that analysis would need.
MODES = {
    "auto": {"label": "Detect from my files", "files": "Upload any mix; the next step recognizes each file."},
    "individual": {"label": "Individual sequences (Sprint 2)",
                   "files": "FASTA sequences + sample metadata."},
    "pooled": {"label": "Pooled sequencing results (Sprint 3)",
               "files": "A processed SNV table (results of upstream variant calling)."},
    "combined": {"label": "Combined (Sprint 2 + Sprint 3)",
                 "files": "Both of the above + one shared reference genome FASTA."},
}
MODE_NEEDS = {
    "individual": [("fasta", "individual sequences (FASTA)"), ("metadata", "sample metadata")],
    "pooled": [("snv", "processed pooled SNV table")],
    "combined": [("fasta", "individual sequences (FASTA)"), ("metadata", "sample metadata"),
                 ("reference", "reference genome FASTA"), ("snv", "processed pooled SNV table")],
}
# Upload fields shown for each chosen mode: "required" or "optional" per mode.
# They follow the backend rules: individual metadata needs only Accession and
# Pathogen; a reference is required only to connect the two sources.
FASTA_ACCEPT = ".fasta,.fa,.fna,.gz"
TABLE_ACCEPT = ".csv,.tsv,.txt,.xlsx,.gz"
UPLOAD_FIELDS = [
    {"role": "fasta", "label": "Individual sequences (FASTA)", "accept": FASTA_ACCEPT, "multiple": True,
     "modes": {"individual": "required", "combined": "required"},
     "summary": "One record per sample",
     "help": "One FASTA record per sample; record names must match the Accession column of the metadata. "
             "Several FASTA files are combined. MAFFT aligns them on the server.",
     "example": "synthetic_individual_sequences.fasta"},
    {"role": "metadata", "label": "Individual sample metadata (CSV, TSV, or XLSX)", "accept": TABLE_ACCEPT,
     "multiple": False, "modes": {"individual": "required", "combined": "required"},
     "summary": "Required columns: Accession and Pathogen",
     "help": "Required columns: Accession (matching the FASTA names) and Pathogen. Optional: Collection_Date "
             "(date comparisons), Location_Name or Latitude/Longitude (distance screening; GeoNames can "
             "estimate coordinates for place names), and Gene_Segment. Missing optional columns only switch "
             "off the analyses that need them.",
     "example": "synthetic_individual_metadata.csv"},
    {"role": "reference", "label": "Reference genome (one-sequence FASTA)", "accept": FASTA_ACCEPT,
     "multiple": False, "modes": {"individual": "optional", "pooled": "optional", "combined": "required"},
     "summary": "Shared genome coordinates",
     "help": "Individual sequences: places them on numbered reference positions (reference-based variants). "
             "Pooled results: checks the reported reference bases and genome coverage. Combined: required, "
             "because both sources must use this same verified reference to be compared.",
     "example": "synthetic_reference.fasta"},
    {"role": "annotation", "label": "Genome annotation (GenBank or GFF3, optional)",
     "accept": ".gb,.gbk,.genbank,.gff,.gff3,.gz", "multiple": False,
     "modes": {"individual": "optional", "pooled": "optional", "combined": "optional"},
     "summary": "Genes and predicted amino-acid changes",
     "help": "The annotation of the same reference used for variant calling (for SARS-CoV-2, NCBI NC_045512.2 as "
             "GenBank or GFF3). It must name the same accession.version; a GenBank file also carries the sequence "
             "used to verify each reference base. Without it, every other analysis still runs.",
     "example": "NC_045512.2.gb"},
    {"role": "snv", "label": "Processed pooled SNV results (table or VCF)",
     "accept": ".csv,.tsv,.txt,.xlsx,.vcf,.gz", "multiple": True,
     "modes": {"pooled": "required", "combined": "required"},
     "summary": "Variant calls already made from sequencing reads",
     "help": "Sample, reference position, reference and alternative base, frequency, and read depth. Accepted: "
             "standardized tables, iVar variant tables, SNV VCF files, and the Fontenele et al. (2021) "
             "supplementary workbook; several iVar or VCF files are combined. Raw FASTQ reads are not accepted.",
     "example": "synthetic_pooled_snvs.csv"},
    {"role": "pooled_metadata", "label": "Pooled-sample metadata", "accept": TABLE_ACCEPT, "multiple": False,
     "modes": {"pooled": "optional", "combined": "optional"},
     "summary": "Dates, sites, coverage, mean depth",
     "help": "Sample_ID plus collection dates, sites, coordinates, genome coverage, or mean depth. Enables the "
             "paper's genome-coverage filter, the depth-versus-SNV regression, and dates or sites missing "
             "from the SNV table."},
    {"role": "positions", "label": "Per-position depth or A/C/G/T counts", "accept": TABLE_ACCEPT,
     "multiple": False, "modes": {"pooled": "optional", "combined": "optional"},
     "summary": "Coverage-aware results and the paper's PCoA",
     "help": "Depth at every listed position lets unreported variants be called 'not detected at adequate "
             "depth' instead of 'unknown'. A/C/G/T counts also enable the paper's Yue–Clayton PCoA "
             "(population structure).",
     "example": "synthetic_pooled_base_counts.csv"},
    {"role": "masked_sites", "label": "Masked problematic sites", "accept": TABLE_ACCEPT, "multiple": False,
     "advanced": True, "modes": {"pooled": "optional", "combined": "optional"},
     "summary": "Positions excluded by quality control",
     "help": "Reference positions to exclude, as in the paper's masked-site quality check."},
    {"role": "signatures", "label": "Documented variant profiles", "accept": TABLE_ACCEPT, "multiple": False,
     "advanced": True, "modes": {"pooled": "optional", "combined": "optional"},
     "summary": "Known variant sets for an overlap table",
     "help": "Lists of variants that define documented profiles; the results show their overlap with the "
             "data (never a lineage assignment).",
     "example": "synthetic_variant_signatures.csv"},
]
# Synthetic example files that visitors may download as format examples.
EXAMPLE_DOWNLOADS = {name: INTEGRATED_FOLDER / name for name in INTEGRATED_EXAMPLE_FILES.values()}
EXAMPLE_DOWNLOADS.update({name: EXAMPLE_FOLDER / name
                          for name in ["synthetic_sequences.fasta", "synthetic_metadata.csv"]})
# Official NCBI RefSeq annotation for SARS-CoV-2 (real, public; matches the published mmc2.xlsx reference).
EXAMPLE_DOWNLOADS["NC_045512.2.gb"] = EXAMPLE_FOLDER / "annotation" / "NC_045512.2.gb"
RESULT_DOWNLOADS = {}  # Short-lived integrated-result files; no database.
ANALYSIS_SLOT = BoundedSemaphore(1)
# Final polish 2B: one build per (results, file) at a time; repeated clicks wait for it.
REPORT_BUILDS = {}
REPORT_WAIT_SECONDS = 150
RESULTS_GONE = ("These results are no longer on the server. Result files are kept in memory for 10 minutes "
                "after an analysis (and are lost if the server restarts); run the analysis again to download them.")
GENETIC_CODES = genetic_code_choices()  # NCBI translation tables for the settings menu
OUT_OF_MEMORY = ("The uploaded dataset exceeds the available memory for this server, so the analysis "
                 "stopped. No results were changed or removed to make it fit. Run it on a local computer "
                 "(see the README), or upload a smaller study.")
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


DEFAULT_RHO_TEXT = ", ".join(f"{value:g}" for value in DEFAULT_SENSITIVITY_RHOS)
DEFAULT_MINIMUM_SHARED_TEXT = str(DEFAULT_MINIMUM_SHARED)


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
        "pooled_date_order": form.get("pooled_date_order", "auto"),
        "overdispersion_rho_values": form.get("overdispersion_rho_values", DEFAULT_RHO_TEXT),
        "grouping_mode": form.get("grouping_mode", PROFILES_MODE),
        "grouping_minimum_shared": form.get("grouping_minimum_shared", DEFAULT_MINIMUM_SHARED_TEXT),
        "genetic_code": str(form.get("genetic_code", "")).strip()[:3],
    }


@app.route("/health")
def health():
    """Let the hosting service check that the website is responding."""

    return "OK", 200


def error_page(message, status):
    """A small page that says what happened and how to go on, even if rendering is failing."""

    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Integrated Genomic '
            f'Surveillance</title></head><body style="font-family:sans-serif;max-width:40rem;margin:3rem auto;'
            f'padding:0 1rem"><h1>Analysis could not be completed</h1><p>{message}</p>'
            '<p><a href="/">Back to the start page</a></p></body></html>', status)


@app.errorhandler(413)
def upload_too_large(error):
    """Show a clear size-limit message instead of a server error page."""

    return error_page("Upload too large. Files in one request must total less than 10 MB.", 413)


@app.errorhandler(500)
def internal_error(error):
    """Keep unexpected internal details out of the public response."""

    return error_page("The analysis could not be completed. Please check your files and try again. "
                      "The server is still running; the details were written to its log.", 500)


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

        try:
            pubtator_status = check_pubtator()
        except Exception:
            app.logger.warning("PubTator diagnostic failed unexpectedly")
            pubtator_status = {"status": "Check failed", "message": "The PubTator check could not finish."}

        results = {"mafft": mafft_status, "geonames": geonames_status,
                   "pubtator": pubtator_status}
        STATUS_CACHE["checked_at"] = now
        STATUS_CACHE["results"] = results
        return results, False


@app.route("/system-status", methods=["POST"])
def system_status():
    """Check external dependencies only when the visitor presses the button."""

    results, cached = get_system_status()
    # The AI line only reads server settings: no request is sent to Claude here.
    results = {**results, "ai": ai_status()}
    values = default_values({})
    return render_template(
        "index.html", values=values, system_status=results, status_cached=cached, step=1
    )


@app.route("/system-status/claude", methods=["POST"])
def claude_status():
    """Check Claude API access only when explicitly requested (token count; no text generated)."""

    results, cached = get_system_status()
    results = {**results, "ai": ai_status(),
               "claude_check": check_claude(request.remote_addr, request.form.get("access_code"))}
    return render_template("index.html", values=default_values({}), system_status=results,
                           status_cached=cached, step=1)


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
                    csv_bytes = safe_csv(organized)
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
                    csv_bytes = safe_csv(display)
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
        if not enough_overlap(pair):
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


def packed(files):
    """Result files compressed for storage: tables shrink about tenfold on a small server."""

    return {name: (zlib.compress(content, 6), mimetype) for name, (content, mimetype) in files.items()}


def store_result_files(files):
    """Keep a few recent result files in memory (compressed) for ten minutes; no database."""

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
        RESULT_DOWNLOADS[download_id] = (now, packed(files))
    return download_id


def build_report_download(download_id, filename):
    """Build the XLSX workbook or PDF report from the stored downloads, once, and cache it.

    Runs only when someone asks for it, one at a time with the analyses, so a
    small server never builds a report and runs an analysis at the same moment.
    Repeated clicks on the same file wait for the first build and then reuse it
    (final polish 2B), instead of building it again or being told to wait.
    """

    with RESULT_LOCK:
        build_lock = REPORT_BUILDS.setdefault((download_id, filename), Lock())
    if not build_lock.acquire(timeout=REPORT_WAIT_SECONDS):
        return error_page("This file is still being prepared for an earlier click. Wait a minute, then click "
                          "the link again; the CSV files are available now.", 503)
    try:
        with RESULT_LOCK:
            item = RESULT_DOWNLOADS.get(download_id)
            cached = item[1].get(filename) if item is not None else None
        if cached is not None:  # an earlier click finished building it
            return send_file(BytesIO(zlib.decompress(cached[0])), mimetype=cached[1], as_attachment=True,
                             download_name=filename)
        return build_report_once(download_id, filename)
    finally:
        build_lock.release()
        with RESULT_LOCK:
            if REPORT_BUILDS.get((download_id, filename)) is build_lock and not build_lock.locked():
                del REPORT_BUILDS[(download_id, filename)]


def build_report_once(download_id, filename):
    """Build one report file from the stored downloads (see build_report_download)."""

    with RESULT_LOCK:
        item = RESULT_DOWNLOADS.get(download_id)
        if item is None or monotonic() - item[0] > RESULT_DOWNLOAD_SECONDS or SUMMARY_NAME not in item[1]:
            return error_page(RESULTS_GONE, 404)
        stored = dict(item[1])
    if not ANALYSIS_SLOT.acquire(blocking=False):
        return error_page("Another analysis is running on this server. Wait a minute, then click the link "
                          "again; the CSV files are available now.", 503)
    try:
        content, _ = build_from_stored(
            filename, lambda name: zlib.decompress(stored[name][0]) if name in stored else None)
    except MemoryError:
        return error_page("The file could not be built because this server ran out of memory. The CSV "
                          "downloads hold the same values.", 503)
    except Exception:
        app.logger.exception("Building %s failed", filename)
        return error_page("The file could not be built. The CSV downloads hold the same values; the details "
                          "were written to the server log.", 500)
    finally:
        ANALYSIS_SLOT.release()
        release_memory()
    mimetype = REPORT_MIMETYPES[filename]
    if len(content) <= 8 * 1024 * 1024:
        with RESULT_LOCK:
            if download_id in RESULT_DOWNLOADS:
                RESULT_DOWNLOADS[download_id][1][filename] = (zlib.compress(content, 6), mimetype)
    return send_file(BytesIO(content), mimetype=mimetype, as_attachment=True, download_name=filename)


@app.route("/results/<download_id>/<filename>")
def result_download(download_id, filename):
    """Download one table or figure from a recent integrated analysis."""

    with RESULT_LOCK:
        item = RESULT_DOWNLOADS.get(download_id)
        # The workbook and the report are built the first time they are asked for (Phase 4A).
        on_demand = item is not None and filename in REPORT_MIMETYPES and filename not in item[1]
    if on_demand:
        return build_report_download(download_id, filename)
    with RESULT_LOCK:
        item = RESULT_DOWNLOADS.get(download_id)
        if (item is None or monotonic() - item[0] > RESULT_DOWNLOAD_SECONDS
                or filename not in item[1]):
            return error_page(RESULTS_GONE, 404)
        content, mimetype = item[1][filename]
    return send_file(BytesIO(zlib.decompress(content)), mimetype=mimetype, as_attachment=True,
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
        related_pairs = read_pair_table(related_path)
        enriched_pairs = read_pair_table(enriched_path)
        membership = pd.read_csv(membership_path, converters={name: str for name in TEXT_LABEL_COLUMNS})
        assignments = pd.read_csv(
            region_assignment_path(metadata_path),
            converters={name: str for name in TEXT_LABEL_COLUMNS},
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
        plain_copy(raw, path, display_name=name)
    finally:
        raw.unlink(missing_ok=True)
    return path, name


def save_inputs(folder, use_example, use_integrated, fasta_upload, metadata_upload, uploads,
                use_literature=False):
    """Copy example files or save uploads; return (paths, display names)."""

    paths = {}
    names = {}
    if use_literature:
        for key, name in LITERATURE_EXAMPLE_FILES.items():
            paths[key] = folder / name
            copyfile(LITERATURE_FOLDER / name, paths[key])
            names[key] = name
        return paths, names
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
        "analysis_plan": None, "declared_pathogen": "", "enrichment": None, "refresh_url": None,
        "is_literature_example": False, "intended_mode": "auto", "mode_check": None,
        "previous_files": [], "info": None, "enrichment_report": None, "enrichment_preview": None,
        "enrichment_note": None, "server_busy": False,
        "enrich": False, "publication_ids": "",
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
    if values["pooled_date_order"] not in DATE_ORDERS:
        return None, None, "Choose how ambiguous day/month collection dates should be read."
    try:
        parse_rho_values(values["overdispersion_rho_values"])
        parse_grouping_mode(values["grouping_mode"])
        parse_minimum_shared(values["grouping_minimum_shared"])
        parse_genetic_code(values["genetic_code"])
    except ValueError as problem:
        return None, None, str(problem)
    return thresholds, date_window, None


def relationship_tables(relationship, thresholds):
    """The Sprint 2 pairwise results that the dashboard turns into figures."""

    return {"pairs": relationship.get("pairs"), "thresholds": thresholds,
            "candidates": relationship["candidates"], "clusters": relationship["clusters"]}


def run_analyses(paths, names, values, thresholds, date_window, temporary_folder, page,
                 pooled_overrides=None, outcome=None):
    """Run every analysis the saved inputs allow, then offer literature enrichment.

    `pooled_overrides` carries approved Sprint 4 values for pooled samples;
    `outcome`, when given, receives the Sprint 2 and integrated results.
    """

    outcome = {} if outcome is None else outcome
    outcome["pooled_overrides"] = pooled_overrides
    # One analysis at a time: two at once could exceed a small server's memory.
    if not ANALYSIS_SLOT.acquire(blocking=False):
        page["server_busy"] = True
        page["error"] = ("Another analysis is running on this server. Your files are still here: "
                         "wait a minute, then try again.")
        return page
    try:
        run_genomic_analyses(paths, names, values, thresholds, date_window, temporary_folder, page, outcome)
        attach_enrichment(page, paths, thresholds, outcome)
    except MemoryError:
        app.logger.exception("Analysis stopped: out of memory")
        page["error"] = OUT_OF_MEMORY
    except Exception as problem:
        app.logger.exception("Analysis stopped at an unexpected error")
        page["error"] = (f"The analysis stopped because of an unexpected server error "
                         f"({type(problem).__name__}). The details were written to the server log.")
    finally:
        ANALYSIS_SLOT.release()
        release_memory()
    return page


def import_details(review, names):
    """Phase 4A: the import record, matching report, and readiness from the review (no calculation)."""

    return {"provenance": names.pop("_provenance", []),
            "readiness": review.get("analysis_readiness"),
            "matching": review.get("matching"),
            "conflicts": review.get("metadata_conflicts")}


def release_memory():
    """Return freed memory to the system after a run (glibc keeps it otherwise)."""

    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except (OSError, AttributeError):
        pass  # not Linux/glibc


def read_raw_table(path):
    """The uploaded table as text (for identifier columns that validation drops), or None."""

    if path is None:
        return None
    try:
        return read_table(path)
    except Exception:
        return None


def attach_enrichment(page, paths, thresholds, outcome):
    """Create an optional enrichment session from the finished analysis (never fails the run)."""

    relationship, result = outcome.get("relationship"), outcome.get("result")
    individual = pairs = pooled = None
    if relationship is not None and relationship["complete"]:
        individual = relationship["analysis_data"].get("samples")
        pairs = relationship.get("pairs")
    if result is not None and result.get("pooled") is not None:
        pooled = result["pooled"]["samples"]
    if individual is None and pooled is None:
        return
    try:
        session = new_session(individual, pooled, pairs, thresholds,
                              synthetic=page["is_literature_example"],
                              label=", ".join(str(name) for name in (page["input_names"] or {}).values()),
                              individual_metadata=read_raw_table(paths.get("metadata")),
                              pooled_metadata=read_raw_table(paths.get("pooled_metadata")))
        page["enrichment"] = view(store_session(session), session)
    except Exception:
        app.logger.exception("Could not prepare literature enrichment")


def run_genomic_analyses(paths, names, values, thresholds, date_window, temporary_folder, page, outcome):
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
        outcome["relationship"] = relationship

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
            annotation_path=paths.get("annotation"),
            settings={"qc_mode": values["qc_mode"], "date_window_days": date_window,
                      "date_order": values["pooled_date_order"],
                      "overdispersion_rhos": parse_rho_values(values["overdispersion_rho_values"]),
                      "grouping_mode": parse_grouping_mode(values["grouping_mode"]),
                      "grouping_minimum_shared": parse_minimum_shared(values["grouping_minimum_shared"]),
                      "genetic_code": parse_genetic_code(values["genetic_code"])},
            input_labels=names,
            relationship=relationship_tables(relationship, thresholds) if completed else None,
            input_notes=page.get("input_notes"),
            pooled_overrides=outcome.get("pooled_overrides"),
            import_details=page.get("import_details"),
        )
    except (ValueError, RuntimeError) as problem:
        page["integrated_error"] = str(problem).replace(str(temporary_folder), "[temporary upload]")
    except MemoryError:
        app.logger.exception("Integrated analysis stopped: out of memory")
        page["integrated_error"] = OUT_OF_MEMORY
    except Exception:
        app.logger.exception("Integrated analysis failed")
        page["integrated_error"] = (
            "The integrated analysis could not be completed. Check the reference and "
            "pooled-sample files."
        )
    else:
        outcome["result"] = result
        files, skipped = download_files(result, relationship_files)
        page["download_id"] = store_result_files(files)
        page["dashboard"] = build_dashboard(result, page["download_id"], set(files))
        page["downloads"] = sorted(files)
        page["skipped_downloads"] = skipped
        return page
    if relationship_files:
        files = {name: (escape_csv_bytes(content), "text/csv") for name, content in relationship_files.items()}
        page["download_id"] = store_result_files(files)
        page["downloads"] = sorted(files)
    return page


def render_page(values, page):
    if page.get("review"):
        toolbar_files = [entry["info"]["name"] for entry in page["review"]["files"]]
    else:
        toolbar_files = list((page.get("input_names") or {}).values()) or page.get("previous_files", [])
    return render_template("index.html", values=values, input_labels=INPUT_LABELS,
                           role_labels=ROLE_LABELS, toolbar_files=toolbar_files[:MAX_GUIDED_FILES],
                           **page)


EXAMPLE_FACTS = {}


@app.context_processor
def shared_page_values():
    """Values every page uses: the glossary, analysis modes, and example facts."""

    if not EXAMPLE_FACTS:
        EXAMPLE_FACTS.update(integrated_example_facts())
    return {"glossary": GLOSSARY, "modes": MODES, "example_facts": EXAMPLE_FACTS or None,
            "upload_fields": UPLOAD_FIELDS,
            # Phase 4A review controls.
            "wide_layout": WIDE_LAYOUT, "base_count_layout": BASE_COUNT_LAYOUT,
            "frequency_scales": FREQUENCY_SCALES, "position_bases": POSITION_BASES,
            "depth_definitions": DEPTH_DEFINITIONS, "resolutions": RESOLUTIONS,
            "genetic_codes": GENETIC_CODES}


def integrated_example_facts():
    """Counts read from the bundled integrated example files (not hard-coded)."""

    try:
        sequences = (INTEGRATED_FOLDER / INTEGRATED_EXAMPLE_FILES["fasta"]).read_text(encoding="utf-8")
        reference = (INTEGRATED_FOLDER / INTEGRATED_EXAMPLE_FILES["reference"]).read_text(encoding="utf-8")
        pooled = pd.read_csv(INTEGRATED_FOLDER / INTEGRATED_EXAMPLE_FILES["snv"])
    except (OSError, ValueError):
        return {}
    bases = "".join(line.strip() for line in reference.splitlines() if not line.startswith(">"))
    return {"individual": sequences.count(">"), "pooled": int(pooled["Sample_ID"].nunique()),
            "sites": int(pooled["Location_Name"].nunique()), "reference_length": len(bases)}


def mode_check(review, intended):
    """Which analysis the reviewed files support, and what the chosen mode would still need."""

    roles = set(review["by_role"])
    individual = {"fasta", "metadata"} <= roles
    if individual and "snv" in roles and "reference" in roles:
        detected = MODES["combined"]["label"]
    elif individual and "snv" in roles:
        detected = "Individual and pooled analyses, run separately (add a shared reference FASTA to connect them)"
    elif individual:
        detected = MODES["individual"]["label"]
    elif "snv" in roles:
        detected = MODES["pooled"]["label"]
    elif "reference" in roles:
        detected = "Reference genome only (add sequences or a pooled SNV table)"
    else:
        detected = "No analysis yet (no sequences or pooled SNV table recognized)"
    if intended not in MODE_NEEDS:
        return {"detected": detected, "intended": None, "missing": []}
    missing = [label for role, label in MODE_NEEDS[intended] if role not in roles]
    return {"detected": detected, "intended": MODES[intended]["label"], "missing": missing}


def intended_mode(form):
    choice = form.get("intended_mode", "auto")
    return choice if choice in MODES else "auto"


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
        use_literature = "run_literature_example" in request.form
        page["is_synthetic_example"] = use_example or use_literature
        page["is_integrated_example"] = use_integrated
        page["is_literature_example"] = use_literature
        uploads = {}
        fasta_upload = None
        metadata_upload = None
        if use_example or use_integrated or use_literature:
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
        elif use_literature and not all(
            (LITERATURE_FOLDER / name).is_file() for name in LITERATURE_EXAMPLE_FILES.values()
        ):
            error = "The bundled synthetic literature example files are missing."
        elif not (use_example or use_integrated or use_literature) and not (has_fasta or has_metadata or "snv" in uploads):
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

        if error is None and use_integrated and request.form.get("enrich_example") == "on":
            # Offline Sprint 4 demonstration: the example files go through the
            # same review and enrichment steps as uploads, before the analysis.
            token = stage_example("combined", INTEGRATED_FOLDER, list(INTEGRATED_EXAMPLE_FILES.values()),
                                  "integrated")
            item = staged_item(token)
            review = review_dataset(item["files"], item["inspections"], {}, {}, {})
            form = {**request.form.to_dict(), "intended_mode": "combined"}
            return start_enrichment(item, token, form, values, {}, review, example="integrated")

        if error is None:
            # Every request gets its own folder, which is deleted automatically.
            with TemporaryDirectory() as temporary_folder:
                folder = Path(temporary_folder)
                try:
                    paths, names = save_inputs(folder, use_example, use_integrated,
                                               fasta_upload, metadata_upload, uploads,
                                               use_literature)
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
    """Delete expired staged uploads and keep at most STAGED_LIMIT sets.

    Folders left behind by an earlier server process (not in STAGED) are
    deleted once they are older than STAGED_SECONDS.
    """

    now = monotonic()
    for token in [key for key, item in STAGED.items() if now - item["created"] > STAGED_SECONDS]:
        remove_staged(token)
    while len(STAGED) >= STAGED_LIMIT:
        remove_staged(min(STAGED, key=lambda key: STAGED[key]["created"]))
    tracked = {Path(item["folder"]).resolve() for item in STAGED.values()}
    for folder in STAGING_ROOT.glob("upload_*"):
        try:
            old = time() - folder.stat().st_mtime > STAGED_SECONDS
        except OSError:
            continue
        if folder.is_dir() and old and folder.resolve() not in tracked:
            rmtree(folder, ignore_errors=True)


PLAIN_EXTENSIONS = {FASTA_KIND: ".fasta", VCF_KIND: ".vcf", EXCEL_KIND: ".xlsx"}


# Role-specific upload fields shown for a chosen analysis mode, in page order.
ROLE_FIELDS = [field["role"] for field in UPLOAD_FIELDS]


def uploads_from_form(files):
    """(upload, declared role or None) for the automatic box and every role field."""

    uploads = [(upload, None) for upload in files.getlist("files") if upload and upload.filename]
    for role in ROLE_FIELDS:
        uploads += [(upload, role) for upload in files.getlist(f"file_{role}") if upload and upload.filename]
    return uploads


def stage_uploads(uploads, enrich=False, publication_ids="", synthetic=None):
    """Save, decompress, and inspect the uploaded files; return a staging token.

    `uploads` holds (file, declared role or None); a declared role is kept
    unless the file's content rules it out.
    """

    if len(uploads) > MAX_GUIDED_FILES:
        raise ValueError(f"Choose at most {MAX_GUIDED_FILES} files at a time.")
    STAGING_ROOT.mkdir(parents=True, exist_ok=True)
    folder = Path(mkdtemp(prefix="upload_", dir=STAGING_ROOT))
    staged = {}
    inspections = []
    # One budget for all files in the request, so many small .gz files cannot
    # together expand to gigabytes on disk.
    remaining = MAX_DECOMPRESSED_BYTES
    try:
        for index, (upload, _) in enumerate(uploads):
            name = Path(upload.filename).name
            raw = folder / f"raw_{index}"
            upload.save(raw)
            plain = folder / f"file_{index}"
            try:
                plain_copy(raw, plain, display_name=name, limit=max(remaining, 0))
                remaining -= plain.stat().st_size
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
    for info, (_, role) in zip(inspections, uploads):
        if role is not None:
            apply_declared_role(info, role)
    with STAGED_LOCK:
        clean_staged()
        token = token_urlsafe(16)
        STAGED[token] = {"folder": folder, "created": monotonic(), "files": staged,
                         "inspections": inspections, "enrich": enrich,
                         "publication_ids": publication_ids[:500], "synthetic": synthetic}
    return token


def apply_sheet_choices(item, form):
    """Read a workbook again when the user selected another worksheet on the review page."""

    changed = False
    for index, info in enumerate(item["inspections"]):
        chosen = form.get(f"sheet_{index}", "")
        sheet_changed = bool(chosen and chosen in info.get("sheets", []) and (
            chosen != info.get("sheet") or info.get("sheet_choice_needed")))
        options = read_file_options(form, index)
        header = options.get("header_row")
        by_name = bool(options.get("sample_from_filename"))
        reading_changed = (header != info.get("header_row_choice")
                           or by_name != bool(info.get("sample_from_filename")))
        has_reading_fields = any(key.startswith(f"opt_{index}_") for key in form)
        if sheet_changed or (has_reading_fields and reading_changed):
            fresh = inspect_file(item["files"][index], info["name"], chosen if sheet_changed else info.get("sheet"),
                                 header if has_reading_fields else info.get("header_row_choice"),
                                 by_name if has_reading_fields else info.get("sample_from_filename"))
            fresh["declared_role"] = info.get("declared_role")
            item["inspections"][index] = fresh
            changed = True
    if changed:
        suggest_roles(item["inspections"])
        for info in item["inspections"]:
            if info.get("declared_role"):
                apply_declared_role(info, info["declared_role"])


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
               "complete_genomes": form.get("genomic_region_mode") == "complete_genomes",
               "intended_mode": intended_mode(form),
               "date_order": (form.get("pooled_date_order") if form.get("pooled_date_order") in DATE_ORDERS
                              else "auto"),
               # Phase 4A: how to read each table, and how to handle metadata conflicts.
               "files": {index: read_file_options(form, index) for index in range(len(inspections))},
               "conflict_resolution": (form.get("conflict_resolution")
                                       if form.get("conflict_resolution") in RESOLUTIONS else ""),
               "allow_unmatched_metadata": form.get("allow_unmatched_metadata") == "on",
               "grouping_mode": form.get("grouping_mode", "profiles"),
               "genetic_code": str(form.get("genetic_code", "")).strip()[:3]}
    return roles, choices, options


def staged_item(token):
    with STAGED_LOCK:
        item = STAGED.get(token or "")
        if item is None or monotonic() - item["created"] > STAGED_SECONDS:
            return None
        return item


def use_profile(item, upload, roles, choices, options):
    """Apply a saved mapping profile to the files it fits; return updated choices and a report."""

    try:
        profile = read_profile(upload.read(64 * 1024 + 1))
    except ValueError as problem:
        return roles, choices, options, [str(problem)]
    # A profile's worksheet is selected first, so its columns can be compared.
    for entry in profile.get("files", []):
        sheet = entry.get("sheet")
        for index, info in enumerate(item["inspections"]):
            if sheet and sheet in info.get("sheets", []) and info.get("sheet") != sheet:
                fresh = inspect_file(item["files"][index], info["name"], sheet)
                fresh["declared_role"] = info.get("declared_role")
                item["inspections"][index] = fresh
    profile_roles, profile_choices, file_options, global_options, report = apply_profile(
        profile, item["inspections"])
    for index, values in file_options.items():
        info = item["inspections"][index]
        if (values.get("header_row") != info.get("header_row_choice")
                or bool(values.get("sample_from_filename")) != bool(info.get("sample_from_filename"))):
            fresh = inspect_file(item["files"][index], info["name"], info.get("sheet"), values.get("header_row"),
                                 values.get("sample_from_filename", False))
            fresh["declared_role"] = info.get("declared_role")
            item["inspections"][index] = fresh
    roles = {**roles, **profile_roles}
    choices = {**choices, **profile_choices}
    options = {**options, **global_options, "files": {**options.get("files", {}), **file_options}}
    return roles, choices, options, report


@app.route("/review/profile", methods=["POST"])
def download_mapping_profile():
    """Download the current review's mapping as a JSON profile (no records, no sample values)."""

    item = staged_item(request.form.get("token"))
    if item is None:
        abort(404)
    roles, choices, options = review_choices(request.form, item["inspections"])
    review = review_dataset(item["files"], item["inspections"], roles, choices, options)
    return send_file(BytesIO(profile_json(review, options)), mimetype="application/json", as_attachment=True,
                     download_name="mapping_profile.json")


@app.route("/review", methods=["POST"])
def review_upload():
    """Step 2: show what each file is and which analyses the dataset allows."""

    values = default_values(request.form)
    page = new_page()
    page["step"] = 2
    page["intended_mode"] = intended_mode(request.form)
    token = request.form.get("token")
    uploads = uploads_from_form(request.files)
    try:
        if uploads:
            token = stage_uploads(uploads, request.form.get("enrich") == "on",
                                  request.form.get("publication_ids", ""))
        item = staged_item(token)
        if item is None:
            raise ValueError("Choose one or more files to upload. Uploaded files are kept for "
                             "30 minutes while you review them.")
        if not uploads and "enrich_shown" in request.form:
            item["enrich"] = request.form.get("enrich") == "on"
            item["publication_ids"] = request.form.get("publication_ids", "")[:500]
        apply_sheet_choices(item, request.form)
        roles, choices, options = review_choices(request.form, item["inspections"])
        profile = request.files.get("mapping_profile")
        if profile is not None and profile.filename and not uploads:
            roles, choices, options, page["profile_report"] = use_profile(item, profile, roles, choices, options)
        page["review"] = review_dataset(item["files"], item["inspections"], roles, choices, options)
        page["mode_check"] = mode_check(page["review"], page["intended_mode"])
        page["token"] = token
        page["conflict_resolution"] = options.get("conflict_resolution", "")
        page["allow_unmatched_metadata"] = options.get("allow_unmatched_metadata")
        page["declared_pathogen"] = options["declared_pathogen"]
        page["enrich"], page["publication_ids"] = item["enrich"], item["publication_ids"]
        if item["enrich"] and page["review"]["ready"]:
            page["enrichment_preview"] = enrichment_preview(item, page["review"], options)
    except ValueError as problem:
        page["error"] = str(problem)
        page["step"] = 1
    except Exception:
        app.logger.exception("Reviewing the uploaded files failed")
        page["error"] = ("The uploaded files could not be reviewed. Check that each file is a "
                         "valid FASTA, CSV/TSV, XLSX, or VCF file and try again.")
        page["step"] = 1
    return render_page(values, page)


# ----------------------------------------------------------------------------
# Optional Sprint 4 before the analysis (upload-page toggle)
# ----------------------------------------------------------------------------

INTEGRATED_CORPUS = INTEGRATED_FOLDER / "synthetic_integrated_articles.json"
INTEGRATED_MOCK = INTEGRATED_FOLDER / "synthetic_integrated_mock_responses.json"
FILES_GONE = ("The uploaded files are no longer on the server (they are kept for at most 30 minutes "
              "and never after an analysis). Please upload them again.")


def enrichment_preview(item, review, options):
    """Metadata gaps and linkable identifiers shown on the review page (nothing is sent anywhere)."""

    try:
        with TemporaryDirectory() as folder:
            paths, _ = prepare_inputs(item["files"], review, folder, options)
            session = preview_session(paths, parse_identifiers(item["publication_ids"]))
    except Exception:
        app.logger.exception("Could not preview metadata gaps")
        return {"error": "The metadata gaps could not be listed; enrichment can still be skipped."}
    counts = {}
    for row in session["gaps"]:
        for gap in row["Gaps"]:
            counts[gap] = counts.get(gap, 0) + 1
    return {"samples": len(session["gaps"]), "eligible": len(session["eligible"]), "counts": counts,
            "rows": [row for row in session["gaps"] if row["Gaps"]][:15],
            "dataset_identifiers": session["dataset_identifiers"]}


def start_enrichment(item, token, form, values, options, review, example=None):
    """Find and check publications at once (free) and show the findings.

    Coming back from the file review reopens the same findings and any Claude
    proposals, unless the files' roles, columns, or identifiers changed.
    """

    signature = repr((sorted(form.items()), item["publication_ids"]))
    previous = item.get("enrichment")
    if previous and previous["signature"] == signature and get_session(previous["token"]) is not None:
        return literature_page(previous["token"], get_session(previous["token"]))
    with TemporaryDirectory() as folder:
        paths, names = prepare_inputs(item["files"], review, folder, options)
        session = preview_session(paths, parse_identifiers(item["publication_ids"]),
                                  synthetic=bool(item.get("synthetic")),
                                  label=", ".join(str(name) for key, name in names.items()
                                                  if not key.startswith("_")))
    session.update(staged_token=token, form=form, example=example)
    if item.get("synthetic") == "integrated":
        session.update(corpus=INTEGRATED_CORPUS, mock=INTEGRATED_MOCK)
    session_token = store_session(session)
    item["enrichment"] = {"signature": signature, "token": session_token}
    try:
        find(session, None)
    except Exception:
        app.logger.exception("Literature search failed")
        return literature_page(session_token, session, "The publication search could not be completed. "
                                                       "You can try again, or continue with the original metadata.")
    return literature_page(session_token, session)


@app.route("/enrich/start", methods=["POST"])
def enrich_start():
    """Sprint 4 before the analysis: list metadata gaps and offer literature discovery."""

    values = default_values(request.form)
    page = new_page()
    token = request.form.get("token")
    item = staged_item(token)
    if item is None:
        page["error"] = FILES_GONE
        return render_page(values, page)
    roles, choices, options = review_choices(request.form, item["inspections"])
    review = review_dataset(item["files"], item["inspections"], roles, choices, options)
    thresholds, date_window, error = check_settings(values)
    if error is None and not review["ready"]:
        error = "Resolve the problems listed in Step 2 before continuing."
    if error is not None:
        page.update(error=error, review=review, token=token, step=2, enrich=item["enrich"],
                    publication_ids=item["publication_ids"], declared_pathogen=options["declared_pathogen"])
        return render_page(values, page)
    try:
        return start_enrichment(item, token, request.form.to_dict(), values, options, review)
    except Exception:
        app.logger.exception("Starting metadata enrichment failed")
        page.update(error="Metadata enrichment could not be prepared; you can run the analysis without it.",
                    review=review, token=token, step=2)
        return render_page(values, page)


def stage_example(intended, folder, files, synthetic):
    """Stage bundled example files as if they had been uploaded (for the enrichment demo)."""

    uploads = []
    handles = [(folder / name).open("rb") for name in files]
    try:
        uploads = [(FileStorage(stream=handle, filename=name), None) for handle, name in zip(handles, files)]
        return stage_uploads(uploads, enrich=True, synthetic=synthetic)
    finally:
        for handle in handles:
            handle.close()


@app.route("/enrich/run", methods=["POST"])
def enrich_run():
    """Apply approved values (or none) and run the analysis on the reviewed files."""

    session_token = request.form.get("token")
    session = get_session(session_token)
    if session is None:
        return expired_page()
    item = staged_item(session.get("staged_token"))
    if item is None:
        page = new_page()
        page["error"] = FILES_GONE
        return render_page(default_values({}), page)
    if session["run"].get("state") == "running":
        return literature_page(session_token, session, "An analysis is running; wait for it or stop it first.")
    skip = request.form.get("skip") == "1"
    derived, log = None, []
    if not skip and session["facts"]:
        record_decisions(session, {key[len("review_"):]: value for key, value in request.form.items()
                                   if key.startswith("review_")})
        derived, log = approved_changes(session)
    form = session["form"]
    values = default_values(form)
    page = new_page()
    page.update(intended_mode=intended_mode(form), is_integrated_example=session.get("example") == "integrated",
                is_synthetic_example=session["synthetic"])
    roles, choices, options = review_choices(form, item["inspections"])
    review = review_dataset(item["files"], item["inspections"], roles, choices, options)
    if "fasta" not in review["by_role"]:
        values["genomic_region_mode"] = "automatic"
    thresholds, date_window, error = check_settings(values)
    outcome = {}
    if error is None:
        with TemporaryDirectory() as temporary_folder:
            try:
                paths, names = prepare_inputs(item["files"], review, temporary_folder, options)
                overrides = None
                if derived is not None:
                    if "metadata" in paths:
                        write_individual_metadata(paths["metadata"], derived)
                    overrides = pooled_overrides(derived) or None
            except ValueError as problem:
                error = str(problem).replace(str(temporary_folder), "[temporary upload]")
            else:
                page["input_notes"] = names.pop("_notes", [])
                page["import_details"] = import_details(review, names)
                run_analyses(paths, names, values, thresholds, date_window, temporary_folder, page,
                             pooled_overrides=overrides, outcome=outcome)
    if error is not None:
        page["error"] = error
    if page.get("server_busy"):  # nothing ran: keep the files and the findings for another try
        return literature_page(session_token, session, page["error"])
    with STAGED_LOCK:
        remove_staged(session["staged_token"])  # uploads are deleted after the analysis, as before
    session["stage"] = "analysis finished"
    if derived is not None and error is None:
        page["enrichment_report"] = enrichment_report(session, derived, log, outcome, thresholds, page)
    elif error is None and session["publications"] is not None:
        page["enrichment_note"] = review_summary(session)["headline"]
    page["analysis_plan"] = review["analyses"]
    return render_page(values, page)


def enrichment_report(session, derived, log, outcome, thresholds, page):
    """What Sprint 4 contributed to this analysis, with evidence downloads."""

    relationship = outcome.get("relationship")
    run_samples = relationship["analysis_data"].get("samples") if relationship and relationship["complete"] else None
    pairs = relationship.get("pairs") if relationship and relationship["complete"] else None
    comparison, screening = before_after(session, derived, run_samples, pairs, thresholds)
    sections = sections_using(derived)
    files = enrichment_downloads(session)
    files["literature_enriched_metadata.csv"] = (safe_csv(derived), "text/csv")
    if comparison is not None:
        files["literature_before_after_pairs.csv"] = (safe_csv(comparison), "text/csv")
    download_id = page.get("download_id")
    with RESULT_LOCK:
        if download_id in RESULT_DOWNLOADS:
            RESULT_DOWNLOADS[download_id][1].update(packed(files))
            page["downloads"] = sorted(RESULT_DOWNLOADS[download_id][1])
    by_id = {fact["Fact_ID"]: fact for fact in session["facts"]}
    changes = [{"sample": item["Sample_ID"], "field": item["Field"].replace("_", " "),
                "original": by_id.get(item["Fact_ID"], {}).get("Original_Value", ""), "value": item["Value"],
                "article": by_id.get(item["Fact_ID"], {}).get("Article_ID", ""),
                "quote": by_id.get(item["Fact_ID"], {}).get("Evidence_Quote", ""), "outcome": item["Outcome"]}
               for item in log]
    changed_pairs = []
    if comparison is not None:
        changed_pairs = comparison[comparison["Change"] != "Unchanged"].fillna("").astype(str).to_dict("records")
    session["applied_count"] = sum(item["Outcome"].startswith("Applied") for item in log)
    return {"summary": contribution(session, derived, log, sections), "changes": changes,
            "per_sample": enrichment_summary(session, log),
            "resolution": gap_resolution(session, derived),
            "screening": screening, "changed_pairs": changed_pairs[:40],
            "changed_pair_total": len(changed_pairs), "files": sorted(files),
            "download_id": download_id, "synthetic": session["synthetic"]}


@app.route("/back", methods=["POST"])
def go_back():
    """Back to the upload or review step, keeping settings but never the deleted files.

    A review can only be shown again while its staged files still exist
    (before the analysis runs, and for at most 30 minutes). After an analysis
    the uploads are gone, so the visitor is asked to select them again.
    """

    values = default_values(request.form)
    page = new_page()
    page["intended_mode"] = intended_mode(request.form)
    page["previous_files"] = [name for name in request.form.getlist("previous_file") if name][:MAX_GUIDED_FILES]
    token = request.form.get("token")
    item = staged_item(token) if token else None
    if request.form.get("target") == "review":
        if item is not None:
            roles, choices, options = review_choices({"intended_mode": page["intended_mode"]}, item["inspections"])
            page["review"] = review_dataset(item["files"], item["inspections"], roles, choices, options)
            page["mode_check"] = mode_check(page["review"], page["intended_mode"])
            page.update(step=2, token=token, previous_files=[])
            return render_page(values, page)
        if request.form.get("example"):
            page["previous_files"] = []
            page["info"] = ("The example used bundled files, so there is no file review to return to. "
                            "Run the example again from the examples card.")
        else:
            page["info"] = ("Your uploaded files are no longer on the server: they are deleted as soon as "
                            "an analysis runs (and after 30 minutes at the latest), and they are never "
                            "kept. Select them again below to return to the file review.")
    elif item is not None:
        with STAGED_LOCK:
            remove_staged(token)  # leaving the review discards its staged files
        page["info"] = "The files you were reviewing were removed from the server. Select files again below."
    return render_page(values, page)


@app.route("/examples/<name>")
def example_file(name):
    """Download one bundled synthetic example file as a format example."""

    path = EXAMPLE_DOWNLOADS.get(name)
    if path is None or not path.is_file():
        abort(404)
    return send_file(path, as_attachment=True, download_name=name)


@app.route("/run", methods=["POST"])
def run_reviewed():
    """Step 3: run the applicable analyses on the reviewed files."""

    values = default_values(request.form)
    page = new_page()
    page["intended_mode"] = intended_mode(request.form)
    token = request.form.get("token")
    item = staged_item(token)
    if item is None:
        page["error"] = "The uploaded files have expired. Please upload them again."
        return render_page(values, page)
    roles, choices, options = review_choices(request.form, item["inspections"])
    try:
        review = review_dataset(item["files"], item["inspections"], roles, choices, options)
    except ServerBusy as busy:
        page.update(error=str(busy), token=token, step=1)
        return render_page(values, page)
    except Exception:
        app.logger.exception("Reviewing the uploaded files failed")
        page["error"] = "The uploaded files could not be reviewed. Please upload them again."
        return render_page(values, page)
    if "fasta" not in review["by_role"]:
        # Genomic-region settings only apply to individual sequences.
        values["genomic_region_mode"] = "automatic"
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
            page["error"] = str(problem).replace(str(temporary_folder), "[temporary upload]")
        except Exception:
            app.logger.exception("Preparing the uploaded files failed")
            page["error"] = "The uploaded files could not be prepared for analysis."
        else:
            page["input_notes"] = names.pop("_notes", [])
            page["import_details"] = import_details(review, names)
            run_analyses(paths, names, values, thresholds, date_window, temporary_folder, page)
    if page.get("server_busy"):  # nothing ran, so the staged files are kept for another try
        page.update(review=review, token=token, step=2, declared_pathogen=options["declared_pathogen"])
        return render_page(values, page)
    failure = page.get("integrated_error") or page.get("error")
    if failure and not page.get("dashboard") and not page.get("analysis_complete"):
        # Nothing was produced, so the uploads are kept and the review is shown again with the
        # reason (audit P1-3); the visitor can change a choice instead of uploading everything again.
        page.pop("integrated_error", None)
        page.update(error=f"The analysis could not be completed: {failure} Your files are still on the server "
                          "for up to 30 minutes; change the choices below or upload corrected files.",
                    review=review, token=token, step=2, declared_pathogen=options["declared_pathogen"])
        return render_page(values, page)
    # The analysis has read the files; the uploads are deleted now (only the
    # results stay, in memory, for the download links).
    with STAGED_LOCK:
        remove_staged(token)
    page["analysis_plan"] = review["analyses"]
    return render_page(values, page)


# ----------------------------------------------------------------------------
# Optional Step 4: literature and metadata enrichment
# ----------------------------------------------------------------------------

EXPIRED_ENRICHMENT = ("This enrichment session has expired (sessions last 30 minutes). "
                      "Run the analysis again to start a new one.")


def literature_page(token, session, error=None, confirm=False):
    """Render the enrichment section on its own results page."""

    page = new_page()
    pre_analysis = session.get("stage") == "pre-analysis"
    page.update(step=4, error=error, is_synthetic_example=session["synthetic"],
                is_literature_example=session["synthetic"] and not session.get("example"))
    if pre_analysis:
        page["token"] = session.get("staged_token")  # Back to file review keeps the staged files
    provider, explanation = (choose_provider(session, run_settings(session)) if session["plan"]
                             else (None, ""))
    state = session["run"].get("state")
    if state not in {None, "not started", "running"} and not session.get("download_id"):
        files = enrichment_downloads(session)
        session["download_id"], session["download_files"] = store_result_files(files), list(files)
    page["enrichment"] = view(token, session, provider, explanation, confirm,
                              session.get("download_id"), session.get("download_files", []))
    if state == "running":
        page["refresh_url"] = url_for("literature_status", token=token)
    return render_page(default_values({}), page)


def expired_page():
    page = new_page()
    page["error"] = EXPIRED_ENRICHMENT
    return render_page(default_values({}), page)


@app.route("/literature/find", methods=["POST"])
def literature_find():
    """Step A: find publications linked to the selected samples (free public services, no AI)."""

    token = request.form.get("token")
    session = get_session(token)
    if session is None:
        return expired_page()
    if session["run"].get("state") == "running":
        return literature_page(token, session, "An analysis is running; wait for it or cancel it first.")
    try:
        find(session, request.form.getlist("sample") or None)
    except Exception:
        app.logger.exception("Literature search failed")
        return literature_page(token, session, "The publication search could not be completed. You can "
                                               "try again, or continue with the original metadata.")
    session["download_id"] = None
    return literature_page(token, session)


@app.route("/literature/analyze", methods=["POST"])
def literature_analyze():
    """Step B: show the cost estimate, then run only after explicit approval."""

    token = request.form.get("token")
    session = get_session(token)
    if session is None:
        return expired_page()
    if not session["plan"]:
        # No usable evidence is an expected result, explained in the section itself, not an error.
        return literature_page(token, session, None if session.get("no_evidence") else
                               "No article passages are ready for analysis. Run 'Find Relevant Literature' first.")
    if session["run"].get("state") == "running":
        return redirect(url_for("literature_status", token=token))
    batch = chosen_batch(session, request.form.get("batch"))
    try:
        provider, explanation, estimate = prepare(session, batch)
    except ProviderError as problem:
        app.logger.warning("Preparing the literature analysis failed: %s", problem)
        return literature_page(token, session, f"The analysis could not be prepared: {problem} Nothing was "
                                               "sent for extraction and nothing was billed.")
    except Exception as problem:
        app.logger.exception("Preparing the literature analysis failed")
        return literature_page(token, session, (
            f"The analysis could not be prepared because of an unexpected server error "
            f"({type(problem).__name__}). Nothing was sent to Claude and nothing was billed; details "
            "are in the server log. You can continue with the original metadata."))
    if provider is None:
        return literature_page(token, session, explanation)
    if request.form.get("confirm") != "1":
        return literature_page(token, session, confirm=True)
    if provider.is_live:
        problem = None
        if request.form.get("acknowledge") != "on":
            problem = "Tick the box to confirm that this analysis uses a paid API before starting it."
        elif not live_request_allowed(request.remote_addr, request.form.get("access_code")):
            problem = ("Paid analysis can only be started on the server computer itself or with the "
                       "access code set in AI_ACCESS_CODE.")
        elif not estimate["fits"]:
            problem = ("The planned analysis could exceed the configured request, token, or spending "
                       "limits. Select fewer samples, or raise the limits on the server.")
        if problem:
            return literature_page(token, session, problem, confirm=True)
    session["download_id"] = None
    start(session, provider, inline=app.config.get("LITERATURE_INLINE_JOBS", False), batch=batch)
    return redirect(url_for("literature_status", token=token))


@app.route("/literature/session/<token>")
def literature_status(token):
    """Progress while the analysis runs; the review table when it is finished."""

    session = get_session(token)
    if session is None:
        return expired_page()
    return literature_page(token, session)


@app.route("/literature/cancel", methods=["POST"])
def literature_cancel():
    """Stop launching new requests; the request already in progress finishes."""

    token = request.form.get("token")
    session = get_session(token)
    if session is None:
        return expired_page()
    if session.get("cancel") is not None:
        session["cancel"].set()
    return redirect(url_for("literature_status", token=token))


@app.route("/literature/apply", methods=["POST"])
def literature_apply():
    """Steps C-D: save review decisions, apply approved values, and recompute affected results."""

    token = request.form.get("token")
    session = get_session(token)
    if session is None:
        return expired_page()
    if session["run"].get("state") == "running" or not session["facts"]:
        return literature_page(token, session, "There are no finished proposals to review.")
    choices = {key[len("review_"):]: value for key, value in request.form.items()
               if key.startswith("review_")}
    try:
        review_and_apply(session, choices)
    except Exception:
        app.logger.exception("Applying literature enrichment failed")
        return literature_page(token, session, "The approved changes could not be applied.")
    session["download_id"] = None
    return literature_page(token, session)


if __name__ == "__main__":
    app.run(debug=False)
