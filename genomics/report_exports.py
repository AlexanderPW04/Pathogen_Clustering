"""Consolidated XLSX workbook and concise PDF report of one analysis (Sprint 5 Phase 4A).

Both are built from the *same* result tables, figures, and findings that the
dashboard shows (no calculation is repeated), so every exported number
matches the website and the CSV downloads.

On the website they are built **on demand**, the first time someone clicks
the link, from the stored CSV and PNG downloads plus a small
``analysis_summary.json`` (findings, notes, readiness, figure titles). This
keeps the analysis request itself as fast as before Phase 4A (the hosted
worker has a 180 s timeout) and guarantees that the workbook equals the CSVs,
because it is read back from them. The command-line runner builds them
directly from the result.

- ``analysis_results.xlsx``: a Summary sheet (findings, table index,
  limitations) and one sheet per key table. Written with openpyxl's
  streaming (write-only) mode. A table longer than WORKBOOK_TABLE_ROWS is
  listed on the Summary sheet as "CSV only" instead of being copied (the CSV
  keeps every row); nothing is cut silently.
- ``analysis_report.pdf``: dataset, main findings, which analyses ran or why
  not, methods, limitations, and the main figures, drawn page by page with
  matplotlib (no new dependency). It describes prototype screening results
  and makes no outbreak, transmission, or lineage claim.
"""

import json
import math
import textwrap
from datetime import date
from io import BytesIO

import matplotlib

matplotlib.use("Agg")
from PIL import Image  # noqa: E402  (installed with matplotlib)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import openpyxl  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

from genomics.result_availability import analysis_availability, distinct_findings  # noqa: E402
from genomics.spreadsheet_safety import escape_text  # noqa: E402
from openpyxl.cell import WriteOnlyCell  # noqa: E402
from openpyxl.styles import Font  # noqa: E402

XLSX_NAME = "analysis_results.xlsx"
PDF_NAME = "analysis_report.pdf"
SUMMARY_NAME = "analysis_summary.json"
MIMETYPES = {XLSX_NAME: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
             PDF_NAME: "application/pdf"}
EXCEL_ROW_LIMIT = 1_048_575  # one header row plus this many data rows
# Detail tables longer than this stay CSV-only (listed on the Summary sheet), which keeps the
# workbook quick to build on a small server; the CSV downloads always hold every row.
WORKBOOK_TABLE_ROWS = 10_000
MAX_FIGURE_PIXELS = 1200  # figures wider or taller are downscaled for the PDF only (PdfPages keeps images until it closes)
SHEETS = [
    ("Readiness", "analysis_readiness.csv"),
    ("Analyses run", "analysis_status.csv"),
    ("Import provenance", "import_provenance.csv"),
    ("Metadata matching", "metadata_matching.csv"),
    ("Samples", "canonical_samples.csv"),
    ("Observations", "canonical_variant_observations.csv"),
    ("Coverage", "canonical_coverage.csv"),
    ("Reference check", "reference_verification.csv"),
    ("Paper QC", "paper_method_reproduction.csv"),
    ("SNVs per sample", "variants_per_sample.csv"),
    ("Variant positions", "variant_positions.csv"),
    ("Temporal trends", "temporal_trends.csv"),
    ("VAF time series", "sprint5_temporal_summary.csv"),
    ("Uncertainty models", "sprint5_uncertainty_models.csv"),
    ("Overdispersion", "sprint5_overdispersion_sensitivity.csv"),
    ("Mutation eligibility", "sprint5_mutation_eligibility.csv"),
    ("Annotated mutations", "sprint5_annotated_mutations.csv"),
    ("Cross-source", "cross_source_comparison.csv"),
    ("Provenance", "provenance.csv"),
]
REPORT_FIGURES = ["variants_per_sample", "qc_outcomes", "depth_distribution", "genome_positions", "vaf_depth",
                  "vaf_time_series", "pcoa", "mutation_annotation_overview", "cross_source", "uncertainty_models"]
MAX_REPORT_FIGURES = 8
NOTE_SECTIONS = ["quality", "temporal", "population", "depth_aware", "uncertainty", "grouping", "annotation",
                 "cross_source", "geography"]
DISCLAIMER = ("Prototype genomic-surveillance screening output. Results describe the uploaded processed data; "
              "they do not confirm outbreaks, transmission, individual strain identity, or lineages. A variant "
              "missing from a pooled sample is unknown unless coverage shows otherwise.")


def cell_value(value):
    """A value Excel can store: NaN/NA become empty cells, numpy scalars become Python numbers."""

    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    try:
        if value != value:  # pandas NA / NaT
            return None
    except (TypeError, ValueError):
        if str(value) in {"<NA>", "NaT"}:
            return None
    text = str(value)
    return text if text not in {"<NA>", "NaT", "nan"} else None


def plain(value):
    """JSON-safe value (numpy scalars become Python numbers)."""

    return value.item() if hasattr(value, "item") else str(value)


def report_summary(result, findings):
    """The small, JSON-serializable part of a result that the workbook and report need."""

    pooled = result.get("pooled")
    return {
        "findings": list(findings),
        "notes": {section: list(result["notes"].get(section, []))[:12] for section in NOTE_SECTIONS
                  if result["notes"].get(section)},
        "readiness": (result.get("import_details") or {}).get("readiness") or [],
        "labels": {key: str(value) for key, value in (result.get("input_labels") or {}).items()
                   if not str(key).startswith("_")},
        "input_format": pooled["input_format"] if pooled is not None else None,
        "canonical": result.get("canonical"),
        "reference": result["summary"].get("reference"),
        "analyses": result.get("analyses", []),
        "figures": {name: {"title": figure["title"], "caption": figure.get("caption", "")}
                    for name, figure in result["figures"].items()},
        # The same availability record the results page uses: title, headline counts, limitations,
        # and which figures are shown first, so the report and the website agree.
        "availability": analysis_availability(result),
    }


def summary_json(result, findings):
    return json.dumps(report_summary(result, findings), default=plain, indent=1).encode("utf-8")


class StoredTables:
    """Result tables read back from stored CSV downloads, only when a sheet needs them."""

    def __init__(self, read):
        self.read = read

    def get(self, name):
        content = self.read(name)
        if content is None:
            return None
        # Only an empty cell is missing: text such as the influenza gene "NA" (neuraminidase) or a
        # sample called "None" stays text, exactly as in the CSV download (final-polish 1B).
        return pd.read_csv(BytesIO(content), float_precision="round_trip", low_memory=False,
                           keep_default_na=False, na_values=[""])


def stored_figures(read, summary):
    figures = {}
    for name, meta in summary.get("figures", {}).items():
        content = read(f"{name}.png")
        if content is not None:
            figures[name] = {"png": content, **meta}
    return figures


def excel_rows(table):
    """Rows as plain Python values (NaN, NA, and infinities become empty cells), converted in bulk."""

    data = table.astype(object)
    data = data.where(data.notna(), None)
    for row in data.itertuples(index=False, name=None):
        # Text that a spreadsheet would run as a formula is stored with a leading apostrophe (audit P1-2).
        yield [escape_text(value) if isinstance(value, str)
               else value if not isinstance(value, float) or math.isfinite(value) else None for value in row]


def build_workbook(tables, summary_values):
    """The consolidated workbook as bytes; returns (bytes, notes about tables left as CSV only).

    `tables` is any mapping with .get(name) (a result's tables or StoredTables).
    """

    workbook = openpyxl.Workbook(write_only=True)
    summary = workbook.create_sheet("Summary")
    bold = Font(bold=True)

    def heading(text):
        cell = WriteOnlyCell(summary, value=escape_text(text))
        cell.font = bold
        summary.append([cell])

    heading("Genomic surveillance analysis results")
    summary.append([f"Created {date.today().isoformat()} by the integrated genomic surveillance application."])
    summary.append([DISCLAIMER])
    summary.append([])
    heading("Main findings (as on the results page)")
    for text in summary_values["findings"]:
        summary.append([escape_text(text)])
    summary.append([])
    heading("Sheets in this workbook")
    summary.append(["Sheet", "Rows", "Source table (CSV download)", "Note"])
    notes, written = [], []
    for sheet_name, table_name in SHEETS:
        table = tables.get(table_name)
        if table is None:
            continue
        rows = len(table)
        if rows > min(WORKBOOK_TABLE_ROWS, EXCEL_ROW_LIMIT):
            summary.append([sheet_name, rows, table_name, "Not copied here (large detail table): use the CSV, "
                                                          "which has every row."])
            notes.append(f"{table_name} ({rows:,} rows) is in the CSV downloads only, not in the workbook.")
            continue
        summary.append([sheet_name, rows, table_name, ""])
        written.append((sheet_name, table))
    summary.append([])
    heading("Limitations recorded during the analysis")
    for section in NOTE_SECTIONS:
        for text in summary_values["notes"].get(section, [])[:12]:
            summary.append([f"[{section}] {text}"])
    for sheet_name, table in written:
        sheet = workbook.create_sheet(sheet_name[:31])
        header = []
        for name in table.columns:
            cell = WriteOnlyCell(sheet, value=escape_text(str(name)))
            cell.font = bold
            header.append(cell)
        sheet.append(header)
        for row in excel_rows(table):
            sheet.append(row)
    buffer = BytesIO()
    workbook.save(buffer)
    return buffer.getvalue(), notes


# ---------------------------------------------------------------------------
# PDF report
# ---------------------------------------------------------------------------

PAGE = (8.27, 11.69)  # A4 portrait, inches


class Page:
    """One text page; `write` wraps lines and opens a new page when this one is full."""

    def __init__(self, pdf, title):
        self.pdf = pdf
        self.title = title
        self.figure = None
        self.y = 0
        self.new()

    def new(self):
        if self.figure is not None:
            self.finish()
        self.figure = plt.figure(figsize=PAGE)
        self.figure.text(0.07, 0.955, self.title, fontsize=14, fontweight="bold", color="#173f5f")
        self.figure.text(0.07, 0.03, "Prototype screening output; see limitations.", fontsize=7, color="#52606d")
        self.y = 0.92

    def write(self, text, size=9, bold=False, indent=0.0, gap=0.004, width=105):
        lines = textwrap.wrap(str(text), width=int(width - indent * 100)) or [""]
        height = size / 72 / PAGE[1] * 1.35
        for line in lines:
            if self.y - height < 0.06:
                self.new()
            self.figure.text(0.07 + indent, self.y, line, fontsize=size,
                             fontweight="bold" if bold else "normal", va="top", color="#1f2937")
            self.y -= height
        self.y -= gap

    def finish(self):
        self.pdf.savefig(self.figure)
        plt.close(self.figure)
        self.figure = None


def dataset_lines(tables, summary):
    lines = [f"{key}: {value}" for key, value in summary["labels"].items()]
    if summary.get("input_format"):
        lines.append(f"Pooled input format: {summary['input_format']}")
    canonical = summary.get("canonical")
    if canonical:
        lines.append(f"Canonical tables: {canonical['samples']} samples, {canonical['observations']} variant "
                     f"observations, {canonical['coverage_records']} coverage records.")
        classes = "; ".join(f"{name}: {count}" for name, count in canonical["measurement_classes"].items())
        lines.append(f"Measurements per observation: {classes}.")
    if summary.get("reference"):
        lines.append(f"Reference: {summary['reference']}")
    references = tables.get("reference_verification.csv")
    if references is not None:
        for row in references.to_dict("records"):
            lines.append(f"Reference {row['Reference_ID']}: {row['Status']} ({row.get('Detail', '')})")
    return lines


def methods_lines(summary):
    ran = {row["Analysis"]: row for row in summary.get("analyses", [])}
    lines = []
    for name, row in ran.items():
        if row["Status"] in {"Completed", "Limited", "Partly completed"}:
            lines.append(f"{name} — {row['Status']}: {row['Detail']}")
    lines.append("Allele frequencies: reported values are kept as published; read-count frequencies are "
                 "alternative reads / total depth. Wilson 95% intervals describe read sampling only; "
                 "beta-binomial intervals use documented technical replicates or labelled hypothetical rho values.")
    lines.append("A missing variant row is 'not reported' (coverage unknown), never 0%, unless a per-position "
                 "table shows adequate coverage.")
    return lines


def report_figure_names(figures, summary):
    """The figures the results page shows first, in report order (older summaries: the fixed list)."""

    availability = summary.get("availability")
    if availability:
        return [name for name in availability["report_figures"] if name in figures]
    return REPORT_FIGURES + [name for name in figures if name not in REPORT_FIGURES]


def build_pdf_report(tables, figures, summary):
    """The concise report as bytes: an executive summary, methods and limitations, then the main figures.

    Long tables are never printed; they are in the CSV downloads and the workbook.
    """

    availability = summary.get("availability") or {}
    buffer = BytesIO()
    with PdfPages(buffer) as pdf:
        page = Page(pdf, "Genomic surveillance analysis report")
        if availability:
            page.write(f"{availability['title']} · {availability['status']}", size=11, bold=True, gap=0.006)
        page.write(f"Created {date.today().isoformat()}.", size=8)
        page.write(DISCLAIMER, size=8.5, gap=0.012)
        if availability:
            page.write("Executive summary", size=11, bold=True)
            for text in availability["key_metrics"]:
                page.write("• " + text, indent=0.02)
            page.y -= 0.01
        page.write("Main findings", size=11, bold=True)
        for text in (distinct_findings(summary["findings"]) if availability else summary["findings"]):
            page.write("• " + text, indent=0.02)
        page.y -= 0.01
        if availability:
            page.write("Important limitations", size=11, bold=True)
            for text in availability["limitations"]:
                page.write("• " + text, indent=0.02, size=8.5)
            page.y -= 0.01
            if availability["unavailable"]:
                page.write("Unavailable or limited analyses", size=11, bold=True)
                for row in availability["unavailable"]:
                    page.write(f"{row['Analysis']} ({row['Status']}): {row['Reason']}", indent=0.02, size=8)
                page.y -= 0.01
            if availability["denominators"]:
                page.write("What each count means", size=11, bold=True)
                for row in availability["denominators"]:
                    page.write(f"{row['label']}: {row['value']:,}. {row['definition']}", indent=0.02, size=8)
                page.y -= 0.01
        page.write("Dataset", size=11, bold=True)
        for line in dataset_lines(tables, summary):
            page.write(line, indent=0.02)
        page.y -= 0.01
        rows = summary.get("readiness") or []
        if rows and not availability:
            page.write("Analyses these data support (decided before the run)", size=11, bold=True)
            for row in rows:
                page.write(f"{row['Status']} — {row['Analysis']}: {row['Reason']}", indent=0.02, size=8.5)
            page.y -= 0.01
        page.write("What ran (methods and results)", size=11, bold=True)
        for line in methods_lines(summary):
            page.write(line, indent=0.02, size=8.5)
        page.y -= 0.01
        provenance = tables.get("import_provenance.csv")
        if provenance is not None and len(provenance):
            changes = provenance[provenance["Decided_By"].isin(["User", "User confirmed",
                                                                 "Importer (user-confirmed minimum)"])]
            if len(changes):
                page.write("Import choices confirmed by the user", size=11, bold=True)
                for row in changes.drop_duplicates(["File", "Item", "Transformation"]).head(30).to_dict("records"):
                    page.write(f"{row['File']}: {row['Item']} — {row['Source']}"
                               + (f" ({row['Transformation']})" if isinstance(row["Transformation"], str)
                                  and row["Transformation"] else ""),
                               indent=0.02, size=8)
                page.y -= 0.01
        page.write("Detailed limitations by analysis", size=11, bold=True)
        for section in NOTE_SECTIONS:
            for text in summary["notes"].get(section, [])[:4 if availability else 6]:
                page.write("• " + text, indent=0.02, size=8)
        page.y -= 0.01
        page.write("Full tables are not printed here: every table is in the CSV downloads, and the key tables are "
                   "also in analysis_results.xlsx (very long detail tables are CSV only).", size=8)
        page.finish()
        shown = 0
        for name in report_figure_names(figures, summary):
            figure = figures.get(name)
            if figure is None or shown >= MAX_REPORT_FIGURES:
                continue
            with Image.open(BytesIO(figure["png"])) as picture:
                picture = picture.convert("RGB")
                picture.thumbnail((MAX_FIGURE_PIXELS, MAX_FIGURE_PIXELS))  # 8-bit, at most MAX_FIGURE_PIXELS
                image = np.asarray(picture)
            height, width = image.shape[:2]
            size = PAGE if width / height < 1.5 else (PAGE[1], PAGE[0])  # wide figures get a landscape page
            canvas = plt.figure(figsize=size)
            canvas.text(0.07, 0.955, figure["title"], fontsize=13, fontweight="bold", color="#173f5f")
            box_width = 0.86
            box_height = min(0.72, box_width * height / width * size[0] / size[1])
            axis = canvas.add_axes([0.07, 0.93 - box_height, box_width, box_height])
            axis.imshow(image, interpolation="none")  # embedded as is, no resampled copy
            axis.axis("off")
            caption = textwrap.wrap(figure.get("caption", ""), 110)
            for number, line in enumerate(caption[:8]):
                canvas.text(0.07, 0.9 - box_height - 0.02 - number * 0.016, line, fontsize=8, va="top",
                            color="#52606d")
            pdf.savefig(canvas)
            plt.close(canvas)
            del image
            shown += 1
    return buffer.getvalue()


def build_report(name, tables, figures, summary):
    """(bytes, notes) for XLSX_NAME or PDF_NAME."""

    if name == XLSX_NAME:
        return build_workbook(tables, summary)
    return build_pdf_report(tables, figures, summary), []


def build_from_stored(name, read):
    """Build XLSX_NAME or PDF_NAME from stored downloads; `read(filename)` returns bytes or None."""

    summary = json.loads(read(SUMMARY_NAME).decode("utf-8"))
    figures = stored_figures(read, summary) if name == PDF_NAME else {}
    return build_report(name, StoredTables(read), figures, summary)


def report_files(result, findings, size_limit=8 * 1024 * 1024):
    """{name: (bytes, mimetype)} for the XLSX and PDF built directly from a result (command line),
    plus messages for anything not produced."""

    files, messages = {}, []
    summary = json.loads(summary_json(result, findings))
    for name in (XLSX_NAME, PDF_NAME):
        mimetype = MIMETYPES[name]

        def build(name=name):
            return build_report(name, result["tables"], result["figures"], summary)

        try:
            content, notes = build()
        except MemoryError:
            messages.append(f"{name} could not be created because this server ran out of memory; the CSV files "
                            "hold the same values.")
            continue
        except Exception as error:  # noqa: BLE001 - reported to the user; the CSVs are unaffected
            messages.append(f"{name} could not be created ({type(error).__name__}); the CSV files hold the same values.")
            continue
        if len(content) > size_limit:
            messages.append(f"{name} is larger than the {size_limit // (1024 * 1024)} MB download limit; use the CSV "
                            "files or run the analysis locally.")
            continue
        files[name] = (content, mimetype)
        messages += notes
    return files, messages
