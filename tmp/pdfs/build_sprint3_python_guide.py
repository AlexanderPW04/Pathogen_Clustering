"""Build a presentation guide from the current, uncommitted Sprint 3 Python files."""

from html import escape
from pathlib import Path
from textwrap import wrap

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    HRFlowable,
    KeepTogether,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
)


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "output" / "pdf" / "BMI461_Sprint3_New_Python_Files_Code_Guide.pdf"
pdfmetrics.registerFont(TTFont("Segoe", r"C:\Windows\Fonts\segoeui.ttf"))
pdfmetrics.registerFont(TTFont("SegoeBold", r"C:\Windows\Fonts\segoeuib.ttf"))
pdfmetrics.registerFont(TTFont("Consolas", r"C:\Windows\Fonts\consola.ttf"))
pdfmetrics.registerFontFamily("Segoe", normal="Segoe", bold="SegoeBold")

NAVY = colors.HexColor("#16334c")
TEAL = colors.HexColor("#167b83")
GRAY = colors.HexColor("#4b5b67")
PALE = colors.HexColor("#e7f1f2")

styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name="TitleX", fontName="SegoeBold", fontSize=24, leading=29, textColor=NAVY, spaceAfter=12))
styles.add(ParagraphStyle(name="SubtitleX", fontName="Segoe", fontSize=11, leading=16, textColor=GRAY, spaceAfter=12))
styles.add(ParagraphStyle(name="H1X", fontName="SegoeBold", fontSize=16, leading=20, textColor=NAVY, spaceBefore=16, spaceAfter=9))
styles.add(ParagraphStyle(name="H2X", fontName="SegoeBold", fontSize=11.5, leading=15, textColor=TEAL, spaceBefore=10, spaceAfter=6))
styles.add(ParagraphStyle(name="BodyX", fontName="Segoe", fontSize=9.3, leading=13.4, textColor=NAVY, spaceAfter=6))
styles.add(ParagraphStyle(name="NoteX", fontName="Segoe", fontSize=8.8, leading=12.5, textColor=NAVY, leftIndent=9, firstLineIndent=-9, spaceAfter=5))
styles.add(ParagraphStyle(name="SmallX", fontName="Segoe", fontSize=8.2, leading=11.3, textColor=GRAY, spaceAfter=4))
styles.add(ParagraphStyle(name="CalloutX", fontName="Segoe", fontSize=9.5, leading=14, textColor=NAVY, backColor=PALE, borderPadding=10, spaceBefore=8, spaceAfter=10))
styles.add(ParagraphStyle(name="CodeX", fontName="Consolas", fontSize=6.65, leading=9.15, textColor=NAVY, spaceAfter=0))
styles.add(ParagraphStyle(name="CoverFootX", fontName="Segoe", fontSize=8.5, leading=12, textColor=GRAY, alignment=TA_CENTER))

story = []


def p(text, style="BodyX"):
    story.append(Paragraph(text, styles[style]))


def h1(text):
    p(escape(text), "H1X")


def h2(text):
    p(escape(text), "H2X")


def notes(items):
    for label, explanation in items:
        p(f"<b>{escape(str(label))}</b>  {escape(explanation)}", "NoteX")


def source_listing(relative_path):
    lines = (ROOT / relative_path).read_text(encoding="utf-8").splitlines()
    h2(f"Full numbered source: {relative_path}")
    p("The numbers are the physical source-file line numbers. A wrapped continuation keeps the same original line number in the file; it is not a new Python statement.", "SmallX")
    chunk = []
    for number, line in enumerate(lines, 1):
        # Consolas is narrow enough to fit about 105 characters in the available width.
        pieces = wrap(line.expandtabs(4), width=107, replace_whitespace=False, drop_whitespace=False, break_long_words=True) or [""]
        for part_number, part in enumerate(pieces):
            label = f"{number:>3} " if part_number == 0 else "    "
            chunk.append(Paragraph(escape(label + part).replace(" ", "&nbsp;"), styles["CodeX"]))
        if len(chunk) >= 38:
            story.extend(chunk)
            chunk = []
    story.extend(chunk)


def page(canvas, doc):
    canvas.saveState()
    width, height = letter
    canvas.setFillColor(NAVY)
    canvas.rect(0, height - 0.19 * inch, width, 0.19 * inch, fill=1, stroke=0)
    canvas.setFont("Segoe", 8)
    canvas.setFillColor(GRAY)
    canvas.drawString(0.72 * inch, 0.45 * inch, "BMI 461 | Sprint 3 code walkthrough | classroom reference")
    canvas.drawRightString(width - 0.72 * inch, 0.45 * inch, f"Page {doc.page}")
    canvas.restoreState()


p("BMI 461 Sprint 3", "SubtitleX")
p("New Python Files:<br/>Line-by-Line Code Review", "TitleX")
p("A professor-ready guide to the processed wastewater SNV mode and the separate environmental-measurement exploration", "SubtitleX")
story.append(HRFlowable(width="100%", thickness=1.5, color=TEAL, spaceAfter=18))
p("Prepared from the current project files on 19 October 2026. This documents the working-tree snapshot, including files not yet committed. It changes no scientific or website code.", "CalloutX")
h2("What this guide covers")
notes([
    ("New production files", "wastewater/__init__.py, wastewater/processed_snv.py, and environmental_analysis.py."),
    ("New test files", "tests/test_wastewater.py and tests/test_environmental_analysis.py."),
    ("Existing file touched for integration", "The new routes and imports in app.py are summarized, but app.py is not called a new file."),
    ("How to read it", "Start with the conceptual flow and function map. For each file, read the annotated line ranges beside the full numbered source listing. Blank lines only separate ideas; indentation and standalone closing brackets end Python blocks or expressions."),
])
h2("One-minute project explanation")
p("Sprint 2 analyzes individual pathogen genomes using FASTA, metadata, MAFFT, genomic distances, dates, geography, and potential-cluster thresholds. The new wastewater mode accepts <b>already processed SNV calls</b> and reports what alleles were observed in pooled wastewater. The environmental mode is a separate descriptive view for already processed measurements. Neither new mode runs the Sprint 2 clustering pipeline, and neither processes raw sequencing reads.")
p("An SNV call means a reported alternative nucleotide at a reference position. Its allele frequency describes the fraction of reads supporting that alternative within the processed sample. A missing SNV row is <b>not</b> proof that a mutation was absent: sequencing coverage at that position may be unknown.", "CalloutX")

story.append(PageBreak())
h1("How the pieces connect")
notes([
    ("Wastewater upload", "Browser POST /wastewater -> app.py saves input to a temporary folder -> load_processed_snvs reads, maps, joins, and validates it -> summarize_snvs counts observed calls -> Flask renders wastewater.html and offers a short-lived CSV download."),
    ("Environmental upload", "Browser POST /environmental -> app.py saves processed measurements -> load_measurements validates them -> summarize_measurements builds tables -> make_charts makes descriptive PNG charts -> Flask renders environmental.html and offers CSV download."),
    ("What stays separate", "Neither path calls run_pipeline, MAFFT, GeoNames, p-distance, geographic outbreak thresholds, or connected-component clustering. The existing Sprint 2 home page and scientific calculations remain conceptually separate."),
])
h2("Python and pandas vocabulary used repeatedly")
notes([
    ("DataFrame", "A pandas table: columns are fields and rows are records. Series is one column-like vector."),
    ("df[\"X\"]", "Select column X; assigning to it replaces or creates that column."),
    (".loc[rows, columns]", "Select by row labels or a Boolean condition; .iloc selects by numeric position."),
    (".isna() / .notna()", "Find missing / present values. Missing is different from numerical zero."),
    (".groupby(...).agg(...) ", "Partition rows by key and calculate summaries such as counts or means."),
    (".nunique() vs .size()", "Count distinct values versus count all rows, including repeated rows."),
    ("errors=\"coerce\"", "Turn unparseable numbers or dates into missing values so validation can count them."),
    ("ValueError", "A readable input-data problem shown to the user; errors are collected before raising so several issues can be reported together."),
    ("with ... as ...", "A context manager: cleans up temporary folders or temporarily replaces a function in a test."),
])
h2("Function map")
notes([
    ("read_table(path)", "Read CSV/XLSX processed SNV calls; identify the published workbook's two-row introduction."),
    ("load_processed_snvs(path, metadata_path=None)", "Convert source fields into a common schema, optionally join sample metadata, validate, retain duplicates with a warning, and return data plus warnings."),
    ("summarize_snvs(data)", "Count observed wastewater samples, rows, unique SNVs, dates, locations, frequencies, depths, top SNVs, and sample-level totals."),
    ("load_measurements(path, metric=\"auto\")", "Read descriptive environmental measurements, standardize field names, validate sample/feature rows, and choose the numeric measurement."),
    ("summarize_measurements(data, metric)", "Build sample, time, location, missing-data, and quality summaries without replacing missing measurements by zero."),
    ("available_value_columns(data)", "Identify present numeric environmental columns."),
    ("table_records(table)", "Convert a pandas table into browser-friendly dictionaries and date strings."),
    ("make_charts(data, metric, tables)", "Make small descriptive plots in memory as PNG bytes, for up to three features."),
])

story.append(PageBreak())
h1("1. wastewater/__init__.py")
p("<b>Purpose:</b> This one-line file makes <i>wastewater</i> a Python package and states its boundary. The period of biological interpretation is important: pooled wastewater observations are kept separate from individual-genome analysis.")
notes([
    ("Line 1", "The triple-quoted module docstring describes the package. It contains no executable analysis function. Python executes this file when importing the wastewater package; here the effect is only documentation/package organization."),
])
source_listing("wastewater/__init__.py")

story.append(PageBreak())
h1("2. wastewater/processed_snv.py")
p("<b>Purpose:</b> This is the actual processed-SNV input and summary module. It does not align sequences or assign lineages. Its input is a row-per-observed-SNV table, either the identified publication supplement or a standardized CSV/XLSX.")
h2("Imports and fixed column names")
notes([
    ("Lines 1-5", "Line 1 is the module description. Path from pathlib (line 3) handles file paths and extensions on Windows/Linux; pandas (line 5) supplies table loading, conversion, validation, grouping, and summarization. Blank lines are visual separators."),
    ("Lines 8-17", "PAPER_COLUMNS maps exact Supplementary Table 2 field names to this module's internal names. Each quoted left side is a source header; each right side is the standardized destination. SAMPLE becomes Sample_ID; COLLECTION_DATE becomes Collection_Date; POSITION, REF and ALT are genomic call fields; ALT_FREQ is allele frequency; TOTAL_DP is read depth. This mapping is a dictionary, not a biological inference."),
    ("Lines 18-22", "REQUIRED_COLUMNS lists the minimum standardized call fields. OUTPUT_COLUMNS adds optional Collection_Date and Location_Name, fixing the order of returned/downloaded columns. Read_Depth must have a column but individual cells may be missing with a warning."),
])
h2("read_table(path): exact line walk")
notes([
    ("25-26", "Define the reader and explain its aim: sample identifiers that look numeric must stay text, because IDs are labels, not quantities."),
    ("28-30", "Convert path to a Path object. A case-insensitive .csv suffix selects pandas.read_csv. dtype=str preserves IDs and original text; keep_default_na=False prevents strings such as NA from being silently interpreted as missing."),
    ("31-35", "For .xlsx, inspect the first two rows without assigning headers. If row 2 starts with SAMPLE, use row 2 as the header; otherwise use row 1. Then read the whole worksheet as strings. This accommodates the paper's introductory prose row."),
    ("36-42", "Read the publication description in first-row column E and record a provenance flag only when the expected table label, 52-sample description, SARS-CoV-2 text, and shifted header all match. The AND conditions are a guard against falsely assigning the paper's reference to an unrelated table."),
    ("43-44", "Return the Excel table. Any extension other than .csv or .xlsx causes a clear ValueError; the function does not guess a format."),
])
h2("load_processed_snvs(path, metadata_path=None): setup and published mapping")
notes([
    ("47-48", "Define the main loader. The optional second path supplies sample metadata when collection date or location is not in the SNV table. Return value is (standardized DataFrame, list of warnings)."),
    ("50-57", "Read the SNV table, trim header whitespace, and decide whether it is the specifically identified published Table 2. The fingerprint flag, expected paper columns, and STATE/LOCATION must all be present. Start an empty warning list."),
    ("59-63", "If the published source is recognized, rename headers, assign the article's documented reference identity NC_045512.2, combine STATE and LOCATION into one readable place label, and parse dates like 30_May_2020. This reference value is deliberately limited to the positively identified publication workbook."),
    ("64-72", "Otherwise require the standardized fields, name missing fields in an error, add empty optional date/location columns, and parse dates generically. An unrelated file must explicitly provide Reference_Genome; it is never guessed from similar-looking columns."),
    ("73-83", "Flag nonempty source dates that could not be parsed; reject an empty table; trim whitespace from IDs, reference names, bases, and location; uppercase Ref/Alt so lowercase source nucleotides validate; initialize a list for errors."),
])
h2("Optional metadata join and conflict detection")
notes([
    ("85-92", "Only execute this block when a metadata file was supplied. Read/trim its columns and Sample_ID values. Require a Sample_ID column with one unique, nonblank ID per metadata row; otherwise the join would be ambiguous."),
    ("93-99", "Compute SNV sample IDs not found in metadata and stop if any are missing. Metadata IDs with no SNV rows produce a warning, not an absent-mutation claim. Set Sample_ID as the metadata lookup index."),
    ("100-103", "For each optional metadata field (date or location), skip it if not supplied. Map each SNV row's Sample_ID to the corresponding metadata value."),
    ("104-112", "For dates, parse metadata text, report bad nonblank dates, report any disagreement between two known dates, then fill only missing SNV dates from metadata. Existing valid SNV dates are not overwritten."),
    ("113-117", "For locations, compare nonblank source and metadata names case-insensitively. Report genuine conflicts. Fill a blank SNV location from metadata but preserve an already supplied location."),
    ("119-125", "Add the earlier bad-source-date count, then count blank IDs, reference names, and Ref/Alt nucleotides. Accumulating errors gives the user more than one useful correction per run."),
])
h2("Biological and numeric validation")
notes([
    ("127-133", "Require Ref and Alt each to be exactly one ordinary A, C, G, or T base; reject ambiguous symbols here. Ref and Alt must differ for an SNV call. A base is a nucleotide, not an amino acid."),
    ("135-138", "Convert Position to a number. It must exist, be greater than zero, equal its rounded value (integer), and not be infinite. Genomic coordinates are positive reference positions."),
    ("139-142", "Convert Alt_Frequency to numeric and require a nonmissing fraction from 0 through 1, inclusive. The code does not automatically change percentages such as 25 into 0.25."),
    ("143-153", "Trim depth text, permit an empty value, and convert other values to numbers. Non-numeric text and non-positive/non-integer/infinite depths are errors; missing depths produce a warning. Crucially, missing depth is not treated as zero coverage."),
    ("155-161", "If any validation errors accumulated, stop with a readable combined message. Otherwise cast position to integer, frequency to float, depth to pandas nullable Int64, and store the parsed dates. Nullable Int64 keeps whole-number depths while allowing missing values."),
])
h2("Within-sample consistency, duplicate warnings, return")
notes([
    ("162-168", "For each sample, check that all known dates agree and all known locations agree. Ignore missing values during the disagreement count. If consistent, propagate a known value to rows of that same sample. Restore blank text for truly unknown locations."),
    ("170-173", "A repeated record is defined by sample ID, reference, position, Ref, and Alt. Count later occurrences as duplicates and warn, but retain all rows because the published workbook itself has repeated call keys (some due to different annotations)."),
    ("174-179", "Count records with missing dates and blank locations and add transparent warnings; no rows are removed for these missing metadata fields."),
    ("181", "Return a copy with only standardized output columns in a fixed order plus the warning list. The copy avoids accidentally mutating the original intermediate DataFrame."),
])
h2("summarize_snvs(data): observed-call summaries")
notes([
    ("184-187", "Define the summary function. The variant identity combines reference genome, genomic position, reference base, and alternative base; position alone is not enough."),
    ("188-193", "Group calls by that variant key. Samples_Observed counts distinct sample IDs, whereas Records counts rows (including repeats). Sort by sample frequency, then row count, then position, and show the first ten."),
    ("194-202", "Group by Sample_ID and take its established first date and location, count call rows, and sort the sample table by ID."),
    ("203-205", "Remove repeated identical calls within each sample for a separate Distinct_SNVs count. Map that count back to each sample; format dates as YYYY-MM-DD for display."),
    ("206-224", "Exclude missing dates/depths from ranges. Return a dictionary of distinct sample count, total records, distinct variant count, date range, listed locations, frequency min/median/max, available depth count/min/median/max, top variants, and per-sample table. The ternary expressions return None when no date or depth exists."),
])
source_listing("wastewater/processed_snv.py")

story.append(PageBreak())
h1("3. environmental_analysis.py")
p("<b>Purpose:</b> This is a separate, descriptive processed-measurement module. It is included because it is also a new Python file in the current Sprint 3 working tree. Do not describe its measurements as wastewater SNV calls or as genomic-distance results.")
h2("Imports, accepted names, and required fields")
notes([
    ("1-6", "The docstring describes the module. BytesIO holds a chart image in memory rather than writing a permanent file; Path handles file paths; pandas reads and summarizes tabular data."),
    ("9", "Every environmental row must identify a sample, collection date field, location field, and feature. Date/location cells may be empty, but the columns must exist."),
    ("10-22", "ALIASES translates reasonable alternative header spellings into a common set. Each dictionary pair is source-style name -> internal field, including measurement value, frequency, quality, and units. This is only header normalization; it does not change scientific units."),
])
h2("load_measurements(path, metric='auto'): file and headers")
notes([
    ("25-26", "Define the loader. metric='auto' chooses a usable numeric column later. One row is intended to represent one sample-feature measurement."),
    ("28-34", "Convert the path to Path. Read UTF-8 CSV with BOM support or XLSX as strings so parsing decisions remain explicit. Reject other extensions rather than guessing."),
    ("36-40", "For every column header, trim, lowercase using casefold, replace hyphens/underscores with spaces, collapse spaces to underscores, then look it up in ALIASES. Set the standardized DataFrame headers."),
    ("41-50", "Reject two source headers that collapse to the same internal field. Require the key columns and at least Frequency or Measurement_Value. Reject empty tables and tables over 10,000 rows to keep this class prototype small."),
    ("52-62", "Trim whitespace from every cell; create blank Quality/Unit columns if not supplied. Start an error list and count blank Sample_ID or Feature_ID cells."),
])
h2("Date, location, and numeric validation")
notes([
    ("64-69", "Remove an Excel midnight time suffix if present, then parse YYYY-MM-DD dates. A nonblank invalid date becomes an error. Store valid dates as pandas datetime values; blank dates remain missing."),
    ("71-78", "For location labels, collapse repeated whitespace and use a case-insensitive dictionary to select a consistent spelling within this upload. This does not geocode places or infer coordinates."),
    ("80-88", "For each available numeric column, turn blank into missing and parse numbers. Reject text/infinity; Frequency specifically must be a fraction in 0-1, not a percentage. Save numeric floating-point values."),
    ("90-98", "Check that a Sample_ID does not have multiple known dates/locations; reject repeated sample-feature pairs so replicates are not silently combined; reject mixed units for one feature because their means would be invalid."),
    ("100-109", "If metric is auto, choose the first available numeric field with a value. Reject a requested missing metric or an all-missing selected metric. Raise accumulated errors, otherwise return rows sorted by sample/feature with a reset index and the chosen metric name."),
])
h2("summarize_measurements(data, metric): descriptive tables")
notes([
    ("112-113", "Define a descriptive summary. The docstring emphasizes that missing measurements do not become zero, since zero and not measured are different claims."),
    ("115-121", "Build one row per Sample_ID: first date and location, number of unique listed features, number of available metric values, and number of listed features lacking a value."),
    ("122-130", "Retain observed numeric values, then rows with dates. Group by feature, location, date to compute a mean and distinct sample count. Sort them, then use shift/diff to show the preceding date and arithmetic change within each feature-location series."),
    ("132-137", "For known locations, calculate feature-by-location mean, minimum, maximum, and number of distinct samples with a value. Blank locations do not become a fabricated place."),
    ("138-142", "Build a missing-data table for required, available numeric, quality, and unit fields. Summarize quality labels, turning blanks into 'Not provided' only for the quality-count display."),
    ("143-158", "Return count totals and the six summary tables. 'features' counts features with an observed numeric value, while 'rows' counts all retained input rows; these describe different denominators."),
])
h2("Small support functions and chart")
notes([
    ("161-164", "available_value_columns returns whichever of Frequency and Measurement_Value actually exist in this DataFrame. It does not choose a scientific threshold."),
    ("167-175", "table_records copies a DataFrame, formats datetime columns as YYYY-MM-DD, changes remaining pandas missing values to empty strings, then makes a list of ordinary dictionaries for the Flask template."),
    ("178-183", "make_charts is descriptive only. It imports Matplotlib inside the function so the plotting package is loaded only for this mode. AutoDateLocator and ConciseDateFormatter format date axes; Figure is the plot container."),
    ("185-191", "Choose up to three observed features. Create a two-column figure with one row per feature. Select time and location summaries for each feature and retrieve that row's two axes."),
    ("192-200", "Plot up to six location time series. If there are no dated values, show an explanatory message; otherwise format the date axis and show the legend."),
    ("201-204", "Likewise, show a message when no location summary exists, or draw a horizontal bar chart for the first six locations."),
    ("205-211", "Find a provided unit for the feature when possible, label frequency as a 0-1 fraction, and set short, readable axis titles. The chart does not convert units or test statistical significance."),
    ("212-214", "Save the figure as PNG into an in-memory BytesIO buffer and return the raw bytes. app.py later turns those bytes into a browser-displayable image."),
])
source_listing("environmental_analysis.py")

story.append(PageBreak())
h1("4. tests/test_wastewater.py")
p("<b>Purpose:</b> These are checks, not analysis code. Each method beginning with <i>test_</i> asserts one expected behavior. A passing test means the particular example behaves as specified, not that every possible biological interpretation is validated.")
h2("Imports and small fixture")
notes([
    ("1-11", "The docstring labels the file. io supplies in-memory upload bytes; tempfile provides disposable directories; unittest supplies TestCase/assertions; Path builds test paths; patch replaces a function temporarily; pandas builds a mock table; the final import selects the two wastewater functions being tested."),
    ("14-18", "STANDARD_ROWS is a tiny standardized CSV string: two samples, the same reference-position C>T call, two frequencies and depths, dates, and sites. Adjacent string literals inside parentheses become one Python string."),
])
h2("WastewaterLoaderTest: what every test proves")
notes([
    ("21-38", "Create a one-row table using actual publication-style headers and a lowercase reference base. Mark it as positively identified Table 2, patch read_table to return it, then assert the ID, reference, parsed date, combined location, uppercase base, depth, and absence of warnings. The patched reader avoids needing the publisher file in the test suite."),
    ("40-49", "Create visually similar columns without the publication provenance flag. Assert that the loader refuses them because Reference_Genome is not provided. This guards against accidentally assigning the paper's reference to unrelated data."),
    ("51-61", "Write a CSV with a duplicate W1 call to a temporary folder. The summary must show 2 unique samples, 3 records, 1 unique SNV, 2 distinct sample observations, and a repeated-record warning. The file disappears when the context exits."),
    ("63-74", "Modify the fixture to include position zero, identical Ref/Alt, frequency greater than one, and nonnumeric depth. Assert one ValueError message contains all four problems, showing that validation reports multiple useful issues together."),
    ("76-88", "Remove dates and sites from the call file, supply them in a second metadata CSV, and assert the join fills them by Sample_ID. Then remove W2 from metadata and assert the missing match is reported."),
    ("90-104", "Supply a call with blank depth/date/location. Assert those remain missing, the summary has no date range/location/depth, and three warnings are returned. This explicitly checks that missing is not silently treated as zero or an absent mutation."),
    ("106-113", "Give W1 conflicting location labels in its call file and metadata. Assert that the loader rejects the contradiction rather than arbitrarily choosing one."),
])
h2("WastewaterFlaskTest: browser-to-backend checks")
notes([
    ("115-140", "Import the existing Flask app inside the test, replace run_pipeline with a mock, and use Flask's test_client to request the home and wastewater pages, upload the two-row CSV, inspect the returned counts/warning text/filename, follow the generated download link, assert run_pipeline was never called, and verify the response is a CSV containing a call row."),
    ("142-151", "Post an incomplete CSV to /wastewater and assert the rendered page contains the readable missing-column error. This checks the user-facing error path, not only the loader exception."),
    ("154-155", "When this file is run directly, unittest.main discovers and runs its test methods. In the normal project test discovery, this guard does not run twice."),
])
source_listing("tests/test_wastewater.py")

story.append(PageBreak())
h1("5. tests/test_environmental_analysis.py")
p("<b>Purpose:</b> This separate test file documents and protects the descriptive environmental-measurement mode. It does not assert wastewater SNV biology or Sprint 2 clustering.")
h2("Imports, reusable rows, and file helper")
notes([
    ("1-11", "The imports have the same test roles as in the wastewater file: in-memory upload, temporary files, unittest assertions, paths, mock patching, pandas workbook construction, and the three functions under test."),
    ("14-21", "CSV_ROWS is a small invented dataset: three sample IDs, two features, dates/sites, a missing frequency, and quality labels. The header deliberately uses aliases such as Sample ID and Feature Identifier to test normalization."),
    ("24-28", "EnvironmentalLoaderTest groups pure analysis tests. write_csv writes CSV_ROWS (or modified rows) into a supplied temporary folder and returns its Path, reducing repeated setup code."),
])
h2("EnvironmentalLoaderTest: each case")
notes([
    ("30-49", "Load the alias-based CSV; assert auto selects Frequency, headers standardize, three samples are counted, one date and one frequency remain missing, time summaries have three groups, missing-data counts are correct, and quality counts cover all five rows."),
    ("51-56", "Change the header to Measurement Value. Assert the alternative numeric metric is selected and numeric values parse. This is a header-support test, not evidence that all measurement values are biologically frequencies."),
    ("58-68", "Create a real temporary XLSX workbook with pandas, load it, and assert its Measurement_Value 12.5 is preserved. This verifies the Excel path without a permanent example workbook."),
    ("70-74", "Replace one numeric value with the text 'unknown' and assert a non-numeric ValueError. The table must not silently discard bad numeric input."),
    ("76-80", "Change a frequency fraction to 25 and assert an error explaining the allowed 0-1 range; the code does not guess that 25 means 25%."),
    ("82-87", "Generate a chart and inspect the first eight PNG signature bytes. This verifies the function returned an image, not its scientific interpretation or visual quality."),
])
h2("EnvironmentalFlaskTest: route behavior")
notes([
    ("90-108", "Import Flask, patch Sprint 2's run_pipeline, GET /environmental, POST the explicitly synthetic example, assert page/result labels and embedded PNG, follow the CSV download link, and ensure the pathogen pipeline was never called."),
    ("110-125", "Upload CSV_ROWS as multipart form data, assert the page returns successfully and shows the filename, three samples, and sample table. The upload is in memory for the test."),
    ("127-136", "Upload an incomplete CSV and assert the user sees a readable missing-column error."),
    ("139-140", "When invoked as a standalone script, unittest.main executes tests; discovery can also import the file without executing this guarded block."),
])
source_listing("tests/test_environmental_analysis.py")

story.append(PageBreak())
h1("Existing Flask integration: app.py additions")
p("The next notes cover only the newly added paths in the existing app.py. They are not a new Python file, but they explain how a browser upload reaches the functions above. The exact line numbers may shift as app.py evolves.")
h2("Wastewater path")
notes([
    ("New imports/constants", "Flask's abort and send_file support downloads; BytesIO wraps in-memory CSV bytes; token_urlsafe makes an unguessable short download token; TemporaryDirectory isolates each upload. The package functions are imported into app.py. A small locked dictionary holds up to four short-lived downloads for ten minutes; this is not a database."),
    ("@app.route('/wastewater')", "GET shows the upload form. POST checks whether the SNV file and optional metadata file have allowed CSV/XLSX extensions, then saves them under safe fixed filenames inside a fresh temporary directory."),
    ("Analysis branch", "The route enforces a 3 MB SNV and 1 MB metadata limit, calls load_processed_snvs, then summarize_snvs. It formats dates for display/CSV, caps the CSV at 4 MB, stores the CSV bytes under a short-lived token, and sends only the first 100 rows to the webpage preview."),
    ("Errors", "ValueError becomes a readable page error. Unexpected exceptions are logged on the server while visitors receive a generic message; raw stack traces are not rendered."),
    ("@app.route('/wastewater/download/<download_id>')", "Looks up the token while holding a lock; absent/expired tokens return 404. Otherwise send_file returns the stored standardized CSV as a download. Uploaded source files have already been removed by TemporaryDirectory."),
])
h2("Environmental path")
notes([
    ("New imports/constants", "environmental_analysis exports loading, summary, chart, and table-to-record helpers. base64 encodes PNG bytes for HTML. A second small locked cache stores recent organized CSV downloads; SYNTHETIC_ENVIRONMENTAL_FILE points to an explicitly invented example."),
    ("@app.route('/environmental')", "GET shows the page. POST either selects the labeled synthetic example or saves a CSV/XLSX upload in a temporary folder; then load_measurements, summarize_measurements, and make_charts perform descriptive analysis."),
    ("Response and download", "Flask converts a few summary tables to browser records, embeds the chart, previews at most 100 measurement rows, and stores a bounded CSV for a short-lived download link. It never calls the individual-genome pipeline."),
])

story.append(PageBreak())
h1("Presentation notes and limits")
h2("Suggested three-minute explanation")
notes([
    ("1. Start with the data", "The wastewater upload is a table of already called SNVs, not raw FASTQ reads. Each row reports sample, reference, position, nucleotide change, allele frequency, and read depth where available."),
    ("2. Explain provenance", "The published supplement has unusual headers and a descriptive first row. read_table recognizes that exact workbook; load_processed_snvs maps it to a documented common schema. An unrelated file must explicitly name its reference genome."),
    ("3. Explain validation", "The loader checks IDs, joins, dates, letters, positions, fractions, depth, and duplicate keys. It reports problems instead of deleting rows; missing coverage is not interpreted as mutation absence."),
    ("4. Explain output", "summarize_snvs counts unique samples separately from call rows and distinct SNVs. Flask displays the summary and table and offers the standardized CSV. No MAFFT or outbreak clustering is applied to pooled wastewater."),
    ("5. Explain tests", "Tiny invented rows test a known mapping, a missing reference, duplicate accounting, bad values, metadata joins, missingness, and browser behavior. External publication data are not required for unit tests."),
])
h2("Questions a professor may ask")
notes([
    ("Why keep Sample_ID as text?", "An identifier such as 003 could lose its leading zeros if parsed as a number. It is a label, not a quantity."),
    ("Why is duplicate_count a warning?", "The publication table contains repeated call keys, sometimes with different annotation context. Retaining source rows is more transparent; summary functions separately use distinct-sample counts where appropriate."),
    ("Why not equate missing rows with absence?", "The input does not prove the position was adequately sequenced in every sample. Without coverage at that position, no-call and true absence cannot be separated."),
    ("Why are Ref and Alt restricted to A/C/G/T?", "This first version handles simple single-nucleotide substitutions. Insertions, deletions, ambiguous calls, and multi-base alleles need an explicit future representation."),
    ("What is not implemented?", "Raw FASTQ processing, BWA mapping, iVar calling, Freyja deconvolution, cross-source matching, individual-genome p-distance applied to mixtures, and outbreak confirmation."),
])
h2("References and where to start in code")
p("Fontenele et al. (2021), <i>Water Research</i>, published full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC8464352/ . The processed per-sample SNV source inspected for this implementation is Supplementary Table 2: https://ars.els-cdn.com/content/image/1-s2.0-S0043135421009040-mmc2.xlsx . Source-code line numbers in this guide refer to the local working-tree snapshot dated 19 October 2026.")
p("For a live class walkthrough, open wastewater/processed_snv.py at lines 25, 47, and 184; then tests/test_wastewater.py at lines 22 and 116; finally app.py at the /wastewater route. Open environmental_analysis.py only if discussing the distinct environmental-measurement mode.", "CalloutX")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
document = SimpleDocTemplate(
    str(OUTPUT), pagesize=letter,
    leftMargin=0.7 * inch, rightMargin=0.7 * inch,
    topMargin=0.62 * inch, bottomMargin=0.68 * inch,
    title="BMI 461 Sprint 3 New Python Files - Code Review Guide",
    author="BMI 461 project documentation",
)
document.build(story, onFirstPage=page, onLaterPages=page)
print(OUTPUT)
