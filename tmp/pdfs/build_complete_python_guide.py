"""Create a source-checked classroom reference without changing project code."""

import ast
import hashlib
import json
from html import escape
from pathlib import Path
from textwrap import wrap

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, HRFlowable, KeepTogether, LongTable,
    PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle,
)
from reportlab.platypus.tableofcontents import TableOfContents


ROOT = Path(__file__).resolve().parents[2]
DESTINATION = ROOT / "output/pdf/BMI461_Complete_New_Python_Line_by_Line_Guide.pdf"
QA_FOLDER = ROOT / "tmp/pdfs/complete_guide_qa"
NAVY = colors.HexColor("#16334c")
TEAL = colors.HexColor("#087c83")
INK = colors.HexColor("#243746")
GRAY = colors.HexColor("#5d6b75")
PALE = colors.HexColor("#edf5f6")
RULE = colors.HexColor("#dce4e8")

pdfmetrics.registerFont(TTFont("GuideSans", r"C:\Windows\Fonts\segoeui.ttf"))
pdfmetrics.registerFont(TTFont("GuideBold", r"C:\Windows\Fonts\segoeuib.ttf"))
pdfmetrics.registerFont(TTFont("GuideCode", r"C:\Windows\Fonts\consola.ttf"))
pdfmetrics.registerFontFamily("GuideSans", normal="GuideSans", bold="GuideBold")
STYLES = getSampleStyleSheet()
for name, font, size, leading, extra in [
    ("GuideTitle", "GuideBold", 27, 33, {"spaceAfter": 17, "textColor": NAVY}),
    ("GuideDeck", "GuideSans", 12, 18, {"spaceAfter": 12, "textColor": GRAY}),
    ("GuideH1", "GuideBold", 18, 23, {"spaceAfter": 11, "textColor": NAVY}),
    ("GuideH2", "GuideBold", 11.5, 15, {"spaceBefore": 13, "spaceAfter": 7, "textColor": TEAL, "keepWithNext": True}),
    ("GuideBody", "GuideSans", 9.4, 13.9, {"spaceAfter": 7, "textColor": INK}),
    ("GuideSmall", "GuideSans", 8.1, 11.4, {"spaceAfter": 5, "textColor": GRAY}),
    ("GuideCell", "GuideSans", 8.15, 11.25, {"textColor": INK}),
    ("GuideCode", "GuideCode", 7.15, 9.2, {"textColor": NAVY}),
    ("GuideLine", "GuideSans", 7.5, 10, {"textColor": GRAY}),
    ("GuideTableHead", "GuideBold", 8, 11, {"textColor": colors.white}),
    ("GuideTOC", "GuideSans", 10, 17, {"textColor": INK, "spaceBefore": 4}),
]:
    STYLES.add(ParagraphStyle(name=name, fontName=font, fontSize=size, leading=leading, **extra))

FILES = [
    "genomic_regions.py", "align_sequences.py", "geocode_locations.py",
    "system_diagnostics.py", "environmental_analysis.py", "wastewater/__init__.py",
    "wastewater/processed_snv.py", "app.py", "tests/test_alignment.py",
    "tests/test_geocoding.py", "tests/test_diagnostics.py",
    "tests/test_environmental_analysis.py", "tests/test_wastewater.py",
]
SOURCES = {name: (ROOT / name).read_text(encoding="utf-8-sig") for name in FILES}
TREES = {name: ast.parse(source, filename=name) for name, source in SOURCES.items()}
FINGERPRINTS = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in FILES}

# Explanations are editorial content. Source is read directly from the project.
FILE_INFO = {}
FUNCTION_INFO = {}
LINE_NOTES = {}
KEY_NOTES = {}
MISSING_NOTES = []


def note_lines(filename, text):
    """Attach a human explanation to an exact stripped source line."""
    target = LINE_NOTES.setdefault(filename, {})
    for line in text.strip().splitlines():
        source, explanation = line.split(" || ", 1)
        target[source.strip()] = explanation.strip()


IMPORTS = {
    "ast": "Python's syntax-tree tools; documentation tooling can inspect source without executing the application.",
    "subprocess": "Starts the external MAFFT program and captures its exit status and output; Python does not perform the alignment itself.",
    "shutil": "Provides which(), which searches the server's PATH for MAFFT.",
    "StringIO": "Treats in-memory text as a readable file, allowing Biopython to parse captured FASTA output or Flask to collect error text.",
    "BytesIO": "Treats in-memory bytes as a file, used for PNG generation and downloadable CSV responses.",
    "Path": "Builds and inspects file paths in a Windows/Linux-friendly way, including extensions, temporary paths, and reads/writes.",
    "TemporaryDirectory": "Creates a disposable directory and removes its files when the with block exits, including when an exception occurs.",
    "SeqIO": "Biopython's sequence-file reader; parses FASTA records into IDs, sequences, and descriptions.",
    "ALLOWED_SEQUENCE_CHARACTERS": "Reuses the validator's accepted nucleotide/gap alphabet rather than defining a competing rule.",
    "clean_sequence": "Reuses the validator's sequence cleanup for MAFFT output.",
    "find_mafft": "Reuses the pipeline's executable discovery, including Windows PATH/PATHEXT behavior.",
    "pandas": "Provides DataFrames and Series for reading, matching, validating, grouping, converting, and summarizing tabular data. pd is its short alias.",
    "requests": "Sends HTTPS requests to GeoNames and supplies network/HTTP exception types.",
    "math": "Supplies isfinite() so infinities and NaN cannot masquerade as usable coordinates.",
    "os": "Reads GEONAMES_USERNAME from the environment; diagnostic responses never include the configured value.",
    "sys": "Provides sys.stderr for warnings and errors, separate from ordinary progress output.",
    "b64encode": "Turns PNG bytes into ASCII base64 for embedding a generated plot in HTML.",
    "redirect_stderr": "Temporarily captures backend diagnostic text for the Flask response. It changes process-level stderr while active.",
    "isfinite": "Rejects infinite and NaN threshold inputs even if float() can parse them.",
    "token_urlsafe": "Creates a random URL-friendly download identifier so recent CSV downloads are not addressed by predictable numbers.",
    "copyfile": "Copies bundled synthetic inputs into the same temporary workflow used for uploads.",
    "Lock": "Protects shared in-memory dictionaries so threads do not update a cache simultaneously. It does not synchronize separate worker processes.",
    "monotonic": "Measures elapsed seconds with a clock that is not shifted by wall-clock corrections; used for cache expiry.",
    "Flask": "Creates the web application and its route/error-handler registry.",
    "abort": "Stops a request with an HTTP error such as 404 when a download token is absent or expired.",
    "render_template": "Renders a Jinja HTML template using ordinary Python values computed by the backend.",
    "request": "Reads the active browser request, including uploaded files and submitted form fields.",
    "send_file": "Builds a downloadable HTTP response from the in-memory CSV file.",
    "check_geonames": "Imports the separate live GeoNames diagnostic; this is not an analysis or geographic-screening calculation.",
    "check_mafft": "Imports the small synthetic MAFFT diagnostic used by the System Status button.",
    "output_paths": "Reuses the runner's names for its genomic, enriched, related-pair, and membership CSVs.",
    "region_assignment_path": "Reuses the runner's per-sample region-label provenance CSV name.",
    "run_pipeline": "Calls the existing six-stage individual-genome workflow directly rather than reproducing its scientific logic in Flask.",
    "FASTA_EXTENSIONS": "Reuses the validator's accepted .fasta/.fa/.fna file extensions.",
    "METADATA_EXTENSIONS": "Reuses the validator's accepted CSV/XLSX extensions.",
    "load_measurements": "Reads and validates the separate processed environmental measurement table.",
    "make_charts": "Creates descriptive environmental PNG plots in memory.",
    "summarize_measurements": "Calculates descriptive sample, date, location, missingness, and quality summaries.",
    "table_records": "Converts pandas values into a list of simple dictionaries for Jinja display.",
    "load_processed_snvs": "Reads, standardizes, joins, and validates processed wastewater SNV observations.",
    "summarize_snvs": "Counts observed wastewater samples, call rows, and distinct reference-position-allele combinations.",
    "AutoDateLocator": "Chooses practical date tick positions for a descriptive Matplotlib time axis.",
    "ConciseDateFormatter": "Formats those date ticks with compact labels.",
    "Figure": "Creates a Matplotlib figure without opening an interactive desktop plotting window.",
    "io": "Supplies BytesIO in tests so uploads can be simulated without a browser or permanent file.",
    "tempfile": "Supplies temporary folders in tests; test fixtures are removed when the context exits.",
    "unittest": "Python's test framework: TestCase, assertions, test discovery, and the standalone test runner.",
    "patch": "Temporarily replaces a function, object, or environment value during a test and restores it afterward.",
    "Mock": "Makes a controllable stand-in such as a GeoNames HTTP response; no live API request is required.",
    "align_samples": "Imports the production alignment function so tests exercise its grouping and output validation.",
    "compare_samples": "Checks that aligned output remains compatible with the inherited direct p-distance comparison.",
    "fill_missing_coordinates": "Imports the actual production geocoder so tests check coordinate retention, cache use, and fallbacks.",
}

KEY_NOTES.update({
    "label": "The standardized region label used downstream for grouping, or None when unresolved.",
    "segment": "A parsed segment identifier; None means this annotation supplied no segment clue.",
    "genes": "A set of explicit gene aliases used when checking for contradictory annotations.",
    "whole": "A Boolean clue indicating explicit whole/complete/full-genome wording.",
    "status": "A short status label for the website; it carries no credential value.",
    "message": "Readable explanatory text; raw external stderr and API credential details are not returned here.",
    "version": "An optional sanitized MAFFT version token; absence does not invalidate a successful alignment.",
    "Sample_ID": "Sample identifier. It is a label, even when it contains only digits.",
    "Reference_Genome": "Reference identity is part of the SNV key so positions from different references are not collapsed together.",
    "Position": "A positive integer reference position for an observed single-nucleotide change.",
    "Ref": "The reported reference nucleotide, restricted to one A/C/G/T in this processed-SNV loader.",
    "Alt": "The reported alternative nucleotide, which must differ from Ref.",
    "Alt_Frequency": "A reported alternative-allele fraction in the inclusive range 0 to 1, not a percentage string.",
    "Read_Depth": "The reported positive integer depth or a nullable missing value; missing depth is not assumed to be zero.",
    "Collection_Date": "Collection date associated with the sample; unavailable dates remain missing.",
    "Location_Name": "Reported or normalized location text. Its presence alone does not prove exact specimen coordinates.",
    "Pathogen": "Pathogen grouping label from the validated input.",
    "Gene_Segment": "Resolved genomic region; comparisons and alignment use this together with Pathogen.",
    "Frequency": "A processed measurement represented as a fraction between 0 and 1.",
    "Measurement_Value": "An alternative finite numeric measurement whose units come from the uploaded table.",
    "Quality": "A provided quality label counted for transparency, not automatically used to exclude rows.",
    "Unit": "An optional unit label; the loader checks consistency but performs no conversion.",
    "samples": "Number of distinct sample IDs, which is different from the number of rows or pairwise relationships.",
    "records": "Number of retained observed-call rows, including repeated keys.",
    "unique_snvs": "Number of distinct reference/position/reference-base/alternative-base combinations.",
    "samples_with_values": "Distinct samples having a nonmissing value for the chosen environmental metric.",
    "features": "Distinct features with an observed value for the selected metric.",
    "rows": "Number of retained measurement rows; missing numeric values do not automatically remove a row.",
    "tables": "A dictionary of DataFrames or display records used by the environmental page.",
    "counts": "Small count summary returned separately from the detailed environmental tables.",
    "measurements": "The full standardized environmental table, including missing values.",
    "temporal_changes": "Descriptive feature/place/date means and differences from the prior observed date.",
    "location_summary": "Descriptive feature/place means, minima, maxima, and distinct sample counts across observed rows.",
    "missing_data": "Per-column missing-row counts, with the input row count as denominator.",
    "quality_summary": "Counts of supplied quality labels; blank labels display as Not provided.",
    "filename": "The original input filename shown to the visitor, or an explicit synthetic-example label.",
    "metric": "The numeric environmental column selected for summaries and plots.",
    "summary": "The backend summary passed to the HTML template.",
    "preview": "A browser preview of at most 100 rows; the downloadable CSV retains the full organized table.",
    "download_id": "An unguessable token referencing a recent CSV held in this worker's memory.",
    "chart": "Base64 PNG text embedded in the environmental webpage.",
    "warnings": "Nonfatal data-quality messages, such as retained duplicates or unavailable metadata.",
    "metadata_filename": "Optional metadata filename, shown only when a metadata upload was used.",
    "sample_1": "First sample ID in a genomic-temporal candidate relationship.",
    "sample_2": "Second sample ID in that relationship; a pair is not another sample.",
    "p_distance": "Existing backend p-distance copied for presentation; this Flask code does not recalculate it.",
    "days_apart": "Existing date difference copied from the enriched CSV.",
    "pathogen": "Pathogen associated with a pair or membership group.",
    "gene_segment": "The pair's or cluster's resolved genomic-region label.",
    "cluster_id": "Existing backend-generated cluster identifier, such as C1.",
    "max_p_distance": "Prototype maximum genomic mismatch fraction read from the form.",
    "max_days_apart": "Prototype maximum absolute collection-date difference.",
    "max_km_apart": "Prototype geographic-distance cutoff in kilometers.",
    "genomic_region_mode": "User's automatic-detection or explicit-complete-genome choice.",
    "confirm_complete_genomes": "Whether the user explicitly confirmed all sequences are complete genomes of one pathogen.",
})

# FILE_INFO, FUNCTION_INFO, and hand-written line notes are added below.

INTRODUCTION = [
    ("What counts as new", "The change log's baseline already includes validation, direct genetic comparison, time/geography enrichment, clustering, the runner, and Flask. This guide covers seven later production modules: genomic_regions, align_sequences, geocode_locations, system_diagnostics, environmental_analysis, and the two wastewater package files. It also covers five corresponding test files. app.py is an existing file, included in full to explain their integration. The inherited validator, runner, and distance/clustering implementations are described as dependencies rather than mislabeled as new files."),
    ("How every line is accounted for", "Every nonblank physical line in the 13 files has an entry. The left column reproduces its source; the right column explains its role. Long lines wrap within a cell without receiving new line numbers. Blank lines are omitted from the tables because their only role is visual separation. Closing brackets and continuation strings are explained too: they complete a statement rather than perform another calculation."),
    ("A function is a reusable operation", "def creates a named function. Parameters are the values it receives. Indented statements form its body. return sends a value back and ends that call. An exception stops the normal path until a matching except block handles it. Creating a function does not immediately execute its body."),
    ("Imports and where work happens", "Standard-library imports provide paths, temporary files, processes, locks, clocks, and tests. Third-party imports provide pandas tables, Biopython FASTA parsing, requests HTTPS, Flask HTML responses, and Matplotlib PNGs. Project imports reuse another local module's function. Importing a package is different from asking MAFFT to align a file or sending a GeoNames request."),
    ("Python punctuation and indentation", "= assigns a value; == compares values. [] can index a list or select a table column; {} can build a dictionary or set. A dictionary maps keys to values; a set keeps distinct items. A colon begins a block or separates a dictionary key from its value. Parentheses can continue an expression over multiple physical lines. A trailing comma separates entries, not a new calculation."),
    ("Missing values and Boolean masks", "None is Python's absent-value object; pandas also uses NA/NaN/NaT for missing cells and dates. Missing is different from zero. isna/notna build masks. With pandas masks, & means row-wise AND, | means OR, and ~ means NOT; parentheses are necessary. Python and/or instead combine scalar truth values and short-circuit evaluation."),
    ("Frequently used pandas operations", "A DataFrame is a table; a Series is a column-like vector. loc selects rows/columns by label or mask; iloc selects by integer position. groupby partitions rows by keys. agg calculates named summaries; nunique counts distinct values; size counts rows; count counts nonmissing values. astype converts type, copy isolates a result, and errors='coerce' makes bad text missing so validation can report it."),
    ("Tests use small controlled examples", "unittest.TestCase classes organize tests even though production modules avoid custom classes. Fixtures are invented data. patch/Mock replace processes or HTTP requests so tests are repeatable. An assertion compares an observed value to an expected one. assertRaisesRegex uses a regular-expression assertion on error text; this is a testing utility, distinct from the production region parser, which uses string operations."),
]

ARCHITECTURE = [
    ("Individual pathogen genomes", "Browser upload -> app.index -> run_pipeline -> validate_inputs (matching and region resolution) -> geocode_locations (optional approximate coordinates) -> align_sequences (separate MAFFT group per Pathogen + Gene_Segment) -> inherited pairwise p-distance -> time/geography -> threshold/connected-component clusters -> CSVs -> Flask summary. genomic_regions supplies annotation clues to the validator. This guide does not replace the original implementations of those inherited scientific stages."),
    ("Server diagnostics", "POST /system-status -> app.get_system_status -> system_diagnostics.check_mafft plus geocode_locations.check_geonames -> short status dictionaries -> index.html. A lock and a 300-second per-process cache reduce repeated external work. GET /health remains a lightweight OK response and never runs these tests."),
    ("Processed environmental measurements", "POST /environmental -> load_measurements -> summarize_measurements -> make_charts -> table_records -> environmental.html. The input is a processed sample/feature table. This route does not call run_pipeline; its means and differences are descriptive measurements, not individual-genome distances."),
    ("Processed wastewater SNV calls", "POST /wastewater -> read_table/load_processed_snvs -> summarize_snvs -> wastewater.html. A row is an observed sample/reference/position/allele call. This mode accepts processed calls, not raw reads. Reference identity stays in variant keys. A missing call row cannot establish that the sample had adequate coverage and no mutation."),
    ("Temporary files and downloads", "Uploads are saved under fixed names in disposable folders. Organized environmental and wastewater CSV bytes then live briefly in separate, bounded dictionaries, each addressed by a random token. The page preview has at most 100 rows, while the download includes the full retained standardized table. Memory caches disappear when the worker restarts and are not shared by multiple worker processes."),
]

WORKED_EXAMPLES = [
    ("1. A numeric Segment column", "metadata_region('4') returns a dictionary with label 'Segment 4', segment '4', an empty gene-alias set, and whole=False. fasta_region('ACC.1 example segment 4') extracts the explicit header clue. regions_conflict returns False for matching segment clues. The validator chooses a resolved label and records its source; this parser itself does not align anything."),
    ("2. Two sequences with an insertion", "A validated group might contain ACGT and ACGTA. align_samples sends both to MAFFT and could receive ACGT- and ACGTA. The inherited comparison ignores the final position because one side has a gap. The four A/C/G/T positions match: Differences=0, Comparable_Positions=4, p-distance=0/4=0, identity=100%. The gap is not silently counted as an ordinary nucleotide mismatch."),
    ("3. A coarse location", "A record labeled 'USA: ND' is recognized by the conservative two-letter US-area test. The geocoder marks its precision as too coarse and skips the state-centroid lookup. The inherited geographic enrichment uses that precision marker to avoid a local-distance claim. A pair that passes genomic/date limits but lacks usable geography can appear in the separate candidate list, rather than the strict cluster result."),
    ("4. Samples, pairs, and clusters have different units", "Ten samples in one comparable group produce 10 x 9 / 2 = 45 unique unordered comparisons. Twenty-two qualifying candidate pairs are 22 relationships, not 22 samples. The website adds both IDs from every candidate to a set to count distinct samples represented. A connected-component cluster can include A and C through A-B and B-C even without a direct A-C passing pair; the backend cluster rule is transitive."),
    ("5. Environmental differences", "If feature_a at North has mean frequency 0.10 on January 1 and 0.30 on January 8, Change_From_Previous is 0.20 in fraction units, equal to 20 percentage points. It is not a 0.20% increase and not a 200% relative increase. The calculation does not divide by elapsed days. A missing measurement remains missing, so an absent value does not create an artificial decline to zero."),
    ("6. Duplicate wastewater call rows", "If W1 has the same call twice and W2 has it once, records=3, samples=2, unique_snvs=1, and Samples_Observed for that variant=2. Repeated keys are warned about and retained. Distinct_SNVs uses deduplicated sample/variant keys. Frequency/depth distributions still use retained rows, so repeated rows can influence their medians; the guide does not imply that every summary is deduplicated."),
    ("7. A diagnostic request", "The first click sends a POST without scientific uploads. The MAFFT test uses two invented short sequences and checks their IDs and equal output lengths; the GeoNames test asks only for fixed Tempe place coordinates. A second click within five minutes returns cached statuses. A valid HTTP 200 from GeoNames can still contain a JSON API status error, so the diagnostic checks both layers."),
]

PROFESSOR_QUESTIONS = [
    ("Which parts did Python implement, and which parts are delegated?", "Python handles validation, grouping, output checks, matching, arithmetic, summaries, and presentation. The external MAFFT executable performs multiple-sequence alignment. GeoNames supplies approximate place coordinates over HTTPS. Flask serves the interface; Gunicorn hosts Flask in Docker. The website does not reimplement genomic analysis in JavaScript."),
    ("Why use both Pathogen and Gene_Segment?", "A shared length does not establish comparable biology. The resolved labels define groups before alignment. Different pathogen/region groups receive separate MAFFT inputs. This depends on metadata/header accuracy; MAFFT is not a gene-identification method."),
    ("Does Operational prove scientific validity?", "No. The diagnostic establishes that a tiny fixed command/request worked under the server's current configuration. It cannot establish correct annotation of a user's dataset, the epidemiological suitability of a threshold, or transmission."),
    ("Are GeoNames coordinates exact specimen locations?", "No. The code labels them as approximate place coordinates and keeps provenance. State/admin-area precision is too coarse for the inherited geographic cutoff. Even supplied coordinates are labeled precision unverified; the code cannot independently verify field collection accuracy."),
    ("Why retain a reference genome in wastewater keys?", "A position has meaning only relative to a reference. Reference/position/Ref/Alt together identify the variant in the table. The special named-workbook branch supplies the reference encoded by that source-specific implementation; arbitrary standardized files must provide it explicitly."),
    ("Why not convert missing calls or measurements to zero?", "No observed record can reflect missing coverage, missing metadata, or no reported measurement. Zero would claim something stronger. The environmental and wastewater summaries preserve that uncertainty and count missingness explicitly."),
    ("Are quality flags used to filter data?", "In the current environmental module, Quality is summarized as provided. It does not automatically remove low-quality rows or weight means. This should be explained as descriptive reporting, not an inferred quality-control pipeline."),
    ("What does a passing test demonstrate?", "It demonstrates the behavior asserted for the controlled fixture and mocked condition. A PNG signature test does not inspect plot layout. A mocked MAFFT success does not prove the installed executable runs. Live server diagnostics and representative real-data tests provide different evidence."),
    ("What happens to downloads after a restart or with more workers?", "They are held only in one worker's memory. Restart clears them; a request reaching another worker may not find the token. The current one-worker deployment supports this simple design. Durable multi-user storage would need a separate future design."),
    ("Why are the gene-parser rules limited?", "They intentionally recognize explicit strings with simple operations so the code is explainable. Numeric segments, complete-genome phrases, gene neighbors, and a name/acronym pattern are supported. This is not a universal annotation ontology or a guarantee that every synonym is equivalent."),
]

VERIFICATION_NOTES = [
    "This document was built by reading the current source text and parsing its Python syntax without executing the scientific pipelines or making new MAFFT/GeoNames calls. No application, analysis, or test source file was edited for this documentation request.",
    "Every nonblank line in the listed files is represented once in the annotated tables. Blank physical lines are layout separators; source line-number gaps represent those blank lines. Complete physical line counts include blanks. Source SHA-256 fingerprints below identify the exact snapshot; code edits will change line numbers and fingerprints.",
    "Historical test counts and publication-workbook results in the project change log are not rerun or promoted to a new live-validation claim here. Existing test implementations are explained exactly as written, including mocked external behavior and the limited scope of image-signature tests.",
    "The single-genome comparison remains disease-agnostic. The separate processed-wastewater loader contains an explicitly source-specific recognition/mapping branch. Describe that branch as a special input adapter, not a universal pathogen/reference inference rule.",
    "Primary source for this code guide: the relative Python paths listed below, read from the working project on 19 October 2026. The earlier PDF and change log helped establish scope; the annotated source text itself controls the explanations. Publication comments are reported as implementation provenance, not independently revalidated literature claims.",
    "Further reading can start with the standard-library documentation for pathlib, tempfile, subprocess, threading, unittest and io, followed by the official pandas, Biopython, Flask, Matplotlib, MAFFT and GeoNames documentation. Consult current official setup documentation before changing deployment behavior.",
]

FILE_INFO.update({
    "genomic_regions.py": {"role": "New annotation helper supporting the inherited FASTA/metadata validator.", "summary": "Parses explicit region wording into a shared dictionary of label, segment, gene aliases, and whole-genome evidence. It uses ordinary string operations and checks contradictions before downstream grouping.", "limits": ["It does not identify a gene from sequence length or sequence content. Custom labels are retained, and only recognized structured clues participate in the conflict checks."]},
    "align_sequences.py": {"role": "New MAFFT adapter in the six-stage backbone.", "summary": "Finds MAFFT, writes one temporary FASTA per pathogen/region group, runs the external aligner, validates its output, and returns a copy of the input table with aligned sequences.", "limits": ["Singleton groups are skipped. The production alignment call has no explicit subprocess timeout in this snapshot; the diagnostic call has its own short timeout. This module does not calculate p-distance."]},
    "geocode_locations.py": {"role": "New optional coordinate adapter plus separate server diagnostic.", "summary": "Keeps complete supplied coordinates, optionally fills missing coordinate pairs from a fixed secure GeoNames endpoint, caches repeated places within a run, and records coordinate provenance/precision. The final function independently tests server connectivity.", "limits": ["The two-letter US-area check is conservative: it tests format, not membership in an official list of state abbreviations. Other administrative results are marked coarse after the API reply.", "The first GeoNames result is used; this is a place estimate, not a validated specimen collection point. A request/API failure ends further optional lookups in that run."]},
    "system_diagnostics.py": {"role": "New on-demand MAFFT execution check.", "summary": "Reuses MAFFT discovery, runs a short invented two-sequence alignment, verifies IDs and equal output lengths, and optionally extracts a simple version token. It returns a public-safe status dictionary.", "limits": ["This test does not assess the correctness of a user's biological grouping or alignment. The optional version command can fail while the functional alignment status stays Operational."]},
    "environmental_analysis.py": {"role": "New separate descriptive processed-measurement module.", "summary": "Loads sample/feature measurements, standardizes headers, validates values and units, summarizes observations by sample/date/place, and draws small descriptive PNG plots.", "limits": ["Quality labels are counted, not used as an automatic exclusion rule. Missing numeric values are retained. The module does not convert units or perform statistical inference.", "Temporal changes are raw differences in the selected metric, not rates per day or relative percentage changes. Charts show up to three alphabetically selected features and six locations per panel."]},
    "wastewater/__init__.py": {"role": "New package marker and documentation.", "summary": "This one-line docstring identifies the processed-wastewater package. It has no analysis function, import, or raw-read processing step."},
    "wastewater/processed_snv.py": {"role": "New separate processed wastewater-SNV loader and summary module.", "summary": "Accepts standardized CSV/XLSX SNV calls or a specifically recognized publication workbook, optionally joins sample metadata, validates observations, retains repeated keys with warnings, and summarizes observed variants.", "limits": ["It does not process FASTQ reads, call variants, assign lineages, or apply individual-genome outbreak clustering to mixtures.", "The reference accession supplied in the named-publication branch is source-specific. Standardized inputs must supply Reference_Genome explicitly. CSV/XLSX reading uses the first/default worksheet and exact standardized field names after trimming.", "Repeated rows remain in frequency/depth distributions, while distinct sample/variant counts use distinct keys. An unlisted call is not established mutation absence."]},
    "app.py": {"role": "Existing Flask file, expanded for the new modules; included in full for context.", "summary": "Connects forms to Python functions, manages temporary uploads and recent in-memory CSV downloads, caches dependency diagnostics, reads original pipeline outputs, and passes display-ready values to Jinja.", "limits": ["Temporary downloads are per process, bounded, and short lived; they vanish at restart. Preview limits do not remove rows from full downloads.", "The inherited threshold reader checks nonnegative finite numbers; it does not enforce a p-distance maximum of 1. Prototype threshold choice is still the user's responsibility.", "The stderr redirection used for the original pipeline is process-wide while active. The one-worker synchronous deployment reduces overlap; it is not a general multi-worker storage/concurrency design."]},
    "tests/test_alignment.py": {"role": "New test file protecting the MAFFT adapter.", "summary": "Exercises executable discovery, separate comparable groups, preserved IDs, downstream gap handling, changed-ID rejection, mixed segments, and Windows batch invocation. MAFFT output is simulated.", "limits": ["A patched which() returning /usr/bin/mafft tests how discovery results are handled; it does not install or execute Linux MAFFT."]},
    "tests/test_geocoding.py": {"role": "New test file protecting coordinate lookup and fallback behavior.", "summary": "Checks supplied-coordinate retention, missing configuration, repeated-place caching, secure HTTPS, state-only skipping, administrative precision, and request failures using mocked HTTP responses."},
    "tests/test_diagnostics.py": {"role": "New test file protecting the on-demand dependency tests.", "summary": "Simulates successful MAFFT/version output and missing, timeout, execution, malformed-output, and Windows conditions. It also simulates GeoNames configuration, valid coordinates, connection/API failures, and out-of-range coordinates."},
    "tests/test_environmental_analysis.py": {"role": "New loader, plot, and Flask tests for processed environmental measurements.", "summary": "Uses invented CSV rows and temporary XLSX files to verify aliases, missingness, alternative metrics, invalid values, frequency limits, PNG output, route isolation, results, and CSV downloads.", "limits": ["The measurement-value test happens to use values from 0 to 1; the production Measurement_Value field is not restricted to that range. The PNG signature test confirms an image format, not visual correctness."]},
    "tests/test_wastewater.py": {"role": "New loader and Flask tests for processed wastewater SNV observations.", "summary": "Verifies the publication mapping guard, standardized reference requirement, duplicate/sample accounting, invalid calls, optional metadata joins, missing values, conflicts, webpage output, and standardized CSV downloads."},
})

FUNCTION_INFO.update({
    "genomic_regions.py": {
        "metadata_region": ("Normalize one explicit metadata region value into a shared clue dictionary.", "Receives value plus an optional source-column label; returns label/segment/genes/whole.", "Numeric segments and explicit gene/whole-genome phrases become consistent grouping labels while custom labels are preserved."),
        "fasta_region": ("Read only explicit annotation words from a Biopython FASTA description.", "Receives the full description, including accession; returns the same four-field clue dictionary.", "The accession token is separated from annotation. Explicit gene/segment evidence, not sequence length, informs grouping."),
        "regions_conflict": ("Compare two annotation dictionaries for recognized contradictions.", "Receives two clue dictionaries; returns a Boolean.", "Different explicit segments/genes or whole-genome versus partial-region clues can be rejected upstream rather than aligned as though compatible."),
    },
    "align_sequences.py": {
        "find_mafft": ("Find the external aligner through the operating system PATH.", "No arguments; returns an executable path or raises RuntimeError.", "Both Linux mafft and Windows batch discovery use the platform-aware shutil.which behavior."),
        "align_samples": ("Align each compatible pathogen/region group and validate the returned FASTA.", "Receives a validated DataFrame and optional executable path; returns an aligned copy, raising RuntimeError on invalid output.", "Preserves identity/group boundaries and makes equal-length aligned strings available to the inherited comparison formula."),
    },
    "geocode_locations.py": {
        "is_state_only_location": ("Recognize conservative US country-code plus two-letter-area labels.", "Receives one location value; returns True/False.", "A coarse label such as USA: ND should not create a local sampling-distance claim."),
        "fill_missing_coordinates": ("Add optional approximate place coordinates with provenance and precision labels.", "Receives a sample DataFrame containing locations/coordinate columns; returns a copy with coordinate/source/precision fields.", "Preserves complete provided coordinates and missingness, while optional external lookup improves usable metadata without silently making a state centroid precise."),
        "check_geonames": ("Perform one live fixed-place diagnostic without returning the configured username.", "No arguments; reads the environment and returns a short status/message dictionary.", "Distinguishes configuration, connectivity, API-response, and coordinate-validation outcomes independently of sample analysis."),
    },
    "system_diagnostics.py": {
        "check_mafft": ("Demonstrate real executable operation using two short synthetic sequences.", "No user input; returns status/message and optional version.", "Executable existence alone is insufficient: the check also requires success, both expected IDs, and equal nonempty aligned length."),
    },
    "environmental_analysis.py": {
        "load_measurements": ("Read and validate one processed sample/feature measurement per row.", "Receives a CSV/XLSX path and optional metric name; returns (sorted standardized DataFrame, selected metric), or ValueError.", "Protects identifiers, numeric ranges, date consistency, duplicate sample/feature keys, and feature units before descriptive summaries."),
        "summarize_measurements": ("Build descriptive counts and sample/time/location/missingness/quality tables.", "Receives validated data and metric; returns a dictionary with counts and DataFrame tables.", "Uses observed values while preserving/reporting missing data. Differences describe the selected measurement, not an outbreak threshold."),
        "available_value_columns": ("List supported numeric columns present in a table.", "Receives a DataFrame; returns a list containing Frequency and/or Measurement_Value.", "Keeps missing-data summaries aligned with the actual uploaded schema."),
        "table_records": ("Prepare pandas values for Jinja presentation.", "Receives a DataFrame; returns a list of dictionaries with date strings and blanks instead of missing objects.", "Presentation formatting does not mutate the analytical DataFrame or redefine missing values as measured zeros."),
        "make_charts": ("Render descriptive time and location plots to PNG bytes.", "Receives validated data, selected metric, and summary tables; returns an in-memory PNG byte string.", "Plots reuse backend summary values and limit displayed features/places for readability; the graphics introduce no statistical inference."),
    },
    "wastewater/processed_snv.py": {
        "read_table": ("Read text-preserving CSV/XLSX and identify the special workbook header/provenance pattern.", "Receives a path; returns a DataFrame with an optional publication-identification attribute.", "Preserves numeric-looking sample IDs and avoids applying source-specific assumptions to arbitrary tables."),
        "load_processed_snvs": ("Standardize and validate processed SNV calls, with an optional metadata join.", "Receives call path and optional metadata path; returns (organized DataFrame, warning list), or ValueError.", "Reference identity, integer positions, single-base alleles, frequencies, depth, sample metadata, and retained duplicates are handled explicitly."),
        "summarize_snvs": ("Summarize observed call rows using distinct sample and reference-aware variant counts.", "Receives validated call rows; returns a summary dictionary including top variants and sample-table records.", "Separates records, unique samples, and unique SNVs; it does not turn unlisted calls into proven absence."),
    },
    "app.py": {
        "health": ("Answer the host's lightweight responsiveness probe.", "No form input; returns text OK with HTTP 200.", "A normal health probe must not repeatedly execute MAFFT or call GeoNames."),
        "upload_too_large": ("Handle Flask's request-size rejection.", "Receives the 413 error; returns readable text plus HTTP 413.", "The global limit is for the complete HTTP request, including multipart overhead."),
        "internal_error": ("Give a generic public response for unexpected server errors.", "Receives the 500 error; returns safe text plus HTTP 500.", "Internal traceback details remain out of the browser response."),
        "get_system_status": ("Run two independent diagnostics with five-minute per-process reuse.", "No inputs; returns (results dictionary, was_cached Boolean).", "A lock prevents simultaneous cache refreshes; one failed dependency does not suppress the other check."),
        "system_status": ("Handle only the explicit diagnostic-button POST.", "No scientific upload; renders the home template with diagnostic results and default form values.", "Normal page loads and /health do not perform dependency checks."),
        "environmental_page": ("Serve the separate processed-measurement upload/example workflow.", "GET shows the form; POST reads the uploaded/synthetic table and renders results plus a recent CSV token.", "Reuses environmental analysis functions, caps files/previews, and cleans temporary uploads."),
        "environmental_download": ("Return one recent organized environmental CSV.", "Receives a random download token from the URL; returns a CSV attachment or HTTP 404.", "The lookup checks expiry while holding the corresponding cache lock."),
        "wastewater_page": ("Serve the separate processed-SNV upload workflow.", "GET shows the form; POST accepts calls and optional sample metadata, then renders summaries/warnings/preview and a CSV token.", "Reuses processed-SNV functions and does not invoke the individual-genome runner."),
        "wastewater_download": ("Return one recent standardized wastewater CSV.", "Receives a URL token; returns the retained CSV bytes or HTTP 404.", "The original uploads have already been deleted; only bounded temporary output bytes remain."),
        "read_thresholds": ("Parse and validate the three user-entered screening limits.", "Receives form values; returns (three numeric limits, None) or (None, readable error).", "Rejects negative/nonfinite values before the runner is called; no scientific threshold is inferred."),
        "make_cluster_list": ("Turn existing membership CSV rows into display dictionaries.", "Receives membership DataFrame; returns cluster/pathogen/region/sample dictionaries.", "It formats the backend membership; it does not find new clusters."),
        "find_geography_unresolved_candidates": ("List existing pairs that pass genomic/date limits but lack usable geography.", "Receives enriched pair DataFrame and threshold list; returns candidate display dictionaries.", "Known distances are excluded even when above the cutoff; named-place strict matches are excluded too. Candidate meaning stays separate from strict related pairs."),
        "index": ("Run the original individual-genome workflow from an upload or bundled synthetic example.", "GET renders the form; POST validates form choices, calls run_pipeline, reads CSVs, and renders strict/candidate results.", "Keeps Python backend outputs as the source of truth and counts original unique samples separately from relationships."),
    },
})

note_lines("genomic_regions.py", '''
GENE_STOP_WORDS = {"complete", "partial", "full", "whole", "coding", || Begin a set of generic annotation words that must not be mistaken for a gene next to the word gene.
text = str(value).strip() || Convert the supplied metadata value to text and remove outer whitespace. The caller is responsible for deciding whether a value is genuinely present.
words = text.split() || Split on whitespace to inspect explicit words without regular expressions.
lower_words = [word.casefold() for word in words] || Make a lowercase-like comparison copy of each token; preserve original text separately for a custom label.
if len(words) == 1 and words[0].isdigit(): || Treat a single all-digit value such as 4 as explicit segment metadata.
segment = words[0] || Keep the numeric token as a text segment identifier, not a measurement.
if len(words) >= 2 and lower_words[0] == "segment": || Recognize an explicit segment prefix, including capitalization variants.
segment = words[1].strip(".,;:()[]{}").upper() || Clean punctuation from the next token and uppercase it for consistent segment comparison.
if lower_words in (["complete", "genome"], ["whole", "genome"], || Start an exact phrase test for complete genome, whole genome, or full genome; a missing region does not satisfy it.
if column in {"Gene", "Gene_Name"} or lower_words[:1] == ["gene"]: || Use explicit gene-column context or a gene-prefixed value as gene evidence.
gene = words[-1].strip(".,;:()[]{}").upper() || Take the final token as the gene label and normalize punctuation/case. This is a simple parser, not a general gene-name ontology.
return {"label": text, "segment": None, "genes": set(), "whole": False} || Preserve an unrecognized but explicit custom region label. No gene/segment/whole-genome clue is invented for conflict checks.
parts = str(description).split(maxsplit=1) || Split the description once: the first item is the FASTA accession token and the second, when present, is annotation text.
annotation = parts[1] if len(parts) == 2 else "" || Inspect only the annotation after the accession. An accession-only header supplies no region evidence.
words = annotation.replace("|", " ").replace(",", " ").replace(";", " ").split() || Replace common header separators with spaces and tokenize the explicit annotation.
clean_words = [word.strip(".,;:()[]{}").casefold() for word in words] || Build punctuation-trimmed, case-insensitive tokens while retaining original tokens for acronym checks.
segment = None || Start with no recognized segment. The scan may fill this when explicit wording is found.
gene = None || Start with no recognized gene. No identity is inferred from nucleotide count.
gene_aliases = set() || Keep recognized explicit gene names/acronyms in a distinct-value set for compatibility checks.
whole = False || Start without whole-genome confirmation; only an explicit phrase can change this clue.
for index, word in enumerate(clean_words): || Scan tokens with their positions so neighboring annotation words can be checked safely.
if word == "segment" and index + 1 < len(words): || Recognize segment only when a following token exists; the bounds check prevents an indexing error.
candidate = clean_words[index + 1] || Read the token immediately after segment.
if candidate.isalnum() and len(candidate) <= 5: || Accept a short alphanumeric segment label. This checks syntax, not membership in a pathogen-specific segment catalog.
segment = candidate.upper() || Store the recognized segment in a consistent uppercase form.
if word in {"complete", "whole", "full"} and index + 1 < len(words): || Look for explicit whole-genome adjectives with a following token.
if clean_words[index + 1] == "genome": || Require genome immediately after the adjective; complete sequence alone does not establish a whole genome.
whole = True || Record the explicit whole-genome clue.
if word == "gene": || Consider tokens immediately before or after the literal word gene.
if (index > 0 and clean_words[index - 1].isalpha() || Start a preceding-token test: it must exist and contain alphabetic characters only.
gene = clean_words[index - 1].upper() || Use the preceding explicit gene token as the display label.
gene_aliases.add(clean_words[index - 1]) || Keep its normalized name as a compatibility clue.
elif (index + 1 < len(words) and clean_words[index + 1].isalpha() || If the preceding token did not qualify, inspect a following alphabetic token within bounds.
gene = clean_words[index + 1].upper() || Use the following explicit gene token as the display label.
gene_aliases.add(clean_words[index + 1]) || Record this normalized gene alias for later conflict checks.
if index > 0 and words[index].startswith("(") and words[index].endswith(")"): || Look for a parenthesized acronym after another word, such as hemagglutinin (HA).
acronym = clean_words[index] || Retrieve the cleaned acronym token.
name = clean_words[index - 1] || Retrieve the immediately preceding putative full gene name.
if (2 <= len(acronym) <= 5 and acronym.isalpha() || Begin the conservative acronym rule: require a short alphabetic acronym, uppercase in the original, after a longer matching-initial name.
gene = acronym.upper() || Use the explicit acronym as the standardized gene label.
gene_aliases.update({acronym, name}) || Retain both the acronym and name as clues, so matching either can avoid a false contradiction.
if whole and (segment or gene): || Detect an internally mixed description claiming both whole genome and a specific segment/gene; return unresolved labeling with clues for upstream rejection.
if segment: || Choose a recognized segment label before a gene or whole-genome label.
label = f"Segment {segment}" || Build the internal display/group label for the explicit segment.
elif whole: || With no segment, use an explicit whole-genome clue when present.
label = "Complete genome" || Standardize the supported complete/whole/full-genome phrases to one label.
elif gene: || With no segment/whole-genome choice, use a recognized explicit gene.
label = f"Gene {gene}" || Build the internal gene display/group label.
label = None || Leave the region unresolved when no supported explicit annotation was found.
if first["whole"] and (second["segment"] or second["genes"]): || Whole-genome evidence conflicts with a second annotation naming a gene or segment.
if second["whole"] and (first["segment"] or first["genes"]): || Check the same contradiction with the two annotations in reversed order.
if first["segment"] and second["segment"]: || Compare segments only when both annotations actually supplied segment clues.
if first["segment"] != second["segment"]: || Two different explicit segment IDs constitute a recognized conflict.
if first["genes"] and second["genes"]: || Compare gene aliases only when both sides contain some gene evidence.
if first["genes"].isdisjoint(second["genes"]): || If the recognized alias sets share no value, treat those explicit gene annotations as contradictory.
return False || No supported contradiction was found. This is not proof that arbitrary unknown/custom annotations are biologically equivalent.
''')

note_lines("align_sequences.py", '''
mafft_path = shutil.which("mafft") || Search PATH with platform rules. On Windows, executable extensions such as .bat can be discovered through PATHEXT; Linux normally returns an executable file.
if mafft_path is None: || Stop discovery when the operating system cannot locate MAFFT; do not silently compare unaligned data.
return mafft_path || Return the located program path to the runner or diagnostic caller.
if mafft_path is None: || If the caller did not already supply the checked path, discover MAFFT here.
mafft_path = find_mafft() || Reuse the discovery function and its readable missing-program exception.
aligned = data.copy() || Create a separate result table so input sequences are not overwritten in the caller's original DataFrame.
groups = data.groupby(["Pathogen", "Gene_Segment"], sort=True) || Partition by both biological grouping keys; unrelated pathogens or regions never enter the same MAFFT group.
for group_names, group in groups: || Process one comparable group at a time.
pathogen, gene_segment = group_names || Unpack the group's two labels for readable error messages.
if len(group) < 2: || A single sample has no within-group partner, so no multi-sequence alignment is needed.
input_path = Path(folder) / "group.fasta" || Create a fixed filename in this group's disposable folder.
with input_path.open("w", encoding="utf-8") as fasta_file: || Open the temporary FASTA for text writing; the file is closed automatically before MAFFT runs.
for _, sample in group.iterrows(): || Visit each sample row. The original row index is unused in this writing step.
fasta_file.write( || Start writing one FASTA record: its original Sample_ID becomes the header and Sequence becomes the next line.
command = [mafft_path, "--auto", "--nuc", str(input_path)] || Build an argument list for MAFFT's automatic nucleotide alignment. No custom alignment algorithm is implemented in Python.
if Path(mafft_path).suffix.lower() in {".bat", ".cmd"}: || Detect the Windows batch-wrapper case without changing Linux invocation.
command = ["cmd.exe", "/c"] + command || Ask Windows cmd to execute the discovered batch file. The command is constructed by server code, not by a browser-supplied shell string.
result = subprocess.run( || Start MAFFT and wait for it to finish, capturing its output. This production call has no explicit timeout argument in the current source.
except (OSError, subprocess.CalledProcessError) as error: || Handle startup/OS failure and a nonzero MAFFT exit code raised because check=True.
details = getattr(error, "stderr", "") or str(error) || Prefer captured program error text when available; otherwise use the exception description.
records = list(SeqIO.parse(StringIO(result.stdout), "fasta")) || Parse captured FASTA text into Biopython records after the command returns successfully.
aligned_by_id = {} || Build a lookup from returned sample ID to cleaned aligned sequence.
for record in records: || Validate every sequence record returned by MAFFT.
sample_id = str(record.id) || Read the exact returned FASTA identifier as a string.
if sample_id in aligned_by_id: || Reject duplicate returned IDs because they would make sample-to-sequence mapping ambiguous.
sequence = clean_sequence(record.seq) || Apply the same sequence normalization used by input validation.
if set(sequence) - ALLOWED_SEQUENCE_CHARACTERS: || Subtract the allowed alphabet; any remaining character means the aligned output is invalid.
aligned_by_id[sample_id] = sequence || Store the validated sequence under its returned ID.
expected_ids = set(group["Sample_ID"].astype(str)) || Collect all IDs that were sent for this group.
if set(aligned_by_id) != expected_ids: || Reject omitted, added, or changed IDs; alignment must preserve sample identity.
aligned_lengths = {len(sequence) for sequence in aligned_by_id.values()} || Calculate the set of output lengths; a valid alignment has exactly one shared length.
if len(aligned_lengths) != 1: || Reject unequal aligned lengths instead of attempting direct comparison of incompatible columns.
for index, sample in group.iterrows(): || Revisit original rows using their indices so aligned strings return to the correct samples.
sample_id = str(sample["Sample_ID"]) || Retrieve the original expected ID for this row.
aligned.at[index, "Sequence"] = aligned_by_id[sample_id] || Replace only the result-copy sequence cell with the matching aligned string.
return aligned || Return all rows with aligned group sequences and the other input columns preserved.
''')

note_lines("system_diagnostics.py", '''
mafft_path = find_mafft() || Reuse the same Linux/Windows discovery as the production pipeline; finding the path is only the first part of this diagnostic.
except RuntimeError: || Convert a missing-MAFFT discovery failure into a short public status rather than exposing an exception traceback.
input_path = Path(folder) / "diagnostic.fasta" || Give the fixed synthetic test input its own disposable file.
input_path.write_text( || Write two invented short nucleotide records with known IDs and different original lengths, so the check exercises an actual alignment.
command = [mafft_path, "--auto", "--nuc", str(input_path)] || Execute the discovered MAFFT with fixed automatic nucleotide-alignment arguments and a server-created input path.
if Path(mafft_path).suffix.lower() in {".bat", ".cmd"}: || Use the Windows wrapper only for discovered batch/command files.
command = ["cmd.exe", "/c"] + command || Preserve Windows batch support; Linux keeps the direct executable argument list.
result = subprocess.run( || Start the real small alignment, capture its output, and apply the diagnostic's short timeout.
except subprocess.TimeoutExpired: || Return a readable timed-out status when the synthetic alignment exceeds its allowed wait.
except OSError: || Handle a discovered program that cannot actually be started.
if result.returncode != 0: || Explicitly inspect command success; existence of an executable alone is insufficient.
records = list(SeqIO.parse(StringIO(result.stdout), "fasta")) || Parse the command's actual captured output as FASTA.
except ValueError: || Treat malformed FASTA parsing as an unsuccessful check.
records = [] || Use an empty parsed result so the following identity/length validation fails clearly.
sequences = {record.id: str(record.seq) for record in records} || Make an ID-to-sequence lookup for checking the two expected synthetic records.
len(records) != 2 || Require exactly two returned records; duplicates or extra records must not pass.
or set(sequences) != {"DIAG_A", "DIAG_B"} || Require the two original IDs. Python or short-circuiting prevents unsafe later ID access when this fails.
or not sequences["DIAG_A"] || Reject an empty aligned sequence after the ID set has been established.
or len(sequences["DIAG_A"]) != len(sequences["DIAG_B"]) || Require equal aligned lengths; the check does not claim to establish biological alignment accuracy for user data.
version = None || Begin with an unavailable version; a successful alignment can remain Operational without a version token.
version_command = [mafft_path, "--version"] || Build a second, fixed optional command for version reporting.
version_command = ["cmd.exe", "/c"] + version_command || Apply the same Windows batch-wrapper handling to the version command.
version_result = subprocess.run( || Ask MAFFT for its version with a separate short wait; failure here does not invalidate the alignment check.
if version_result.returncode == 0: || Only inspect version output from a successful version command.
for word in (version_result.stdout + " " + version_result.stderr).split(): || MAFFT may print its version to stdout or stderr, so inspect both as whitespace-separated tokens.
candidate = word.strip(" ,;()") || Trim surrounding punctuation before testing a possible version token.
candidate.startswith("v") || Require the simple supported version prefix.
and any(character.isdigit() for character in candidate) || Require at least one digit rather than accepting a lone letter v.
and all(character in "v0123456789." for character in candidate) || Restrict displayed version characters to a small safe version-token alphabet; raw error output is not shown.
version = candidate || Store the first clearly recognized version token.
except (OSError, subprocess.TimeoutExpired): || An optional version-start or timeout failure is handled independently of the successful alignment.
''')

note_lines("geocode_locations.py", '''
GEONAMES_URL = "https://secure.geonames.org/searchJSON" || Keep one fixed secure endpoint shared by sample lookup and the diagnostic. Visitors cannot supply an arbitrary API URL.
if pd.isna(location): || Missing location text is not a recognizable state-only label.
country, separator, area = " ".join(str(location).split()).partition(":") || Collapse whitespace and split once at the colon into country, separator, and area text.
bool(separator) || Require an actual colon in the country:area format.
and country.casefold() in {"usa", "us", "united states"} || Limit this conservative short-code rule to the listed US country labels, ignoring capitalization.
and len(area.strip()) == 2 || Require a two-character area code, as in ND; this is a format check, not a state-code registry lookup.
and area.strip().isalpha() || Require alphabetic area characters rather than treating a number as a state label.
samples = data.copy() || Preserve the caller's table while enriching a result copy.
samples["Coordinate_Source"] = pd.NA || Initialize provenance as missing until a provided or GeoNames coordinate pair is established.
samples["Geographic_Precision"] = "Named location only (precision unverified)" || Start with an explicit uncertainty label for text-only places.
provided = samples["Latitude"].notna() & samples["Longitude"].notna() || A usable supplied pair requires both coordinate cells; one nonmissing component is insufficient.
samples.loc[provided, "Coordinate_Source"] = "Provided" || Label complete original coordinate pairs without overwriting them.
samples.loc[provided, "Geographic_Precision"] = ( || Label supplied coordinates as precision unverified; numerical validity does not establish specimen-level accuracy.
state_only = samples["Location_Name"].map(is_state_only_location) || Apply the conservative state-format test to every location label.
samples.loc[state_only, "Geographic_Precision"] = ( || Mark state-only locations too coarse even if some coordinates were supplied; the downstream screen reads this precision field.
if state_only.any(): || Report when at least one sample has a state-only label.
missing_count = int((~provided & ~state_only).sum()) || Count only incomplete coordinate pairs that are not already excluded as state-only places.
if missing_count == 0: || If no eligible missing coordinates remain, return without checking credentials or sending API requests.
username = os.getenv("GEONAMES_USERNAME", "").strip() || Read the username from the server environment and remove outer whitespace. Do not hard-code or return the value.
if not username: || Without configuration, skip optional external lookup and retain the inherited specific-place fallback.
cache = {} || Start a per-analysis lookup cache keyed by normalized place text; it is not a permanent database.
geocoded_count = 0 || Count samples filled from GeoNames, not distinct place strings.
coarse_count = 0 || Count filled results identified as administrative/coarse places.
for index, sample in samples.loc[~provided & ~state_only].iterrows(): || Visit only samples missing a complete supplied coordinate pair and not flagged state-only.
location = sample["Location_Name"] || Retrieve this sample's reported place text.
query = " ".join(str(location).split()) || Collapse repeated whitespace for a stable search query and cache key.
if not query: || Skip blank place names rather than sending an empty search.
cache_key = query.casefold() || Merge capitalization variants of the same normalized query in the per-run cache.
if cache_key not in cache: || Request an uncached place once; repeated place rows reuse the result, including a cached no-result.
response = requests.get( || Send the fixed HTTPS request using the environment-provided username and explicit certificate verification.
response.raise_for_status() || Check transport-level HTTP success; GeoNames can still report an API problem inside successful HTTP JSON.
result = response.json() || Decode the API response body into Python dictionaries/lists.
except (requests.RequestException, ValueError): || Handle connection/HTTP/JSON problems with a warning, then end further optional lookups for this analysis.
if not isinstance(result, dict) or "status" in result: || Treat a non-object response or a GeoNames JSON status message as unavailable lookup service.
cache[cache_key] = None || Cache the default no-usable-result state to avoid requesting this place repeatedly.
places = result.get("geonames", []) || Read the result list, defaulting to empty when no place was found.
if places: || Only try coordinate extraction when at least one result exists; maxRows=1 means the first match is used.
latitude = float(places[0]["lat"]) || Convert the first place's latitude string to a numeric value.
longitude = float(places[0]["lng"]) || Convert its longitude string to a numeric value.
except (KeyError, TypeError, ValueError, IndexError): || Treat missing or malformed coordinate fields as no usable result for this place.
valid_latitude = math.isfinite(latitude) and -90 <= latitude <= 90 || Require a finite latitude inside its geographic numerical range.
valid_longitude = math.isfinite(longitude) and -180 <= longitude <= 180 || Require a finite longitude inside its geographic numerical range.
if valid_latitude and valid_longitude: || Store a result only when both coordinate components pass validation.
feature_code = str(places[0].get("fcode", "")).upper() || Read the GeoNames feature code, if supplied, for precision classification.
feature_class = str(places[0].get("fcl", "")).upper() || Read the broad feature class for the same conservative precision safeguard.
if feature_class == "A" or feature_code.startswith( || Administrative/country-style results are too coarse; the following prefixes also identify administrative or political units.
precision = "Administrative area (too coarse for distance cutoff)" || Retain the estimate but flag it as unusable for a local kilometer relationship.
precision = "GeoNames place (approximate)" || A nonadministrative match is still an approximate place coordinate, not an exact collection site.
cache[cache_key] = (latitude, longitude, precision) || Cache both numeric coordinates and the precision description together.
coordinates = cache[cache_key] || Retrieve the new or reused lookup outcome for this normalized place.
if coordinates is None: || Leave this sample unresolved when no valid coordinate pair was found.
samples.at[index, "Latitude"] = coordinates[0] || Fill the result-copy latitude from the resolved place estimate.
samples.at[index, "Longitude"] = coordinates[1] || Fill its paired longitude too; a partial original pair is replaced together.
samples.at[index, "Coordinate_Source"] = "GeoNames" || Preserve the provenance of these estimated coordinates.
samples.at[index, "Geographic_Precision"] = coordinates[2] || Preserve the classification needed by downstream geographic screening and display.
geocoded_count += 1 || Increase the number of sample rows filled; repeated places can fill several samples from one request.
if "too coarse" in coordinates[2]: || Count administrative estimates separately from more specific place estimates.
coarse_count += 1 || Increase the number of filled-but-coarse rows reported to the caller.
unresolved_count = int( || Begin counting rows still missing either coordinate; a coarse but numerically filled pair is not included in this numerical-missing count.
if coarse_count: || Warn explicitly if some GeoNames estimates cannot satisfy the local-distance criterion.
return samples || Return the enriched copy with source/precision fields; downstream code decides which evidence can satisfy geographic screening.
response = requests.get( || The same fixed secure endpoint is used here; the exact query/timeout differs between the sample-enrichment function and diagnostic.
except requests.HTTPError: || Distinguish an HTTP API rejection from a failure to connect at all.
except requests.RequestException: || Report connectivity/timeout/certificate failures without returning exception details or the username.
if not isinstance(result, dict) or "status" in result: || A GeoNames API status object, even with HTTP 200, is a diagnostic failure; do not display its raw message.
place = result["geonames"][0] || Select the fixed Tempe test lookup's first returned place; no arbitrary browser-supplied place is accepted here.
latitude = float(place["lat"]) || Parse the diagnostic's returned latitude for explicit numerical validation.
longitude = float(place["lng"]) || Parse its paired longitude for explicit numerical validation.
except (ValueError, TypeError, KeyError, IndexError): || Convert invalid JSON structure or coordinate values into a readable Check failed result.
if not (valid_latitude and valid_longitude): || A response with out-of-range, NaN, or infinite coordinates must not be labeled Operational.
''')

note_lines("environmental_analysis.py", '''
REQUIRED_COLUMNS = ["Sample_ID", "Collection_Date", "Location_Name", "Feature_ID"] || Declare required field names. A required column can still contain a missing date/location cell; blank sample/feature IDs are errors.
ALIASES = { || Start the dictionary mapping common normalized header spellings into this module's standardized schema.
path = Path(path) || Normalize the supplied file-path object before extension and file operations.
if path.suffix.lower() == ".csv": || Choose CSV parsing based on a case-insensitive extension.
data = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig") || Load source cells as text to preserve IDs; accept a UTF-8 BOM and preserve default-NA-like text for explicit validation.
elif path.suffix.lower() == ".xlsx": || Otherwise accept an XLSX workbook, using pandas' default worksheet.
data = pd.read_excel(path, dtype=str, keep_default_na=False) || Read worksheet cells as text before controlled number/date conversion.
names = [] || Collect one standardized destination name per original header.
for name in data.columns: || Process every supplied header in its existing order.
key = "_".join(str(name).strip().casefold().replace("-", " ").replace("_", " ").split()) || Trim/case-normalize headers, treat hyphens/underscores as spaces, collapse whitespace, and join tokens with underscores; no regex is used.
names.append(ALIASES.get(key, str(name).strip())) || Map a known alias or preserve a trimmed unknown header.
data.columns = names || Replace the DataFrame's header list with the standardized names.
if data.columns.duplicated().any(): || Reject two source columns collapsing onto the same internal field; silently choosing one would lose information.
missing = [name for name in REQUIRED_COLUMNS if name not in data] || List the mandatory standardized columns absent from the file.
if missing: || Report the missing fields before attempting calculations that need them.
available_metrics = [name for name in ["Frequency", "Measurement_Value"] if name in data] || Discover supported numeric columns, with Frequency first in automatic-selection priority.
if not available_metrics: || Require at least one supported measurement column.
if data.empty or len(data) > 10000: || Enforce the prototype's input bound of 1 through 10,000 measurement rows.
for column in data: || Visit each source column, including retained optional columns.
data[column] = data[column].fillna("").astype(str).str.strip() || Standardize missing textual cells to blank and remove outer whitespace before field-specific conversion.
for column in ["Quality", "Unit"]: || Ensure optional quality/unit fields exist even when omitted from the input.
if column not in data: || Add an optional field only when no such column was supplied.
data[column] = "" || Create the optional column as blank rather than inventing a quality assessment or unit.
errors = [] || Collect multiple validation problems before reporting them together.
for column in ["Sample_ID", "Feature_ID"]: || Check the two identifiers necessary to define a measurement row.
count = int(data[column].eq("").sum()) || Count blank IDs explicitly rather than deleting those rows.
if count: || Add a readable error only when at least one blank ID exists.
raw_dates = data["Collection_Date"].str.removesuffix(" 00:00:00") || Accept Excel date cells represented as midnight timestamps by removing that exact trailing time.
dates = pd.to_datetime(raw_dates.replace("", pd.NA), format="%Y-%m-%d", errors="coerce") || Parse explicit YYYY-MM-DD dates; blanks and bad values become missing for validation, while bad nonblank strings are reported next.
if (raw_dates.ne("") & dates.isna()).any(): || Reject a supplied nonblank date that could not be parsed; an honestly blank date remains allowed.
data["Collection_Date"] = dates || Replace raw strings with pandas datetime values for chronological grouping.
locations = {} || Build a display-label lookup that merges capitalization/spacing variants of a location.
for name in data["Location_Name"]: || Inspect every reported location label, including blanks.
label = " ".join(name.split()) || Collapse repeated internal whitespace and remove leading/trailing whitespace.
locations.setdefault(label.casefold(), label) || Keep the first display spelling for each case-insensitive normalized place key.
data["Location_Name"] = data["Location_Name"].map( || Replace each original label with its consistent display spelling from the lookup.
lambda name: locations[" ".join(name.split()).casefold()] || For each place string, normalize it the same way and retrieve the retained display label.
for column in available_metrics: || Validate every supported numeric column present, even when only one will be selected for summaries.
raw = data[column] || Keep source text to distinguish blank values from nonblank invalid numeric text.
values = pd.to_numeric(raw.replace("", pd.NA), errors="coerce") || Convert numeric text; preserve blanks as missing and mark unparseable strings missing for the next check.
invalid = raw.ne("") & (values.isna() | values.isin([float("inf"), float("-inf")])) || Flag nonblank text that became missing or infinite; a genuinely blank measurement is allowed.
if invalid.any(): || Accumulate a numeric-validation error rather than silently excluding invalid rows.
if column == "Frequency" and (values.notna() & ~values.between(0, 1)).any(): || Frequencies must be fractions in 0-1. Measurement_Value does not receive this specific range restriction.
data[column] = values.astype(float) || Store the parsed numerical series; missing values remain missing float cells.
for column in ["Collection_Date", "Location_Name"]: || Check whether one sample ID has contradictory known date or place metadata.
present = data[column].replace("", pd.NA) || Treat blank location strings as missing during consistency checks.
if present.groupby(data["Sample_ID"]).nunique().gt(1).any(): || More than one distinct known value within a sample indicates conflicting metadata.
if data.duplicated(["Sample_ID", "Feature_ID"]).any(): || Reject repeated sample/feature keys; biological or technical replicates must be resolved or separately identified upstream.
units = data["Unit"].replace("", pd.NA).groupby(data["Feature_ID"]).nunique() || Count distinct nonblank unit labels for each feature across the file.
if units.gt(1).any(): || Reject a feature represented in multiple units because no unit conversion is implemented.
if metric == "auto": || Select a numeric measurement only when the caller did not explicitly choose one.
metric = next((name for name in available_metrics if data[name].notna().any()), available_metrics[0]) || Choose the first supported column with observed values, preferring Frequency; if all are missing, retain a column so the next validation can explain the problem.
if metric not in available_metrics: || Reject an explicitly selected metric column that is absent.
elif not data[metric].notna().any(): || Reject a selected metric with no usable observed values.
if errors: || Stop only after the collected validation problems can be reported together.
raise ValueError(" ".join(errors)) || Join all validation messages into one exception for the Flask page or caller.
return data.sort_values(["Sample_ID", "Feature_ID"]).reset_index(drop=True), metric || Return a predictable sample/feature order and fresh row index together with the chosen metric.
samples = data.groupby("Sample_ID", as_index=False).agg( || Build one descriptive sample-summary row per distinct Sample_ID.
Collection_Date=("Collection_Date", "first"), || Use the first nonmissing collection date in the already consistency-checked sample group.
Location_Name=("Location_Name", "first"), || Take the group's first location label. Missing metadata on individual measurement rows is not globally imputed here.
Features_Listed=("Feature_ID", "nunique"), || Count distinct features listed for the sample, whether or not their measurement is missing.
Values_Available=(metric, "count"), || Count only nonmissing values in the selected metric.
samples["Missing_Values"] = samples["Features_Listed"] - samples["Values_Available"] || Calculate the number of listed features lacking a value; duplicate sample/feature rows were rejected earlier.
observed = data.loc[data[metric].notna()] || Select observed numeric rows for means, ranges, and observed-feature counts; retain the full table separately.
dated = observed.loc[observed["Collection_Date"].notna()] || Temporal summaries require both an observed value and a known date.
temporal = dated.groupby( || Start grouping valid dated observations by feature, place, and collection date.
).agg(Mean_Value=(metric, "mean"), Samples_With_Value=("Sample_ID", "nunique")) || Calculate an ordinary mean and a distinct sample count in each feature/place/date group; no depth weighting is applied.
temporal = temporal.sort_values(["Feature_ID", "Location_Name", "Collection_Date"]) || Put dates in chronological order within each feature/place before previous-date calculations.
groups = temporal.groupby(["Feature_ID", "Location_Name"]) || Keep previous-date comparisons within one feature and location.
temporal["Previous_Date"] = groups["Collection_Date"].shift() || Shift dates by one row within each group to identify the preceding observed date; the first is missing.
temporal["Change_From_Previous"] = groups["Mean_Value"].diff() || Subtract the preceding mean from the current mean within the group. This is not normalized by days or expressed as a relative percentage change.
locations = observed.loc[observed["Location_Name"].ne("")].groupby( || Exclude blank location labels from location summaries and group observed values by feature/place across dates.
Mean_Value=(metric, "mean"), Minimum=(metric, "min"), Maximum=(metric, "max"), || Summarize observed values with an ordinary mean and range; these are descriptive statistics.
Samples_With_Value=("Sample_ID", "nunique"), || Count distinct samples contributing observed values to the group.
missing = [] || Collect one missingness-summary record per relevant field.
for column in REQUIRED_COLUMNS + available_value_columns(data) + ["Quality", "Unit"]: || Include required metadata, actual numeric measurement columns, and optional quality/unit fields in the missing-data report.
count = int((data[column].isna() | data[column].eq("")).sum()) || Count both pandas missing values and blank text cells.
missing.append({"Field": column, "Missing_Rows": count, "Total_Rows": len(data)}) || Save the field name, missing count, and input-row denominator for transparent reporting.
quality = data["Quality"].replace("", "Not provided").value_counts().rename_axis("Quality").reset_index(name="Rows") || Count supplied quality labels, labeling blanks only in this display summary; quality is not used to filter or weight measurements.
return [name for name in ["Frequency", "Measurement_Value"] if name in data] || Return supported numeric columns that actually exist, in a fixed order.
display = table.copy() || Format a copy so the underlying analytical table is not changed for display.
for column in display: || Inspect each result column for display formatting.
if pd.api.types.is_datetime64_any_dtype(display[column]): || Detect pandas datetime columns before applying date-string formatting.
display[column] = display[column].dt.strftime("%Y-%m-%d") || Present known dates consistently as YYYY-MM-DD strings.
display = display.astype(object).where(pd.notna(display), "") || Allow mixed display values and replace missing objects with empty strings for HTML.
return display.to_dict("records") || Return one ordinary Python dictionary per table row for Jinja iteration.
features = sorted(data.loc[data[metric].notna(), "Feature_ID"].unique())[:3] || Choose up to three alphabetically ordered features with observed values; this is a display limit, not significance selection.
figure = Figure(figsize=(10, 3.6 * len(features)), layout="constrained") || Make a noninteractive figure with enough height for each selected feature and constrained layout spacing.
axes = figure.subplots(len(features), 2, squeeze=False) || Create two plot columns per feature and keep a two-dimensional axes array even for one feature.
for index, feature in enumerate(features): || Draw one row of time/location panels for each selected feature.
dates = tables["temporal_changes"].loc[tables["temporal_changes"]["Feature_ID"].eq(feature)] || Reuse this feature's already calculated temporal summary rather than recalculate analysis in HTML.
locations = tables["location_summary"].loc[tables["location_summary"]["Feature_ID"].eq(feature)].head(6) || Reuse up to the first six location-summary rows for readable bars; they are not necessarily the largest values.
time_axis, location_axis = axes[index] || Assign the left and right axes for this feature.
for location, group in list(dates.groupby("Location_Name"))[:6]: || Plot at most six place-specific series for readability.
if dates.empty: || Display an explanation when this feature lacks dated observations.
locator = AutoDateLocator() || Select readable date tick locations based on the plotted range.
time_axis.xaxis.set_major_locator(locator) || Apply the date tick positions to the time axis.
time_axis.xaxis.set_major_formatter(ConciseDateFormatter(locator)) || Format the ticks without repetitive long date strings.
time_axis.legend(fontsize=7) || Identify the plotted location series with a small legend.
if locations.empty: || Display an explanation when no location-summary values are available.
units = data.loc[data["Feature_ID"].eq(feature), "Unit"] || Collect unit labels associated with the current feature.
units = units.loc[units.ne("")] || Ignore missing unit labels when choosing an axis label.
label = "Frequency (0–1)" if metric == "Frequency" else (units.iloc[0] if not units.empty else "Measurement value") || Use a fraction label for Frequency, otherwise a provided consistent unit or a generic measurement label. No conversion is performed.
time_axis.set_title(f"{str(feature)[:48]}: over time") || Label the feature's time panel while capping long feature identifiers for layout.
time_axis.set_ylabel(label) || Label the time panel's numeric axis with its measurement interpretation.
location_axis.set_title(f"{str(feature)[:48]}: location means") || Identify the descriptive location-mean panel.
location_axis.set_xlabel(label) || Label the bar lengths with the same metric/unit interpretation.
image = BytesIO() || Allocate an in-memory binary output buffer instead of writing a plot file into the application folder.
figure.savefig(image, format="png", dpi=110) || Render a display-resolution PNG into the buffer.
return image.getvalue() || Return PNG bytes; Flask later base64-encodes them for the webpage.
''')

note_lines("wastewater/processed_snv.py", '''
PAPER_COLUMNS = { || Define an exact source-header adapter for the named publication workbook. This is input translation, not a general reference-genome inference rule.
REQUIRED_COLUMNS = [ || Declare the standardized call fields that arbitrary inputs must provide, including explicit reference identity and a depth column.
OUTPUT_COLUMNS = REQUIRED_COLUMNS + ["Collection_Date", "Location_Name"] || Fix the organized/download column order while retaining optional date and location metadata.
path = Path(path) || Normalize the path object for extension selection and pandas input.
if path.suffix.lower() == ".csv": || Accept a case-insensitive CSV extension using the standardized table-reader path.
return pd.read_csv(path, dtype=str, keep_default_na=False) || Preserve numeric-looking sample IDs as strings and keep NA-like text for explicit validation.
if path.suffix.lower() == ".xlsx": || Inspect an XLSX workbook, using the default sheet.
first_rows = pd.read_excel(path, header=None, nrows=2, dtype=str) || Inspect the first two worksheet rows without assuming a header; the special source has one introduction row.
header_row = 1 if len(first_rows) > 1 and first_rows.shape[1] > 0 and str(first_rows.iloc[1, 0]).strip() == "SAMPLE" else 0 || Use row index 1 as the header only when its first cell is SAMPLE; otherwise use the ordinary first row.
data = pd.read_excel(path, header=header_row, dtype=str, keep_default_na=False) || Read the complete selected sheet with the detected header while preserving identifier text.
description = str(first_rows.iloc[0, 4]) if first_rows.shape[1] > 4 else "" || Inspect the expected introductory description cell only when that column exists.
data.attrs["fontenele_2021_table2"] = ( || Store a table-provenance flag in DataFrame attributes; all following source-fingerprint conditions must match.
return data || Return the table, including its provenance attribute when the special workbook pattern was recognized.
data = read_table(path) || Reuse the format/header reader for processed calls.
data.columns = [str(name).strip() for name in data.columns] || Trim header edges; this module uses exact standardized names after trimming, unlike the broad environmental alias map.
published = ( || Require the provenance attribute plus expected paper fields and location fields before applying the special source adapter.
warnings = [] || Accumulate nonfatal problems while retaining rows transparently.
if published: || Enter the source-specific translation path only for the recognized workbook fingerprint.
data = data.rename(columns=PAPER_COLUMNS) || Translate paper field names into the loader's common schema.
data["Reference_Genome"] = "NC_045512.2" || Assign the reference encoded by this source-specific adapter. Arbitrary files do not receive this default.
data["Location_Name"] = data["STATE"].astype(str).str.strip() + ": " + data["LOCATION"].astype(str).str.strip() || Combine the paper's state and site label into one readable location string.
dates = pd.to_datetime(data["Collection_Date"], format="%d_%b_%Y", errors="coerce") || Parse the source's explicit day_month-name_year convention, such as 30_May_2020.
missing = sorted(set(REQUIRED_COLUMNS) - set(data.columns)) || For general standardized files, list required call fields that are absent.
if missing: || Reject an incomplete standardized schema with a named-field error.
if "Collection_Date" not in data: || Add a blank optional date column when it was omitted.
data["Collection_Date"] = "" || Represent unknown dates honestly rather than infer them.
if "Location_Name" not in data: || Add a blank optional location column when omitted.
data["Location_Name"] = "" || Represent unknown places as blank text.
dates = pd.to_datetime(data["Collection_Date"].replace("", pd.NA), errors="coerce") || Parse general source dates; blank cells become missing, and invalid nonblank dates are recorded next.
invalid_source_dates = data["Collection_Date"].fillna("").astype(str).str.strip().ne("") & dates.isna() || Distinguish a supplied but invalid date from an allowed missing date.
if data.empty: || Reject a table with no call records.
for column in ["Sample_ID", "Reference_Genome", "Ref", "Alt", "Location_Name"]: || Clean the text fields needed for matching, reference-aware keys, nucleotides, and place labels.
data[column] = data[column].fillna("").astype(str).str.strip() || Convert null text to blank and remove outer whitespace without dropping rows.
data["Ref"] = data["Ref"].str.upper() || Normalize reference nucleotides to uppercase before alphabet checking.
data["Alt"] = data["Alt"].str.upper() || Normalize alternative nucleotides the same way.
errors = [] || Start an aggregated validation-error list distinct from nonfatal warnings.
if metadata_path: || Perform the optional sample metadata join only when a file was supplied.
metadata = read_table(metadata_path) || Reuse the same CSV/XLSX reader for sample metadata.
metadata.columns = [str(name).strip() for name in metadata.columns] || Trim metadata headers before checking the exact Sample_ID field.
if "Sample_ID" not in metadata: || Require the join key explicitly.
metadata["Sample_ID"] = metadata["Sample_ID"].fillna("").astype(str).str.strip() || Preserve text IDs while removing missing/outer-whitespace inconsistencies.
if metadata["Sample_ID"].eq("").any() or metadata["Sample_ID"].duplicated().any(): || Reject blank or repeated metadata IDs; one-to-many ambiguous metadata matches are not acceptable.
unknown = sorted(set(data["Sample_ID"]) - set(metadata["Sample_ID"])) || Find call-table sample IDs without any metadata row. Matching is exact text after trimming, not FASTA accession-version normalization.
if unknown: || Stop when supplied metadata cannot cover one or more observed call samples.
unused = set(metadata["Sample_ID"]) - set(data["Sample_ID"]) || Find metadata-only samples with no call-table row.
if unused: || Warn about these unmatched metadata samples without interpreting their missing calls as mutation absence.
metadata = metadata.set_index("Sample_ID") || Make the unique sample ID the lookup index for efficient map-based joining.
for column in ["Collection_Date", "Location_Name"]: || Join only these supported sample metadata fields.
if column not in metadata: || Skip an optional metadata field that was not provided.
matched = data["Sample_ID"].map(metadata[column]).fillna("").astype(str).str.strip() || Look up the metadata value for every call row's sample ID, preserving the row order.
if column == "Collection_Date": || Parse/check date joins separately from textual location joins.
matched_dates = pd.to_datetime(matched.replace("", pd.NA), errors="coerce") || Convert joined date strings and preserve blanks as missing.
bad = matched.ne("") & matched_dates.isna() || Detect nonblank metadata date values that failed parsing.
if bad.any(): || Add metadata date errors to the combined validation report.
conflict = dates.notna() & matched_dates.notna() & dates.ne(matched_dates) || A date conflict requires two known dates that disagree; missing metadata is not itself a contradiction.
if conflict.any(): || Report recognized source-versus-metadata disagreements instead of arbitrarily choosing one.
dates = dates.fillna(matched_dates) || Fill only unavailable call-table dates from metadata; keep existing known source dates.
conflict = data[column].ne("") & matched.ne("") & data[column].str.casefold().ne(matched.str.casefold()) || Compare two known location strings case-insensitively. Internal-spacing variants are not normalized by this expression.
data[column] = data[column].where(data[column].ne(""), matched) || Keep nonblank source locations and fill blanks from joined metadata.
if invalid_source_dates.any(): || Report invalid nonblank source dates even if metadata happened to supply a usable date later.
for column in ["Sample_ID", "Reference_Genome", "Ref", "Alt"]: || Check mandatory nonblank identifiers and allele text fields.
count = int(data[column].eq("").sum()) || Count empty text values without silently removing their call rows.
if count: || Add a missing-value error when at least one required text field is blank.
for column in ["Ref", "Alt"]: || Validate each side of the proposed single-nucleotide substitution.
invalid = ~data[column].isin(["A", "C", "G", "T"]) || Require exactly one ordinary DNA base; ambiguous bases, gaps, and multi-base alleles are outside this loader's SNV representation.
if invalid.any(): || Record alphabet errors for the full combined report.
same_base = data["Ref"].eq(data["Alt"]) || Find records that propose no nucleotide change.
if same_base.any(): || Reject identical Ref and Alt because the row would not describe a substitution.
positions = pd.to_numeric(data["Position"], errors="coerce") || Convert reference-position text to numbers for validity checking.
invalid_position = positions.isna() | positions.le(0) | positions.ne(positions.round()) | positions.isin([float("inf"), float("-inf")]) || Reject missing, nonpositive, fractional, or infinite positions; a valid reference coordinate is a positive integer.
if invalid_position.any(): || Add the total invalid-position count to the errors.
frequencies = pd.to_numeric(data["Alt_Frequency"], errors="coerce") || Convert reported allele-frequency text to numerical fractions.
invalid_frequency = frequencies.isna() | ~frequencies.between(0, 1) || Reject unavailable or out-of-range fractions. The code does not automatically divide percentage inputs by 100.
if invalid_frequency.any(): || Add invalid frequency counts to the errors.
raw_depths = data["Read_Depth"].fillna("").astype(str).str.strip() || Preserve the distinction between blank depth and a nonblank malformed depth.
depths = pd.to_numeric(raw_depths.replace("", pd.NA), errors="coerce") || Parse depths as numbers while leaving blank depth missing.
invalid_depth_text = raw_depths.ne("") & depths.isna() || Identify nonblank depth text that could not be parsed.
if invalid_depth_text.any(): || Reject invalid text rather than relabel it as a permitted missing depth.
bad_depth = depths.notna() & (depths.lt(1) | depths.ne(depths.round()) | depths.isin([float("inf"), float("-inf")])) || A present depth must be a finite positive integer; missing depths are handled by a warning instead.
if bad_depth.any(): || Add invalid present-depth counts to the errors.
depth_missing = int(depths.isna().sum()) || Count unavailable depths explicitly.
if depth_missing: || Warn about missing coverage information rather than substituting zero.
if errors: || Halt after collecting all recognized validation issues.
data["Position"] = positions.astype(int) || Store checked reference coordinates as ordinary integers.
data["Alt_Frequency"] = frequencies.astype(float) || Store checked allele fractions as floats.
data["Read_Depth"] = depths.astype("Int64") || Use pandas' nullable integer dtype so valid whole-number depths coexist with missing values.
data["Collection_Date"] = dates || Store the parsed/joined datetime series.
values = data[column].replace("", pd.NA) if column == "Location_Name" else data[column] || Treat blank location labels as missing for within-sample consistency; dates already use datetime missing values.
conflicting_samples = values.groupby(data["Sample_ID"]).nunique(dropna=True).gt(1) || Flag a sample with more than one distinct known date or place across its SNV rows.
if conflicting_samples.any(): || Reject inconsistent sample-level metadata rather than choosing an arbitrary known value.
data[column] = values.groupby(data["Sample_ID"]).transform("first") || Broadcast the sample's first consistent known metadata value back to all its call rows.
data["Location_Name"] = data["Location_Name"].fillna("") || Restore blank display text when no location is known for the sample.
keys = ["Sample_ID", "Reference_Genome", "Position", "Ref", "Alt"] || Define repeated-call identity including sample and reference, not position alone.
duplicate_count = int(data.duplicated(keys).sum()) || Count later repeated occurrences; keep the first plus all duplicates in the retained table.
if duplicate_count: || Emit a warning without deleting repeated records.
missing_dates = int(data["Collection_Date"].isna().sum()) || Count retained call rows whose sample date remains missing.
missing_locations = int(data["Location_Name"].eq("").sum()) || Count retained call rows with no place label.
if missing_dates: || Report missing date records as a warning, not a false inferred date.
if missing_locations: || Report missing location records as a warning.
return data[OUTPUT_COLUMNS].copy(), warnings || Return the fixed standardized column set and warning list. Extra source annotations are not included in this organized output schema.
variant = ["Reference_Genome", "Position", "Ref", "Alt"] || Define a reference-aware variant independently of which sample reported it.
top = ( || Start the descriptive ranking of up to ten variants by how many samples and records report them.
data.groupby(variant, as_index=False) || Group identical reference/position/Ref/Alt combinations without merging different reference genomes.
.agg(Samples_Observed=("Sample_ID", "nunique"), Records=("Sample_ID", "size")) || Count distinct observed samples separately from retained call-row count, including repeats.
.sort_values(["Samples_Observed", "Records", "Position"], ascending=[False, False, True]) || Rank by descending distinct sample count, then records, then ascending reference position; this is not a significance test.
.head(10) || Limit this displayed top-variant table to ten rows without changing the underlying retained calls.
sample_table = ( || Build one summary row per observed wastewater sample.
data.groupby("Sample_ID", as_index=False) || Group call rows by their exact sample ID.
Collection_Date=("Collection_Date", "first"), || Use the sample's first consistent known date after loader checks/propagation.
Location_Name=("Location_Name", "first"), || Use the sample's consistent retained place label.
SNV_Records=("Position", "size"), || Count all retained call records for the sample, including repeated keys.
.sort_values("Sample_ID") || Put sample-summary rows in a predictable identifier order.
observed = data.drop_duplicates(["Sample_ID"] + variant) || Create a deduplicated view only for distinct sample/variant accounting; the original data is retained unchanged.
sample_table["Distinct_SNVs"] = sample_table["Sample_ID"].map(observed.groupby("Sample_ID").size()) || Map each sample to its number of distinct observed reference-aware variants.
sample_table["Collection_Date"] = sample_table["Collection_Date"].dt.strftime("%Y-%m-%d").fillna("") || Format summary dates for HTML, leaving unknown dates blank.
dates = data["Collection_Date"].dropna() || Use only known dates for the displayed date range.
depths = data["Read_Depth"].dropna() || Use only reported depths for depth distribution summaries; missing coverage is not treated as zero.
''')


def expression_text(node):
    text = ast.unparse(node)
    return text if len(text) <= 180 else text[:177] + "..."


def func_name(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    return "expression"


def automatic_statement(node, filename):
    """Explain Python mechanics precisely where a custom note is unnecessary."""
    if isinstance(node, (ast.Import, ast.ImportFrom)):
        descriptions = []
        for name in node.names:
            explanation = IMPORTS.get(name.name, "Imports the named helper from the indicated module so it can be called here.")
            alias = f" as {name.asname}" if name.asname else ""
            descriptions.append(f"{name.name}{alias}: {explanation}")
        return " ".join(descriptions)
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        info = FUNCTION_INFO.get(filename, {}).get(node.name)
        if info:
            return "Define " + node.name + "(): " + info[0]
        if node.name.startswith("test_"):
            return "Define this independent test case: " + node.name[5:].replace("_", " ") + ". unittest discovers methods beginning with test_."
        if node.name.startswith("fake_"):
            return "Define a temporary simulated external command for this test. It returns fixed outputs instead of starting real MAFFT."
        return "Define " + node.name + "(); its indented body runs when called, not merely because Python reads the definition."
    if isinstance(node, ast.ClassDef):
        return "Group related test methods in " + node.name + ". Inheriting unittest.TestCase provides assertions and automatic test discovery; production analysis modules do not introduce these test classes."
    if isinstance(node, ast.Assign):
        target = ", ".join(expression_text(item) for item in node.targets)
        if isinstance(node.value, ast.Dict):
            return "Create the dictionary assigned to " + target + "; following key/value entries define its fields."
        if isinstance(node.value, (ast.List, ast.Set, ast.Tuple)):
            return "Assign an ordered collection or set to " + target + ". Its following entries are data/configuration values rather than extra analysis steps."
        return "Evaluate the right-hand expression and store its result in " + target + ". " + call_meaning(node.value)
    if isinstance(node, ast.AugAssign):
        return "Update " + expression_text(node.target) + " using " + expression_text(node.value) + "; += increases the existing counter instead of replacing it with a fresh zero."
    if isinstance(node, ast.Return):
        if isinstance(node.value, ast.Dict):
            return "Return the following dictionary to the caller and finish this function; later key/value lines define the output fields."
        if isinstance(node.value, ast.Constant):
            return "Return " + repr(node.value.value) + " immediately. In this branch the remaining function body does not run."
        return "Return " + (expression_text(node.value) if node.value else "no explicit value") + " to the caller and finish this function."
    if isinstance(node, ast.Raise):
        return "Raise the named exception with a readable explanation. The normal return path stops; Flask or the caller can catch the error."
    if isinstance(node, ast.If):
        if expression_text(node.test) == "__name__ == '__main__'":
            return "Run the following startup block only when this file is executed directly; importing it does not launch the server/test runner again."
        return "Run the indented branch only if this condition is true: " + expression_text(node.test) + "."
    if isinstance(node, ast.For):
        return "Loop over " + expression_text(node.iter) + ", placing each next item in " + expression_text(node.target) + ". The indented body runs once per item."
    if isinstance(node, ast.With):
        names = [func_name(item.context_expr.func) if isinstance(item.context_expr, ast.Call) else expression_text(item.context_expr) for item in node.items]
        if any("patch" in expression_text(item.context_expr) for item in node.items):
            return "Temporarily replace the named dependency/environment/object for this test. The original is restored when this with block exits, even if the assertion fails."
        if any(name.startswith("assertRaises") for name in names):
            return "Expect the indented call to raise the specified exception; the test fails if it does not. assertRaisesRegex additionally matches text in the exception message."
        if "test_client" in names:
            return "Open Flask's simulated browser client so GET/POST requests can exercise routes without starting a public server."
        if "TemporaryDirectory" in names:
            return "Create a temporary folder for this block. Its files are removed when the block exits, including on errors."
        if any("LOCK" in name for name in names):
            return "Hold this thread lock while accessing shared process memory; the lock is released automatically afterward."
        return "Enter the context manager for " + ", ".join(names) + ", with automatic cleanup/restoration when the block ends."
    if isinstance(node, ast.Try):
        return "Start a protected block. If a listed exception occurs, control moves to the matching except branch instead of following the normal path."
    if isinstance(node, ast.ExceptHandler):
        return "Handle " + (expression_text(node.type) if node.type else "an exception") + ". This branch explains or recovers from that specific failure."
    if isinstance(node, ast.Continue):
        return "Skip the rest of this loop iteration and move to its next item; this does not stop the whole function."
    if isinstance(node, ast.Break):
        return "Stop the nearest loop immediately; the function continues with code after that loop."
    if isinstance(node, ast.Pass):
        return "Deliberately do nothing in this branch. An optional version lookup may fail while the successful alignment status is retained."
    if isinstance(node, ast.Delete):
        return "Remove this cache entry from process memory. This expires/evicts a temporary download, not an uploaded scientific record."
    if isinstance(node, ast.Expr):
        if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            return "Documentation string: records the module/function purpose in __doc__; it is not a scientific calculation."
        return call_meaning(node.value)
    return "Continue this Python statement; indentation determines which enclosing block it belongs to."


def call_meaning(value):
    if not isinstance(value, ast.Call):
        return "The expression combines the selected values using the operators visible in the code."
    name = func_name(value.func)
    if name.startswith("assert"):
        args = [expression_text(arg) for arg in value.args]
        if name == "assertEqual" and len(args) >= 2:
            return "Check that " + args[0] + " equals " + args[1] + ". A mismatch fails the test."
        if name in {"assertIn", "assertNotIn"} and len(args) >= 2:
            relation = "is present in" if name == "assertIn" else "is absent from"
            return "Check that " + args[0] + " " + relation + " " + args[1] + "."
        if name == "assertTrue":
            return "Require this condition to be true: " + (args[0] if args else "the supplied expression") + "."
        if name == "assertIsNone":
            return "Require the result to be None, representing an unavailable optional value."
        if name == "assertIs" and len(args) >= 2:
            return "Check object identity: " + args[0] + " must be the same object as " + args[1] + ", not merely an equal-looking value."
        return "Run unittest assertion " + name + " to check the stated expected behavior; failure marks the test as failed."
    meanings = {
        "assert_not_called": "Verify the mocked dependency was never invoked; this protects the separation from the Sprint 2 pipeline or an unnecessary external call.",
        "assert_called_once": "Verify exactly one call occurred, which checks request reuse or an independent dependency check.",
        "write_text": "Write this fixed fixture/synthetic text to the specified temporary file using the stated encoding.",
        "append": "Append this item to the existing list; warnings/errors are accumulated so the caller can report them together.",
        "add": "Add this value to a set; repeated occurrences do not create extra distinct items.",
        "update": "Add these values to the existing set/dictionary while retaining previous entries.",
        "raise_for_status": "Turn an HTTP error status into a requests.HTTPError; a successful HTTP response can still contain a separate GeoNames API error in its JSON.",
        "save": "Save the uploaded file bytes under the fixed temporary filename selected by the server.",
        "exception": "Log the unexpected problem on the server; the route returns a shorter public-facing error instead of the traceback.",
        "warning": "Log a generic diagnostic warning without placing credential-bearing exception details in the response.",
        "main": "Run unittest's test discovery/runner when this test file is executed directly.",
        "run": "Start the requested external process or local Flask server with the explicit arguments shown here.",
        "plot": "Draw the selected observed date/value points and connecting line; this is a descriptive plot, not a statistical trend test.",
        "text": "Put an explanatory label on the empty plot rather than inventing measurements.",
        "barh": "Draw horizontal bars for the selected location means; no new measurements or thresholds are calculated by the chart.",
        "savefig": "Encode the figure into the in-memory PNG buffer at the chosen display resolution.",
        "abort": "Return the HTTP error status and stop this request, typically because a temporary download is unavailable.",
    }
    return meanings.get(name, "Call " + expression_text(value.func) + " with the shown arguments. The surrounding assignment/return determines where its result is used.")


def syntax_note(text, owner_line):
    stripped = text.strip()
    if stripped in {"else:", "elif"} or stripped.startswith("elif "):
        return "Alternative branch: evaluate this only when the preceding if/elif condition did not choose a branch."
    if stripped == "else:":
        return "Run this alternative block when the preceding condition is false."
    if stripped.startswith("except "):
        return "Handle the listed exception from the preceding try block."
    if stripped.startswith("#"):
        return "Comment for the reader: " + stripped.lstrip("# ") + ". Python does not execute the comment."
    if stripped.startswith("@"):
        return "Decorator: registers the following function as a Flask route/error handler with the stated path, methods, or status code."
    if stripped in {"}", "]", ")", "},", "],", "),", "):", "),", "}", "])"}:
        return "Close the multi-line collection, call, or condition begun above. This punctuation completes the expression; it adds no separate operation."
    if stripped.startswith(("and ", "or ")):
        relation = "Also require this condition" if stripped.startswith("and ") else "Accept this additional alternative"
        return relation + "; it continues the logical test begun above. Python short-circuits and/or tests when the result is already determined."
    if stripped.startswith('"') or stripped.startswith("'") or stripped.startswith('f"'):
        if ":" not in stripped.split("=", 1)[0]:
            return "Continue the string/fixture/message begun above. Adjacent quoted literals are concatenated; \\n encodes a newline and is not a new source-code line."
    return f"Continue the multi-line statement begun at line {owner_line}. This line supplies part of that same expression, not a separate Python step."


def annotate_file(filename):
    source = SOURCES[filename]
    lines = source.splitlines()
    tree = TREES[filename]
    starts = {}
    owners = {}
    extra = {}
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.stmt, ast.ExceptHandler)) and hasattr(node, "lineno"):
            starts.setdefault(node.lineno, []).append(node)
            for line in range(node.lineno, getattr(node, "end_lineno", node.lineno) + 1):
                previous = owners.get(line)
                if previous is None or getattr(node, "end_lineno", node.lineno) - node.lineno < getattr(previous, "end_lineno", previous.lineno) - previous.lineno:
                    owners[line] = node
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
            docstrings.update(range(node.lineno, node.end_lineno + 1))
        if isinstance(node, ast.Dict):
            for key, value in zip(node.keys, node.values):
                if isinstance(key, ast.Constant) and isinstance(key.value, str):
                    extra.setdefault(key.lineno, []).append(KEY_NOTES.get(key.value, "Define output/configuration field '" + key.value + "' using the expression on its right."))
        if isinstance(node, ast.keyword):
            keyword_notes = {
                "capture_output": "Capture standard output and standard error in the CompletedProcess object rather than printing them directly.",
                "text": "Decode process output as text, so StringIO/FASTA parsing receives strings instead of bytes.",
                "check": "If True, a nonzero external exit status raises CalledProcessError instead of being accepted as valid output.",
                "timeout": "Set the diagnostic/request timeout; this is operational protection, not a biological threshold.",
                "verify": "Keep HTTPS certificate verification enabled for GeoNames.",
                "params": "Supply fixed query fields and the environment-provided username to the secure GeoNames endpoint; the username is not returned to the webpage.",
                "encoding": "Choose the file/byte encoding explicitly; UTF-8 preserves text identifiers consistently.",
                "dtype": "Choose pandas column types intentionally, often strings to preserve IDs or nullable integers for missing depths.",
                "errors": "With errors='coerce', unparseable numeric/date values become missing so validation can find and report them.",
                "format": "Choose the required date representation or PNG output format for the surrounding call.",
                "index": "With index=False, omit pandas' extra row-number index from the written CSV/XLSX file.",
                "as_index": "With as_index=False, grouping keys stay ordinary result columns rather than becoming a pandas row index.",
                "mimetype": "Label the HTTP download as CSV for the browser.",
                "as_attachment": "Ask the browser to download the response as a file.",
                "download_name": "Set the fixed filename suggested to the browser for this temporary CSV.",
                "sort": "Control group order for predictable output; sorting changes presentation order, not genomic similarity.",
                "keep_default_na": "False preserves textual labels such as 'NA' instead of letting pandas treat them automatically as missing.",
            }
            if node.arg in keyword_notes:
                extra.setdefault(node.value.lineno, []).append(keyword_notes[node.arg])
    result = {}
    custom = LINE_NOTES.get(filename, {})
    for number, line in enumerate(lines, 1):
        if not line.strip():
            continue
        stripped = line.strip()
        if stripped in custom:
            explanation = custom[stripped]
        elif number in docstrings:
            explanation = "Documentation string: describes the purpose of this module/function. It is retained as __doc__, not executed as an analysis instruction."
        elif number in starts:
            nodes = starts[number]
            node = min(nodes, key=lambda item: getattr(item, "end_lineno", item.lineno) - item.lineno)
            explanation = automatic_statement(node, filename)
            if explanation.startswith("Evaluate the right-hand") or explanation.startswith("Run the indented branch") or explanation.startswith("Call "):
                MISSING_NOTES.append((filename, number, stripped, explanation))
        elif number in extra:
            explanation = " ".join(dict.fromkeys(extra[number]))
        else:
            owner = owners.get(number)
            explanation = syntax_note(line, owner.lineno if owner else number)
        result[number] = explanation
    return result


def paragraph(text, style="GuideBody"):
    return Paragraph(escape(text), STYLES[style])


class GuideDocument(BaseDocTemplate):
    def __init__(self, path):
        super().__init__(str(path), pagesize=letter, leftMargin=38, rightMargin=38, topMargin=45, bottomMargin=43, title="BMI 461 - Complete New Python Code Guide", author="BMI 461 project documentation")
        self.chapter = "Classroom reference"
        frame = Frame(self.leftMargin, self.bottomMargin, self.width, self.height, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        self.addPageTemplates(PageTemplate(id="main", frames=frame, onPage=self.draw_page))

    def draw_page(self, canvas, doc):
        canvas.saveState()
        canvas.setFillColor(NAVY)
        canvas.rect(0, 780, 612, 12, fill=1, stroke=0)
        canvas.setFont("GuideSans", 7.7)
        canvas.setFillColor(GRAY)
        canvas.drawString(38, 765, "BMI 461  |  PYTHON CODE WALKTHROUGH  |  19 OCTOBER 2026")
        canvas.setStrokeColor(RULE)
        canvas.line(38, 35, 574, 35)
        canvas.drawString(38, 23, self.chapter[:88])
        canvas.drawRightString(574, 23, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph) and hasattr(flowable, "bookmark_key"):
            key = flowable.bookmark_key
            text = flowable.getPlainText()
            self.canv.bookmarkPage(key)
            level = getattr(flowable, "outline_level", 0)
            self.canv.addOutlineEntry(text, key, level=level, closed=True)
            if level == 0:
                self.chapter = text
                self.notify("TOCEntry", (0, text, self.page, key))


STORY = []


def add(text, style="GuideBody"):
    STORY.append(paragraph(text, style))


def heading(text, key, level=0):
    item = paragraph(text, "GuideH1" if level == 0 else "GuideH2")
    item.bookmark_key = key
    item.outline_level = level
    STORY.append(item)


def callout(title, text):
    table = Table([[paragraph(title, "GuideH2")], [paragraph(text)]], colWidths=[536])
    table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), PALE), ("LEFTPADDING", (0, 0), (-1, -1), 11), ("RIGHTPADDING", (0, 0), (-1, -1), 11), ("TOPPADDING", (0, 0), (-1, -1), 6), ("BOTTOMPADDING", (0, 0), (-1, -1), 6)]))
    STORY.extend([table, Spacer(1, 10)])


def code_cell(text):
    pieces = wrap(text.expandtabs(4), width=56, replace_whitespace=False, drop_whitespace=False, break_long_words=True, break_on_hyphens=False) or [""]
    formatted = "<br/>".join(escape(piece).replace(" ", "&nbsp;") for piece in pieces)
    return Paragraph(formatted, STYLES["GuideCode"])


def line_table(filename, start, end, annotations):
    lines = SOURCES[filename].splitlines()
    rows = [[paragraph("Line", "GuideTableHead"), paragraph("Exact source", "GuideTableHead"), paragraph("Meaning and importance", "GuideTableHead")]]
    for number in range(start, end + 1):
        if number not in annotations:
            continue
        rows.append([paragraph(str(number), "GuideLine"), code_cell(lines[number - 1]), paragraph(annotations[number], "GuideCell")])
    if len(rows) == 1:
        return
    table = LongTable(rows, colWidths=[29, 229, 278], repeatRows=1, hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, TEAL),
        ("LINEBELOW", (0, 1), (-1, -1), 0.25, RULE),
        ("BACKGROUND", (0, 1), (1, -1), colors.HexColor("#f7f9fb")),
    ]))
    STORY.append(table)


def sections(filename):
    tree = TREES[filename]
    size = len(SOURCES[filename].splitlines())
    groups = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            start = min([node.lineno] + [decorator.lineno for decorator in node.decorator_list])
            groups.append((start, node.end_lineno, node.name, node))
        elif isinstance(node, ast.ClassDef):
            methods = [item for item in node.body if isinstance(item, ast.FunctionDef)]
            if methods:
                groups.append((node.lineno, methods[0].lineno - 1, node.name, node))
            for method in methods:
                groups.append((method.lineno, method.end_lineno, method.name, method))
    result = []
    cursor = 1
    for start, end, name, node in groups:
        if start > cursor:
            result.append((cursor, start - 1, "Module setup / continued source", None))
        result.append((start, end, name, node))
        cursor = end + 1
    if cursor <= size:
        result.append((cursor, size, "Standalone entry point / continued source", None))
    return result or [(1, size, "Package/module setup", None)]


def build_guide():
    add("BMI 461  /  CURRENT PROJECT SNAPSHOT", "GuideDeck")
    add("Understand the New Python Code", "GuideTitle")
    add("A line-by-line classroom reference for the updated genomic-surveillance backbone, server diagnostics, processed environmental measurements, and processed wastewater SNVs.", "GuideDeck")
    STORY.append(HRFlowable(width="100%", thickness=2, color=TEAL, spaceAfter=18))
    total = sum(len(source.splitlines()) for source in SOURCES.values())
    nonblank = sum(sum(bool(line.strip()) for line in source.splitlines()) for source in SOURCES.values())
    callout("Prepared for a live code review", f"13 Python files | {total:,} physical source lines | {nonblank:,} individually annotated nonblank lines. The current app.py is included in full because it connects the additions. Source files are unchanged by this documentation.")
    add("Prepared: 19 October 2026", "GuideBody")
    add("This guide describes the files actually present in the working folder, including uncommitted environmental and wastewater additions. It preserves the distinction between new modules and updates to an existing module. It is an explanation of a code snapshot, not a claim that the current deployment or every possible dataset has been verified.", "GuideBody")
    add("Classroom route: read the architecture and vocabulary first; use function summaries to explain intent; then follow the numbered rows while opening the same file in your editor. The PDF has clickable contents and bookmarks.", "GuideBody")
    STORY.append(PageBreak())
    heading("Contents", "contents")
    toc = TableOfContents()
    toc.levelStyles = [STYLES["GuideTOC"]]
    toc.dotsMinLevel = 0
    STORY.append(toc)
    STORY.append(PageBreak())
    heading("Reading the guide and the project", "reading")
    for title, text in INTRODUCTION:
        add(title, "GuideH2")
        add(text)
    STORY.append(PageBreak())
    heading("How the modules connect", "architecture")
    for title, text in ARCHITECTURE:
        add(title, "GuideH2")
        add(text)
    add("Files treated as new since the recorded baseline", "GuideH2")
    for filename in FILES:
        info = FILE_INFO[filename]
        add(filename + " - " + info["role"], "GuideSmall")
    for index, filename in enumerate(FILES, 1):
        STORY.append(PageBreak())
        info = FILE_INFO[filename]
        heading(f"{index}. {filename}", "file_" + str(index))
        add(info["summary"])
        add("Role: " + info["role"], "GuideSmall")
        add(f"Snapshot: {len(SOURCES[filename].splitlines())} physical lines; SHA-256 prefix {FINGERPRINTS[filename][:16]}. Exact code and comments are shown in the left column. Wrapped text keeps its original line number.", "GuideSmall")
        for label, text in info.get("examples", []):
            add(label, "GuideH2")
            add(text)
        for text in info.get("limits", []):
            add("Interpretation note: " + text, "GuideSmall")
        annotations = annotate_file(filename)
        for start, end, name, node in sections(filename):
            if not any(number in annotations for number in range(start, end + 1)):
                continue
            add(name + ("()" if isinstance(node, ast.FunctionDef) else ""), "GuideH2")
            if isinstance(node, ast.FunctionDef):
                info_fn = FUNCTION_INFO.get(filename, {}).get(name)
                if info_fn:
                    for label, text in zip(("Purpose", "Inputs and output", "Why it matters"), info_fn):
                        add(label + ": " + text, "GuideSmall")
                elif name.startswith("test_"):
                    add("Test aim: " + name[5:].replace("_", " ") + ". Read the fixture first, then the action, then the assertions. Mocked successes verify program behavior; they do not establish live server availability.", "GuideSmall")
            line_table(filename, start, end, annotations)
    STORY.append(PageBreak())
    heading("Worked examples for a professor discussion", "examples")
    for title, text in WORKED_EXAMPLES:
        add(title, "GuideH2")
        add(text)
    STORY.append(PageBreak())
    heading("Questions to prepare for class", "questions")
    for question, answer in PROFESSOR_QUESTIONS:
        add(question, "GuideH2")
        add(answer)
    STORY.append(PageBreak())
    heading("Verification, source identity, and further reading", "verification")
    for text in VERIFICATION_NOTES:
        add(text)
    for filename in FILES:
        add(filename, "GuideH2")
        add(FINGERPRINTS[filename], "GuideSmall")
    DESTINATION.parent.mkdir(parents=True, exist_ok=True)
    QA_FOLDER.mkdir(parents=True, exist_ok=True)
    document = GuideDocument(DESTINATION)
    document.multiBuild(STORY)
    for filename in FILES:
        assert hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == FINGERPRINTS[filename], filename + " changed during documentation"
    (QA_FOLDER / "annotation_audit.json").write_text(json.dumps({"files": FILES, "source_lines": total, "annotated_nonblank_lines": nonblank, "generic_statement_notes": MISSING_NOTES, "fingerprints": FINGERPRINTS}, indent=2), encoding="utf-8")
    print(DESTINATION)
    print(f"Physical lines: {total}; annotated nonblank lines: {nonblank}; generic statement notes to review: {len(MISSING_NOTES)}")


if __name__ == "__main__":
    build_guide()
