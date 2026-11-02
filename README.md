# BMI 461 Integrated Genomic Surveillance

A classroom Flask application that analyzes **two kinds of genomic evidence
in one workflow**: individual (consensus) sequences and pooled-sample
sequencing results such as wastewater SNV tables. Both are placed on the same
reference-genome coordinates when the reference is verified, so the
dashboard can show sequencing quality, genomic variation, temporal trends,
population structure (PCoA), geographic context, cross-source variant
overlap, and the original Sprint 2 sequence-relationship screening.

The application can run on a personal computer or in the supplied Docker
container. This README describes the code on the `sprint4-integration`
branch as of November 1, 2026. A potential outbreak cluster, a shared variant,
or a variant-profile overlap is a screening result; none of them confirms
transmission, an outbreak, or the lineage present in a mixed sample.

> **Integration summary:** see [Integrated analysis (Sprint 2 and Sprint 3)](#integrated-analysis-sprint-2-and-sprint-3)
> below, `INTEGRATION_IMPLEMENTATION_LOG.md` for what was implemented, and
> `ARTICLE_METHOD_VALIDATION.md` for which methods of Fontenele et al. (2021)
> are reproduced, partially reproduced, or exploratory, with a comparison
> against the paper's published supplementary data.

## Start locally on Windows

These instructions use **Windows Command Prompt**. Replace the example project
path with the folder that contains `app.py` and `requirements.txt`.

### First-time setup

1. Install [Python 3.11](https://www.python.org/downloads/) if it is not already
   available. Confirm the installation:

   ```bat
   py -3.11 --version
   ```

2. Download and extract the complete project, or clone its GitHub repository.
   Keep the Python modules, `templates`, `examples`, and `wastewater` together.
   Open Command Prompt in that folder.

3. Create a virtual environment, install the packages, and start Flask:

   ```bat
   cd /d "C:\path\to\project"
   py -3.11 -m venv .venv
   .venv\Scripts\activate.bat
   python -m pip install --upgrade pip
   python -m pip install -r requirements.txt
   python app.py
   ```

4. Open a browser at:

   ```text
   http://127.0.0.1:5000
   ```

The server stays running in that terminal. Press `Ctrl+C` to stop it.
`python app.py` does not open the browser automatically.

The website and processed wastewater pages can start without MAFFT or a
GeoNames account. **Install MAFFT before running the individual-genome
analysis.** GeoNames is optional. Both are explained below.

If the `py` launcher is unavailable but `python --version` reports Python 3.11,
use `python -m venv .venv` instead. Create the environment on each computer;
do not copy another computer's `.venv` folder.

### Starting the project again

After the first installation, use:

```bat
cd /d "C:\path\to\project"
.venv\Scripts\activate.bat
python app.py
```

If using a session-only MAFFT path or GeoNames environment variable, set it
again in this same terminal before starting Flask. After downloading a new
version of the project, reinstall any updated requirements:

```bat
python -m pip install -r requirements.txt
```

### Windows PowerShell

PowerShell can run the environment's Python directly without activating a
script or changing the execution policy:

```powershell
Set-Location "C:\path\to\project"
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe app.py
```

For later runs, only the folder change and final command are needed.

### macOS or Linux

Use an installed Python 3.11 interpreter:

```sh
cd "/path/to/project"
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python app.py
```

If `python3` is already Python 3.11, it can replace `python3.11`.
The local browser address is the same on all platforms.

Later `python` commands assume the virtual environment is active. If you used
the no-activation approach, replace `python` with
`.\.venv\Scripts\python.exe` on Windows or `.venv/bin/python` on macOS/Linux.

## Pages and a first demonstration

| Page | Local address | Input and purpose |
|---|---|---|
| **Integrated dashboard (primary)** | `http://127.0.0.1:5000/` | Guided workflow: upload any mix of files, review what each one was recognized as, then run every applicable analysis; individual sequences need MAFFT |
| Wastewater Data Analysis (Sprint 3, legacy) | `http://127.0.0.1:5000/wastewater-data` | One processed feature-measurement CSV/XLSX |
| Wastewater Surveillance (Sprint 3, legacy) | `http://127.0.0.1:5000/wastewater` | Processed SNV CSV/XLSX, with optional sample metadata |
| Earlier environmental page | `http://127.0.0.1:5000/environmental` | Same measurement analysis as Wastewater Data Analysis |
| Health check | `http://127.0.0.1:5000/health` | Returns `OK`; performs no external dependency tests |

The main page also has **Run Integrated Example** (needs MAFFT), which runs
the entirely synthetic paired demonstration described in
[Integrated analysis (Sprint 2 and Sprint 3)](#integrated-analysis-sprint-2-and-sprint-3). Pooled-only
uploads (an SNV table without FASTA files) do not need MAFFT.

For a first local demonstration of a Sprint 3 page without installing MAFFT:

1. Open `/wastewater-data`.
2. Click **Run Synthetic Example**.
3. Inspect sample counts, features, tables, and charts.
4. Select **Download organized measurements (CSV)**.

The example uses `examples/synthetic_environmental_measurements.csv` and is
clearly labeled synthetic. It contains invented samples and frequencies;
coverage was not supplied, so the page explains that it is unavailable.

After installing MAFFT, the main page also has **Run Example Analysis**.
That uses `examples/synthetic_sequences.fasta` and
`examples/synthetic_metadata.csv` through the six-stage genomic pipeline.
One example record has a real place name with synthetic sample information
so optional GeoNames lookup can be demonstrated. All example records are
classroom data rather than actual outbreaks.

## Integrated analysis (Sprint 2 and Sprint 3)

This section describes the integration of the Sprint 2 individual-sequence
pipeline with the Sprint 3 pooled-sequencing work. It is an update to those
sprints, not Sprint 4 work. The optional Sprint 4 AI literature and metadata
enrichment (dashboard section 9, disabled for paid requests by default) is
documented in `SPRINT4_AI_NLP_INTEGRATION_REPORT.md`; this README will be
updated for it after manual validation.

### Architecture

```text
Individual FASTA + metadata          Reference FASTA          Pooled SNV table (+ optional metadata,
        |                                  |                   per-position depth/A-C-G-T counts,
        v                                  |                   masked sites, variant profiles)
 Sprint 2 pipeline (unchanged):            |                              |
 validate -> GeoNames -> MAFFT ->          |                              v
 p-distance -> time/place -> clusters      |            wastewater/processed_snv.py (Sprint 3 importer,
        |                                  |            now keeps published annotations, reads iVar TSV)
        |   validated, geocoded samples    |                              |
        +--------------------------------> genomics/  <-------------------+
                                           references.py         verified reference identity (SHA-256)
                                           sequence_variants.py  MAFFT to reference -> 1-based variants
                                           pooled_observations.py shared variant/sample tables, coverage
                                           quality_control.py    paper preset or exploratory flags
                                           temporal_analysis.py  per-site frequency trends, regression
                                           diversity_analysis.py Yue-Clayton PCoA (paper) / Jaccard (exploratory)
                                           variant_matching.py   cross-source categories, profiles
                                           figures.py            PNG figures (no Flask code)
                                           integrated_analysis.py orchestration + provenance
                                                      |
                              app.py + integrated_view.py + templates/index.html
                              (one dashboard, CSV/PNG downloads)   run_integrated_analysis.py (command line)
```

Both evidence types use one shared variant table (`genomics/unified_schema.py`):
`Variant_Key` (= `Reference_ID:Position:Ref>Alt`), `Reference_ID`, 1-based
`Position`, `Ref`, `Alt`, `Sample_ID`, `Sample_Type` (individual sequence or
pooled sample), `Observation_Type` (consensus base call or allele frequency),
`Alt_Frequency`, `Read_Depth`, `Alt_Depth`, `Callable_Status`,
`Quality_Status`, date, location, and `Data_Source`. Source-provided extras are
kept as `Published_*` (paper annotations) or `Source_*` (variant-caller output)
columns, separate from values this application calculates.

Three modules support the guided upload: `genomics/file_formats.py` (gzip,
FASTA/table/VCF recognition, VCF conversion), `genomics/column_detection.py`
(flexible column names with explicit ambiguity), and
`genomics/input_inspection.py` (file roles, dataset review, analysis plan, and
standardized copies for the existing loaders).

### Guided workflow: upload, review, run

1. **Upload files.** Select every file at once (up to 20 files, 10 MB in total).
   Gzip-compressed files are decompressed with a 60 MB size limit.
2. **Review dataset.** Each file is inspected by content and shown with its
   detected format, record count, sample IDs, reference IDs, and a suggested
   role with the evidence for it (for example "5 of 5 values in column
   'Sample' match the FASTA records", or "its ID SYNREF.1 is the reference
   named in pooled.csv"). When the evidence is not decisive (a one-sequence
   FASTA, metadata whose IDs match no other file) the role is left blank and
   must be chosen. Every table's detected columns are listed; when two columns
   could be the same field (for example `Date` and `Collection_Date`) the run is
   blocked until the user picks one. The review also shows reference
   compatibility, metadata completeness, missing optional information, required
   information that is still missing, and which analyses will run, run with
   limits, or be skipped (and why).
3. **Run analysis.** Advanced thresholds and the QC mode are in a collapsible
   section. Only applicable analyses run; the results start with a list of what
   ran and what was skipped.

Uploaded files wait in a temporary folder for at most 30 minutes between
review and run (at most six upload sets per server process) and are never
stored permanently. **Run Integrated Example** remains on the first page.
The named-field form fields (`fasta_file`, `metadata_file`, `reference_file`,
`snv_file`, ...) are still accepted by `POST /` for scripts and tests.

### Supported file formats

| Format | Use | Support |
|---|---|---|
| FASTA `.fasta`, `.fa`, `.fna` (and `.gz`) | Individual sequences; a one-record FASTA as the reference | Full. Content is checked, so FASTA in a `.txt` file is also recognized (and noted). |
| CSV, TSV/TXT (tab), XLSX first sheet; CSV/TSV/TXT may be `.gz` | Metadata, pooled SNV tables, per-position tables, masked sites, profiles | Full. `.xlsx.gz` is not needed (XLSX is already compressed). |
| Fontenele et al. (2021) Supplementary Table 2 (`mmc2.xlsx`) | Pooled SNV table with published annotations | Full; verified against the real workbook (see `ARTICLE_METHOD_VALIDATION.md`). |
| iVar `variants` TSV | Pooled SNV table | Full for SNVs; insertions/deletions are set aside and counted. Without a `SAMPLE` column, one file = one sample named after the file (shown on the review page). Several iVar files are combined. |
| VCF `.vcf`, `.vcf.gz` | Pooled SNV calls | Partial by design: single-nucleotide alleles only. Frequencies come only from allele-specific values: `FORMAT/AD` (one value per allele), `FORMAT/AO`+`RO`, `FORMAT/AF`, iVar `FORMAT/ALT_FREQ` (single ALT), or `INFO/AF` only for a sites-only file (LoFreq style; with sample columns INFO/AF can be genotype-based). With allele counts, depth is their total, so frequency = alternate reads / depth. Multi-allelic records are split per allele only when such values exist. Genotype-only (`GT`) entries, indels, symbolic alleles, and entries with no alternate reads are set aside and counted; they are never turned into 0% frequencies. `FILTER` is recorded, not applied. VCF has no dates or sites, so add pooled metadata for those analyses. A sites-only VCF uses its file name as the sample ID (shown). |
| FASTQ (`.fastq`, `.fq`), SAM/BAM/CRAM | Raw reads / alignments | **Not accepted** as results. They are recognized and refused with a pointer to `offline_workflow/`. |
| Other files | — | Refused with a message. |

### Required versus optional information

Column headers are matched without regard to case, spaces, underscores,
hyphens, or dots, using role-specific names, for example: `Sample_ID` /
`Sample` / `SampleName` / `Sample Name`; `Accession`; `Collection_Date` /
`Date` / `Sample_Date`; `Location_Name` / `Location` / `Geo_Location` / `Site`;
`Latitude` / `Lat`; `Longitude` / `Lon` / `Long`; `Reference_Genome` /
`Reference` / `Reference_ID` / `CHROM`; `Read_Depth` / `Depth` / `DP`;
`Alt_Frequency` / `AF` / `Frequency`. Units are never rescaled: coverage
written as a fraction (all values ≤ 1) or frequencies written as percentages
(values > 1) block the run with an explanation. Explicit FASTA description
annotations such as `[collection_date=2020-06-01]`, `[organism=...]`,
`[country=...]`, and `[host=...]` fill a metadata field only when the metadata
file has no such column; the source is recorded. Nothing is inferred from file
names (except the visible fallback sample ID above) or from sequence length.

| Analysis | Needs | When it is missing |
|---|---|---|
| Sequence comparison (MAFFT, p-distance) | FASTA IDs, valid sequences, a pathogen (metadata column, FASTA annotation, or stated once on the review page), a genomic region (metadata, FASTA description, or confirmed complete genomes) | Not run; the review page says what is missing |
| Temporal screening of pairs | Collection dates (full dates; a year-only or year-month date is not given a day) | Pairs without two dates are compared genetically only |
| Geographic-distance screening | Coordinates, or specific place names GeoNames can locate | State-level or missing locations stay "unresolved", never 0 km |
| Reference-coordinate variants and cross-source comparison | A reference FASTA; pooled records on the same verified reference | Skipped, or records kept apart as incomparable |
| Variant-frequency analysis | Interpretable 0–1 frequencies | Records without one are rejected by the importer |
| Coverage-aware non-detection | Per-position depth | Unreported variants stay "coverage unknown", never 0% |
| Paper Yue–Clayton PCoA | Per-position A/C/G/T counts | A labeled exploratory ordination is shown instead |

**Metadata rule change (documented and tested).** Individual metadata now
requires only `Accession` and `Pathogen`. A missing `Collection_Date`,
`Location_Name`, or `Host_Source` column, or a blank value, is a visible note
instead of an error, so the genetic comparison still runs. Unparseable dates
remain errors. Year-only or year-month dates (for example `2020` or
`2020-06`) are no longer read as 1 January / day 1; they are treated as
unknown for day-level comparisons and listed in the notes.

### Dashboard sections

Integrated sections come first, then source-specific ones. Every section has a
purpose statement, summary statistics ("Not available" is shown instead of a
misleading 0), figures with captions and PNG downloads, tables with CSV
downloads, and an "Uncertainty and limitations" list.

| Section | Figures |
|---|---|
| 1. Overview | Metadata completeness by source; how the two sources were connected (reference checksum, mapping rule, label and Ref-base checks, date overlap); analyses run or skipped |
| 2. Cross-source evidence | Category counts; variant-overlap matrix (variant × sample, individual and pooled side by side, distinguishing detection, non-detection with adequate depth, unknown coverage, and QC exclusion); reference positions by category; when each variant was seen in each source |
| 3. Genomic variation | Variants along the reference (Fig 3 analogue): samples per position, alternative base by color, source by marker shape, multi-allele positions, published novelty flags, annotation bands |
| 4. Temporal patterns | Sampling dates and sites for both sources at their real dates; frequency trajectories; first-versus-last frequency changes with day spans; missing and unknown observations per sample; sample-by-variant heatmap |
| 5. Geographic context | Location precision by source (source-provided coordinates, GeoNames estimates, broad administrative areas, named sites, none) and a plain coordinate plot (not a map) |
| 6. Sequencing quality | Coverage by sample; which QC checks could be applied (pass / fail / not verifiable); read-depth distributions; callable and uncallable positions in individual sequences; SNVs per sample (Fig 2A); depth versus SNVs (Fig 2B) |
| 7. Population structure | PCoA colored by date with site (or state) markers and an eigenvalue panel; paper method or labeled exploratory method |
| 8. Individual sequence relationships | The original Sprint 2 results, plus a genetic-distance distribution with genetic-versus-temporal distance and a relationship network (solid = all three cutoffs, dashed = geography unresolved; never "transmission") |

### Rules that keep the comparison honest

- Positions are compared only on the same reference. Two different accession
  labels count as the same reference only when both sequences are supplied
  and identical (checked by SHA-256). A pooled record whose `Ref` base
  disagrees with the supplied reference is marked incomparable.
- Individual sequences are aligned **to the reference**, so alignment columns
  are never used as positions. Partial sequences keep their real span; `N`,
  IUPAC ambiguity codes, deletions, and positions outside the span are not
  callable; insertions/deletions are listed but not compared. A sequence that
  does not map (fewer than 20 callable positions or below 90% identity — a
  prototype safeguard, not a paper rule) keeps its Sprint 2 distances but gets
  no variant profile.
- An unreported pooled variant is never set to 0%. With a per-position depth
  table it becomes *adequate depth, not reported*, *low coverage*, or *not
  covered*; without one it stays *coverage unknown*. A reported call that
  failed QC is shown as *excluded by QC*, not as a non-detection.
- Cross-source categories: **Shared**, **Individual only** (not detected in
  pooled samples with adequate depth), **Pooled only** (callable individual
  sequences show the reference base), **Unresolved** (the other source lacks
  coverage), **Incomparable** (reference not verified).
- Pooled site labels are not geocoded and are never treated as precise points;
  state-only individual locations keep the existing coarse-location rules.
- Variant-profile overlap is reported as overlap only, never as a lineage call.

### Fontenele et al. (2021) methods: reproduced, applied, or unavailable

The full method-to-implementation table (what the study did, what is
implemented, faithful / partial / exploratory status, required inputs, tests,
and limits for Methods 2.1–2.6 and Figures 1–5) is in
`ARTICLE_METHOD_VALIDATION.md`. Summary:

| Paper method | In this project |
|---|---|
| >90% genome coverage at 10x for detailed analysis | Applied when genome coverage is supplied in metadata or calculable from a complete per-position depth table |
| Variant depth >= 20x | Applied to supplied `Read_Depth` |
| Fisher's exact test p < 0.05 (iVar) | Applied to supplied p-values or iVar PASS flags; never recalculated or invented |
| No minimum allele frequency | Matched |
| Minimum base quality 20 | Not verifiable from processed records; applied when converting `samtools mpileup` text offline |
| Problematic-site masking (De Maio et al. 2020) | Applied only to an uploaded site list; the historical list is not bundled |
| SNVs per sample (Fig 2A) | Reproduced as a per-site bar chart |
| SNV count vs mean depth with 95% CI (Fig 2B) | Least-squares line, seeded bootstrap 95% band, permutation p-value (model details not stated in the paper); needs mean depth |
| SNVs by genome position (Fig 3) | Reproduced; published novelty flags shown only when present in the source table |
| Yue–Clayton PCoA (Fig 5) | Reproduced exactly when per-position A/C/G/T counts are supplied (per-position dissimilarity summed over positions; zero-depth positions contribute 0, as stated in the paper). Otherwise a separately labeled exploratory Jaccard ordination of reported SNV sets |
| BWA-MEM / iVar read processing | Not run here; command template and tested converters in `offline_workflow/` |
| GISAID comparison and PANGOLIN lineage support | Not reproduced (requires GISAID access); optional documented variant profiles are compared as overlap only |

The quality-control mode can be **paper** (failures excluded from analyses,
still listed in downloads) or **exploratory** (failures flagged only). The
dashboard shows a method-reproduction table for every run.

### Run the synthetic integrated example

Every file in `examples/integrated/` is invented; it was generated by
`examples/integrated/make_integrated_example.py`, which also lists the
expected category for every planted variant. With MAFFT installed:

```console
python app.py
```

Open `http://127.0.0.1:5000`, click **Run Integrated Example**, and use the
section links (Overview through Individual sequence relationships). Or, from
the command line:

```console
python run_integrated_analysis.py --example --output-folder results/integrated_example
```

Expected highlights: 2 shared variants, 1 individual-only, 1 pooled-only,
1 unresolved for each source, 2 incomparable records (another reference and a
reference-base mismatch), WW_S4 excluded by the paper preset (78% genome
coverage), and the paper-method PCoA on 7 pooled samples.

### Real data

- **Pooled:** upload Fontenele et al. (2021) Supplementary Table 2
  (`mmc2.xlsx`, from the article's supplementary materials) as the SNV table.
  It has SNV calls only, so PCoA falls back to the exploratory ordination and
  unreported positions stay *coverage unknown*. Genome coverage per sample is
  not in that table; supply it as pooled metadata (from the paper's Table 1)
  to apply the >90% rule. The real workbook was checked on November 1, 2026:
  7,973 records, 52 samples, 10 states, 24–793 SNVs per sample, and D614G in
  all 52 samples match the paper exactly; the full comparison (including small
  differences) is in `ARTICLE_METHOD_VALIDATION.md`, and
  `python validation/check_fontenele_2021.py path\to\mmc2.xlsx` repeats it.
  The whole table runs in about 8 seconds.
- **Individual + cross-source:** upload SARS-CoV-2 consensus sequences from a
  comparable place and period, the NC_045512.2 reference FASTA (one record),
  and the pooled table. The bundled H5N1 influenza Sprint 2 data and the
  SARS-CoV-2 wastewater data are different pathogens and cannot be compared.
- Data larger than the 10 MB web limit (for example a genome-wide A/C/G/T
  table for 52 samples) should be analyzed with `run_integrated_analysis.py`.

## Requirements

| Requirement | Used for |
|---|---|
| Python 3.11 and pip | Runs the application and installs packages |
| A web browser | Displays the local website |
| MAFFT on `PATH` | Individual-genome alignment and the MAFFT diagnostic |
| Optional GeoNames account and internet access | Missing-coordinate lookup and the GeoNames diagnostic |
| Optional Docker Desktop/Engine | Runs the supplied Linux container |

MAFFT and GeoNames are not needed for the processed-measurement or
processed-SNV pages. Those pages summarize the uploaded tables locally.
All pages share one Flask application, so install the complete
`requirements.txt` even when using only one page.

| Python dependency | Purpose |
|---|---|
| `pandas>=2.1,<3.0` | Tables, cleaning, grouping, and CSV output |
| `numpy>=1.26,<3.0` | Yue–Clayton distances, PCoA eigen-decomposition, regression (already a pandas dependency; now listed because the code imports it directly) |
| `openpyxl>=3.1,<4.0` | Reading XLSX workbooks |
| `Flask>=3.0,<4.0` | Web forms and rendered results |
| `biopython>=1.83,<2.0` | FASTA parsing |
| `requests>=2.31,<3.0` | GeoNames HTTPS requests |
| `gunicorn>=22.0,<24.0` | Web server in the Linux Docker image |
| `matplotlib>=3.8,<4.0` | Dashboard figures and Sprint 3 measurement charts |

Local Windows use starts with `python app.py`. Gunicorn is used in Linux
hosting. MAFFT is an external program and is not installed by pip.

## Install MAFFT for individual-genome analysis

MAFFT aligns comparable nucleotide sequences before direct position-by-position
comparison. The code finds it with `shutil.which("mafft")` and runs
`mafft --auto --nuc` separately for each pathogen and genomic region.
If Flask is running in that terminal, stop it with `Ctrl+C` before changing
the session settings, then restart it with `python app.py`.

### Windows

1. Download the all-in-one Windows package from the
   [official MAFFT Windows page](https://mafft.ddbj.nig.ac.jp/alignment/software/windows_without_cygwin.html).
2. Extract the entire package and keep its folders together. Find `mafft.bat`.
3. Add its containing folder to `PATH`. For the current Command Prompt session:

   ```bat
   set "PATH=C:\tools\mafft-win;%PATH%"
   where mafft
   mafft --version
   ```

Replace `C:\tools\mafft-win` with your actual folder. In PowerShell:

```powershell
$env:Path = "C:\tools\mafft-win;" + $env:Path
Get-Command mafft
mafft --version
```

For a permanent setting, add that folder to your Windows user `Path` through
the environment-variable settings, then open a new terminal. The application
supports `mafft.bat` through `cmd.exe /c`.

### macOS and Linux

Install MAFFT using the [official software instructions](https://mafft.ddbj.nig.ac.jp/alignment/software/)
or your operating system's package manager. On Debian/Ubuntu:

```sh
sudo apt-get update
sudo apt-get install mafft
command -v mafft
mafft --version
```

Start Flask from a terminal that can find MAFFT. If it cannot be found, the
genomic pipeline stops with an installation message. Uploaded sequences can
have different original lengths; MAFFT must return equal aligned lengths
within each comparable group.

## Configure optional GeoNames lookup

1. Create your own [GeoNames account](https://www.geonames.org/login).
2. Confirm the registration email.
3. Sign in and enable **Free Web Services** on
   [Manage Account](https://www.geonames.org/manageaccount).
4. Set `GEONAMES_USERNAME` in the terminal that will start Flask.

Windows Command Prompt:

```bat
set "GEONAMES_USERNAME=your_username"
python app.py
```

PowerShell:

```powershell
$env:GEONAMES_USERNAME = "your_username"
.\.venv\Scripts\python.exe app.py
```

macOS/Linux, with the virtual environment active:

```sh
export GEONAMES_USERNAME="your_username"
python app.py
```

Replace the placeholder with your account username. The application needs no
GeoNames password. Keep the value in the environment rather than the source.
A `.env` file alone is not automatically loaded by `python app.py`.
Restart an already-running Flask process after changing its environment.

The application uses `https://secure.geonames.org/searchJSON` with SSL
verification enabled. If configuration or lookup is unavailable, it continues
with missing coordinates and the existing specific-location-name fallback.

GeoNames is used only when coordinates are incomplete. Both provided coordinates
are retained when present. Results are approximate place coordinates and may
be ambiguous; repeated location names are cached within an analysis run.
See the [GeoNames web-service documentation](https://www.geonames.org/export/web-services.html).

State-only labels such as `USA: ND` and administrative-area estimates cannot
establish the kilometer cutoff or the named-location fallback. Their precision
is recorded. The processed-data pages group supplied location labels without
calling GeoNames.

## Check the local or hosted system status

On the main page, click **Check MAFFT & GeoNames**. No uploads are required.

| Status | Meaning |
|---|---|
| Operational | The real test alignment or live location lookup succeeded |
| Unavailable | An executable could not be found/started, or a connection failed |
| Not configured | `GEONAMES_USERNAME` is absent |
| Check failed | Execution, timeout, API, or returned-data checks failed |

The checks are independent. MAFFT aligns two short synthetic sequences and
verifies both IDs and equal output lengths. A version is shown when it can
be read reliably. GeoNames makes one fixed Tempe, Arizona lookup and checks
the returned latitude and longitude. The username is never displayed.

Results are cached for five minutes per application process. Checks run on
button presses rather than page loads. `/health` remains a separate lightweight
endpoint. Automated tests mock external calls; clicking this button verifies
the dependencies on the actual computer or cloud server.

## Sprint 2: FASTA and metadata input

### FASTA sequences

Accepted extensions: `.fasta`, `.fa`, and `.fna`. Each record contains one
nucleotide sequence, with its accession/sample identifier first in the header:

```text
>ACC001.1 segment 4
ACGTACGTACGT
>ACC002.1 segment 4
ACGTACGTACGA
```

Biopython reads the ID and full description. Sequences are uppercased and
whitespace is removed. Allowed characters are `A C G T N -` and, since the
Sprint 2 + Sprint 3 integration update, the IUPAC ambiguity codes `R Y S W K M B D H V`, which real
consensus genomes often contain. Every stage treats an ambiguity code like
`N`: it is never counted as a nucleotide difference or a variant, so p-distance
results for `A C G T N -` input are unchanged. Empty sequences, other
characters, and duplicate normalized accessions are reported.

### Metadata CSV or XLSX

Use a header row and one metadata row per accession. The loader accepts NCBI
Virus nucleotide FASTA plus a Results Table CSV when the required fields below
are included in the export.

| Standard column | Accepted aliases | Requirement |
|---|---|---|
| `Accession` | `Accession` | Required |
| `Pathogen` | `Organism_Name` | Required (or stated once on the review page) |
| `Collection_Date` | `Collection_Date` | Optional context; values must be valid dates when given |
| `Location_Name` | `Geo_Location` | Optional context |
| `Host_Source` | `Host` | Optional context |
| `Gene_Segment` | `Segment`, `Gene`, `Gene_Name`, `Region` | Optional if the region can be established another way |
| `Latitude` | `Latitude` | Optional; valid range -90 to 90 |
| `Longitude` | `Longitude` | Optional; valid range -180 to 180 |

Prefer ISO dates such as `2026-01-08`. Metadata headers tolerate inconsistent
capitalization and spacing. The guided upload recognizes more header spellings
(see [Required versus optional information](#required-versus-optional-information));
a missing optional column is a visible note, and year-only dates are treated as
unknown for day-level comparisons. A segment value `4` becomes `Segment 4`.
A text gene label such as `NA` is preserved as text.

Accession matching trims surrounding spaces and ignores one final numeric
version suffix: metadata `QB061804` matches FASTA `QB061804.1`. Matching is
case-sensitive. Full original values are retained where useful. Duplicate
base accessions and unmatched records on either side are errors.

### Resolve the genomic region

The website's **Automatic detection** uses explicit metadata first, then
explicit FASTA descriptions. **Complete genomes** additionally requires a
confirmation checkbox that all sequences are complete genomes of one pathogen.
Conflicting annotations stop the run.

The program does not determine a gene from length or use MAFFT to identify it.
Region-label sources are reported, and only equal pathogen/region groups are
aligned and compared.

For backward compatibility, the command line also accepts
`--default-gene-segment "Region A"` or `--confirm-same-region`.
The latter assigns `Unspecified_Region` only with explicit confirmation and
one pathogen. These options do not override contradictory annotations.

### Analysis and results

The six stages are:

1. Validate and match FASTA records and metadata.
2. Fill missing coordinates through GeoNames when available.
3. Align each pathogen/region group with MAFFT.
4. Compare aligned sequences.
5. Add collection-date and geographic relationships.
6. Apply thresholds and find connected potential clusters.

Positions containing `N` or `-` in either sequence are excluded. The formulas
remain:

```text
P_Distance = Differences / Comparable_Positions
Percent_Identity = (Comparable_Positions - Differences) / Comparable_Positions * 100
```

If there are no comparable positions, p-distance and percent identity remain
missing rather than dividing by zero.

The webpage takes three non-negative thresholds: maximum p-distance, days,
and kilometers. Strict related pairs need genetic and date criteria plus usable
geographic evidence. A specific matching location name can supply the existing
fallback when kilometer distance is unavailable. State-level locations cannot.

The results keep strict potential clusters separate from genomic-temporal
candidate pairs with unresolved geography. Candidates pass genetic and date
criteria but lack usable geographic evidence; pairs with known distances over
the cutoff are not categorized as unresolved.

The summary distinguishes samples from relationships. The current NCBI
regression fixture has 10 comparable samples, 45 unique comparisons, and
22 unresolved-geography candidate pairs involving 9 samples under its
demonstration settings. Those relationships do not create additional samples.
P-distances are displayed to six decimals while calculations and CSV values
retain their precision. Candidate lists start collapsed under **Show candidate
pairs**.

Thresholds are prototype settings rather than universal biological definitions.
Connected clusters can include A and C through B even without a direct A-C link.

## Run the genomic pipeline from the command line

Activate your environment and make MAFFT available first:

```console
python run_pipeline.py "path/to/sequences.fasta" "path/to/metadata.csv" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50
```

Replace both paths with your inputs; `metadata.xlsx` works too. Relative paths
start from the current folder. If the data are confirmed complete genomes:

```console
python run_pipeline.py "path/to/sequences.fasta" "path/to/metadata.xlsx" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50 --complete-genomes
```

The numbers above are demonstration settings. All five output CSVs are written
beside the metadata, so use a writable working folder. To try the bundled
synthetic inputs without creating outputs among the source examples, use this
Command Prompt example:

```bat
mkdir outputs\local_demo
copy examples\synthetic_sequences.fasta outputs\local_demo\synthetic_sequences.fasta
copy examples\synthetic_metadata.csv outputs\local_demo\synthetic_metadata.csv
python run_pipeline.py "outputs/local_demo/synthetic_sequences.fasta" "outputs/local_demo/synthetic_metadata.csv" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50
```

If `outputs\local_demo` already exists, reuse it and omit `mkdir`. That directory
is ignored by Git.

For `metadata.csv`, the outputs are:

| Filename | Contents |
|---|---|
| `metadata_pairwise_comparisons.csv` | Valid within-group genetic comparisons |
| `metadata_pairwise_comparisons_temporal_geographic.csv` | Genetic pairs plus days, distance, locations, coordinate sources, and precision |
| `metadata_pairwise_comparisons_temporal_geographic_related_pairs.csv` | Strict related pairs |
| `metadata_pairwise_comparisons_temporal_geographic_cluster_membership.csv` | Cluster IDs and sample membership |
| `metadata_region_assignments.csv` | Resolved genomic region and its source per sample |

The main Flask page uses temporary uploads and outputs, which are removed after
the request. Use the command-line workflow when you want to keep these genomic
CSV files on disk. The separate candidate list is a webpage summary, not an
addition to the strict related-pairs CSV.

## Wastewater Data Analysis: processed measurement tables

Open `http://127.0.0.1:5000/wastewater-data`. The earlier `/environmental` route
uses the same analysis with its earlier page title.

Use one row per sample and feature. CSV files should be UTF-8; for XLSX use a
simple header row on the first worksheet.

| Column | Meaning |
|---|---|
| `Sample_ID` | Required; collected sample identifier |
| `Feature_ID` | Required; genomic feature identifier |
| `Collection_Date` | Optional ISO date or Excel date cell |
| `Location_Name` | Optional sampling location |
| `Frequency` | Optional fraction from 0 to 1 |
| `Measurement_Value` | Optional numeric measurement |
| `Coverage` | Optional non-negative supplied coverage value |
| `Read_Depth` | Optional non-negative supplied read depth |
| `Quality` | Optional supplied quality label or information |
| `Unit` | Optional measurement unit |

Aliases include `Sample`, `Sample ID`, `Feature`, `Feature Identifier`, `Date`,
`Site`, `Sampling Location`, `Feature Frequency`, `Alt_Frequency`, `ALT_FREQ`,
`Measurement`, `Value`, `Coverage Depth`, `Depth`, and `Quality Flag`.
Header spacing and capitalization are normalized.

Example:

```csv
Sample_ID,Feature_ID,Collection_Date,Location_Name,Frequency,Coverage
W1,feature_alpha,2026-01-01,North site,0.10,100
W2,feature_alpha,2026-01-08,North site,0.15,
W3,feature_beta,2026-01-08,South site,,
```

Only `Sample_ID` and `Feature_ID` are required. A feature-only table still
produces counts. Feature IDs should distinguish references, organisms, or assays
when those affect comparability.

The page displays:

- Unique samples and listed features, with feature counts per sample.
- Available frequencies and explicit frequency/coverage availability messages.
- Sample/feature counts by date and location.
- Means and ranges for supplied measurements, with changes from the previous
  date within a feature/location group.
- Missing-data and supplied-quality summaries.
- Basic charts and the first 100 organized rows.
- A download of every validated, organized row.

Frequency is summarized when it has usable values; otherwise available
`Measurement_Value` values are summarized. Coverage and read depth stay
separate. Missing numerical values remain blank and are excluded from means;
a supplied zero remains zero. Explicit missing dated measurements produce gaps,
and the first date has no previous-date change. Frequency or coverage is never
inferred from feature counts.

Feature counts count entries listed in the table. Unlisted features do not
establish absence. Dates and locations can be missing, while invalid supplied
values are reported. Repeated sample/feature pairs, conflicting sample metadata,
and conflicting supplied units for a feature are rejected.

Charts show the first 20 samples and numerical trends for up to three features
and six locations per feature. Tables retain the other features; full organized
rows remain in the CSV. These are descriptive comparisons, not statistical
significance tests. No external alignment, geocoding, or raw-read software is
run by this page.

## Existing processed-SNV page

Open `http://127.0.0.1:5000/wastewater` for **Wastewater Surveillance**.

Its standardized CSV/XLSX input requires these columns:

```text
Sample_ID,Reference_Genome,Position,Ref,Alt,Alt_Frequency,Read_Depth
```

`Collection_Date` and `Location_Name` are optional. `Position` is a positive
integer; `Ref` and `Alt` are distinct single `A/C/G/T` bases.
`Alt_Frequency` is required and uses fractions from 0 to 1. The `Read_Depth`
column is required, but individual depth values can be blank; supplied depths
must be positive integers.

An optional metadata CSV/XLSX can provide `Sample_ID`, `Collection_Date`, and
`Location_Name`. Metadata IDs must be unique and cover every SNV sample.
Conflicting dates or locations are errors. Repeated SNV keys are reported and
retained; distinct-sample counts avoid counting a sample twice.

This page also recognizes the existing supported format of Fontenele et al.
(2021), Supplementary Table 2 (`mmc2.xlsx`), available from the
[article's supplementary materials](https://pmc.ncbi.nlm.nih.gov/articles/PMC8464352/).
For that recognized workbook only, the loader maps `SAMPLE`, `COLLECTION_DATE`,
`POSITION`, `REF`, `ALT`, `ALT_FREQ`, and `TOTAL_DP`, combines `STATE` and
`LOCATION`, and supplies its configured reference label `NC_045512.2`.
Other files must explicitly provide `Reference_Genome`. The published workbook
is not bundled in GitHub; the generic measurement page uses a different schema.

Results include counts, date/location information, frequency/depth summaries,
the most frequently reported SNVs, and a CSV download. A missing call is unknown
when coverage is unavailable. The page does not run the individual-genome
pipeline or infer lineages.

## Uploads, temporary storage, and downloads

| Limit | Current behavior |
|---|---|
| Whole web request | 10 MB, including multipart overhead |
| Measurement upload (`/wastewater-data` or `/environmental`) | Under 3 MB; 1 to 10,000 rows |
| Processed-SNV upload (`/wastewater`) | Under 3 MB |
| Optional SNV metadata upload | Under 1 MB |
| Organized processed-data CSV | At most 4 MB |
| Detailed processed-data preview | First 100 rows |
| Processed-data download cache | Up to 10 minutes; at most four recent results per mode and application process |
| Guided upload (review step) | Up to 20 files per request; gzip files expand to at most 60 MB; kept in a temporary folder for at most 30 minutes, at most six upload sets |

Uploads are handled in temporary directories. The processed-data download
cache holds CSV bytes in memory, so a restart clears links and a later run may
evict an older link before ten minutes. Download promptly and rerun the upload
if a link expires. The application has no database.

## Tests and local verification

With the environment active:

```console
python -m unittest discover -s tests -p "test_environmental_analysis.py" -v
```

This portable measurement-page test file covers CSV/XLSX loading, supplied zeros
versus missing values, counts, date/location summaries, chart generation,
uploads, downloads, and the bundled synthetic example.

For the full suite in the development folder:

```console
python -m unittest discover -s tests -v
```

The latest full run (November 1, 2026, Python 3.11.9, with real MAFFT v7 on
`PATH`) ran **204 tests: 202 passed and 2 were skipped** (the two optional
real-data tests below). The suite is self-contained:

- The invented DEMO/NCBI-style example files live in `tests/fixtures/`.
- The state-only Flask upload check generates its own NCBI-style records.
- `tests/test_optional_real_ncbi_data.py` keeps the original real-NCBI
  (H5N1 segment 4) check as an **optional integration test**. It is skipped
  with a stated reason unless `sequences.fasta` and `sequences.csv` (untracked)
  are placed in the project root; it passed when those files were present.
- `RealMafftReferenceAlignmentTest` runs only when real MAFFT is on `PATH`
  (it ran and passed with `C:\Tools\mafft-win\mafft.bat`).
- `tests/test_optional_fontenele_table.py` checks the real Supplementary
  Table 2 against the paper's counts when `FONTENELE_MMC2` points to
  `mmc2.xlsx` (not stored in Git); it passed when the workbook was supplied.
- `tests/test_input_detection.py` and `tests/test_guided_workflow.py` cover
  gzip, VCF, FASTQ refusal, header detection, ambiguous columns, optional
  metadata, missing dates/geography/reference, compatible and incompatible
  cross-source comparisons, the review/run pages, and the new figures.

**Mocked versus real external checks.** MAFFT is mocked in the automated
tests: older tests use a pass-through, and the integrated tests use
`tests/fake_mafft.py`, a small Needleman–Wunsch aligner, so partial
sequences and indels can be exercised deterministically. GeoNames network
calls are always mocked. Before deploying, confirm the real tools with the
System Status button and by clicking **Run Integrated Example** on a machine
or container where MAFFT is installed.

Real-tool checks run outside the automated mocks on November 1, 2026 (Windows,
MAFFT from `C:\Tools\mafft-win`): the synthetic integrated example gave the
same cross-source categories as the mocked tests (2 shared, 1 individual-only,
1 pooled-only, 1 + 1 unresolved, 2 incomparable), and the real NCBI H5N1
segment 4 data through the guided workflow gave 10 sequences, 45 comparisons,
22 genomic-temporal candidate pairs among 9 samples, and 0 strict clusters
(state-level locations). GeoNames was **not** checked live (no
`GEONAMES_USERNAME` in this environment).

## Troubleshooting locally

| Problem | What to check |
|---|---|
| `python` or `py` is not recognized | Install Python 3.11 and reopen the terminal. Use the available launcher to create the environment. |
| `ModuleNotFoundError` | Use the intended environment and run `python -m pip install -r requirements.txt` with that same Python. |
| PowerShell blocks activation | Run `.\.venv\Scripts\python.exe` directly; no execution-policy change is needed. |
| A copied environment does not work | Create a fresh virtual environment on the current computer and reinstall requirements. |
| Browser cannot connect | Keep the terminal running and use `http://127.0.0.1:5000`, not a filesystem path. |
| Port 5000 is busy | Stop another local copy, or run `python -m flask --app app run --host 127.0.0.1 --port 5001` and visit port 5001. |
| New page is missing or code changes are not visible | Stop Flask with `Ctrl+C` and restart it from the current project folder. Debug auto-reload is disabled. |
| MAFFT is unavailable | Check `where mafft` (Windows) or `command -v mafft` (macOS/Linux) in the same terminal used to start Flask. |
| GeoNames is not configured | Set `GEONAMES_USERNAME` and restart Flask. This is optional for the analysis. |
| GeoNames is unavailable or rejected | Check the account's web-service activation and network access; use System Status for the readable result. Keep SSL verification enabled. |
| Unknown genomic region | Supply explicit metadata/header annotation, or explicitly confirm genuinely complete genomes through the form. |
| Accessions do not match | Compare FASTA first-header IDs with metadata accessions; only a final numeric version suffix is ignored. |
| Missing measurements | Leave unavailable cells blank. Supply numeric frequencies as fractions, not percentage text. |
| Zero strict clusters but candidate pairs exist | Genetic/date relationships exist, while usable geographic evidence is insufficient for the strict result. |
| Download returns 404 | Its temporary cache entry expired, was evicted, or the server restarted. Rerun the analysis and download again. |
| Full-suite `FileNotFoundError` after cloning | Check the local regression-fixture list above. |

The local Flask development server uses `debug=False`. Each computer's
`127.0.0.1` points to that computer; running locally does not automatically
create a public website.

## Optional local Docker workflow

Docker provides a separate way to run the project with Python and Linux MAFFT
inside the container. It is not required for the Python setup above.
Start Docker Desktop/Engine, then from the project folder run:

```console
docker build -t bmi461-sprint2 .
docker run --rm bmi461-sprint2 mafft --version
docker run --rm -p 10000:10000 bmi461-sprint2
```

Open `http://127.0.0.1:10000` or
`http://127.0.0.1:10000/wastewater-data`. Stop the foreground container with
`Ctrl+C`.

To use GeoNames in Docker, first set `GEONAMES_USERNAME` in your host terminal
as shown earlier, then pass the variable by name:

```console
docker run --rm -p 10000:10000 --env GEONAMES_USERNAME bmi461-sprint2
```

The Dockerfile installs MAFFT through `apt-get`, installs `requirements.txt`,
copies runtime modules/templates/examples, and starts one Gunicorn worker with
a 180-second timeout. It listens on `0.0.0.0` using `PORT`, defaulting to
`10000`. This README does not claim a new Docker build or hosted verification
was performed during its documentation update.

## GitHub and Render hosting

Local use does not require GitHub or Render. To publish the existing code:

1. Keep the root Python modules, `wastewater/`, `templates/`, `examples/`,
   `requirements.txt`, `Dockerfile`, `.dockerignore`, and `.gitignore` in GitHub.
   Include tests and this README for maintenance.
2. Review `git status` before committing. Keep credentials, uploaded sample
   data, caches, and generated output files out of commits. The ignore files
   do not remove anything already tracked.
3. In Render, create or update a **Web Service** connected to the repository
   and branch. Select **Docker**, use the repository's Dockerfile, and choose
   **Free** if suitable for the class demonstration.
4. Use the Dockerfile's default start command. Set the health-check path to
   `/health`. If desired, set `GEONAMES_USERNAME` in Render's Environment page.
5. Deploy the current commit and open the service's assigned
   `https://...onrender.com` URL. Use the same route suffixes listed above.
6. Verify **Run Synthetic Example** on `/wastewater-data` and
   **Check MAFFT & GeoNames** on the main page. Then test a permitted NCBI
   FASTA/results-table upload with the chosen demonstration thresholds.

For an existing Git-backed service, push updates to its linked branch and
check Render's deploy result. If automatic deployment is disabled, trigger
a deployment from the dashboard. Visitors to the hosted service need only
a browser; MAFFT runs on the cloud server.

Official setup references:
[Render Docker](https://render.com/docs/docker),
[Web Services](https://render.com/docs/web-services),
[Environment Variables](https://render.com/docs/configure-environment-variables),
and [Free hosting limitations](https://render.com/docs/free).
Free services can sleep when idle and have an ephemeral filesystem. The app's
temporary uploads/downloads are designed around that behavior. Use permitted
classroom or appropriately de-identified data on this unauthenticated prototype.

## Main files

| File or folder | Role |
|---|---|
| `app.py` | Flask forms, results, diagnostics, and temporary downloads |
| `validate_inputs.py` | FASTA/metadata validation and accession matching |
| `genomic_regions.py` | Explicit genomic-region labels and consistency checks |
| `geocode_locations.py` | Optional GeoNames lookup and its diagnostic |
| `align_sequences.py` | MAFFT alignment within comparable groups |
| `compare_genetic_distances.py` | Pairwise nucleotide differences and identity |
| `add_temporal_geographic.py` | Date differences, Haversine distances, and location checks |
| `detect_potential_clusters.py` | Strict related pairs and connected components |
| `run_pipeline.py` | Six-stage individual-genome runner |
| `system_diagnostics.py` | Small MAFFT test alignment |
| `environmental_analysis.py` | Processed measurement loading, summaries, and charts (Sprint 3 legacy page) |
| `wastewater/processed_snv.py` | Processed-SNV importer: paper Table 2, iVar TSV, SNV VCF, standardized tables |
| `genomics/` | Shared reference-based model, QC, analyses, matching, figures, orchestration (Sprint 2 + 3 integration) |
| `genomics/file_formats.py`, `genomics/column_detection.py`, `genomics/input_inspection.py` | Guided upload: formats and VCF, flexible columns, file roles, dataset review |
| `integrated_view.py` | Arranges integrated results into the eight dashboard sections |
| `literature/`, `literature_view.py` | Optional Sprint 4 AI literature and metadata enrichment (see `SPRINT4_AI_NLP_INTEGRATION_REPORT.md`) |
| `validation/check_fontenele_2021.py` | Compares the import of the paper's real Supplementary Table 2 with its reported counts |
| `ARTICLE_METHOD_VALIDATION.md` | Method-to-implementation table for Fontenele et al. (2021) |
| `run_integrated_analysis.py` | Command-line integrated analysis (larger processed datasets) |
| `offline_workflow/` | Optional BWA/iVar command template and tested output converters (not used by the website) |
| `templates/` | Flask HTML templates (`index.html` is the unified dashboard) |
| `examples/` | Bundled synthetic inputs; `examples/integrated/` is the paired demonstration |
| `tests/` | Unit and integration tests; `tests/fixtures/` holds small invented regression files |
| `INTEGRATION_IMPLEMENTATION_LOG.md` | Integration research facts, progress, and limitations |
| `requirements.txt` | Python package requirements |
| `Dockerfile` | Linux runtime, MAFFT, and Gunicorn |
| `SPRINT_3_CHANGELOG.md` | Record of changes after the Sprint 2 baseline |

The local historical `output/pdf/` guides and `PYTHON_LINE_BY_LINE_GUIDE.md`
describe earlier code. Some `output/` and `outputs/` files are not included in
GitHub. Use this README and the current source for the current setup.

The website performs no raw FASTQ processing, read mapping, variant calling,
lineage deconvolution, or advanced statistical modeling; an optional offline
template for read processing is in `offline_workflow/`. Genomic similarity,
shared variants, time/location proximity, and processed-table patterns need
scientific context before interpretation.

## AI tools acknowledgment

ChatGPT and Claude were used as supporting references during development for
brainstorming, code explanation, debugging assistance, and documentation review.
Their suggestions were reviewed by the project team. They are not scientific
evidence or substitutes for peer-reviewed genomic epidemiology guidance.

- OpenAI. *ChatGPT*. <https://chatgpt.com/>. Accessed September 28, 2026.
- Anthropic. *Claude*. <https://claude.ai/>. Accessed September 28, 2026.
