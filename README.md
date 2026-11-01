# BMI 461 Genomic Surveillance and Wastewater Data Analysis

A classroom Flask application with an individual-genome surveillance pipeline
and pages for summarizing already processed environmental measurements.

The current application can run on a personal computer or in the supplied
Docker container. This README describes the current code as of October 19,
2026. A potential outbreak cluster is a screening result; it does not confirm
transmission or an outbreak.

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
| Individual-genome cluster analysis | `http://127.0.0.1:5000/` | FASTA + CSV/XLSX metadata; requires MAFFT for analysis |
| Wastewater Data Analysis | `http://127.0.0.1:5000/wastewater-data` | One processed feature-measurement CSV/XLSX |
| Wastewater Surveillance | `http://127.0.0.1:5000/wastewater` | Processed SNV CSV/XLSX, with optional sample metadata |
| Earlier environmental page | `http://127.0.0.1:5000/environmental` | Same measurement analysis as Wastewater Data Analysis |
| Health check | `http://127.0.0.1:5000/health` | Returns `OK`; performs no external dependency tests |

For a first local demonstration without installing MAFFT:

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
| `openpyxl>=3.1,<4.0` | Reading XLSX workbooks |
| `Flask>=3.0,<4.0` | Web forms and rendered results |
| `biopython>=1.83,<2.0` | FASTA parsing |
| `requests>=2.31,<3.0` | GeoNames HTTPS requests |
| `gunicorn>=22.0,<24.0` | Web server in the Linux Docker image |
| `matplotlib>=3.8,<4.0` | Basic measurement charts |

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
whitespace is removed. Allowed characters are `A C G T N -`.
Empty sequences, invalid characters, and duplicate normalized accessions
are reported.

### Metadata CSV or XLSX

Use a header row and one metadata row per accession. The loader accepts NCBI
Virus nucleotide FASTA plus a Results Table CSV when the required fields below
are included in the export.

| Standard column | Accepted aliases | Requirement |
|---|---|---|
| `Accession` | `Accession` | Required |
| `Pathogen` | `Organism_Name` | Required |
| `Collection_Date` | `Collection_Date` | Required, with valid values |
| `Location_Name` | `Geo_Location` | Required, with values |
| `Host_Source` | `Host` | Required, with values |
| `Gene_Segment` | `Segment`, `Gene`, `Gene_Name`, `Region` | Optional if the region can be established another way |
| `Latitude` | `Latitude` | Optional; valid range -90 to 90 |
| `Longitude` | `Longitude` | Optional; valid range -180 to 180 |

Prefer ISO dates such as `2026-01-08`. Metadata headers tolerate inconsistent
capitalization and spacing. A segment value `4` becomes `Segment 4`.
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

The latest full run passed **110 tests on Python 3.11.9** on October 19, 2026.
External MAFFT and GeoNames behavior is mocked in automated tests; also use
the real System Status button and a real example analysis.

Two full-suite regression methods use local data that are excluded by
`.gitignore` and may be absent from a fresh GitHub download:

```text
sequences.fasta
sequences.csv
outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/example_sequences.fasta
outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/example_metadata.csv
outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/ncbi_virus_style_sequences.fasta
outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/ncbi_virus_style_metadata.csv
```

These are used by
`test_state_only_ncbi_upload_keeps_strict_and_candidate_results_separate` in
`test_app.py` and `test_both_existing_examples_still_produce_clusters` in
`test_pipeline.py`. Obtain the original matching fixtures from the project
maintainer to reproduce those exact tests. Their absence does not prevent
running the application or the portable measurement tests.

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
| `environmental_analysis.py` | Processed measurement loading, summaries, and charts |
| `wastewater/processed_snv.py` | Existing processed-SNV loader and summaries |
| `templates/` | Flask HTML templates |
| `examples/` | Three bundled synthetic input files |
| `tests/` | Unit and integration tests |
| `requirements.txt` | Python package requirements |
| `Dockerfile` | Linux runtime, MAFFT, and Gunicorn |
| `SPRINT_3_CHANGELOG.md` | Record of changes after the Sprint 2 baseline |

The local historical `output/pdf/` guides and `PYTHON_LINE_BY_LINE_GUIDE.md`
describe earlier code. Some `output/` and `outputs/` files are not included in
GitHub. Use this README and the current source for the current setup.

The prototype performs no raw FASTQ processing, read mapping, variant calling,
lineage deconvolution, or advanced statistical modeling. Genomic similarity,
time/location proximity, and processed-table patterns need scientific context
before interpretation.

## AI tools acknowledgment

ChatGPT and Claude were used as supporting references during development for
brainstorming, code explanation, debugging assistance, and documentation review.
Their suggestions were reviewed by the project team. They are not scientific
evidence or substitutes for peer-reviewed genomic epidemiology guidance.

- OpenAI. *ChatGPT*. <https://chatgpt.com/>. Accessed September 28, 2026.
- Anthropic. *Claude*. <https://claude.ai/>. Accessed September 28, 2026.
