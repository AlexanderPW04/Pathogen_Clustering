# Disease-Agnostic Genomic Epidemiology Prototype

This BMI 461 Sprint project is a Python prototype for screening pathogen
samples for **potential outbreak clusters**. It combines nucleotide-sequence
similarity with collection dates and locations, then presents cluster membership
through either command-line CSV outputs or a small Flask webpage. The webpage
can run locally or on a Render Docker web service.

> **Important:** A potential cluster is a screening signal, not a confirmed
> outbreak or proof of direct transmission. Thresholds must be selected for the
> pathogen, genomic region, sampling design, and public-health question.

## Current input workflow

The program uses two source files, similar to a sequence-plus-results-table
download from NCBI Virus:

```text
FASTA nucleotide sequences        XLSX or CSV metadata
             \                         /
              \                       /
        accession matching, region resolution, and validation
                              |
                              v
          matched internal DataFrame used by all analysis stages
                              |
                              v
              GeoNames fills missing coordinates when available
                              |
                              v
             MAFFT aligns each pathogen/gene group
                              |
                              v
       genetic comparison -> time/location -> potential clusters
                              |
                              v
                    CSV files or Flask results
```

The FASTA and metadata files are matched through a normalized `Accession`.
A final numeric FASTA version such as `.1` is ignored for matching, so metadata
`QB061804` matches FASTA `QB061804.1`. The full metadata accession is preserved
in `Accession`, while `Sample_ID` uses the full matched FASTA ID. The existing
genomic, temporal/geographic, and clustering calculations keep their established
internal column names and scientific rules.

## Project structure

```text
app.py                          Local Flask interface
Dockerfile                      Python 3.11 + Linux MAFFT + Gunicorn hosting
.dockerignore / .gitignore      Keep local data and caches out of deployments
run_pipeline.py                 Six-step command-line runner
validate_inputs.py              FASTA/metadata loading, matching, and validation
genomic_regions.py              Simple explicit-region parsing and checks
geocode_locations.py            Optional GeoNames coordinate lookup
align_sequences.py              MAFFT alignment within comparable groups
system_diagnostics.py           On-demand synthetic MAFFT check
compare_genetic_distances.py    Pairwise genomic comparison
add_temporal_geographic.py      Date and Haversine-distance enrichment
detect_potential_clusters.py    Threshold filtering and DFS clustering
requirements.txt                Third-party Python dependencies
PYTHON_LINE_BY_LINE_GUIDE.md     Historical pre-MAFFT line-by-line explanation
SPRINT_3_CHANGELOG.md           Changes made after the Sprint 2 baseline
templates/index.html            Flask/Jinja webpage template
examples/                       Bundled, explicitly synthetic demo inputs
tests/                          Unit and integration tests
outputs/                        Demonstration data and generated results
```

The earlier conceptual walkthrough is available as
`output/pdf/BMI461_Sprint2_Code_Review_Guide.pdf`, with editable Markdown beside
it. That PDF and `PYTHON_LINE_BY_LINE_GUIDE.md` describe the pre-MAFFT/GeoNames
baseline; their line numbers and pipeline order are not current. This README
describes the updated workflow.

### Beginner-friendly code style

The main Python files are intentionally written for an undergraduate
presentation. They use ordinary functions, named intermediate variables, basic
`for` loops, `if` statements, lists, dictionaries, and sets. The project avoids
classes, workflow frameworks, graph libraries, regular expressions, and
JavaScript scientific calculations. Comments separate the major steps, and the
scientific stages remain in individual files so each calculation can be
explained and tested separately.

## Requirements and installation

The project can run locally on Windows, macOS, or Linux, or in its Docker
container. It does not require a database. For local use, MAFFT must be
installed separately; Docker installs it automatically. GeoNames needs internet
access only when optional coordinate lookup is used.

Required software:

- Python 3.11 is the project target. The current update was also tested locally
  with Python 3.12 because Python 3.11 was unavailable in this test environment.
- `pip`, which is normally included with Python.
- MAFFT, the external nucleotide alignment program, available on `PATH`.
- A web browser if the optional Flask interface will be used.

The packages in `requirements.txt` are:

| Package requirement | Purpose |
|---|---|
| `pandas>=2.1,<3.0` | Table loading, conversion, grouping, and CSV output |
| `openpyxl>=3.1,<4.0` | XLSX metadata files |
| `biopython>=1.83,<2.0` | FASTA parsing |
| `Flask>=3.0,<4.0` | Local webpage |
| `requests>=2.31,<3.0` | Optional GeoNames web-service calls |
| `gunicorn>=22.0,<24.0` | Production web server inside Docker (not used by Windows `python app.py`) |

MAFFT is **not** a pip package in this project. Python starts the separately
installed program using `subprocess`.

### 1. Download or copy the project

Place the entire project folder on the computer. Keep `requirements.txt`, the
Python files, `templates/`, and `tests/` together. Open Command Prompt,
PowerShell, Terminal, or another shell in the project folder—the folder that
contains `run_pipeline.py`.

For example, the general form is:

```console
cd "path/to/sprint_2"
```

Replace `path/to/sprint_2` with the folder's actual location. Quotation marks are
recommended when a path contains spaces.

### 2. Confirm Python is available

On Windows, either of these commands may be available:

```console
py -3.11 --version
python --version
```

On macOS or Linux, use:

```console
python3 --version
```

The displayed version should begin with `Python 3.11`. If none of these commands
works, install Python 3.11 from <https://www.python.org/downloads/>. On Windows,
the installer can add Python to `PATH`; the `py` launcher can also select version
3.11 directly.

### 3. Create an isolated virtual environment

A virtual environment keeps this project's packages separate from other Python
projects.

Windows Command Prompt:

```console
py -3.11 -m venv .venv
.venv\Scripts\activate.bat
```

If the Windows `py` launcher is unavailable but `python --version` reports Python
3.11, use `python -m venv .venv` instead.

Windows PowerShell:

```console
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
```

macOS or Linux:

```console
python3 -m venv .venv
source .venv/bin/activate
```

If `python3` points to another Python version but `python3.11` is installed, use
`python3.11 -m venv .venv` instead.

After activation, the command prompt normally begins with `(.venv)`. If
activation is unavailable, the environment's Python executable can still be
called directly as `.venv\Scripts\python.exe` on Windows or
`.venv/bin/python` on macOS/Linux.

If PowerShell blocks its activation script, use the Windows Command Prompt
instructions or call `.venv\Scripts\python.exe` directly. Changing the computer's
execution policy is not required.

### 4. Install the required packages

With the virtual environment active, run:

```console
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Package installation requires an internet connection. Later analysis runs can
stay local if the metadata already has coordinates, or if GeoNames lookup is
skipped; MAFFT still runs locally on the computer.

### 5. Install MAFFT and put it on PATH

MAFFT aligns sequences that may have different lengths or gaps before genomic
distance is calculated. The code uses `shutil.which("mafft")` to check that the
executable is discoverable, then runs `mafft --auto --nuc` once per pathogen and
gene/segment group. It does not align unrelated groups together. All aligned
sequences in a group must have the same length before direct position-by-position
comparison. The existing formulas still ignore positions containing `N` or `-`.

Download MAFFT from the [official MAFFT site](https://mafft.ddbj.nig.ac.jp/alignment/software/).
For Windows, the [official Windows instructions](https://mafft.ddbj.nig.ac.jp/alignment/software/windows_without_cygwin.html)
describe the portable `mafft-win` download, which includes `mafft.bat`. Extract
it and add the folder containing `mafft.bat` to the system or user `PATH`. For
the current **Command Prompt window only**, an example is:

```bat
set "PATH=C:\path\to\mafft-win;%PATH%"
where mafft
```

Replace `C:\path\to\mafft-win` with the actual folder containing `mafft.bat`.
The `where mafft` command should print its location. On macOS or Linux, install
MAFFT using the official instructions or a suitable package manager and check
with `command -v mafft`. Restart the terminal after a permanent PATH change.
If MAFFT is missing, the pipeline stops before analysis with a clear error; it
does not silently compare unaligned sequences.

### 6. Optional: create a GeoNames account

If metadata has missing latitude or longitude, the pipeline can ask
[GeoNames](https://www.geonames.org/) for coordinates based on `Location_Name`.
Create a free GeoNames account, confirm the email address, and enable the free
web service in your account settings, as described in the
[GeoNames web-services guide](https://www.geonames.org/export/web-services.html).
Do not put the username in the source code. In the current Windows Command
Prompt window, set:

```bat
set "GEONAMES_USERNAME=your_username"
```

In PowerShell, use `$env:GEONAMES_USERNAME = "your_username"`; on macOS/Linux,
use `export GEONAMES_USERNAME=your_username`. Set the variable in the same
terminal from which you run the pipeline or Flask site. No password is needed
in the project.

If the username is absent, GeoNames lookup is skipped with a readable warning.
If its service or network is unavailable, unresolved coordinates remain blank.
In either case, the existing normalized `Same_Location` fallback still applies
for sufficiently specific place names. Lookups use the secure
`https://secure.geonames.org/searchJSON` endpoint with normal SSL certificate
verification.
The first GeoNames search result is used, so ambiguous place names should be
checked by a person. GeoNames gives approximate *place* coordinates, not an
exact specimen collection point. Provided coordinates are never replaced when
both latitude and longitude are present; if only one is present, both values
are replaced together when GeoNames resolves the place. Repeated location names
are looked up only once per run. When enabled, location names are sent to the
GeoNames service, so avoid including private specimen-level details in them.

A state-only name such as `USA: ND` is not an exact collection site. It is not
sent to GeoNames for a state centroid. GeoNames results identified as
administrative areas are also marked too coarse. For any pair containing such a
sample, `Km_Apart` and `Same_Location` are left blank, so the state name cannot
pass a 50 km cutoff or the name-based fallback. Existing coordinates are kept
and their source remains visible, but they are not used for that pair's
geographic comparison when the only location description is state-level. The
enriched CSV includes `Geographic_Precision_1` and
`Geographic_Precision_2`; the program also warns about state-only records.
Specific sites or appropriate sample coordinates are needed for a meaningful
local-distance assessment.

To confirm a **real** GeoNames request independently of MAFFT, set
`GEONAMES_USERNAME` as above, then run this from the project folder in Windows
Command Prompt:

```bat
python -c "import pandas as pd; from geocode_locations import fill_missing_coordinates; print(fill_missing_coordinates(pd.DataFrame({'Location_Name':['Tempe, Arizona'],'Latitude':[None],'Longitude':[None]})).to_string(index=False))"
```

The output should show `Coordinate_Source` as `GeoNames`, approximate latitude
and longitude, and `Geographic_Precision` as a GeoNames place estimate. The
automated tests mock the web request and confirm the secure URL; they do not
contact the live service.

### 7. Confirm the installation

Run the automated tests from the project folder:

```console
python -m unittest discover -s tests -v
```

A successful setup ends with `OK`. Tests
mock the external MAFFT program and GeoNames service, so also verify a real
MAFFT installation with `where mafft` or `command -v mafft` before a real run.

All remaining commands in this README use `python`, assuming the virtual
environment is active. On a system where it is not active, replace `python` with
the appropriate environment executable described above.

## FASTA format

Accepted extensions are `.fasta`, `.fa`, and `.fna`. Each FASTA record must
contain one nucleotide sequence. Biopython uses the first identifier in the
header as the record ID. For example, the ID below is `ACC001`:

```text
>ACC001 optional description
ACGTACGTNN--
```

Sequences are converted to uppercase. Whitespace is removed with ordinary
string splitting and joining; the implementation does not use regular
expressions. The accepted sequence characters are `A`, `C`, `G`, `T`, `N`, and
`-`. The program reports empty sequences, invalid characters, and duplicate
FASTA accessions. Biopython also reads each record's full description so an
explicit genomic-region label can be used when metadata does not provide one.

## Metadata format

Metadata can be `.xlsx` or `.csv`. The standardized internal columns and accepted
NCBI Virus aliases are:

| Standard column | Accepted NCBI alias | Required? | Purpose |
|---|---|---|---|
| `Accession` | `Accession` | Yes | Identifier used to match a FASTA record |
| `Pathogen` | `Organism_Name` | Yes | Prevents cross-pathogen comparison |
| `Gene_Segment` | `Segment`, `Gene`, `Gene_Name`, `Region` | No | Separates different genomic regions when supplied |
| `Collection_Date` | `Collection_Date` | Yes | Calculates absolute days apart |
| `Location_Name` | `Geo_Location` | Yes | Provides location comparison and fallback |
| `Host_Source` | `Host` | Yes | Descriptive epidemiological metadata |
| `Latitude` | `Latitude` | No | Haversine input when available |
| `Longitude` | `Longitude` | No | Haversine input when available |

Headers are matched without regard to capitalization or extra spaces; for
example, `  sEgMeNt  ` and `Gene Name` are recognized. An NCBI Virus Results
Table CSV using `Segment = 4` can be uploaded directly and becomes internal
`Gene_Segment = Segment 4`. Text values such as the gene name `NA` are kept as
text instead of being mistaken for blank spreadsheet cells. Latitude and longitude
may be absent or blank. If a GeoNames username is configured, the pipeline tries
to fill missing coordinates from the location name. Otherwise `Km_Apart` remains
blank and the existing normalized location-name fallback is used for specific
place names, but not state-level names.

Accession values are read as text and surrounding spaces are removed. For
matching only, one final dot followed entirely by digits is removed. Thus
`ABC123`, `ABC123.1`, and `ABC123.2` share the normalized base `ABC123`.
Original accessions are retained for reporting. Matching remains case-sensitive.

### How the genomic region is chosen

For each matched sample, the input stage uses the first available source:

1. An explicit metadata value in `Gene_Segment`, `Segment`, `Gene_Name`, `Gene`,
   or `Region`.
2. An explicit FASTA description, such as `segment 4`, `hemagglutinin (HA)`,
   or `complete genome`.
3. For the website's **Complete genomes** option only, `Complete genome` is
   assigned to still-unannotated records after the user separately confirms
   that all sequences are complete genomes of one pathogen.

The website defaults to **Automatic detection**. If metadata and FASTA are
both silent about a sample's region, automatic detection stops with a clear
error; it does not guess from sequence length. **Complete genomes** requires
the confirmation checkbox, one pathogen, and no conflicting gene or segment
annotation. Explicit partial-sequence or CDS descriptions are also treated as
contradictions. Neither
sequence length nor MAFFT is used to identify a gene. Conflicting metadata and
FASTA annotations cause validation to stop. If no region can be established,
the error asks for an annotation or an explicit complete-genome confirmation.
A generic metadata label is the user's assertion of comparability;
review it before interpreting results. MAFFT aligns only within the resolved
`Pathogen` and `Gene_Segment` groups.

The older command-line `--default-gene-segment` and `--confirm-same-region`
options remain for backwards compatibility, but are not controls on the
website. `--confirm-same-region` assigns `Unspecified_Region`, not a complete
genome. Do not use either option without independently checking sequence
comparability.

The validator shows region-label sources, and the pipeline writes a
`*_region_assignments.csv` file listing every sample's label and source:
`Metadata`, `FASTA header`, `User override`, or
`User-confirmed unspecified region`, or `User-confirmed complete genome`. The Flask results also show the source
for each clustered sample.

Validation clearly reports:

- missing required columns or values
- duplicate metadata accessions
- duplicate FASTA accessions
- metadata accessions without a FASTA sequence
- FASTA accessions without metadata
- invalid sequence characters
- invalid dates or coordinates
- unknown or conflicting genomic-region annotations

Rows are not silently deleted. The complete pipeline stops before genetic
comparison when validation or matching fails.

## Run the complete command-line pipeline

```console
python run_pipeline.py "path/to/sequences.fasta" "path/to/metadata.xlsx" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50
```

CSV metadata works the same way:

```console
python run_pipeline.py "path/to/sequences.fna" "path/to/metadata.csv" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50
```

Replace both example paths with the actual input-file locations. Relative paths
are interpreted from the current project folder. Forward slashes work in Python
paths on Windows, macOS, and Linux, which is why they are used in these examples.
The five command-line output CSVs are written into the same folder as the
metadata file, so that folder must allow the current user to create files.

If every sequence is truly a complete genome of the same pathogen, explicitly
confirm that on the command line:

```console
python run_pipeline.py "path/to/sequences.fasta" "path/to/metadata.csv" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50 --complete-genomes
```

The flag is a user assertion, not a method for inferring genome completeness.
It cannot override contradictory metadata or FASTA headers.

An NCBI Virus nucleotide FASTA and Results Table CSV can be supplied directly
when the CSV contains the accepted fields above. Extra NCBI columns are ignored
by the internal analysis contract.

Threshold meanings:

- `--max-p-distance`: largest allowed mismatch proportion among comparable bases.
- `--max-days-apart`: largest allowed absolute collection-date difference.
- `--max-km-apart`: largest allowed Haversine distance when coordinates exist.

These are prototype parameters, not universal biological outbreak definitions.
If coordinates remain missing after optional geocoding, a pair can satisfy the
geographic rule only when its normalized location names match. Missing
p-distance or date information prevents a pair from being related.

To run only the new input validation and matching step:

```console
python validate_inputs.py "path/to/sequences.fasta" "path/to/metadata.csv"
```

To check the 10-record NCBI files currently in this project folder:

```console
python validate_inputs.py "sequences.fasta" "sequences.csv"
```

After installing MAFFT, run those same files through the complete pipeline
with demonstration-only thresholds:

```console
python run_pipeline.py "sequences.fasta" "sequences.csv" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50
```

## Run the local website

```console
python app.py
```

Open:

```text
http://127.0.0.1:5000
```

The page has separate upload areas for FASTA sequences and XLSX/CSV metadata,
two genomic-region choices, and three thresholds. Automatic detection is the
default; Complete genomes needs its own confirmation checkbox. Flask saves
uploads in a unique temporary directory, calls `run_pipeline()` directly,
then reads the existing CSV outputs. It displays strict potential clusters
separately from **genomic-temporal candidate pairs with unresolved geography**.
Those candidate pairs pass the genetic and date limits, have no usable
kilometer comparison, and cannot pass the named-location fallback. Pairs with
a known distance above the limit are *not* called unresolved. Candidate pairs
are not added to the strict related-pairs or cluster-membership CSVs and are
not confirmed outbreaks. The page also shows MAFFT, GeoNames, comparison, and
location-precision counts plus the input filenames. Its main summary separates
unique samples from pairwise relationships: for example, 22 candidate pairs
can be relationships among the original 10 samples, not 22 extra samples.
The candidate-pair list starts collapsed; select **Show candidate pairs** to
inspect it. P-distances are rounded to six decimal places on the webpage only;
the backend calculations and CSV values retain their precision. The HTML
performs no scientific calculations. Keep the terminal open; press `Ctrl+C`
to stop the site.

### Check MAFFT and GeoNames on the server

Near the top of the website, click **Check MAFFT & GeoNames**. No FASTA or
metadata upload is needed. This button is separate from an analysis run and
from Render's lightweight `/health` endpoint. The checks do **not** run when
someone merely opens the page.

- MAFFT is located through the same cross-platform lookup used by the pipeline.
  The diagnostic runs a real alignment of two short, invented sequences and
  checks that both return at equal aligned lengths. It shows a version only
  when MAFFT returns one clearly. Linux `mafft` and Windows `mafft.bat` are
  supported.
- GeoNames makes one secure, certificate-verified lookup of Tempe, Arizona
  through the existing GeoNames integration and checks that valid coordinates
  come back. The page never displays `GEONAMES_USERNAME`. A missing username,
  connection problem, or API error gets a distinct readable status.

Statuses are **Operational**, **Unavailable**, **Not configured**, or
**Check failed**. The server reuses a recent result for five minutes to limit
repeated MAFFT work and GeoNames requests. A successful status confirms these
dependencies responded to this small test; it does not validate a user's
dataset or prove an outbreak.

To test locally, start the website with `python app.py`, open
`http://127.0.0.1:5000`, and click the button. On Render, push the changed
files to GitHub, wait for the Docker Web Service to redeploy, open its public
`onrender.com` page, and click the same button. If GeoNames says **Not
configured**, set `GEONAMES_USERNAME` in Render's Environment settings. If
MAFFT is **Unavailable**, inspect the Docker build log for the Linux `mafft`
package. `/health` should still simply return `OK` without running either
diagnostic.

**Run Example Analysis** on the webpage uses `examples/synthetic_sequences.fasta`
and `examples/synthetic_metadata.csv`. These are invented records, not NCBI
accessions, real specimens, or actual outbreaks. One synthetic sample has a
Tempe place name but no coordinates; with a configured GeoNames account its
approximate place coordinates can be looked up. Without GeoNames, the existing
specific-place name fallback remains. The example includes sequences of
different original lengths so MAFFT must align them. The button uses the same
threshold settings and the same six-stage backend as an upload.
With GeoNames disabled and MAFFT mocked in the automated test, the synthetic
example produces 28 genomic pairs and two strict potential clusters. A real
MAFFT run may change the exact pairwise values.

To try the synthetic data from the command line instead:

```console
python run_pipeline.py "examples/synthetic_sequences.fasta" "examples/synthetic_metadata.csv" --max-p-distance 0.05 --max-days-apart 14 --max-km-apart 50
```

## Publish the website with GitHub and Render Docker

Visitors to the hosted website need only a browser. The Docker image installs
Python 3.11, the requirements, and the Linux `mafft` executable on the cloud
server. Windows users still run locally with `python app.py`; Gunicorn is for
the Linux container. The container runs one Gunicorn worker with a 180-second
analysis timeout. The combined FASTA-and-metadata upload request is limited to
10 MB. Large datasets may need a different hosting plan or workflow.

### 1. Put the source on GitHub

Create an [empty repository in GitHub](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository)
(do not add a GitHub-generated README or `.gitignore`). In Windows Command
Prompt, from this project folder:

```bat
git init
git add .
git status --short
git commit -m "Prepare Sprint 2 Docker website"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/YOUR-REPOSITORY.git
git push -u origin main
```

Replace the example repository URL with your own. Before `git commit`, review
`git status --short` and confirm **no real sample data, credentials, `.env`
files, generated outputs, or virtual environment** will be published. The
provided `.gitignore` omits root-level FASTA/CSV/XLSX data and `outputs/`, but
it cannot remove files that were already committed in an earlier Git history.
The files needed by Render are the root Python modules, `requirements.txt`,
`Dockerfile`, `templates/index.html`, and `examples/` synthetic inputs. Include
the README and tests in GitHub for explanation and verification; Docker only
copies the runtime files it needs.

### 2. Create a Render Free Docker Web Service

1. Sign in to [Render](https://dashboard.render.com/) and connect your GitHub
   account/repository.
2. Choose **New > Web Service**, select this repository and its `main` branch.
3. Choose **Docker** as the runtime and **Free** as the instance type. Keep the
   repository root as the service root and use the root `Dockerfile`. No separate
   build or start command is needed because the Dockerfile contains both.
4. Under advanced settings, set the health-check path to `/health` if offered.
   The container listens on `0.0.0.0:${PORT}` (default `10000`).
5. In the service's **Environment** page, add the key `GEONAMES_USERNAME` with
   your GeoNames web-service username as its value. Do not put that value in
   GitHub, Dockerfile, or the README. Save and deploy the change. If this
   variable is omitted or GeoNames is unavailable, analysis continues with
   missing coordinates and the existing specific-place fallback.
6. Wait for Render to report a successful deploy, then open its assigned
   `https://...onrender.com` address. Check `/health` for `OK` and open `/` for
   the upload form. **This repository has not been deployed to Render yet.**

Render's [Docker guide](https://render.com/docs/docker),
[Web Service setup](https://render.com/docs/web-services), and
[environment-variable instructions](https://render.com/docs/configure-environment-variables)
describe the current dashboard. The [Free plan](https://render.com/docs/free)
may spin down after idle time, so the first visit can be slow. Its filesystem
is temporary; this project intentionally deletes uploaded files and results
after each request. Do not upload identifiable patient data to this public,
unauthenticated class prototype. GeoNames queries also send location names to
the external GeoNames service.

### 3. Check the hosted analysis

- Click **Run Example Analysis** with the displayed demonstration thresholds.
  Its sequences and metadata are explicitly synthetic. Confirm the page shows
  an aligned-sequence count and results; it must not claim real outbreaks.
- For the real NCBI exercise, use your own downloaded NCBI Virus FASTA and
  Results Table CSV in the two upload fields. With demonstration thresholds
  `0.05`, `14`, and `50`, the existing 10-sample, state-only dataset produces
  45 unique comparisons and (in the current tested run) 22 genomic-temporal
  candidate *pairs* among 9 samples, but no strict geographic clusters. These
  are screening settings, not universal outbreak definitions. Real MAFFT
  alignment or different source data may change the exact pair counts.
- If a hosted upload fails, read the message on the page and the Render deploy
  or runtime logs. Confirm MAFFT installed during the Docker build and that
  the optional GeoNames variable was saved in Render. `GEONAMES_USERNAME` is
  never needed in a visitor's browser.

If Docker Desktop is available locally, an optional pre-deploy check is:

```bat
docker build -t bmi461-sprint2 .
docker run --rm bmi461-sprint2 mafft --version
docker run --rm -p 10000:10000 bmi461-sprint2
```

Then visit `http://127.0.0.1:10000`. Stop the container with `Ctrl+C`.

## Output files

Output names are based on the metadata filename. For `metadata.csv`, the runner
creates:

| File | Role |
|---|---|
| `metadata_pairwise_comparisons.csv` | All valid within-group genomic pairs |
| `metadata_pairwise_comparisons_temporal_geographic.csv` | Genomic pairs plus days, kilometers, same-location status, coordinates used for each sample, their sources (`Provided`, `GeoNames`, or blank), and geographic-precision labels |
| `metadata_pairwise_comparisons_temporal_geographic_related_pairs.csv` | Only pairs passing all prototype thresholds |
| `metadata_pairwise_comparisons_temporal_geographic_cluster_membership.csv` | `C1`, `C2`, and later cluster-membership rows |
| `metadata_region_assignments.csv` | Each sample's resolved `Gene_Segment` and its source |

## Included dual-input demonstration

The `outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/` folder contains
`example_sequences.fasta` and `example_metadata.csv`. These historical output
files are local-only and excluded from GitHub and the
Docker image. The bundled `examples/` synthetic files above are the portable
hosted demonstration. Use your own NCBI downloads for the hosted real-data test.
If you have the local historical files, run them with:

```console
python run_pipeline.py "outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/example_sequences.fasta" "outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/example_metadata.csv" --max-p-distance 0.05 --max-days-apart 7 --max-km-apart 10
```

With GeoNames lookup disabled, these demonstration-only settings produced 36
genomic pairs, 4 related pairs, and 2 potential clusters containing 5 samples
in the automated end-to-end test with a mocked MAFFT pass-through. A real MAFFT
run is still required to verify final results on this laptop. These settings
are not biological cutoff recommendations.

The same folder also contains a 10-record NCBI-style compatibility example with
versioned FASTA IDs, NCBI column aliases, and no latitude/longitude columns:

```console
python run_pipeline.py "outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/ncbi_virus_style_sequences.fasta" "outputs/01a0e5b9-470b-7032-a8c4-c8465c89040b/ncbi_virus_style_metadata.csv" --max-p-distance 0.05 --max-days-apart 7 --max-km-apart 10
```

The automated end-to-end test matched all 10 metadata rows to all 10 FASTA
records with zero validation errors. With GeoNames lookup disabled and a mocked
MAFFT pass-through, it produced 45 genomic pairs, 4 related pairs, and 2
potential clusters. Kilometer values remained blank, and matching normalized
location names supplied the existing geographic fallback. Actual MAFFT output
or optional GeoNames place coordinates may change pairwise results.

## Run tests

```console
python -m unittest discover -s tests -v
```

The tests cover FASTA/XLSX matching, FASTA/CSV
matching, NCBI aliases, versioned accessions, optional coordinates, leading-zero
accessions, true unmatched and duplicate accessions, invalid sequences, metadata
validation, optional genomic-region labels and contradictions, MAFFT grouping
and output checks, GeoNames caching and fallback,
genomic comparison, temporal/geographic enrichment, threshold logic, transitive
clustering, complete backend runs, a Flask dual-upload run, and the on-demand
MAFFT/GeoNames status checks. External program/network behavior is mocked;
these tests do not install or run MAFFT.

## Scientific and technical limitations

- Sequences within a pathogen/gene group must represent the same genomic region;
  MAFFT now aligns them, but it cannot make unrelated regions comparable.
- FASTA descriptions must state a region explicitly; unknown descriptions are
  not interpreted from length or presumed to be a whole genome.
- MAFFT must be available on `PATH` (installed locally by the user or in the
  hosted Docker image).
- GeoNames may choose an ambiguous place and supplies approximate place
  coordinates, not specimen-level collection coordinates.
- Genomic similarity does not prove direct transmission.
- Temporal and geographic proximity do not prove epidemiological linkage.
- Appropriate thresholds and mutation rates differ by pathogen and region.
- Incomplete sampling and missing metadata can change connected components.
- `N` and `-` positions are excluded from direct nucleotide comparison.
- Transitive clustering can connect A and C through B without an A-C link.
- The prototype builds no phylogenetic tree or phylodynamic model.
- It uses no machine learning and produces no confirmed outbreak labels.
- NCBI exports using other column names or complex identifier conventions may
  still require preprocessing.
- `python app.py` starts Flask's local development server; the hosted Docker
  image uses Gunicorn. Public hosting has no authentication, so use only
  synthetic or appropriately de-identified, permitted data.

## AI tools acknowledgment

ChatGPT and Claude were used as supporting references during project development,
including brainstorming, code explanation, debugging assistance, and documentation
review. Their suggestions were reviewed by the project team, and the Python
backend was verified with the automated tests described above. These tools were
not treated as scientific evidence or as substitutes for peer-reviewed genomic
epidemiology guidance.

- OpenAI. *ChatGPT*. <https://chatgpt.com/>. Accessed September 28, 2026.
- Anthropic. *Claude*. <https://claude.ai/>. Accessed September 28, 2026.
