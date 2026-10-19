# Beginner-Friendly Python Code Guide

This guide documents the **pre-MAFFT/GeoNames Sprint 2 baseline**. It explains
the six Python files that existed then. The October 19 backbone update added
`align_sequences.py`, `geocode_locations.py`, and `genomic_regions.py`, and
changed the runner and validator, so some
line numbers and the pipeline order below are historical. For current setup and
workflow, use `README.md`; for changes since the baseline, use
`SPRINT_3_CHANGELOG.md`. Multiline statements are described together because
the opening line, indented values, and closing parenthesis are one Python
instruction. Blank lines only separate sections and do not affect the program.

## Historical program flow (before the MAFFT/GeoNames update)

```text
FASTA + CSV/XLSX metadata
          |
          v
validate_inputs.py
          |
          v
compare_genetic_distances.py
          |
          v
add_temporal_geographic.py
          |
          v
detect_potential_clusters.py
          |
          v
CSV results

run_pipeline.py calls all four stages.
app.py lets a browser call run_pipeline.py.
```

The refactor deliberately uses:

- ordinary functions instead of classes
- basic `for` loops instead of dense comprehensions
- named variables instead of long nested expressions
- dictionaries and sets instead of a graph library
- comments before each major step
- no regular expressions

---

# `validate_inputs.py` — lines 1–549

## Why this file matters

This is the pipeline's input gatekeeper. It creates the single standardized
DataFrame used by every later stage. It reports errors without silently deleting
rows.

## Imports — lines 3–8

| Line | Import | What it does and why it matters |
|---|---|---|
| 3 | `argparse` | Reads FASTA and metadata paths from a command-line command. |
| 4 | `sys` | Provides `sys.stderr`, the separate stream used for error messages. |
| 5 | `Path` | Handles file paths, extensions, existence checks, and absolute paths on different operating systems. |
| 7 | `pandas as pd` | Loads CSV/XLSX tables, handles missing values, converts dates/numbers, and stores the internal DataFrame. |
| 8 | `SeqIO` | Biopython's tested FASTA parser. It avoids writing a fragile custom parser. |

## Constants — lines 12–60

- Lines 12–19 list the required standardized metadata fields.
- Line 22 documents that latitude and longitude are optional.
- Lines 25–34 give the full metadata order used in summaries.
- Lines 37–42 map NCBI Virus names to internal names.
- Lines 45–56 define the exact downstream DataFrame column order.
- Lines 58–59 list accepted FASTA and metadata extensions.
- Line 60 creates the allowed nucleotide character set `A`, `C`, `G`, `T`, `N`, and `-`.

These constants keep important rules visible near the top instead of scattering
them through many functions.

## `blank_mask` — lines 63–68

1. Converts a pandas column to nullable strings.
2. Strips surrounding whitespace and checks for an empty string.
3. Combines that check with pandas' true missing-value check.
4. Returns one Boolean per row.

Importance: an empty Excel cell, a CSV null, and a cell containing only spaces
are all treated consistently.

## `metadata_rows` — lines 71–79

1. Starts an empty row-number list.
2. Loops through the Boolean error mask with `enumerate`.
3. Adds 2 to each zero-based pandas position: one for zero-based counting and one
   for the metadata header.
4. Joins the row numbers for a readable error message.

## `clean_sequence` — lines 82–86

1. Converts the sequence to text.
2. Uses ordinary `split()` to separate all whitespace-delimited pieces.
3. Joins the pieces without spaces.
4. Converts bases to uppercase.

No regular expression is used, and nucleotide order is not changed.

## `normalize_accession` — lines 89–100

1. Returns an empty matching key for a missing accession.
2. Converts the accession to stripped text.
3. Splits at the final period using `rpartition`.
4. Removes the final part only when it contains digits.
5. Otherwise returns the original stripped accession.

Example: `QB061804.1` becomes `QB061804` for matching, but the full original ID
is still preserved for output.

## `standardize_metadata_columns` — lines 103–115

1. Starts an empty rename dictionary.
2. Loops through the four NCBI-to-standard mappings.
3. Renames an NCBI alias only when its standardized column is not already there.
4. Returns the renamed DataFrame.

The explicit loop is longer than a comprehension but easier to explain one step
at a time.

## `load_metadata` — lines 118–129

- Converts the input to a `Path` and reads its lowercase extension.
- Uses `pd.read_excel(..., engine="openpyxl")` for XLSX.
- Uses `pd.read_csv` for CSV.
- Reads initial values as strings so identifiers such as `001` keep leading zeros.
- Raises a clear error for unsupported file types.

## `validate_metadata` — lines 132–249

This is the longest function because it checks every metadata requirement.

### Lines 135–147: setup and required columns

- Copies the input so the caller's DataFrame is not unexpectedly changed.
- Applies NCBI aliases.
- Starts an error list.
- Loops through required columns and reports all missing headers together.

### Lines 149–153: missing-value counts

- Creates a dictionary of missing counts for columns that exist.
- Optional absent coordinates are not errors.

### Lines 155–175: text cleaning

- Lists the text fields.
- Skips a field if its column is absent.
- Strips spaces.
- Uses `blank_mask` and reports exact metadata rows.

### Lines 177–199: duplicate metadata accessions

- Creates normalized base accessions.
- Marks every repeated nonempty base with `duplicated(keep=False)`.
- Loops through duplicate values.
- Reports row numbers and full original accessions.

This prevents ambiguous matching such as two metadata rows that normalize to the
same accession.

### Lines 201–222: collection dates

- Saves the original date values.
- Separates blank dates from invalid date text.
- Converts dates with `pd.to_datetime`.
- Stores real datetime values for later subtraction.
- Reports missing and invalid rows separately.

### Lines 224–249: coordinates and return values

- Defines valid latitude and longitude ranges.
- Skips a coordinate column when it is absent.
- Converts present values to numbers.
- Reports nonnumeric and out-of-range values.
- Returns the validated table, complete error list, and missing counts.

Blank latitude/longitude values are allowed; bad nonblank coordinate values are
still reported.

## `load_and_validate_metadata` — lines 252–256

This small wrapper first calls `load_metadata`, then sends the resulting table to
`validate_metadata`. Other modules can reuse both operations with one call.

## `load_fasta` — lines 259–282

1. Opens the FASTA as UTF-8 inside a `with` block so it closes automatically.
2. Uses `SeqIO.parse` to read FASTA records.
3. Numbers records starting at 1 for understandable errors.
4. Saves full accession, normalized accession, cleaned sequence, and record number.
5. Returns a DataFrame with a predictable column order, even when no records exist.

## `match_fasta_to_metadata` — lines 285–411

### Lines 288–321: FASTA identifier checks

- Reports an empty FASTA.
- Finds blank FASTA IDs.
- Detects duplicates after version normalization.
- Includes record numbers and original full IDs in errors.

### Lines 323–348: sequence checks

- Loops through every record.
- Reports empty sequences.
- Uses set subtraction to find unsupported nucleotide symbols.
- Reports accession, record number, invalid symbols, and the allowed alphabet.

### Lines 350–400: matching dictionaries

- Builds normalized accession dictionaries for metadata and FASTA.
- Uses ordinary loops to find metadata-only and FASTA-only accessions.
- Sorts unmatched lists for stable messages.
- Builds lookup dictionaries for sequence and full FASTA ID.

### Lines 402–411: internal DataFrame

- Copies metadata.
- Maps normalized accessions to full FASTA IDs and sequences.
- Creates missing columns when necessary.
- Reorders to `INTERNAL_COLUMNS`.
- Stores FASTA record count in DataFrame attributes.
- Returns matched data plus errors.

## `load_and_validate_inputs` — lines 414–427

This is the main function imported by `run_pipeline.py`:

1. Loads and validates metadata.
2. Loads FASTA records.
3. Matches both sources.
4. Counts missing internal IDs and sequences.
5. Combines metadata and matching errors.
6. Returns the internal DataFrame, errors, and counts.

## `category_summary` — lines 430–452

Creates text such as `Human (5), Animal (2)`. It removes blanks, counts values,
sorts labels, uses a simple loop to create entries, and limits very long output
to ten categories plus a “more” message.

## `print_summary` — lines 455–490

Prints FASTA count, metadata count, matched rows, pathogen/gene groups, date
range, hosts, locations, missing values, and every validation error. It prints
`column absent` for optional metadata columns not included in a file.

## Command-line functions — lines 493–549

- `parse_arguments` creates the two required file arguments.
- `main` checks extensions and file existence, runs the complete validator,
  prints the summary, and returns a status code.
- Lines 548–549 call `main` only when this file is executed directly.

---

# `compare_genetic_distances.py` — lines 1–209

## Imports — lines 3–13

- `argparse`: command-line arguments.
- `sys`: stderr errors.
- `Path`: path handling.
- `pandas`: grouping and CSV output.
- Three imports from `validate_inputs`: supported suffixes and the shared input
  function.

Lines 16–25 define the exact genomic output columns.

## `compare_sequences` — lines 28–60

1. Rejects unequal lengths because this prototype does not align sequences.
2. Starts mismatch and comparable-position counters at zero.
3. Uses `zip` to visit corresponding bases.
4. Skips a position when either base is `N` or `-`.
5. Counts all remaining positions and A/C/G/T mismatches.
6. Returns missing rates when the denominator is zero.
7. Otherwise calculates:

```text
P_Distance = Differences / Comparable_Positions
Percent_Identity = Matching_Positions / Comparable_Positions * 100
```

8. Returns counts and rates in a dictionary.

## `compare_samples` — lines 63–108

- Groups rows by `Pathogen` and `Gene_Segment`.
- Converts each group into simple record dictionaries.
- Uses two nested index loops to create every unique pair exactly once.
- Calls `compare_sequences` for each pair.
- Reports unequal lengths without attempting alignment.
- Creates one output dictionary per valid pair.
- Returns the result DataFrame and all comparison errors.

The nested loops replaced the less-familiar `itertools.combinations` import.

## `run_comparison` — lines 111–148

- Calls `compare_samples`.
- Writes the genomic CSV.
- Counts and describes pathogen/gene groups with an ordinary loop.
- Prints three summary lines.
- Prints comparison errors to stderr.
- Returns 1 if any pair failed or 0 when all pairs succeeded.

## `main` and entry point — lines 151–209

- Reads FASTA, metadata, and optional output paths.
- Creates the default output name beside metadata.
- Checks extensions and file existence in separate beginner-readable steps.
- Calls the shared validator.
- Stops on validation problems.
- Calls `run_comparison` only with valid matched data.

---

# `add_temporal_geographic.py` — lines 1–242

## Imports — lines 3–10

- `argparse`: standalone command.
- `math`: Haversine trigonometry.
- `sys`: stderr errors.
- `Path`: paths and output names.
- `pandas`: CSV, missing values, dates, and DataFrames.
- Shared metadata extension rules and validation.

## `haversine_km` — lines 13–43

1. Uses an approximate Earth radius of 6,371 km.
2. Converts all coordinates from degrees to radians.
3. Calculates latitude and longitude changes.
4. Separates the formula into `first_part`, `second_part`, and
   `haversine_value` for presentation.
5. Clamps small floating-point drift to 0–1.
6. Calculates the central angle.
7. Multiplies by Earth radius to return kilometers.

## `enrich_pairs` — lines 46–121

- Copies input tables.
- Cleans sample IDs and makes them the metadata index.
- Starts three result lists.
- Looks up both samples for every pair.
- Raises a clear error when metadata is missing.
- Calculates absolute date difference when both dates exist.
- Reads all four coordinate values into named variables.
- Calls Haversine only when all coordinates exist.
- Leaves kilometers blank otherwise.
- Normalizes location text by whitespace and capitalization.
- Adds nullable day, kilometer, and same-location columns.

## `run_stage` — lines 124–161

- Reads the genomic pairwise CSV.
- Checks the two sample-ID columns with a simple loop.
- Calls `enrich_pairs` and writes a different CSV.
- Prints pair count, valid date count, and valid distance count.
- Returns status 0, 1, or 2.

## `main` and entry point — lines 164–242

- Reads metadata, pairwise CSV, and optional output path.
- Creates a default enriched filename.
- validates paths separately.
- Refuses to overwrite the genomic CSV.
- Reuses metadata validation.
- Allows unavailable dates or coordinates in this enrichment-only command.
- Copies `Accession` into `Sample_ID` only for standalone use.

The recommended full workflow uses `run_pipeline.py`, which passes the actual
FASTA-derived `Sample_ID` directly into `run_stage`.

---

# `detect_potential_clusters.py` — lines 1–243

## Imports and schemas — lines 3–28

- `argparse`: enriched CSV and thresholds.
- `math`: rejects NaN and infinite thresholds.
- `sys`: stderr.
- `Path`: paths.
- `pandas`: CSVs, missing values, grouping, and outputs.
- `PAIR_COLUMNS`: required related-pair output.
- `CLUSTER_COLUMNS`: cluster-membership output.

## `find_related_pairs` — lines 31–73

For each pair:

1. Read p-distance, days, and kilometers.
2. Reject missing genetic or date evidence.
3. Calculate named `genetic_passes` and `time_passes` Booleans.
4. Reject the pair if either is false.
5. Convert `Same_Location` to a reliable Boolean.
6. When kilometers exist, compare them with the kilometer threshold.
7. When kilometers are missing, use same-location as the fallback.
8. Copy passing output fields with a basic loop.
9. Return all passing pairs as a DataFrame.

## `build_clusters` — lines 76–145

1. Group pairs by pathogen and gene/segment.
2. Build a dictionary where each sample has a set of neighbors.
3. Add every passing connection in both directions.
4. Keep a `visited` set.
5. Start a list called `stack` for each unvisited sample.
6. Pop samples until all reachable neighbors have been visited.
7. Save components with at least two samples.
8. Give them IDs `C1`, `C2`, and later values.
9. Return one membership row per clustered sample.

This is basic depth-first search. If A-B and B-C pass, all three samples are in
one connected component.

## `run_stage` — lines 148–205

- Reads enriched pairs.
- Checks required columns with a loop.
- Converts decision fields to numeric values.
- Calls pair filtering and cluster building.
- Creates two output paths.
- Writes related-pair and membership CSVs.
- Prints four summary counts.

## `main` and entry point — lines 208–243

- Reads the enriched CSV and all three required thresholds.
- Requires an existing CSV.
- Uses a loop to reject negative, NaN, or infinite thresholds.
- Calls `run_stage` and publishes its status code.

---

# `run_pipeline.py` — lines 1–146

## Imports — lines 3–15

- `argparse`: command-line interface.
- `math`: finite threshold checks.
- `sys`: stderr.
- `Path`: path handling.
- Existing stage functions are imported and renamed clearly:
  `add_time_and_geography` and `detect_clusters`.
- The shared validator and supported extensions are imported from
  `validate_inputs.py`.

## `output_paths` — lines 18–33

Creates each filename in four explicit steps:

1. genomic comparisons
2. temporal/geographic enrichment
3. related pairs
4. cluster membership

It returns all four paths in that order. Flask uses the same helper, so the web
interface and command line cannot disagree about filenames.

## `run_pipeline` — lines 36–93

1. Gets all four output paths.
2. Prints validation progress.
3. Calls the shared validator and stops on errors.
4. Calls `run_comparison` and checks its status.
5. Calls `add_time_and_geography` and checks its status.
6. Calls `detect_clusters` with all thresholds and checks its status.
7. Prints final related-pair and membership paths.
8. Returns 0 only when every stage succeeds.

Each stage status has a named variable, making the stop conditions easy to
present.

## `main` and entry point — lines 96–146

- Reads both source files and three thresholds.
- Checks FASTA and metadata types and existence separately.
- Loops through thresholds and rejects negative/non-finite values.
- Calls the full reusable `run_pipeline` function.

---

# `app.py` — lines 1–146

## Imports — lines 3–13

| Import | Purpose |
|---|---|
| `redirect_stderr` | Captures backend validation errors for the webpage. |
| `StringIO` | Stores captured error text in memory. |
| `isfinite` | Rejects NaN/infinite form thresholds. |
| `Path` | Checks upload suffixes and creates temporary paths. |
| `TemporaryDirectory` | Isolates each run and deletes temporary files automatically. |
| `pandas` | Reads final result CSVs. |
| Flask objects | Creates the app, reads requests, and renders HTML. |
| Pipeline imports | Calls the real backend and uses its output names. |
| Extension constants | Keeps browser file rules consistent with validation. |

Line 16 creates the Flask app.

## `read_thresholds` — lines 19–35

- Converts all three form strings to floats.
- Returns a simple message when conversion fails.
- Loops through values and rejects negatives, NaN, and infinity.
- Returns the three-number list plus no error when valid.

## `make_cluster_list` — lines 38–57

- Groups membership rows by cluster, pathogen, and gene/segment.
- Converts sample IDs to a normal Python list.
- Creates one basic dictionary per cluster for the HTML template.

This helper keeps pandas grouping details out of the main webpage function.

## `index` route — lines 61–142

### Initial values

- Starts error, completion, cluster, and pair-count values.
- Reads the submitted thresholds or uses demonstration defaults.

### Upload validation

- Handles only POST submissions.
- Reads the two upload fields.
- Checks missing filenames and extensions in a clear `if`/`elif` chain.
- Calls `read_thresholds` only after file checks pass.

### Temporary pipeline run

- Creates one automatic temporary directory.
- Gives uploads fixed safe internal filenames.
- Saves both files.
- Captures stderr while calling `run_pipeline` directly.
- Shows captured errors if the pipeline returns a nonzero status.

### Successful result display

- Gets final paths from `output_paths`.
- Reads related pairs and cluster membership.
- Counts related pairs.
- Calls `make_cluster_list`.
- Marks analysis complete, including a successful result with zero clusters.
- Sends all values to `templates/index.html`.

Lines 145–146 start Flask only when `python app.py` is run directly.

---

# Function summary for presentation

| File | Most important functions | One-sentence explanation |
|---|---|---|
| `validate_inputs.py` | `load_and_validate_inputs` | Creates one valid matched sample table from FASTA and metadata. |
| `compare_genetic_distances.py` | `compare_sequences`, `compare_samples` | Counts direct aligned A/C/G/T differences within biological groups. |
| `add_temporal_geographic.py` | `haversine_km`, `enrich_pairs` | Adds days, kilometers, and normalized-location equality. |
| `detect_potential_clusters.py` | `find_related_pairs`, `build_clusters` | Applies thresholds and finds connected sample groups. |
| `run_pipeline.py` | `run_pipeline` | Runs the existing four backend stages in order. |
| `app.py` | `index` | Moves browser uploads into the Python pipeline and displays final CSV results. |

# Scientific interpretation reminder

The program produces screening signals under user-selected prototype thresholds.
It does not prove direct transmission or confirm an outbreak. Confirmation needs
pathogen-specific evidence, sequence-quality review, sampling context,
epidemiological investigation, and public-health judgment.
