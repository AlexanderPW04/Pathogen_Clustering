# Final polish: progress record

Read this file first if the work is interrupted. Continue from the last milestone marked **Done**.

## Milestone 0: recovery checkpoint — Done (2026-11-16)

| Item | Value |
|---|---|
| Baseline application commit | `5a20fb073aa828a5dfd411309a08f83b2e039baf` (validated on Render staging) |
| Backup branch (GitHub) | `backup/sprint5-staging-validated-20261115` |
| Backup commit | `24f81fe3d3a5d1b87e5bbe08930315eae4f7f693` |
| Backup verification | Verified with `git ls-remote`, `git fetch`, and the GitHub API branch endpoint. The backup has 231 files, identical to the local commit |
| Local tag | `checkpoint/sprint5-staging-validated-20261115` (annotated, on `24f81fe`; not pushed). No existing tag was changed |
| Working branch | `fix/final-polish-20261115`, created from `24f81fe`. All improvements happen here |
| Staging (unchanged) | https://pathogen-clustering-2.onrender.com, built from `render-staging` = `5a20fb0` |
| Production (unchanged) | Builds from `sprint4-integration` = `6e122ff` |
| `main` (unchanged) | `eb72013` |
| Baseline tests | `Ran 685 tests in 616.7 s — OK (skipped=4)` on `fix/final-polish-20261115` with real MAFFT and `FONTENELE_MMC2`, `CRITS_CHRISTOPH_DIR`, `SCHUMANN_DIR`, `ANNOTATION_MPOX_GENBANK`. The 4 skips: 3 SciPy cross-checks, and the NCBI Virus files test |

### What the backup contains

`24f81fe` is commit `5a20fb0` plus one more commit. That commit adds only:
- `docs/RENDER_HOSTED_VALIDATION_REPORT.md`;
- `docs/RENDER_HOSTED_VALIDATION_PROGRESS.md` (one account name generalised);
- `validation/hosted_smoke_test.py`.

The backup holds all source code, tests, the Dockerfile, the examples, and the documentation.

**Excluded:** `.env` (git-ignored, never tracked), credentials, and the published datasets (kept
outside the repository). Before the push, the checks found no key pattern in the tree and no
`.env`, `.pem`, `.key` or secret files.

### How to restore the validated version safely

Nothing is restored automatically. Choose the case that applies.

1. **Inspect or run the validated code locally**, without changing any branch:
   `git fetch origin` and then
   `git switch --detach origin/backup/sprint5-staging-validated-20261115`.
2. **Start new work from the checkpoint:**
   `git switch -c <new-branch> origin/backup/sprint5-staging-validated-20261115`.
3. **Put staging back on the validated code.** `render-staging` already points at `5a20fb0`, so
   there is nothing to do as long as that branch is not pushed to. If it ever changes, ask before
   pushing, and do not force-push. Then point the Render staging service's branch at
   `backup/sprint5-staging-validated-20261115` in the dashboard, and deploy.
4. **Production:**
   - Production stays on `sprint4-integration`. Promotion is a separate, explicitly approved step:
     a normal (non-force) push of the chosen validated commit to that branch, followed by
     re-running `validation/hosted_smoke_test.py` against production.
   - Never delete or rewrite `backup/sprint5-staging-validated-20261115`.

## Milestones 1-6

Not started.

## Milestone 1: scientific labelling and biological correctness — Done (2026-11-16)

**1A: UTR genome band.**
- `genomics/integrated_analysis.py`:
  - `annotation_file_bands` draws the gene track from the uploaded GenBank/GFF3 CDS coordinates.
    Overlapping CDS of one gene (ORF1a/ORF1ab) are merged into one band.
  - `annotation_ranges` is used without an annotation file. It now draws a source label only when its
    records form one uninterrupted run, so a recurring label such as Fontenele's intergenic "UTR" is
    left unannotated.
  - `gene_track` chooses between the two and writes the figure caption, which states the track's
    source.
  - The annotation is now read once, before the figure, and reused by `run_mutation_annotation`.
- Verified on NC_045512.2:
  - bands equal the NCBI CDS: S 21563–25384, ORF3a 25393–26220, N 28274–29533, ORF1ab 266–21555;
  - 25385–25392 stays unannotated;
  - without the file, the Fontenele source labels give 13 bands and no "UTR".

**1B: influenza "NA".**
- `genomics/report_exports.StoredTables.get` reads with `keep_default_na=False, na_values=[""]`. The
  XLSX and PDF are rebuilt from the CSVs, so "NA" stays text in both; only empty cells are missing.
- Sprint 2 readers also kept "NA" as text after this change:
  - `compare_genetic_distances.read_pair_table`;
  - `app.run_relationship_analysis`, using the new `TEXT_LABEL_COLUMNS`: Sample IDs, Pathogen,
    Gene_Segment.

**1C: genetic code.**
- `genomics/genome_annotation.py` records `Translation_Table_Source` for every CDS: stated / verified
  (table 1 reproduces the GenBank `/translation`) / confirmed by the user.
- Otherwise the CDS is unsupported with `TABLE_UNVERIFIED_DETAIL`, so no amino-acid consequence is
  predicted for it.
- A new optional setting, "Genetic code for coding regions that state none":
  - website: Advanced settings `genetic_code` (validated by `parse_genetic_code`);
  - CLI: `--genetic-code`.
  - The review page warns when an annotation needs it.
- A stated table is never overridden.
- The new `Translation_Table_Source` column appears in `sprint5_annotation_sources.csv`.

**Intentional test changes (documented in each test).** GFF3-based tests now confirm the standard
code explicitly (`genetic_code=1`), because GFF3 states no code and the silent default was the defect:
- `tests/test_genome_annotation.py`: 9 calls;
- `tests/test_snv_consequences.py`: 2 calls;
- `tests/test_sprint5_phase4_web.py`: 4 uploads.

**New tests.** `tests/test_final_polish_annotation.py` (12): code 1 versus code 4 (TGG>TGA is W2* versus
W2W), stated / verified / confirmed / unverified codes, NCBI bands, recurring labels, NA / HA / M1; M2 /
NS versus a blank, and Sprint 2 "NA" segments.

**Results.**
- Annotation, website and integrated modules: 64 tests OK.
- Import, web, analyses, app and Sprint 2 modules: 90 tests OK.
- The full suite runs at Milestone 5.

**Remaining (documented, not changed):** a start-codon change to another start codon (for example
ATG>CTG, or GTG>ATG under table 11) is still reported as "Start lost".

## Milestone 2: long-run and download feedback — Done (2026-11-16)

**2A: processing feedback** (`templates/index.html`, one script). Any form that waits for the server
immediately shows a fixed banner with:
- a spinner (static when the visitor prefers reduced motion);
- a message for the step, for example "Processing genomic data. Large sequencing datasets may take
  several minutes.";
- the **real elapsed time**. No percentage is shown, because the server reports none.

Repeated submissions are ignored, and all buttons are marked busy only after the browser has
collected the form, so button values are still sent. The state resets on Back (`pageshow`). The
profile download (`data-no-busy`) leaves the page usable. The old "up to two minutes" text is gone.

**2B: report downloads.**
- The XLSX and PDF links show "Preparing Excel workbook / PDF report… the download starts by itself",
  then "ready", or the server's error text.
- Repeated clicks are ignored on the page.
- On the server (`app.build_report_download`), a per-file lock makes repeated or concurrent requests
  wait for the first build and reuse its cached copy. The lock gives up after 150 s with a clear 503;
  the 8 MB cache rule and the one-at-a-time slot are unchanged.
- Expired or unknown results now return a clear 404 page ("kept in memory for 10 minutes … lost if the
  server restarts").

**Browser check** (local server through the browser pane):
- The banner appears at once with the message and an elapsed clock (0:01, then 0:02). All 5 buttons
  are marked busy, the 2nd and 3rd clicks are blocked, and Back resets the state.
- The XLSX and PDF links show "Preparing…", ignore repeated clicks, then show "ready".
- An expired link shows "could not be prepared: …".

**Tests.** `tests/test_final_polish_feedback.py` (4): three concurrent requests give one build and
identical bytes; a failed build gives a clear 500 and can be retried; expired links give 404 with an
explanation; the page has the feedback markup.

No background-job system was added: the request architecture finished every hosted run, so it was
not needed.

## Milestone 3: website organisation and navigation — Done (2026-11-16)

Every result and figure is kept; no page was redesigned.

**3A: homepage** (`templates/_upload.html`, `_how_it_works.html`, `index.html` CSS).
- "How this works" starts collapsed.
- A compact **"Which files do I need?"** box leads Step 1. It says you never need every type, and
  names individual FASTA + metadata, pooled SNV data, complete A/C/G/T counts, a reference FASTA, a
  gene annotation, and supplementary metadata.
- The file fields now come before the optional Sprint 4 opt-in.
- The detailed help ("More detail on which files to use") follows the Review button.
- The upload card takes about two-thirds of the width; the mode choices are in two columns and show
  only the selected mode's description (`:has`, with a fallback that shows all of them).
- Automatic recognition is unchanged.

**3B: review page.** A status box at the top shows "Ready to run …" with a link to Run, or
**"Required information to fix before running"** with the problems. The problems are now listed once,
at the top; before, they were at the bottom. Sections are numbered (1. Files recognized … 4. Run
analysis), and the Run bar stays sticky.

**3C: results page.**
- The jump bar adds **Downloads** and the Sprint 5 parts: Uncertainty, Mutation annotations, and
  Experimental grouping.
- Sections with nothing to show are marked "(n/a)".
- On phones the bar is a sticky, sideways-scrolling single line.
- The downloads bar has `id="downloads"`.

**3D: Sprint 4.** Already states, in order: publications found, location information, whether it
applies, possible improvement, whether Claude ran, and whether anything changed. This was verified on
staging (hosted report §7). No change was needed.

**Browser measurements** (built-in browser, local server; scripted, not a human test):

| Element | Before (hosted, 1024×768) | After |
|---|---|---|
| File picker | 1,778 px | **956 px** |
| Review button | 1,939 px | **1,332 px** |
| File guide | after the button (2,010 px) | **374 px**, first screen |
| Homepage height | 2,737 px | 2,130 px |

- **Phone (375×812):** guide at 595 px, file picker at 1,490 px, examples after the form, no sideways
  scrolling.
- **Review page (1280×800):** the status box is at 504 px (first screen), and the sticky Run button is
  visible at the bottom of the window.
- **Results page:** Downloads at 696 px and Major findings at 1,165 px. Every jump link works, and on a
  phone the bar stays at the top (40 px tall).

**Remaining for a human test:**
- At 1024×768 the file picker is still one scroll down, because the mode chooser comes first.
- The review and results pages are still long (6.7k px and 15.8k px), though navigable.

**Tests.** `tests/test_final_polish_navigation.py` (3). The UI, web, literature, and Sprint 5 web
modules all pass (77 tests).

## Milestone 4: robustness — Done (2026-11-16)

**4A: mapping profiles** (`genomics/mapping_profiles.validate_profile`, called by `read_profile`).
- **Checked:**
  - only known top-level and per-file keys;
  - a supported format and version;
  - 1–20 file entries, each a JSON object;
  - roles from the review page (never "ignore" or paths);
  - field names that exist for that role;
  - plain-text values of at most 200 characters, and no empty column names;
  - a valid creation date;
  - `date_order` in `DATE_ORDERS`.
- **Per-file options** must be exactly what the review form accepts: the same `read_file_options`
  rules (header row 1–50, numbering 0/1, fraction/percent, base-count minimum 1–1000, known depth
  definitions).
- **Worksheets:** a profile naming a worksheet the workbook lacks does not fit.
- **On failure,** the whole profile is refused with a specific message ("'header_row' is -5; it must
  be a whole number from 1 to 50"); the review stays open and nothing is applied.
- **Compatibility:** all 12 profiles saved earlier (local and staging, old code) still validate.

**4B: result cache** (`app.make_room_for_results`, `store_result_files`, `result_download`).
- **Bounds:** a memory budget of 48 MB compressed by default (`RESULT_CACHE_MB`), plus at most 12
  result sets (was 4, count only).
- **Eviction order:** expired sets first, then the least recently used. Creating a set and
  downloading from it both count as use.
- A first design (three priority groups) was replaced after the full suite found that it could evict
  the newest, not-yet-downloaded result first (test `test_upload_enrichment`).
- **Use tracking:** each download refreshes a set's last use. A cached XLSX/PDF is kept only while the
  budget allows it.
- **Unchanged:** the 10-minute lifetime, the 128-bit random IDs (results stay separate), memory-only
  storage (lost on restart, explained on the 404 page), and the 8 MB per-file rule.

**4C: Schumann categorical storage — not implemented (documented).**
- The observation matrix (533 → 50 MB as categories, identical CSV bytes) feeds at least 9 consumers.
- Several of them group or join on those columns: `temporal_analysis.frequency_trends`,
  `vaf_trends` merges, `figures` groupby, `mutation_profiles`.
- With pandas 2.x, grouping by a categorical column defaults to `observed=False`, which creates every
  category combination. That could change outputs and raise memory.
- Proving equivalence would need every consumer audited, so it stays remaining work.
  `MAX_OBSERVATION_CELLS` keeps refusing Berlin35 on the free plan with a clear message.

**Tests.** `tests/test_final_polish_robustness.py` (8): every invalid option, key, role, field, and
structure is refused; a web upload of a bad profile shows the reason; an active visitor keeps results
while another runs many; the memory budget is never exceeded; expiry and separation hold; 30
concurrent stores stay bounded.

## Milestone 5: regression testing — Done (2026-11-16)

- **Full suite:**
  - First run after Milestones 1–4: 712 tests, **1 error**. The first cache design evicted the newest,
    not-yet-downloaded result (`test_upload_enrichment`). Fixed with least-recently-used eviction plus
    a regression test.
  - Final: **`Ran 713 tests in 634.5 s — OK (skipped=4)`** (685 baseline + 28 new).
- **Published datasets:** 33 PASS, 0 FAILED (Schumann VEP 1,830/1,830; files unchanged).
- **Output equivalence** with the validated version, browser-like test, same inputs: **129 of 130 CSVs
  identical**. The only difference is the new `Translation_Table_Source` column. Main findings,
  readiness, Sprint 4, the formula check, unsupported input, and the busy message are identical.
- **Real H5N1:** gene "NA" kept in 10/10 workbook cells (blank before), with code 1 confirmed.
- **Free-plan imitation** (512 MB / 10% CPU), `mmc2`: run 127.1 s (123.5 s before; suite running in
  parallel), peak 300 MB (302), XLSX 41.0 s, PDF 5.9 s. No regression found.

## Milestone 6: final report — Done (2026-11-16)

`docs/FINAL_POLISH_IMPLEMENTATION_REPORT.md`: **FINAL POLISH COMPLETE — READY FOR STAGING
REVALIDATION**. Nothing has been pushed or deployed.
