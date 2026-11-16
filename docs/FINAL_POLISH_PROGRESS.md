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
