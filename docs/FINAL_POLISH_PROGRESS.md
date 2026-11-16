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
