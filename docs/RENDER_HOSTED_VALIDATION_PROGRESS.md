# Render staging hosted validation: progress record

This file records the progress of the Render staging validation, so the work can be resumed if Render
or browser access is interrupted.

## Decisions (from the user, 2026-11-15)

**Production branch.** The production Render website deploys from `sprint4-integration`. Pushing that
branch would therefore update production.

**Staging branch.** Commit `5a20fb0` is pushed to a **new** GitHub branch, `render-staging`, and the
staging service builds from that branch. Nothing else changes:
- `sprint4-integration` on GitHub keeps its old commit `6e122ff`;
- `main` is not touched;
- the production service is not touched.

This differs from the goal text ("push sprint4-integration"), and the user approved it to keep
production unchanged.

**Who creates the service.** The user signs in to Render in the browser pane. Claude creates only the
new staging service, on the Free plan.

## Stage 1: preflight (local)

| Check | Result |
|---|---|
| Branch / HEAD / tree | `sprint4-integration`, `5a20fb0`, clean |
| Remote | `origin` = github.com/AlexanderPW04/Pathogen_Clustering (also `class` = BMI461-ASU repo, not used) |
| GitHub before push | `origin/sprint4-integration` = `6e122ff` (production source); `origin/main` = `eb72013` |
| Secrets in Git | `.env` ignored and never tracked; `ai_enrichment.example.env` holds placeholders only (model name, commented URL) |
| Dockerfile modules | Every module imported by `app.py` is copied: Sprint 2 scripts, `genomics/`, `literature/`, `wastewater/`, `templates/`, `examples/`. No `static/` folder exists or is needed |
| Gunicorn | `app:app`, `0.0.0.0:${PORT:-10000}`, 1 worker, gthread 4 threads, timeout 180 |
| Health route | `/health` returns `OK` (`app.py:256`) |
| Temporary storage | `tempfile.gettempdir()` (`/tmp`), writable for `appuser`; `MPLCONFIGDIR=/tmp/matplotlib` |
| Limits | 10 MB per request, `MAX_OBSERVATION_CELLS=700000`, `MAX_WORKSHEET_CELLS=1500000`; one review at a time; one analysis at a time |
| Windows paths | None in runtime code |
| Flask secret | Not used: the app sets no cookies and has no `SECRET_KEY` setting |

## Staging service settings (no secrets)

| Setting | Value |
|---|---|
| Type / runtime | Web Service, Docker, existing `Dockerfile` |
| Repository / branch | AlexanderPW04/Pathogen_Clustering, `render-staging` |
| Instance | Free |
| Health check path | `/health` |
| Auto-deploy | On (staging only) |

**Environment variables**

| Variable | Staging value | Required? |
|---|---|---|
| `AI_ENRICHMENT_ENABLED` | `false` | Keep paid Claude off |
| `ANTHROPIC_API_KEY` | not set | Optional, and must stay unset |
| `AI_ACCESS_CODE` | not set | Only needed if AI is ever enabled |
| `GEONAMES_USERNAME` | the user's own value, entered by the user in Render; without it, GeoNames lookups are skipped with a notice | Optional |
| `NCBI_EMAIL`, `NCBI_API_KEY` | optional; free NCBI E-utilities work without them at a lower rate | Optional |
| `MAX_WORKSHEET_CELLS`, `MAX_OBSERVATION_CELLS` | already set in the Dockerfile | — |
| `PORT` | provided by Render | — |

## Log
- 2026-11-15: preflight done. Production branch confirmed as `sprint4-integration`; staging branch
  decided as `render-staging`. Render is not signed in yet in the browser pane. Local Docker engine is
  starting.
- 2026-11-15, blockers:
  - **Docker:** Docker Desktop launched, but its Linux engine did not start (the `docker-desktop` WSL
    distribution stays stopped). The local image test was not run.
  - **Push refused:** `git push origin 5a20fb0:refs/heads/render-staging` failed with
    "Repository not found". The only credential on this machine was a different GitHub account, which
    could not see `AlexanderPW04/Pathogen_Clustering` (presumably private). Nothing was pushed.
  - **Render:** after the user's sign-in, the built-in browser pane still shows Render's login form.
- 2026-11-15: **the user will push `render-staging` and create the staging service** (dashboard steps
  given). Claude waits for the staging URL.
- 2026-11-15: **local baseline recorded.**
  - Script: `validation/hosted_smoke_test.py`.
  - Server: Flask on 127.0.0.1:5055 with AI off, GeoNames unset, `MAX_WORKSHEET_CELLS=1500000` and
    `MAX_OBSERVATION_CELLS=700000`.
  - Output: `<scratchpad>/hosted/local/record.json`.
  - All scenarios behaved as expected:
    - **Sprint 2 synthetic:** 4 strict pairs, 2 clusters.
    - **NCBI H5N1 ten-record set:** 0 strict pairs; 2 same-segment pairs blocked only by state-level
      geography.
    - **Integrated example:** run 14.6 s.
    - **mmc2 with FASTA and GenBank:** run 12.5 s.
    - **mmc2 alone:** runs.
    - **Formula check:** 0 formula cells.
    - **Unsupported input:** FASTQ refused; 1.6 M-cell CSV refused; an upload over 10 MB gets a 413.
    - **Busy check:** the second run is told the server is busy.
    - **Sprint 4 H5N1:** the 2022 article is rejected for the 2026 samples; no change; Claude not run.
  - Next: run the same script with `--base <staging URL>`.
- 2026-11-15: **hosted run done.**
  - Staging: https://pathogen-clustering-2.onrender.com, branch `render-staging` = `5a20fb0`
    (verified by `git ls-remote`). Production branch unchanged.
  - Output: `<scratchpad>/hosted/staging/record.json`; content comparison in
    `<scratchpad>/hosted/compare.txt`.
  - Report: `docs/RENDER_HOSTED_VALIDATION_REPORT.md`, verdict **PASSED WITH LIMITATIONS**.
  - Still requested from the Render dashboard: deploy event, build log, runtime log, memory metrics.
- 2026-11-15: the user supplied the build and runtime logs. Commit `5a20fb073aa8…` was checked out from
  `render-staging`; the build succeeded; gunicorn gthread started; one worker boot; no errors. The
  report was updated. Only the memory graph is still missing.
