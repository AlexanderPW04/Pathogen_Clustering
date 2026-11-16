# Homepage + Sprint 4 UI — progress log

Branch `feature/homepage-sprint4-ui-20261116`, created from `145dc4c` (local `fix/final-polish-20261115`).
Its code is identical to the staging-validated `8690561`: the two later commits only add
`docs/FINAL_POLISH_STAGING_VALIDATION_REPORT.md`. Nothing was pushed, merged, or deployed. The recovery
branch `backup/sprint5-staging-validated-20261115` (`24f81fe`) was not touched.

## 2026-11-16

### 1. Inventory (no code changed)
- **Homepage, before:**
  - order: header (long paragraph and five step pills), "How this works", then a two-column grid;
  - left column: "Step 1 · Choose your files", a quick guide, then the mode radios (Detect / Individual
    / Pooled / Combined) about 450 px down at 1024×768, the upload fields, the Sprint 4 checkbox and
    the file guide;
  - right column: the examples.
- **Mode switching:** `MODES` and `UPLOAD_FIELDS` in `app.py`; the script in `_upload.html` shows,
  hides, and disables fields per mode. Reused unchanged.
- **Sprint 4, before:**
  - already three steps in `_literature.html` (find, review, continue) with a sticky action bar;
  - findings in a 6-column table;
  - headlines and statuses written by `literature/review.py`;
  - search started from the file review (`/enrich/start`) or from the Sprint 4 page
    (`/literature/find`);
  - the existing busy banner already blocks repeated submissions and shows the elapsed time, with no
    percentage.
- **Gaps against the brief:**
  - the mode choice was not near the top;
  - no "What is…?" help or workflow strips per mode;
  - the combined fields were not grouped;
  - there were no "Changes applied?" answers per publication;
  - "Could not check" instead of "Not checked";
  - the H5N1 headline was not the requested wording;
  - the Claude sentence was buried in a list;
  - the step 3 button was not named "Continue to Genomic Analysis";
  - there was no automatic scroll to the conclusion;
  - there was no way to preview the Sprint 4 states locally.

### 2. Baseline
- **First baseline run:** started on the working tree, then stopped, because files were being edited
  during the run.
- **Clean baseline:** in a separate worktree of `145dc4c`:
  **Ran 714 tests — OK (skipped=16)**.
  - The earlier record was 714 OK with 4 skipped. Those runs had the optional real-data variables (`FONTENELE_MMC2` and others) and MAFFT on `PATH`. Here they were not set, so 16 optional tests skip on both the baseline and this branch.

### 3. Homepage
- **Header:** the requested one-sentence description, then the screening-signal label and a link to
  the examples.
- **Choosing an analysis:**
  - three analysis cards: Individual Sequences (Sprint 2), Pooled Sequencing (Sprint 3 + 5), Combined
    Analysis (Sprint 2 + 3 + 5);
  - each card has its one-line purpose and a small HTML/CSS workflow strip;
  - "Detect from my files" stays as a fourth, smaller choice and remains the default.
- **Upload fields:**
  - "2 · Add your files", with a one-line "Needed / Optional" hint for the chosen mode;
  - fields grouped into Individual sequences, Pooled sequencing results, and Reference genome and gene
    annotation;
  - every field has a native `<details>` "What is…?" answer, with the earlier text under "Technical
    details".
- **Sprint 4 box:** labelled "3 · Optional: publication-based metadata (Sprint 4)".
- **Moved below the upload:** the file guide (one collapsed section), the examples ("No files yet? Try
  a synthetic example"), and "How this works".
- **Mode labels in the review banner:** they now name Sprint 5:
  - "Pooled sequencing results (Sprint 3 + 5)";
  - "Combined (Sprint 2 + 3 + 5)".
- **Phones:** only the current step pill is shown.

### 4. Sprint 4
- **`literature/review.py` (summaries only):**
  - "Not checked";
  - per publication: `changes` (Yes only when one of its proposals was applied) and `text_checked`
    (retrieved with text and assessed);
  - "reviewed" appears in a headline only when every publication was read;
  - the H5N1-shaped case uses the requested sentences;
  - one Claude sentence (`claude_line`);
  - `searched` is separate from the headline;
  - an offline demonstration reports "No (offline demonstration)", never a Claude analysis.
- **`_literature.html`:**
  - the conclusion box first, with the Claude sentence and the next step;
  - one card per publication with the five answers, and its details in an expandable section;
  - the per-sample table is collapsed;
  - the sticky bar has "← Back to file review" and "Continue to Genomic Analysis".
- **`app.py`:**
  - the results-page note is prefixed with the search scope;
  - the applied proposal IDs are kept, so a revisit shows "Changes applied? Yes".
- **`index.html`:**
  - "Find and Check Publications";
  - plainer loading messages;
  - after a search, the conclusion scrolls into view and receives focus.

### 5. Local preview
- **`validation/local_ui_preview.py`:** port 5058, with `/preview` linking every homepage mode and
  every Sprint 4 state.
- **Fixtures:** synthetic NCBI-shaped test records and in-memory articles; the offline mock provider
  is used.
- **Labels:** a "LOCAL UI PREVIEW" banner, plus the app's own SYNTHETIC labels.

### 6. Browser inspection
Done in the built-in browser and headless Chrome; found and fixed:
- **Retrieval failure:** the wording said "0 associated publications were retrieved" and "No usable …
  evidence was established", which overstates an incomplete review.
- **Claude step:** the step 3 hint said "no metadata changes are available" while the Claude step was
  still open.
- **Headline:** the pending case read "One publication was found. One publication may hold …".
- **Wrong place note:** "general mention only" was shown for a link that was still uncertain.
- **Offline demonstration:** "Claude analysis performed? Yes" was shown for an offline demonstration.
- **Header height:** the step pills used 150 px on phones; there they now show only the current step.

### 7. Tests
- **Changed tests:** updated for the intended wording.
- **New file:** `tests/test_homepage_sprint4_ui.py`.
- **Full suite (end):** Ran 720 tests — OK (skipped=16). That is 714 + 6 new tests, with the same skips as the baseline.
- **Example outputs:** the integrated and Sprint 2 example CSV downloads (52 files) are byte-for-byte identical to the baseline.
