# Things To Improve — Wetlab Scientist UX Review

## Overview
I tested the system as a wetlab biologist with zero coding experience. My goal was to design TaqMan primers for a pathogen detection assay. Below are the friction points I encountered.

---

## Critical Issues

### 1. No file upload — must know server paths
- The only way to use local data is to type a filesystem path (e.g. `test_data/target`).
- A wetlab scientist expects an **Upload** button to select FASTA files from their computer.
- The backend actually has a `/api/jobs/upload` endpoint that accepts file uploads, but the frontend has **no upload UI** — it only offers "browse" which opens a server-side file picker.

### 2. Hardcoded default paths that don't exist
- Defaults are `test_data/target` and `test_data/background` — these don't exist for a new user.
- Running with defaults fails with a confusing error ("Could not find valid primers for any length") instead of a clear message like "Selected target directory is empty or contains no FASTA files."

### 3. No onboarding or welcome screen
- Landing page immediately redirects to an empty dashboard with no explanation.
- No "Getting Started" guide, tooltips, or workflow wizard.
- The `About` page has useful info but is hidden among the navigation tabs.

### 4. NCBI email field is unexplained
- The form requires an NCBI email but never says *why*.
- No link to NCBI's policy or guidance on what email to use.

---

## Usability Issues

### 5. Terminal logs are developer-oriented
- Job progress is shown as raw CLI-style logs in a `<pre>` block.
- A biologist wants to see: *"Step 2/4: Designing primers... 45% complete"*, not raw text like `[Phase 1] Mining Target Pairs (1 strains) with Step=3...`
- No ETA or progress percentage anywhere.

### 6. AI chat shows raw JSON proposals
- When the AI proposes a configuration, it shows raw JSON in a `<pre>` block.
- A wetlab scientist cannot read JSON — they need a visual summary or auto-filled form.

### 7. Validation page needs bulk primer import
- Primers must be typed one pair at a time.
- No CSV/Excel import for bulk validation of hundreds of primers.

### 8. No template or saved configurations
- All parameters, paths, and queries must be re-entered each session.
- No way to save a "SARS-CoV-2 assay design" template for reuse.

### 9. Design forms don't show live parameter values
- Core parameters (Tm, GC%, sensitivity) are in the shared sidebar, not on the form.
- A user may change a parameter on the sidebar and forget they changed it.

---

## UI Polish Issues

### 10. Dashboard shows no results for failed jobs
- A failed job shows empty tables with no error explanation.
- Need to click the job then scroll to the terminal to see what went wrong.

### 11. No result comparison across jobs
- Cannot compare primer sets from two different runs side by side.

### 12. Primers displayed as raw tables with no visualization
- No genome browser, primer alignment map, or melting temperature plot.
- Tables with 8+ columns are hard to read at a glance.

### 13. Multiplex mode selector is confusing
- Three modes (local / auto / analyze) with no explanation of when to use each.

### 14. Download options limited to ZIP
- Can't download an individual CSV result table.
- Cannot copy a formatted results table to clipboard.

---

## Technical Debt

### 15. Monolithic 2400+ line Dashboard component
- Single file holds all views, state, API calls, and sub-components.
- Makes future UI improvements difficult and risky.

### 16. State bleed between views
- Changing a run name on one page affects all other pages.
- Language toggle only partially translated (dataset names stay in English).

### 17. Status panel shows stale data on initial load
- On first visit, all status fields show "-" or "0 B" until the user clicks refresh.
