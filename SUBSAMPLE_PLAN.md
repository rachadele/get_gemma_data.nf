# Per-study cell subsampling plan

## Problem

Current state: the eval pipeline subsamples **100 cells per sample**. Because studies vary widely in sample count (Ling-2024 has hundreds of samples, GSE180670 has a handful), aggregate metrics are dominated by high-powered studies — e.g. Ling-2024 contributes ~19k cells to evaluation while GSE180670 contributes ~400.

Caveat: downstream modeling operates at the **sample level**, not the cell level. Cell subsampling alone cannot equalize sample counts across studies — some studies only have ~10 samples and there is no way to manufacture more. So:

- **Cell-level subsampling** (this plan) fixes per-sample noise variance: each sample contributes roughly comparable amounts of data to its own per-sample estimate.
- **Sample-count imbalance** must be handled at the reporting stage by macro-averaging per-sample metrics within each study before averaging across studies (rather than pooling).

## Scope

Implement **per-study cell subsampling** inside `get_gemma_data.nf`, before the eval pipeline consumes the data. Downstream modeling still needs **per-sample h5ads**, so the subsampled study h5ad must be split back out by `sample_id`.

## Parameters

- `params.subsample_n = 10000` — target cell count per study.
- Studies with fewer than `subsample_n` cells **pass through unchanged** (do not drop underpowered studies).
- Random uniform sampling across the whole study, **not stratified by cell type / subclass** — preserves natural within-study composition.
- Fixed seed = 42 (consistent with eval pipeline `SEED`).

## Expected totals at N=10,000

**Human — 17 active studies:** ~169,651 cells (GSE211870 at 9,651 passes through; rest capped at 10,000).

**Mouse — 7 active studies** (GSE231868 intentionally excluded; GSE212068 relabel deleted as orphan): ~59,280 cells. GSE199460.2 (283 cells, Endothelial-only) and GSE181021.2 (8,997 cells) pass through; rest capped at 10,000.

## Architecture

Insert a new process after the per-study h5ad is built. Both `process_samples=true` and `process_samples=false` paths should converge to the same output contract.

### New process: `subsample_study`
- **Location:** `modules/processes/subsample_study.nf`
- **Input:** per-study h5ad (`{study_name}.h5ad`)
- **Output:** subsampled per-study h5ad (`{study_name}.subsampled.h5ad`), published to `${params.outdir}/h5ad_subsampled/`
- **Logic:**
  - Load h5ad.
  - If `n_obs <= subsample_n`: pass through (write the same cells out).
  - Else: random uniform sample of `subsample_n` cell indices with `np.random.default_rng(42)`.
  - Write subsampled h5ad.
  - Append a row to `subsampling_log.tsv` with `study_name, n_before, n_after, seed`.

### New process: `split_to_samples`
- **Location:** `modules/processes/split_to_samples.nf`
- **Input:** subsampled per-study h5ad
- **Output:** one h5ad per `sample_id` in `${params.outdir}/h5ad_sample_subsampled/{study_name}/`
- **Logic:**
  - Load subsampled study h5ad.
  - Group by `obs.sample_id`.
  - For each group with `n_obs >= 50` (existing `check_size` threshold), write `{study_name}_{sample_id}.h5ad`.
  - Smaller groups go to `small_samples/` (same convention as `process_query_samples.py`).

### Wiring in `main.nf`
- Append `subsample_study` after `PROCESS_QUERY_COMBINED` (and after the merge step in the `process_samples=true` branch — both paths feed `subsample_study` with a single per-study h5ad).
- `split_to_samples` consumes `subsample_study.out`.
- Both new processes publish to distinct directories; the existing `h5ad/` outputs remain so we can A/B compare if needed.

## Eval pipeline changes (separate PR, after `get_gemma_data.nf` lands)

In `nextflow_eval_pipeline`:
- Point query inputs at `h5ad_sample_subsampled/{study}/` instead of the existing per-sample h5ads.
- **Remove per-sample subsampling**, or change its default to `null` so it is off by default. Once the upstream subsampling is the source of truth, leaving the per-sample cap on would double-subsample and bias results.
  - Identify the param (likely something like `params.subsample_query` or similar) and the process that consumes it.
  - Default to `null`; keep the code path available behind an explicit flag for legacy runs if needed.

## Reproducibility

- `subsampling_log.tsv` published at `${params.outdir}/h5ad_subsampled/subsampling_log.tsv` with one row per study.
- Seed is hardcoded in the process script (not surfaced as a param) to ensure consistency across runs.

## Decided

- N = 10,000 cells per study (both organisms).
- Pass through studies smaller than N — do not drop.
- No subclass stratification — preserve natural distribution.
- GSE231868 intentionally excluded (developmental data, doesn't fit adult cortex benchmark).
- GSE212068 orphan relabel file already deleted.

## Open

- Sample-count imbalance is **not** addressed by this plan. Some studies have ~10 samples while others have hundreds; cell subsampling cannot change this. Partially mitigated downstream — sample-level metrics, macro-averaging across studies depending on the analysis, and modeling with study as a random effect — but residual concern remains. Sample-level subsampling (capping to M samples per study) is rejected because the floor is too low (~10) and would throw away most of the data for marginal bias reduction.
- GSE231868 still listed in `study_names_mouse.txt` despite being intentionally dropped — optional cleanup to stop wasting download time.
- Exact location of the existing per-sample subsampling param/process in the eval pipeline needs to be identified before that side can be implemented. Default should be flipped to `null` (off) once upstream subsampling is in place.

## Status

Not implemented yet. Waiting for current eval pipeline run to finish.
