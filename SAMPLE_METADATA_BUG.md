# Sample metadata ID bug — fixed in `bin/get_gemma_meta.py`, but existing downloads need reprocessing

## The bug

`bin/get_gemma_meta.py` used to build each study's sample metadata table from
**two separate API calls**, then zipped the results together by list
position, assuming both calls returned samples in the same order:

```python
samples_raw = client.raw.get_dataset_samples(study_name)                          # call #1
samples = client.get_dataset_samples(study_name, use_processed_quantitation_type=False)  # call #2
sample_names = [x for x in samples["sample_name"]]
sample_ids = [x.id for x in samples_raw.data]
sample_meta = [df for df in samples["sample_characteristics"]]
sample_meta_updated1 = [df.assign(sample_id=sample_ids[i]) for i, df in enumerate(sample_meta)]        # zipped by position i
sample_meta_updated2 = [df.assign(sample_name=sample_names[i]) for i, df in enumerate(sample_meta_updated1)]  # zipped by position i
```

That assumption is false. Verified live on the CMC study:

```python
ids_from_raw = [x.id for x in samples_raw.data]                # BioAssay IDs (call #1)
ids_from_processed = samples['sample_ID'].tolist()               # call #2's own "sample_ID" column
set(ids_from_raw) & set(ids_from_processed)   # -> set() -- ZERO overlap
```

`call #1`'s IDs are **BioAssay IDs** (the same IDs used everywhere else —
filenames, MEX export, cell type assignments). `call #2`'s `sample_ID`
column is actually the nested **BioMaterial ID** (`BioAssay.sample.id`,
confirmed directly) — a different Gemma entity/table entirely. Two
unrelated ID spaces, zipped by position, with rows that don't even
correspond to the same sample at the same index. Concrete proof: the
positional pairing assigns BioAssay `1203011` the name `CMC_MSSM_035`; the
real name (read directly off `1203011`'s own nested BioMaterial) is
`CMC_MSSM_192`.

## The fix

`bin/get_gemma_meta.py` now builds everything from `call #1` alone —
`sample_id`, `sample_name`, and every characteristic all come off the same
`BioAssay` object (`s.id`, `s.name`, `s.sample.characteristics`), so there's
nothing left to misalign:

```python
samples_raw = client.raw.get_dataset_samples(study_name)
rows = []
for s in samples_raw.data:
    row = {
        "sample_id": s.id,
        "sample_name": s.name,
        "organism": s.array_design.taxon.scientific_name.lower().replace(" ", "_"),
    }
    for c in s.sample.characteristics:
        row[c.category] = c.value
    rows.append(row)
sample_meta_df = pd.DataFrame(rows)
```

## Why this corrupts downstream data

`process_query_samples.py`'s `add_sample_meta()` merges this sample-metadata
table into each h5ad's `obs` **by `sample_id`** (a real, correct join key —
that part is fine). But if the metadata table itself has the wrong
name/characteristics attached to a given `sample_id` (the bug above), that
wrong data gets faithfully merged into the h5ad. So for affected studies:

- `sample_id` in the h5ad/filename is correct (always was).
- `sample_name`, disease/diagnosis, age, sex, PMI, and every other
  characteristic baked into `obs` for that sample can be **wrong** —
  silently attached from a different, unrelated sample.
- The per-sample MEX/h5ad filenames (built independently, directly from
  `BioAssay.name` at export time) are **not** affected by this bug and can
  be trusted.

## What still needs to happen

Every study already downloaded through this pipeline has
`metadata/*_sample_meta.tsv` and `metadata_standardized/*_sample_meta_std.tsv`
files built by the **old, buggy** script, and every per-sample h5ad's `obs`
(everything except `sample_id` itself) was merged in from those same buggy
files.

1. Re-run the pipeline with `-resume` for every already-downloaded study
   (both this directory's study list and the mouse one). Only
   `bin/get_gemma_meta.py` changed, so Nextflow's cache should correctly
   re-run just `getGemmaMeta` → `standardizeMetadata` → `PROCESS_QUERY_SAMPLE`
   (or `PROCESS_QUERY_COMBINED`, whichever mode was used) for every study,
   while `DOWNLOAD_STUDIES`/`downloadCelltypes`/MEX download stay cached
   (unaffected by this bug).
2. Don't assume which studies were actually affected — in the human set,
   CMC/HBCC_Cohort/Ling-2024/SZBDMulti-Seq/Batiuk-2022 showed real
   mismatches (up to 100% of samples) but MultiomeBrain/GSE254569 showed
   zero; **the mouse studies haven't been checked at all**. Verify each
   study rather than assuming.
3. Anything downstream that was built from the old (wrong) metadata needs
   regenerating too once this is fixed — e.g.
   `aim2-reannotation/meta/diagnosis_normalized/*.tsv`.

## How to verify a given study's old metadata was actually wrong

```python
import gemmapy, os
client = gemmapy.GemmaPy(auth=[os.environ['GEMMA_USERNAME'], os.environ['GEMMA_PASSWORD']], path='staging')
samples_raw = client.raw.get_dataset_samples(study_name)
samples = client.get_dataset_samples(study_name, use_processed_quantitation_type=False)
overlap = set(x.id for x in samples_raw.data) & set(samples['sample_ID'].tolist())
print(study_name, "affected:" , len(overlap) == 0)
```
