# Design: Integrate `andes` as a search engine in quantms (→ default)

**Date:** 2026-06-20
**Branch:** `andes-integration` (from `dev`)
**Status:** Approved design — ready for implementation plan

## 1. Goal

Add `andes` (the pure-Rust, data-driven peptide search engine; formerly msgf-rust/cimas)
as a selectable DDA search engine in quantms, validate it on the EMBL-EBI Codon cluster
(Singularity) against the existing engines, then make it the **default** search engine.

### Locked-in decisions

- **quantms owns rescoring/FDR (engine-agnostic).** andes emits plain OpenMS `.idparquet`;
  the existing downstream (`PSM_CLEAN → MSRESCORE_FEATURES → PERCOLATOR → PSM_FDR_CONTROL`)
  handles rescoring exactly as it does for comet/msgf/sage. andes does **not** run its own
  `--rescore` inside the pipeline.
- **Default flip is gated** on benchmark results — it lands in a later commit on the same branch,
  not as part of the initial wiring.
- **Benchmark scope = baseline + andes sweep.** First establish a baseline by running current
  quantms (comet/msgf/sage) on small datasets to learn the pipeline and produce figures, then
  run the andes sweep head-to-head.

## 2. Why this integration is clean

quantms is **idparquet-native**: every search engine emits an OpenMS QPX `.idparquet` bundle and
the entire downstream is engine-agnostic on that format. andes already emits a byte-compatible
QPX `.idparquet` via `--output-parquet`. So andes occupies the exact slot comet/sage occupy with
**no new downstream code**. andes's own levers (`--chimeric`, `--refine`, `--score strong`) become
*search-time* parameters; "with/without rescoring" is just toggling quantms's
`skip_rescoring` / `ms2features_enable`, which keeps cross-engine comparison fair.

## 3. Pipeline integration

### 3.1 New module — `modules/local/andes/main.nf`
- Process `ANDES`, mirroring `modules/local/openms/comet/main.nf` in shape.
- Inputs: `tuple val(meta), path(mzml_file), path(database)`.
- Command: `andes --spectrum <mzml> --database <fasta> --output-parquet <out>_andes.idparquet`
  plus mapped search params (see 3.4). Emit a `*.log` and `versions.yml`.
- Output: `id_files_andes` = `${mzml_file.baseName}_andes.idparquet`.
- Lives under `modules/local/andes/` (not `openms/`) because andes is a standalone binary, not
  an OpenMS adapter.

### 3.2 Subworkflow wiring — `subworkflows/local/peptide_database_search/main.nf`
- Add, alongside the existing comet/msgf/sage `if` blocks:
  ```groovy
  if (params.search_engines.contains("andes")) {
      ANDES(ch_mzmls_search.combine(ch_searchengine_in_db))
      ch_versions = ch_versions.mix(ANDES.out.versions)
      ch_id_andes = ch_id_andes.mix(ANDES.out.id_files_andes)
  }
  ```
- Mix `ch_id_andes` into the merged `ch_id_files_out`. No other downstream change required.

### 3.3 Schema / config
- Add `andes` to the `search_engines` valid list and help text in `nextflow_schema.json`.
- Keep `search_engines = 'comet'` as the default **until Phase 4**.

### 3.4 Parameter mapping (reuse existing shared params)
| quantms param | andes flag |
|---|---|
| `precursor_mass_tolerance` (+ `_unit`) | `--precursor-tol-ppm` / `--precursor-tol-da` |
| `fragment_mass_tolerance` (+ `_unit`) | `--fragment-tol-ppm` / `--fragment-tol-da` (mzML uses model default; flags are MGF-only — see note) |
| `enzyme` (from SDRF/meta) | `--enzyme` |
| `num_enzyme_termini` (`fully`/`semi`/`none`) | `--enzyme-specificity` (`--ntt`) |
| `allowed_missed_cleavages` | `--max-missed-cleavages` |
| `variable_mods` / fixed mods | `--mods` (generated mods file) |
| `min_precursor_charge` / `max_precursor_charge` | `--charge-min` / `--charge-max` |
| `min_peptide_length` / `max_peptide_length` | `--min-length` / `--max-length` |
| `max_mods` | `--max-mods` |
| `isotope_error_range` | `--isotope-error-min` / `--isotope-error-max` |
| `num_hits` | `--top-n` |

New **andes-only** params, all **default-off** so they don't change baseline behaviour:
- `andes_score` — `rank` (default) | `strong`
- `andes_chimeric` — boolean (default false) → `--chimeric`
- `andes_refine` — boolean (default false) → `--refine`

**Note on fragment tolerance:** for mzML, andes auto-selects a per-model fragment tolerance;
the `--fragment-tol-*` flags are MGF-only. The mapping must not force a tolerance that overrides
the model default on mzML input (this previously broke andes's own models — see msgf-rust memory
`gap-corpus-execution`). Confirm behaviour during Phase 1 smoke test.

## 4. Containerization (phased)

### Phase 1 — local binary (now)
- Build on Codon: `cargo build --release -p andes --features thermo` (rust toolchain module).
- Wrap the built binary + `resources/` into a local Singularity `.sif`.
- Module gains an override param (`andes_container` / local image path) so we can iterate before
  any image is published.

### Phase 2 — released image (later)
- Add `andes/Dockerfile` to the existing **`quantms-containers`** repo
  (`/Users/yperez/work/quantms-workspace/quantms-containers`, which already builds diann/relink via
  `.github/workflows/quantms-containers.yml`).
- Publish `ghcr.io/bigbio/andes:<version>` + `oras://ghcr.io/bigbio/andes-sif:<version>`.
- Module switches to the standard dual-URI `ghcr.io/bigbio` pattern used by every other engine:
  ```groovy
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'oras://ghcr.io/bigbio/andes-sif:<ver>' : 'ghcr.io/bigbio/andes:<ver>' }"
  ```

## 5. Benchmark harness — `quantms_benchmark/`

A self-contained area holding the Nextflow runner (configs + a Codon Singularity profile), the
`.sif` containers, dataset references, results, and a small `figures.py` (matplotlib/pmultiqc).
Version-controlled scripts/configs live in the repo/workspace; heavy artifacts (`.sif`, raw data,
nf-work, results) live on Codon under `/hps/nobackup/juan/pride/reanalysis/`.

### Datasets (Codon)
- Test datasets are staged at `/hps/nobackup/juan/pride/reanalysis/quantms-test-datasets`.
- Raw spectra to be downloaded go into `/hps/nobackup/juan/pride/reanalysis/quantms-test-datasets/raw`.

### Phase 0 — baseline + learn the pipeline (in-repo CI data, fast)
Run current quantms (comet/msgf/sage) end-to-end on the small in-repo datasets to learn the pipeline
and produce the comparison baseline + figures:
- **BSA LFQ** (`lfq_ci/BSA`, 6 × ~5 MB mzML)
- **ProteoBench HYE DDA** (`dda_ci`, mixed-species quant benchmark)
- **PXD000001 TMT** (`tmt_ci`, Erwinia)

Figures: PSMs / peptides / proteins @ 1% FDR + wall time per engine.

### Phase 2 — andes sweep (UPS1 / Astral)
⚠️ UPS1/Astral DDA are **not** in `quantms-test-datasets` (only a DIA UPS1 FASTA exists). These are
the andes-campaign PRIDE datasets and must be downloaded into the Codon `raw/` folder above.
- Head-to-head: `{comet, msgf, sage, andes}`.
- andes variants: `{plain, chimeric, refine, score=strong}` × rescoring `{off, percolator}`.
- Metrics: PSMs/peptides/proteins @ 1% FDR, entrapment-FDP where the FASTA supports it, wall time.

### Execution
Via the **codon-cluster** skill: `pst_prd` service account, `/hps/nobackup` storage, `sbatch`,
Singularity. Use the quantms `singularity` profile.

## 6. Phasing (milestone commits on `andes-integration`, single closing PR)
0. Baseline runs + figures (understand the pipeline).
1. `ANDES` module + subworkflow wiring + param mapping + local binary; smoke-test on BSA.
2. andes sweep benchmarks on Codon → figures.
3. andes Dockerfile in `quantms-containers` → published image; module switches to it.
4. **Flip default** to `search_engines = 'andes'` — gated on Phase 2 showing a win at honest FDP.

## 7. Testing
- nf-test for the `ANDES` module (mirror `modules/local/openms/comet/tests`).
- `conf/tests/test_andes.config` profile running andes on BSA LFQ.
- All cluster validation runs through Singularity on Codon.

## 8. Out of scope (for this round)
- andes-internal `--rescore` path inside quantms (quantms owns rescoring).
- DIA / timsTOF benchmarking (DDA first).
- Multi-engine ConsensusID tuning involving andes.
- Training/retraining andes models (uses bundled `models.parquet`).
