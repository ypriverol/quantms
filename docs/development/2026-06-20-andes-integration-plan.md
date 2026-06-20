# Andes Search Engine Integration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `andes` as a selectable DDA search engine in quantms (emitting OpenMS `.idparquet`), validate it on Codon via Singularity against comet/msgf/sage, then flip it to the default.

**Architecture:** andes is a standalone Rust binary that emits an OpenMS-compatible QPX `.idparquet` (`--output-parquet`). A new `ANDES` Nextflow process drops into the existing engine-agnostic search subworkflow exactly where comet/sage sit; all downstream rescoring/FDR is unchanged (quantms-owned). andes-specific levers (`--chimeric`, `--refine`, `--score strong`) are exposed as opt-in, default-off params.

**Tech Stack:** Nextflow DSL2 (Groovy), nf-test, OpenMS `.idparquet` (QPX), Rust (andes build), Singularity/Apptainer, EMBL-EBI Codon SLURM.

## Global Constraints

- Branch: `andes-integration` (off `dev`). Milestone commits on this one branch; single closing PR.
- Public repo: commit messages and PRs MUST NOT contain Claude/AI/"superpowers" attribution.
- Default `search_engines` stays `'comet'` until Task 10 (gated on benchmark).
- andes reads per-run params from `meta` (SDRF-derived: `meta.enzyme`, `meta.fixedmodifications`, `meta.variablemodifications`, `meta.precursormasstolerance`, `meta.precursormasstoleranceunit`, `meta.fragmentmasstolerance`, `meta.fragmentmasstoleranceunit`) and pipeline-wide params from `params.*` — mirroring comet/sage. Never override mzML fragment tolerance with a forced value (andes auto-selects per-model; forcing it degrades its models).
- andes `mods.txt` line format: `<formula-or-mass>,<residue|*>,<fix|opt>,<location>,<name>`, preceded by `NumMods=<max_mods>`. Locations: `any`, `N-term`, `C-term`, `Prot-N-term`, `Prot-C-term`.
- New andes-only params default OFF: `andes_score='rank'`, `andes_chimeric=false`, `andes_refine=false`.
- andes binary name: `andes`. Search invocation: `andes --spectrum <mzml> --database <fasta> --output-parquet <out>.idparquet [flags]`.
- Codon paths: test data `/hps/nobackup/juan/pride/reanalysis/quantms-test-datasets`; downloads `/hps/nobackup/juan/pride/reanalysis/quantms-test-datasets/raw`; benchmark workspace `/hps/nobackup/juan/pride/reanalysis/quantms_benchmark`. Cluster access via the `codon-cluster` skill (pst_prd service account).

---

## File Structure

| Path | Responsibility |
|---|---|
| `modules/local/andes/main.nf` (create) | `ANDES` process: build mods.txt from meta, run andes, emit `*_andes.idparquet` |
| `modules/local/andes/tests/main.nf.test` (create) | nf-test for the module on BSA mzML + FASTA |
| `modules/local/andes/tests/main.nf.test.snap` (create) | nf-test snapshot |
| `modules/local/andes/tests/nextflow.config` (create) | per-test config (params the module reads) |
| `subworkflows/local/peptide_database_search/main.nf` (modify) | include + invoke ANDES; thread `ch_id_andes` through every mix point |
| `nextflow.config` (modify) | add `andes_score`, `andes_chimeric`, `andes_refine` defaults |
| `nextflow_schema.json` (modify) | add `andes` to `search_engines`; document new params |
| `conf/tests/test_andes.config` (create) | test profile running andes on BSA LFQ |
| `quantms_benchmark/` (create, workspace-level) | Codon profile, sif build script, dataset download, run driver, `figures.py` |

---

## Task 1: andes→OpenMS mods.txt mapping (Groovy helper, unit-validated by nf-test)

**Files:**
- Create: `modules/local/andes/main.nf` (helper closure `andesMods` + process skeleton)
- Test: `modules/local/andes/tests/main.nf.test` (asserts generated mods.txt content)

**Interfaces:**
- Produces: a Groovy closure that converts an OpenMS mod-name string (e.g. `"Oxidation (M)"`, `"Carbamidomethyl (C)"`, `"TMT6plex (K)"`, `"Acetyl (Protein N-term)"`) to an andes mods.txt line `<mass>,<residue|*>,<fix|opt>,<location>,<name>`.

The OpenMS mod string encodes residue/position in parentheses. Build a static name→(mono-mass, name) table for the mods quantms uses, and parse residue/location from the parenthetical.

- [ ] **Step 1: Write the failing test** — `modules/local/andes/tests/main.nf.test`

```groovy
nextflow_process {
    name "Test Process ANDES"
    script "../main.nf"
    process "ANDES"
    config "./nextflow.config"

    test("BSA mzML produces non-empty idparquet and correct mods.txt") {
        when {
            process {
                """
                input[0] = [
                    [ id:'BSA1_F1', mzml_id:'BSA1_F1', enzyme:'Trypsin',
                      fixedmodifications:'Carbamidomethyl (C)',
                      variablemodifications:'Oxidation (M)',
                      precursormasstolerance:10, precursormasstoleranceunit:'ppm',
                      fragmentmasstolerance:0.05, fragmentmasstoleranceunit:'Da' ],
                    file(params.test_mzml, checkIfExists: true),
                    file(params.test_fasta, checkIfExists: true)
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.id_files_andes.get(0).get(1)).exists()
            def mods = path("${process.out.get(0)}").parent.resolve('mods.txt')
            assert workflow.success
        }
    }
}
```

- [ ] **Step 2: Run test to verify it fails** — `nf-test test modules/local/andes/tests/main.nf.test` → FAIL ("ANDES" process not found / script missing).

- [ ] **Step 3: Write the helper + minimal process** in `modules/local/andes/main.nf`:

```groovy
// OpenMS Unimod name -> monoisotopic delta mass (Da)
def ANDES_MOD_MASS = [
    'Carbamidomethyl': 57.02146,
    'Oxidation'      : 15.99491,
    'Acetyl'         : 42.01057,
    'Phospho'        : 79.96633,
    'TMT6plex'       : 229.16293,
    'TMT10plex'      : 229.16293,
    'TMT16plex'      : 304.20715,
    'TMT18plex'      : 304.20715,
    'iTRAQ4plex'     : 144.10206,
    'iTRAQ8plex'     : 304.20536,
    'Deamidated'     : 0.98402,
]

// "Oxidation (M)" / "Acetyl (Protein N-term)" / "TMT6plex (N-term)" -> andes line
def andesModLine = { String openmsMod, String kind ->   // kind = 'fix' | 'opt'
    def m = (openmsMod =~ /^(.+?)\s*\((.+)\)\s*$/)
    assert m.matches() : "Unrecognised modification: ${openmsMod}"
    def name = m.group(1).trim()
    def target = m.group(2).trim()
    def mass = ANDES_MOD_MASS[name]
    assert mass != null : "No mass mapping for modification '${name}'. Add it to ANDES_MOD_MASS."
    def residue = '*'
    def location = 'any'
    switch (target.toLowerCase()) {
        case 'protein n-term': residue='*'; location='Prot-N-term'; break
        case 'protein c-term': residue='*'; location='Prot-C-term'; break
        case 'n-term':         residue='*'; location='N-term';      break
        case 'c-term':         residue='*'; location='C-term';      break
        default:               residue=target; location='any';      break  // single residue e.g. "M","C","K"
    }
    return "${mass},${residue},${kind},${location},${name}"
}
```

- [ ] **Step 4: Run test to verify it passes** — `nf-test test modules/local/andes/tests/main.nf.test` → PASS (after Task 2 fills the process body; if running Task 1 alone, assert the helper via a tiny Groovy `assert` harness instead).

- [ ] **Step 5: Commit**

```bash
git add modules/local/andes/main.nf modules/local/andes/tests/main.nf.test
git commit -m "feat(andes): mods.txt mapping helper for andes module"
```

---

## Task 2: ANDES process body

**Files:**
- Modify: `modules/local/andes/main.nf`
- Test: `modules/local/andes/tests/main.nf.test`, `modules/local/andes/tests/nextflow.config`

**Interfaces:**
- Consumes: `andesModLine` closure (Task 1).
- Produces: process `ANDES` with input `tuple val(meta), path(mzml_file), path(database)`; outputs `id_files_andes = tuple(meta, "${mzml_file.baseName}_andes.idparquet")`, `versions`, `log`.

- [ ] **Step 1: Write `tests/nextflow.config`** (params the module reads):

```groovy
params {
    allowed_missed_cleavages = 2
    min_peptide_length = 6
    max_peptide_length = 40
    num_hits = 1
    max_mods = 3
    min_precursor_charge = 2
    max_precursor_charge = 4
    num_enzyme_termini = 'fully'
    isotope_error_range = '0,1'
    andes_score = 'rank'
    andes_chimeric = false
    andes_refine = false
    test_mzml  = "${projectDir}/../../../../../testdata/lfq_ci/BSA/BSA1_F1.mzML"
    test_fasta = "${projectDir}/../../../../../testdata/lfq_ci/BSA/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta"
}
process {
    withName: 'ANDES' {
        container = 'andes:local'   // local sif/docker tag built in Task 7
    }
}
```

- [ ] **Step 2: Run test to verify it fails** — `nf-test test modules/local/andes/tests/main.nf.test` → FAIL (process body empty / no command).

- [ ] **Step 3: Write the process body** in `modules/local/andes/main.nf` (after the helper from Task 1):

```groovy
process ANDES {
    tag "$meta.mzml_id"
    label 'process_medium'
    label 'andes'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        ( params.andes_container ?: 'oras://ghcr.io/bigbio/andes-sif:0.1.0' ) :
        ( params.andes_container ?: 'ghcr.io/bigbio/andes:0.1.0' ) }"

    input:
    tuple val(meta), path(mzml_file), path(database)

    output:
    tuple val(meta), path("${mzml_file.baseName}_andes.idparquet"), emit: id_files_andes
    path "versions.yml", emit: versions
    path "*.log",        emit: log

    script:
    def args = task.ext.args ?: ''

    // enzyme name -> andes enzyme slug
    def enzymeMap = [ 'Trypsin':'trypsin', 'Trypsin/P':'trypsin', 'Arg-C':'argc',
                      'Asp-N':'aspn', 'Chymotrypsin':'chymotrypsin', 'Lys-C':'lysc',
                      'Lys-N':'lysn', 'Glu-C':'gluc', 'unspecific cleavage':'nonspecific' ]
    def andesEnzyme = enzymeMap[meta.enzyme] ?: 'trypsin'

    def ntt = (meta.enzyme == 'unspecific cleavage') ? 'non-specific' :
              (params.num_enzyme_termini == 'fully') ? 'fully' : 'semi'

    def iso = params.isotope_error_range.split(',')
    def isoMin = iso[0].trim(); def isoMax = iso[1].trim()

    // precursor tolerance: ppm vs Da
    def precFlag = (meta.precursormasstoleranceunit == 'ppm') ?
        "--precursor-tol-ppm ${meta.precursormasstolerance}" :
        "--precursor-tol-da ${meta.precursormasstolerance}"

    // build mods.txt from meta (fixed + variable)
    def fixedLines = meta.fixedmodifications?.trim() ?
        meta.fixedmodifications.tokenize(',').collect { andesModLine(it.trim(), 'fix') } : []
    def varLines = meta.variablemodifications?.trim() ?
        meta.variablemodifications.tokenize(',').collect { andesModLine(it.trim(), 'opt') } : []
    def modsContent = (["NumMods=${params.max_mods}"] + fixedLines + varLines).join('\n')

    def scoreFlag    = params.andes_score == 'strong' ? '--score strong' : '--score rank'
    def chimericFlag = params.andes_chimeric ? '--chimeric' : ''
    def refineFlag   = params.andes_refine   ? '--refine'   : ''

    """
    cat > mods.txt <<'EOF'
    ${modsContent}
    EOF

    andes \\
        --spectrum ${mzml_file} \\
        --database "${database}" \\
        --output-parquet ${mzml_file.baseName}_andes.idparquet \\
        --threads $task.cpus \\
        --enzyme ${andesEnzyme} \\
        --enzyme-specificity ${ntt} \\
        --max-missed-cleavages $params.allowed_missed_cleavages \\
        --min-length $params.min_peptide_length \\
        --max-length $params.max_peptide_length \\
        --top-n $params.num_hits \\
        --max-mods $params.max_mods \\
        --charge-min $params.min_precursor_charge \\
        --charge-max $params.max_precursor_charge \\
        --isotope-error-min ${isoMin} \\
        --isotope-error-max ${isoMax} \\
        ${precFlag} \\
        --mods mods.txt \\
        ${scoreFlag} ${chimericFlag} ${refineFlag} \\
        $args \\
        2>&1 | tee ${mzml_file.baseName}_andes.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        andes: \$(andes --version 2>&1 | sed 's/andes //g')
    END_VERSIONS
    """
}
```

- [ ] **Step 4: Run test (needs the local container from Task 7)** — `nf-test test modules/local/andes/tests/main.nf.test` → PASS once `andes:local` exists. Assert `*_andes.idparquet` exists and is a non-empty directory.

- [ ] **Step 5: Snapshot + commit**

```bash
git add modules/local/andes/
git commit -m "feat(andes): ANDES process emitting OpenMS idparquet"
```

---

## Task 3: New params in nextflow.config

**Files:**
- Modify: `nextflow.config` (search-engine params block, near `search_engines = 'comet'`)

- [ ] **Step 1: Add the three andes params + container override** after the `search_engines`/`sage_processes` lines:

```groovy
    // andes search engine (opt-in levers; default off)
    andes_score              = 'rank'   // 'rank' | 'strong'
    andes_chimeric           = false    // two-pass co-isolated peptide cascade
    andes_refine             = false    // PTM-discovery pass-2
    andes_container          = null     // override image/local sif path for dev
```

- [ ] **Step 2: Verify config parses** — `nextflow config -profile test_lfq . > /dev/null && echo OK` → prints `OK`.

- [ ] **Step 3: Commit**

```bash
git add nextflow.config
git commit -m "feat(andes): add andes_score/chimeric/refine params"
```

---

## Task 4: Schema update for `search_engines` + new params

**Files:**
- Modify: `nextflow_schema.json` (lines ~237–243 `search_engines`, and the search-engine params group)

- [ ] **Step 1: Edit the `search_engines` block** — change description/help/enum-ish text to include `andes`:

```json
                "search_engines": {
                    "type": "string",
                    "description": "A comma separated list of search engines to use (and combine). Valid: comet, msgf, sage, andes",
                    "default": "comet",
                    "fa_icon": "fas fa-tasks",
                    "help_text": "A comma-separated list of search engines to run in parallel on each mzML file. Currently supported: comet, msgf, sage and andes (default: comet). andes is a Rust engine emitting OpenMS idparquet; rescoring/FDR are handled by quantms as for the other engines."
                },
```

- [ ] **Step 2: Add the three andes params** to the same `properties` group (after `sage_processes`):

```json
                "andes_score": { "type": "string", "default": "rank", "enum": ["rank", "strong"], "description": "andes PIN RawScore source: rank (default) or strong (fused intensity+competition)." },
                "andes_chimeric": { "type": "boolean", "default": false, "description": "Enable andes two-pass co-isolated-peptide cascade (mzML/.raw only)." },
                "andes_refine": { "type": "boolean", "default": false, "description": "Enable andes PTM-discovery refinement pass." },
                "andes_container": { "type": "string", "description": "Override andes container image or local sif path (dev)." },
```

- [ ] **Step 3: Validate schema** — `nextflow schema validate` is not standard; instead run `python -c "import json,sys; json.load(open('nextflow_schema.json'))" && echo JSON_OK` → `JSON_OK`.

- [ ] **Step 4: Commit**

```bash
git add nextflow_schema.json
git commit -m "feat(andes): register andes in search_engines schema"
```

---

## Task 5: Wire ANDES into the search subworkflow

**Files:**
- Modify: `subworkflows/local/peptide_database_search/main.nf`

**Interfaces:**
- Consumes: `ANDES.out.id_files_andes` (tuple meta, idparquet), `ANDES.out.versions`.
- Produces: `ch_id_andes` threaded into every place the existing engines are mixed.

- [ ] **Step 1: Add the include** (after line 5, the SAGE include):

```groovy
include { ANDES } from '../../../modules/local/andes/main'
```

- [ ] **Step 2: Add `ch_id_andes` to the init tuple** (line 18):

```groovy
    (ch_id_msgf, ch_id_comet, ch_id_sage, ch_id_andes, ch_versions) = [ channel.empty(), channel.empty(), channel.empty(), channel.empty(), channel.empty() ]
```

- [ ] **Step 3: Add the ANDES invocation** (after the comet `if` block, ~line 33):

```groovy
    if (params.search_engines.contains("andes")) {
        ANDES(ch_mzmls_search.combine(ch_searchengine_in_db))
        ch_versions = ch_versions.mix(ANDES.out.versions)
        ch_id_andes = ch_id_andes.mix(ANDES.out.id_files_andes)
    }
```

- [ ] **Step 4: Thread `ch_id_andes` into every mix point.** Replace each `ch_id_msgf.mix(ch_id_comet).mix(ch_id_sage)` occurrence (lines 115, 119, 130, 134, 152, 155, 159, 164) with `ch_id_msgf.mix(ch_id_comet).mix(ch_id_sage).mix(ch_id_andes)`. Also add andes to the fine-tuning engine pool (after line 97):

```groovy
                    if (params.search_engines.contains("andes")) engine_opts.add("andes")
```
and the selector ternary (line 100–102):

```groovy
                    ch_selected_engine = (selected_engine == "sage")  ? ch_id_sage :
                                        (selected_engine == "msgf")  ? ch_id_msgf :
                                        (selected_engine == "andes") ? ch_id_andes :
                                        ch_id_comet
```

- [ ] **Step 5: Run the existing search-subworkflow nf-test to confirm no regression** (comet path untouched):

Run: `nf-test test tests/default.nf.test` (or the LFQ pipeline test). Expected: existing comet/sage tests still PASS.

- [ ] **Step 6: Commit**

```bash
git add subworkflows/local/peptide_database_search/main.nf
git commit -m "feat(andes): wire ANDES into peptide_database_search subworkflow"
```

---

## Task 6: test_andes.config profile

**Files:**
- Create: `conf/tests/test_andes.config`
- Modify: `nextflow.config` (register the profile in the `profiles` block alongside other `test_*`)

- [ ] **Step 1: Create `conf/tests/test_andes.config`** (mirror `test_lfq.config`, BSA data, andes engine):

```groovy
process {
    resourceLimits = [ cpus: 4, memory: '12.GB', time: '48.h' ]
}
params {
    config_profile_name        = 'Test profile for DDA LFQ with andes'
    config_profile_description = 'Minimal BSA test exercising the andes search engine'
    outdir = "./results_andes"
    input    = 'https://raw.githubusercontent.com/bigbio/quantms-test-datasets/refs/heads/quantms/testdata/lfq_ci/BSA/BSA_design.sdrf.tsv'
    database = 'https://raw.githubusercontent.com/bigbio/quantms-test-datasets/quantms/testdata/lfq_ci/BSA/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta'
    search_engines = "andes"
    decoy_string = "rev"
    protein_level_fdr_cutoff = 1.0
    psm_level_fdr_cutoff = 1.0
    quantify_decoys = true
    mzml_features = true
}
```

- [ ] **Step 2: Register the profile** in `nextflow.config` `profiles {}` (next to `test_lfq`):

```groovy
    test_andes      { includeConfig 'conf/tests/test_andes.config' }
```

- [ ] **Step 3: Verify it resolves** — `nextflow config -profile test_andes,docker . > /dev/null && echo OK` → `OK`.

- [ ] **Step 4: Commit**

```bash
git add conf/tests/test_andes.config nextflow.config
git commit -m "test(andes): add test_andes profile on BSA LFQ"
```

---

## Task 7: Build the local andes container (Codon)

**Files:**
- Create: `quantms_benchmark/containers/build_andes_sif.sh`

**Deliverable:** a `andes_<ver>.sif` (and/or `docker:andes:local`) usable by the module via `--andes_container`.

- [ ] **Step 1: Write `quantms_benchmark/containers/build_andes_sif.sh`** (run on a Codon compute node via codon-cluster skill):

```bash
#!/usr/bin/env bash
set -euo pipefail
# 1. Build andes release binary from the msgf-rust checkout
SRC=${1:?path to msgf-rust checkout}
OUT=${2:?output dir on /hps}
( cd "$SRC" && RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features thermo )
# 2. Stage binary + resources into a rootfs
STAGE="$OUT/andes_root"; rm -rf "$STAGE"; mkdir -p "$STAGE/opt/andes"
cp "$SRC/target/release/andes" "$STAGE/opt/andes/andes"
cp -r "$SRC/resources" "$STAGE/opt/andes/resources"
# 3. Build a Singularity image from a minimal def
cat > "$OUT/andes.def" <<'DEF'
Bootstrap: docker
From: debian:bookworm-slim
%files
    ANDES_ROOT/opt/andes /opt/andes
%post
    ln -s /opt/andes/andes /usr/local/bin/andes
%environment
    export PATH=/opt/andes:$PATH
DEF
sed -i "s#ANDES_ROOT#$STAGE#" "$OUT/andes.def"
singularity build --fakeroot "$OUT/andes_local.sif" "$OUT/andes.def"
echo "Built $OUT/andes_local.sif"
```

- [ ] **Step 2: Run it on Codon** (codon-cluster skill, compute node):

Run: `bash quantms_benchmark/containers/build_andes_sif.sh /hps/.../msgf-rust /hps/nobackup/juan/pride/reanalysis/quantms_benchmark/containers`
Expected: `Built …/andes_local.sif`, and `singularity run …/andes_local.sif andes --version` prints `andes 0.1.0`.

- [ ] **Step 3: Smoke-run the module on BSA** with the local sif:

Run: `nextflow run . -profile test_andes,singularity --andes_container '/hps/.../andes_local.sif' --outdir results_andes_smoke -resume`
Expected: pipeline completes; `results_andes_smoke` contains percolator-filtered PSMs; andes log shows non-zero rank-1 PSMs.

- [ ] **Step 4: Commit the build script**

```bash
git add quantms_benchmark/containers/build_andes_sif.sh
git commit -m "build(andes): local Singularity image build script"
```

---

## Task 8: Benchmark harness — Codon profile, dataset staging, run driver, figures

**Files:**
- Create: `quantms_benchmark/conf/codon_singularity.config` (SLURM executor + singularity + resource limits)
- Create: `quantms_benchmark/datasets/download.sh` (fetch UPS1/Astral raw → convert to mzML → `…/quantms-test-datasets/raw`)
- Create: `quantms_benchmark/run/run_benchmark.sh` (loop over {dataset × engine × variant}, launch nextflow)
- Create: `quantms_benchmark/figures/figures.py` (parse pmultiqc/PSM tables → barplots + FDP curves)
- Create: `quantms_benchmark/README.md`

- [ ] **Step 1: Write `codon_singularity.config`** — SLURM process executor, `singularity.enabled = true`, `singularity.autoMounts = true`, queue/account for pst_prd, resourceLimits matching the cluster.

- [ ] **Step 2: Write `datasets/download.sh`** — download the chosen UPS1 + Astral PRIDE raw files into `/hps/nobackup/juan/pride/reanalysis/quantms-test-datasets/raw`, convert `.raw`→`.mzML` (ThermoRawFileParser container) if needed, and write a small SDRF for each. (Accessions filled in when the user confirms them.)

- [ ] **Step 3: Write `run/run_benchmark.sh`** — for each dataset, run quantms with `-profile <dataset>,singularity -c quantms_benchmark/conf/codon_singularity.config`, sweeping `--search_engines` over `comet`, `msgf`, `sage`, `andes`, and for andes additionally toggling `--andes_chimeric`, `--andes_refine`, `--andes_score strong`, and `--skip_rescoring` true/false. Tag each run's `--outdir` uniquely.

- [ ] **Step 4: Write `figures/figures.py`** — collect per-run PSM/peptide/protein counts @1% FDR + wall time from the pmultiqc/openms output, emit grouped barplots (engine on x-axis, dataset facets) and an entrapment-FDP curve where the FASTA supports it; save PNGs under `quantms_benchmark/figures/out/`.

- [ ] **Step 5: Commit the harness**

```bash
git add quantms_benchmark/
git commit -m "bench(andes): Codon benchmark harness (profile, datasets, runner, figures)"
```

> Note: Steps 2–4 are operational scripts; their "test" is a successful dry-run on the BSA dataset before the real UPS1/Astral runs.

---

## Task 9 (ops milestone): Phase 0 baseline + Phase 2 andes sweep on Codon

This task produces results/figures, not code. Deliverable = committed figures + a short results note.

- [ ] **Step 1: Phase 0 baseline** — run current quantms (comet/msgf/sage) on BSA LFQ, ProteoBench HYE DDA, PXD000001 TMT via `run_benchmark.sh` on Codon. Produce baseline figures. This validates the harness and teaches the pipeline.
- [ ] **Step 2: Stage UPS1/Astral** — run `datasets/download.sh` (accessions confirmed with user) into the Codon `raw/` folder.
- [ ] **Step 3: andes sweep** — run the full {engine × andes-variant × rescoring} matrix on UPS1/Astral.
- [ ] **Step 4: Figures + note** — generate comparison figures; write `quantms_benchmark/RESULTS.md` summarising PSMs/peptides/proteins @1% FDR (+ FDP) and wall time.
- [ ] **Step 5: Commit** results note + figures.

```bash
git add quantms_benchmark/RESULTS.md quantms_benchmark/figures/out/
git commit -m "bench(andes): baseline + andes sweep results on UPS1/Astral"
```

---

## Task 10 (gated): publish container + flip default

**Precondition:** Task 9 shows andes winning (or matching at better FDP) on the benchmark.

- [ ] **Step 1: Add `andes/Dockerfile` to the `quantms-containers` repo** (`/Users/yperez/work/quantms-workspace/quantms-containers`), modelled on the existing relink/diann Dockerfiles: build the andes release binary, copy binary + `resources/`, set `ENTRYPOINT ["andes"]`. Wire it into `.github/workflows/quantms-containers.yml`. Tag `ghcr.io/bigbio/andes:0.1.0` + `oras://ghcr.io/bigbio/andes-sif:0.1.0`.
- [ ] **Step 2: Point the module at the published image** — confirm `modules/local/andes/main.nf` default container tags match the published version; drop the `--andes_container` override from benchmark runs.
- [ ] **Step 3: Flip the default** — in `nextflow.config` set `search_engines = 'andes'`; update the schema `default` to `'andes'`.
- [ ] **Step 4: Run the full test matrix** — `nf-test test` (module) + `nextflow run -profile test_andes,docker` + a `test_lfq` run to confirm the default change is green.
- [ ] **Step 5: Commit**

```bash
git add nextflow.config nextflow_schema.json modules/local/andes/main.nf
git commit -m "feat(andes): make andes the default search engine"
```

---

## Self-Review notes

- **Spec coverage:** Section 3 (module+wiring+schema+params) → Tasks 1–6; Section 4 (containerization) → Tasks 7 & 10; Section 5 (benchmark harness + datasets) → Tasks 8 & 9; Section 6 phasing → Task ordering; Section 7 (testing) → Tasks 2,5,6. All covered.
- **mods mapping** is the one genuinely novel unit and is isolated in Task 1 with an explicit "add to ANDES_MOD_MASS" failure path for unmapped mods.
- **Open inputs to confirm at execution time:** UPS1/Astral PRIDE accessions (Task 8/9); the published andes version tag (Task 10); the exact Codon SLURM queue/account flags for `codon_singularity.config` (Task 8, via codon-cluster skill).
