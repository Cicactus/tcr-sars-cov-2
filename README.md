# covid-tcr-markers

Comparison of published computational methods for finding SARS-CoV-2-associated
TCR clonotypes ("markers") in T-cell repertoire sequencing data. Master's thesis
project.

## Data

[Zenodo 10.5281/zenodo.17780804](https://zenodo.org/records/17780804) (v5, open
access), Genomics of Adaptive Immunity Laboratory — "Robust detection of
SARS-CoV-2 exposure in population using T-cell repertoire profiling". TCRβ
repertoires, 9 sequencing batches, COVID-positive vs healthy donors.

## Methods compared

| Method | Tool | Tests |
| --- | --- | --- |
| Exact Fisher's test | `vdjtools.biomarker.association` | individual clonotype |
| ALICE | `vdjtools.overlap.alice` | individual clonotype |
| TCRnet | `vdjtools.overlap.tcrnet` | individual clonotype |
| tcrdist3 meta-clonotypes | `tcrdist.repertoire.TCRrep` | individual clonotype |
| REDCEA | official `antigenomics/redcea` package | cluster |

## Repository layout

```
scripts/
  01_prepare_data.py   run ONCE: download, batch-correct, select shared candidates
  02_run_fisher.py     \
  03_run_alice.py       } each reads scripts/01's output, no re-download/re-correction
  04_run_tcrnet.py      }
  05_run_tcrdist3.py   /
  run_all.sh            orchestrates 01-05 on one machine (venv, checkpointed, idempotent)
redcea/
  prepare_redcea_input.py   thin adapter: shared prepared data -> REDCEA's TSV schema
  run_redcea.sh              REDCEA needs conda, so it runs separately (own machine/env)
notebooks/
  01_batch_effect.ipynb        diagnoses + corrects the sequencing-batch effect
  02_additional_figures.ipynb   the comparison figures for the thesis
```

## The core design decision: one data-prep step, everyone reads from it

Every method needs the same three things: the batch-corrected case/control
repertoires, and a list of candidate clonotypes to test. Originally each
method's script downloaded and batch-corrected the cohort itself — wasteful,
and each one ended up selecting a *slightly* different candidate set (their
own internal incidence filters disagreed with each other), which silently
distorted the FDR comparison between methods.

`01_prepare_data.py` now does this exactly once and writes it to disk:

- `case_repertoire_corrected.tsv` / `control_repertoire_corrected.tsv` — the
  full train cohort, batch-corrected.
- `test_case_repertoire_corrected.tsv` / `test_control_repertoire_corrected.tsv`
  — the held-out donors, run through the *same* fitted correction (fit-on-train,
  transform-on-test — the correction itself never sees the test donors).
- `candidates.tsv` — clonotypes seen in ≥3 case donors (train only). This is
  the single, shared list of hypotheses every method is evaluated against.

Scripts `02`-`05` (and `redcea/prepare_redcea_input.py`) only *read* these
files. Re-running a method is fast; re-running `01_prepare_data.py` is the
only step that touches the network or the raw per-donor files.

### Why the shared candidate list is applied differently per method

- **Fisher** takes the list through `vdjtools.biomarker.association`'s native
  `candidates=` parameter: the full cohort is still used for the incidence
  statistics, only the set of *tested* hypotheses is restricted — which is
  exactly what that parameter is for.
- **ALICE, TCRnet, tcrdist3** have no equivalent parameter. Restricting their
  *input* to the candidate list ahead of time would distort neighbour
  counting (a candidate's neighbour that itself falls short of the donor
  threshold would silently stop being counted, deflating the signal). So
  these three run on the **full** batch-corrected repertoire, and the
  candidate restriction + a fresh Benjamini–Hochberg correction are applied
  **afterwards**, to the raw p-values each method reports. The result column
  `fdr_shared` in every method's output file is therefore computed over the
  exact same number of hypotheses across all four methods.
- **REDCEA** tests clusters, not individual clonotypes — a genuinely
  different unit of hypothesis, kept as a separate comparison axis rather
  than forced into the same FDR pool (see `notebooks/02_final_report_figures.ipynb`).


## Running

```bash
# on the main machine (Fisher / ALICE / TCRnet / tcrdist3):
scp -r scripts/ user@server:~/
ssh user@server
cd scripts && ./run_all.sh

# on the REDCEA machine (needs conda; can be the same server or a separate one):
scp -r redcea/ user@redcea-server:~/
rsync -avz user@server:~/covid_tcr_markers_run/prepared/ ~/prepared_data/
ssh user@redcea-server
PREPARED_DIR=~/prepared_data ./redcea/run_redcea.sh
```

Then run `notebooks/01_batch_effect.ipynb` and
`notebooks/02_final_report_figures.ipynb` locally, pointing them at the
downloaded `prepared/` and `results/` directories.