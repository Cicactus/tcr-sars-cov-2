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