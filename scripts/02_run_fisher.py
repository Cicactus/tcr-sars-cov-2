#!/usr/bin/env python3
"""
Method: exact Fisher's test (vdjtools.biomarker.association).

Reads the already batch-corrected repertoires and the shared candidate list
produced by 01_prepare_data.py. No downloading, no batch correction happens
here -- this script is fast and can be re-run freely.
"""
import argparse
from pathlib import Path

import pandas as pd
import polars as pl
from statsmodels.stats.multitest import multipletests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="output-dir of 01_prepare_data.py")
    parser.add_argument("--output", required=True)
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    # sample_id/batch must be read as strings -- some donor ids look purely
    # numeric ("02000020") while others are alphanumeric ("p17_17_DNA"); if
    # polars infers the column from the numeric-looking rows it sees first,
    # it picks int64 and then crashes on the first alphanumeric id it meets
    schema_overrides = {"sample_id": pl.Utf8, "batch": pl.Utf8}
    case_bc = pl.read_csv(str(input_dir / "case_repertoire_corrected.tsv"), separator="\t", schema_overrides=schema_overrides)
    control_bc = pl.read_csv(str(input_dir / "control_repertoire_corrected.tsv"), separator="\t", schema_overrides=schema_overrides)
    candidates_pl = pl.read_csv(str(input_dir / "candidates.tsv"), separator="\t")
    metadata = pd.read_csv(input_dir / "metadata_split.csv")

    from vdjtools import biomarker

    cohort_train = pl.concat([case_bc, control_bc]).with_columns(pl.col("sample_id").cast(pl.Utf8))
    case_ids = case_bc["sample_id"].unique().to_list()
    control_ids = control_bc["sample_id"].unique().to_list()
    train_meta = metadata[metadata["id"].isin(case_ids + control_ids)][["id", "COVID_status"]]
    meta_pl = pl.from_pandas(train_meta.rename(columns={"id": "sample_id"})).with_columns(
        pl.col("sample_id").cast(pl.Utf8),
        (pl.col("COVID_status") == "COVID").alias("is_covid"),
    )

    # candidates= restricts which hypotheses are tested while the full cohort
    # is still used for the incidence statistics themselves -- see README
    # "Why the shared candidate list is applied differently per method".
    result = biomarker.association(
        cohort_train, meta_pl,
        pheno_col="is_covid", test=["fisher"],
        candidates=candidates_pl.select(["junction_aa", "v_call", "j_call"]),
        min_incidence=1, min_incidence_frac=0.0,
        threads=args.cpus,
    )
    result_pd = result.to_pandas().sort_values("p_value")
    result_pd["fdr_shared"] = multipletests(result_pd["p_value"], method="fdr_bh")[1]
    result_pd.to_csv(args.output, sep="\t", index=False)

    n_sig = int((result_pd["fdr_shared"] < 0.05).sum())
    print(f"Fisher: {len(result_pd)} candidates tested, {n_sig} significant (shared FDR < 0.05)")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
