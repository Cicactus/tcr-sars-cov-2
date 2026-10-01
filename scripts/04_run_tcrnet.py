#!/usr/bin/env python3
"""
Method: TCRnet (vdjtools.overlap.tcrnet) -- neighbourhood enrichment against
an empirical control repertoire (Pogorelyy et al., 2019, PLoS Biol).

Reads the already batch-corrected case/control repertoires produced by
01_prepare_data.py. Runs on the FULL case cohort against the FULL control
background; the shared candidate list is applied AFTER, to the raw
p-values, with a fresh FDR correction restricted to that shared set -- see
README, "Why the shared candidate list is applied differently per method".
"""
import argparse
from pathlib import Path

import pandas as pd
import polars as pl
from statsmodels.stats.multitest import multipletests


def restrict_and_refdr(df: pd.DataFrame, candidate_keys: set, p_col: str) -> pd.DataFrame:
    keys = list(zip(df["junction_aa"], df["v_call"], df["j_call"]))
    mask = pd.Series(keys, index=df.index).isin(candidate_keys)
    restricted = df.loc[mask].copy()
    if len(restricted) == 0:
        print("[warning] empty intersection with the shared candidate list")
        return restricted
    restricted["fdr_shared"] = multipletests(restricted[p_col], method="fdr_bh")[1]
    return restricted.sort_values(p_col)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="output-dir of 01_prepare_data.py")
    parser.add_argument("--output", required=True)
    parser.add_argument("--cpus", type=int, default=1)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    schema_overrides = {"sample_id": pl.Utf8, "batch": pl.Utf8}
    case_bc = pl.read_csv(str(input_dir / "case_repertoire_corrected.tsv"), separator="\t", schema_overrides=schema_overrides)
    control_bc = pl.read_csv(str(input_dir / "control_repertoire_corrected.tsv"), separator="\t", schema_overrides=schema_overrides)
    candidates_pl = pl.read_csv(str(input_dir / "candidates.tsv"), separator="\t")
    candidate_keys = set(zip(
        candidates_pl["junction_aa"].to_list(),
        candidates_pl["v_call"].to_list(),
        candidates_pl["j_call"].to_list(),
    ))

    import seqtree
    from vdjtools import overlap

    control_cdr3s = control_bc["junction_aa"].unique().to_list()
    control_index = seqtree.Index.build(control_cdr3s, alphabet="aa")

    tcrnet_result = overlap.tcrnet(case_bc, control=control_index, threads=args.cpus).to_pandas()
    result = restrict_and_refdr(tcrnet_result, candidate_keys, "p_enrichment")
    result.to_csv(args.output, sep="\t", index=False)

    n_sig = int((result["fdr_shared"] < 0.05).sum()) if len(result) else 0
    print(f"TCRnet: {len(result)} shared candidates tested, {n_sig} significant (shared FDR < 0.05)")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
