#!/usr/bin/env python3
"""
Method: ALICE (vdjtools.overlap.alice) -- neighbourhood enrichment against a
generative V(D)J-recombination null model (Pogorelyy et al., 2019, PLoS Biol).

Reads the already batch-corrected case repertoire produced by
01_prepare_data.py. Runs on the FULL case cohort (not pre-filtered to
candidates) so that neighbour counting sees the true local context; the
shared candidate list is applied AFTER, to the raw p-values, together with a
fresh FDR correction restricted to that shared set -- see README, "Why the
shared candidate list is applied differently per method".
"""
import argparse
from pathlib import Path

import pandas as pd
import polars as pl
from statsmodels.stats.multitest import multipletests

CHAIN = "TRB"


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
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    schema_overrides = {"sample_id": pl.Utf8, "batch": pl.Utf8}
    case_bc = pl.read_csv(str(input_dir / "case_repertoire_corrected.tsv"), separator="\t", schema_overrides=schema_overrides)
    candidates_pl = pl.read_csv(str(input_dir / "candidates.tsv"), separator="\t")
    candidate_keys = set(zip(
        candidates_pl["junction_aa"].to_list(),
        candidates_pl["v_call"].to_list(),
        candidates_pl["j_call"].to_list(),
    ))

    from vdjtools.model import load_bundled
    from vdjtools import overlap

    model = load_bundled(CHAIN, source="olga")
    valid_v = set(model.genomic["genes_v"]["gene"].to_list())
    valid_j = set(model.genomic["genes_j"]["gene"].to_list())
    filtered = case_bc.filter(pl.col("v_call").is_in(list(valid_v)) & pl.col("j_call").is_in(list(valid_j)))
    dropped = case_bc.height - filtered.height
    if dropped:
        print(f"[ALICE] dropped {dropped} clonotypes with genes outside the model "
              f"({dropped / case_bc.height:.1%} of {case_bc.height})")

    alice_result = overlap.alice(filtered).to_pandas()
    result = restrict_and_refdr(alice_result, candidate_keys, "p_enrichment")
    result.to_csv(args.output, sep="\t", index=False)

    n_sig = int((result["fdr_shared"] < 0.05).sum()) if len(result) else 0
    print(f"ALICE: {len(result)} shared candidates tested, {n_sig} significant (shared FDR < 0.05)")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
