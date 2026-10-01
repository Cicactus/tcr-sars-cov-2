#!/usr/bin/env python3
"""
Method: tcrdist3 meta-clonotypes -- biochemically-informed CDR3 distance
(Dash et al., 2017, Nature; Mayer-Blackwell et al., 2021, eLife).

Reads the already batch-corrected case/control repertoires produced by
01_prepare_data.py. Neighbour counts (case-vs-case, case-vs-control) are
computed against the FULL repertoires so distance context is not distorted;
the shared candidate list is then applied to decide which clonotypes are
individually Fisher-tested, with FDR restricted to that shared set -- see
README, "Why the shared candidate list is applied differently per method".

Uses scipy sparse output from compute_sparse_rect_distances(): a distance
strictly greater than `radius` is stored as a literal 0 (i.e. "not stored"),
while a true distance of 0 (identical CDR3) is stored as -1 to keep it from
being confused with "no entry". Neighbour counts must therefore be read via
.getnnz(axis=1) (count of non-zero sparse entries per row), NOT via a
"<= radius" comparison on the dense values (which would count every
non-neighbour pair too, since 0 <= radius is always true).
"""
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from scipy import stats
from statsmodels.stats.multitest import multipletests

RADIUS = 10  # tcrdist units; lowered from the tutorial default of 12 after
             # diagnosing excess false negatives -- REDCEA's own vDBSCAN
             # independently estimated eps ~9.6-9.8 on the same cohort
GENE_SUFFIX = "*01"


def add_allele_suffix(series: pd.Series) -> pd.Series:
    return series.apply(lambda g: g if "*" in g else f"{g}{GENE_SUFFIX}")


def to_tcrdist_df(df_pl: pl.DataFrame) -> pd.DataFrame:
    pdf = df_pl.to_pandas()
    unique_df = pdf.drop_duplicates(subset=["junction_aa", "v_call", "j_call"])
    out = pd.DataFrame({
        "cdr3_b_aa": unique_df["junction_aa"],
        "v_b_gene": add_allele_suffix(unique_df["v_call"]),
        "j_b_gene": add_allele_suffix(unique_df["j_call"]),
        "count": unique_df["duplicate_count"] if "duplicate_count" in unique_df.columns else 1,
    })
    return out.reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", required=True, help="output-dir of 01_prepare_data.py")
    parser.add_argument("--output", required=True)
    parser.add_argument("--cpus", type=int, default=1)
    parser.add_argument("--radius", type=int, default=RADIUS)
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

    from tcrdist.repertoire import TCRrep
    from tcrdist.breadth import get_safe_chunk

    case_df = to_tcrdist_df(case_bc)
    control_df = to_tcrdist_df(control_bc)

    tr_case = TCRrep(cell_df=case_df, organism="human", chains=["beta"], compute_distances=False)
    tr_control = TCRrep(cell_df=control_df, organism="human", chains=["beta"], compute_distances=False)
    tr_case.cpus = args.cpus
    tr_control.cpus = args.cpus

    n_case, n_ctrl = len(tr_case.clone_df), len(tr_control.clone_df)
    print(f"tcrdist3: case={n_case} clonotypes, control={n_ctrl} clonotypes, radius={args.radius}")

    t0 = time.time()
    safe_chunk_self = get_safe_chunk(n_case, n_case, target=10**7)
    tr_case.compute_sparse_rect_distances(df=tr_case.clone_df, df2=tr_case.clone_df,
                                           radius=args.radius, chunk_size=safe_chunk_self)
    case_self_rw = tr_case.rw_beta.copy()  # copy before the next call overwrites tr_case.rw_beta
    case_neighbors = np.asarray(case_self_rw.getnnz(axis=1)).flatten() - 1  # -1: exclude self-match
    print(f"  self-distances (case vs case): {time.time() - t0:.1f} s")

    t0 = time.time()
    safe_chunk_cross = get_safe_chunk(n_case, n_ctrl, target=10**7)
    tr_case.compute_sparse_rect_distances(df=tr_case.clone_df, df2=tr_control.clone_df,
                                           radius=args.radius, chunk_size=safe_chunk_cross)
    control_neighbors = np.asarray(tr_case.rw_beta.getnnz(axis=1)).flatten()
    print(f"  cross-distances (case vs control): {time.time() - t0:.1f} s")

    clone_df = tr_case.clone_df
    keys = list(zip(
        clone_df["cdr3_b_aa"],
        clone_df["v_b_gene"].str.split("*").str[0],
        clone_df["j_b_gene"].str.split("*").str[0],
    ))
    candidate_mask = pd.Series(keys, index=clone_df.index).isin(candidate_keys).to_numpy()
    n_candidates = int(candidate_mask.sum())
    print(f"tcrdist3: {n_candidates} of {len(clone_df)} clonotypes match the shared candidate list")

    rows = []
    for i in np.where(candidate_mask)[0]:
        a, c = int(case_neighbors[i]), int(control_neighbors[i])
        b, d = n_case - a, n_ctrl - c
        _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "cdr3_b_aa": clone_df["cdr3_b_aa"].iloc[i],
            "v_b_gene": clone_df["v_b_gene"].iloc[i],
            "j_b_gene": clone_df["j_b_gene"].iloc[i],
            "n_case_neighbors": a,
            "n_control_neighbors": c,
            "p_value": p,
        })
    result = pd.DataFrame(rows).sort_values("p_value")
    if len(result):
        result["fdr_shared"] = multipletests(result["p_value"], method="fdr_bh")[1]
    result.to_csv(args.output, sep="\t", index=False)

    n_sig = int((result["fdr_shared"] < 0.05).sum()) if len(result) else 0
    print(f"tcrdist3: {len(result)} candidates tested, {n_sig} significant (shared FDR < 0.05)")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
