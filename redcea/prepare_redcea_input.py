#!/usr/bin/env python3
"""
Converts the shared, already batch-corrected repertoires (written by
scripts/01_prepare_data.py) into the TSV schema REDCEA expects
(v_call, j_call, junction_aa, locus). Case is restricted to the shared
candidate list; control is passed through in full as the background.

Run this on the machine/venv that has pandas installed -- it does NOT need
vdjtools, since the input files are already in AIRR schema.
"""
import argparse
from pathlib import Path

import pandas as pd

CHAIN = "TRB"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-dir", required=True, help="output-dir of scripts/01_prepare_data.py")
    parser.add_argument("--out-case", required=True)
    parser.add_argument("--out-control", required=True)
    args = parser.parse_args()

    prepared_dir = Path(args.prepared_dir)
    case_df = pd.read_csv(prepared_dir / "case_repertoire_corrected.tsv", sep="\t")
    control_df = pd.read_csv(prepared_dir / "control_repertoire_corrected.tsv", sep="\t")
    candidates = pd.read_csv(prepared_dir / "candidates.tsv", sep="\t")

    candidate_keys = set(zip(candidates["junction_aa"], candidates["v_call"], candidates["j_call"]))
    case_keys = list(zip(case_df["junction_aa"], case_df["v_call"], case_df["j_call"]))
    case_df = case_df[pd.Series(case_keys, index=case_df.index).isin(candidate_keys)]

    def to_redcea_tsv(df, out_path):
        unique_df = df.drop_duplicates(subset=["junction_aa", "v_call", "j_call"])
        out = pd.DataFrame({
            "v_call": unique_df["v_call"],
            "j_call": unique_df["j_call"],
            "junction_aa": unique_df["junction_aa"],
            "locus": CHAIN,
        })
        out.to_csv(out_path, sep="\t", index=False)
        print(f"{out_path}: {len(out)} unique clonotypes")

    to_redcea_tsv(case_df, args.out_case)       # candidates only
    to_redcea_tsv(control_df, args.out_control)  # full background, not filtered


if __name__ == "__main__":
    main()
