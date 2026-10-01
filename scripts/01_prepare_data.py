#!/usr/bin/env python3
"""
Step 1 of the pipeline -- run this ONCE.

Downloads the cohort from Zenodo, corrects the V/J-usage batch effect
(9 sequencing batches, see notebooks/01_batch_effect.ipynb for why this is
necessary), and selects the shared candidate clonotype list (clonotypes seen
in >= MIN_DONORS distinct COVID train-donors).

Every method script under scripts/ (02-05) reads its inputs from the files
this script writes -- none of them re-download or re-correct anything. This
is the fix for the "each method recomputes its own slightly-different
candidate set" problem: there is exactly one batch correction, computed
exactly once, and exactly one candidate list, shared by every method.

Outputs (all written to --output-dir):
    metadata_split.csv                 donor id, COVID status, batch, train/test split
    case_repertoire_corrected.tsv      full COVID train cohort, batch-corrected
    control_repertoire_corrected.tsv   full healthy train cohort, batch-corrected
    candidates.tsv                     clonotypes seen in >= MIN_DONORS case donors
    test_case_repertoire_corrected.tsv held-out COVID donors, same correction applied
    test_control_repertoire_corrected.tsv  held-out healthy donors, same correction applied

Why the held-out (test) donors go through the SAME correction, not their own:
the correction (preprocess.correct_vj_usage) is fit ONCE on the train pool
only. Applying that already-fitted correction to the test donors (rather than
re-fitting on test, or leaving test uncorrected) is the standard
fit-on-train/transform-on-test practice -- it evaluates methods under the
same technical adjustment they were discovered under, without ever letting
the test donors influence the correction itself.
"""
import argparse
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import requests
from tqdm import tqdm

ZENODO_RECORD = "17780804"
CHAIN = "TRB"
SEED = 42
MIN_DONORS = 3


def zenodo_download(record: str, filename: str, out_path: Path) -> None:
    url = f"https://zenodo.org/api/records/{record}/files/{filename}/content"
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1 << 16):
            f.write(chunk)


def load_group(ids_list, all_files, id_to_batch, desc):
    """Read a group of donors through vdjtools.io (schema auto-detection +
    non-productive-clonotype filtering), tagging each row with its donor id
    and sequencing batch."""
    from vdjtools import io as vio, preprocess

    frames = []
    for donor_id in tqdm(ids_list, desc=desc):
        matches = [f for f in all_files if donor_id in f.name]
        if not matches:
            continue
        path = matches[0]
        # Zenodo's per-donor files are comma-separated but use the classic
        # VDJtools column names (count, freq, cdr3nt, cdr3aa, v, d, j, ...);
        # vio.read()'s format sniffer only recognises that schema over a tab
        # delimiter, so we convert once and cache the .tsv next to the .csv.
        tsv_path = path.with_suffix(".tsv")
        if not tsv_path.exists():
            pd.read_csv(path).to_csv(tsv_path, sep="\t", index=False)
        df = vio.read(str(tsv_path))
        df = preprocess.filter_productive(df)  # drops "_" (out-of-frame) / "*" (stop codon) junctions
        df = df.with_columns(
            pl.lit(donor_id).alias("sample_id"),
            pl.lit(id_to_batch.get(donor_id, "")).alias("batch"),
        )
        frames.append(df)
    if not frames:
        raise RuntimeError(f"No repertoire files found for group '{desc}'")
    return pl.concat(frames)


def fit_correction(case_train_pl, control_train_pl):
    """Fit the V/J-usage batch correction on the train pool ONLY."""
    from vdjtools import preprocess

    train_pool = pl.concat([case_train_pl, control_train_pl])
    return preprocess.correct_vj_usage(train_pool, batch_col="batch", transform="sigmoid")


def apply_correction(df, usage_correction):
    """Apply an already-fitted correction to any set of donors (train or test)."""
    from vdjtools import preprocess

    out = []
    for sid in df["sample_id"].unique():
        sub = df.filter(pl.col("sample_id") == sid)
        try:
            corrected = preprocess.apply_vj_correction(sub, usage_correction, sample_id=sid)
        except Exception as e:
            print(f"  [warning] correction failed for donor {sid}, keeping uncorrected rows: {e}")
            corrected = sub
        out.append(corrected)
    return pl.concat(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", required=True, help="raw download cache (Zenodo files, extracted archive)")
    parser.add_argument("--output-dir", required=True, help="where the prepared, reusable files are written")
    parser.add_argument("--min-donors", type=int, default=MIN_DONORS,
                         help="a clonotype is a 'candidate' if seen in at least this many case donors")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- download -----------------------------------------------------
    for fname in ["patient_metadata.csv", "TRB.zip", "covid_clonotypes.csv"]:
        out = data_dir / fname
        if not out.exists():
            print(f"[download] {fname} ...")
            zenodo_download(ZENODO_RECORD, fname, out)

    # --- train/test split (fixed seed -- reproduced identically by every
    #     script and notebook that needs it) --------------------------
    metadata = pd.read_csv(data_dir / "patient_metadata.csv")
    rng = np.random.default_rng(SEED)
    ids = metadata["id"].unique()
    rng.shuffle(ids)
    n_test = max(1, len(ids) // 5)
    test_ids = set(ids[:n_test])
    metadata["split"] = np.where(metadata["id"].isin(test_ids), "test", "train")
    metadata.to_csv(output_dir / "metadata_split.csv", index=False)

    case_train_ids = metadata.loc[(metadata.COVID_status == "COVID") & (metadata.split == "train"), "id"].tolist()
    control_train_ids = metadata.loc[(metadata.COVID_status == "healthy") & (metadata.split == "train"), "id"].tolist()
    case_test_ids = metadata.loc[(metadata.COVID_status == "COVID") & (metadata.split == "test"), "id"].tolist()
    control_test_ids = metadata.loc[(metadata.COVID_status == "healthy") & (metadata.split == "test"), "id"].tolist()
    print(f"[split] train: {len(case_train_ids)} COVID / {len(control_train_ids)} healthy")
    print(f"[split] test:  {len(case_test_ids)} COVID / {len(control_test_ids)} healthy")

    # --- extract repertoire archive -------------------------------------
    extract_dir = data_dir / f"{CHAIN}_extracted"
    if not extract_dir.exists():
        print("[extract] unzipping repertoire archive ...")
        with zipfile.ZipFile(data_dir / f"{CHAIN}.zip") as z:
            z.extractall(extract_dir)
    all_files = list(extract_dir.rglob("*.csv"))
    id_to_batch = metadata.set_index("id")["batch"].to_dict()

    # --- load + batch-correct -------------------------------------------
    case_train_pl = load_group(case_train_ids, all_files, id_to_batch, "case (train)")
    control_train_pl = load_group(control_train_ids, all_files, id_to_batch, "control (train)")
    case_test_pl = load_group(case_test_ids, all_files, id_to_batch, "case (test)")
    control_test_pl = load_group(control_test_ids, all_files, id_to_batch, "control (test)")

    print("[batch-correction] fitting on the train pool ...")
    usage_correction = fit_correction(case_train_pl, control_train_pl)

    case_train_bc = apply_correction(case_train_pl, usage_correction)
    control_train_bc = apply_correction(control_train_pl, usage_correction)
    case_test_bc = apply_correction(case_test_pl, usage_correction)
    control_test_bc = apply_correction(control_test_pl, usage_correction)

    case_train_bc.write_csv(str(output_dir / "case_repertoire_corrected.tsv"), separator="\t")
    control_train_bc.write_csv(str(output_dir / "control_repertoire_corrected.tsv"), separator="\t")
    case_test_bc.write_csv(str(output_dir / "test_case_repertoire_corrected.tsv"), separator="\t")
    control_test_bc.write_csv(str(output_dir / "test_control_repertoire_corrected.tsv"), separator="\t")
    print(f"[save] case_repertoire_corrected.tsv: {case_train_bc.height} rows")
    print(f"[save] control_repertoire_corrected.tsv: {control_train_bc.height} rows")
    print(f"[save] test_case_repertoire_corrected.tsv: {case_test_bc.height} rows")
    print(f"[save] test_control_repertoire_corrected.tsv: {control_test_bc.height} rows")

    # --- shared candidate list -------------------------------------------
    counts = case_train_bc.group_by(["junction_aa", "v_call", "j_call"]).agg(
        pl.col("sample_id").n_unique().alias("n_donors")
    )
    candidates = counts.filter(pl.col("n_donors") >= args.min_donors)
    candidates.write_csv(str(output_dir / "candidates.tsv"), separator="\t")
    print(f"[candidates] {candidates.height} candidate clonotypes (seen in >= {args.min_donors} "
          f"case donors) out of {counts.height} unique clonotypes total")

    print(f"\nDone -- every file used by scripts/02-05 now lives in {output_dir}")


if __name__ == "__main__":
    main()
