#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Runs the whole pipeline in order: prepare data once (01), then
# Fisher / ALICE / TCRnet / tcrdist3 (02-05), each reading the files
# 01 already wrote. REDCEA lives separately under redcea/ (its own
# conda environment) -- see redcea/README.md.
# ============================================================

WORKDIR="${WORKDIR:-$HOME/covid_tcr_markers_run}"
DATA_DIR="$WORKDIR/raw_data"
PREPARED_DIR="$WORKDIR/prepared"
RESULTS_DIR="$WORKDIR/results"
VENV_DIR="$WORKDIR/venv"
LOG_FILE="$WORKDIR/run.log"
N_CORES="$(nproc)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$WORKDIR" "$DATA_DIR" "$PREPARED_DIR" "$RESULTS_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

LOCK_FILE="$WORKDIR/.running.lock"
if [ -f "$LOCK_FILE" ] && kill -0 "$(cat "$LOCK_FILE")" 2>/dev/null; then
    echo "[error] already running (PID $(cat "$LOCK_FILE"))."
    exit 1
fi
echo $$ > "$LOCK_FILE"
trap 'rm -f "$LOCK_FILE"' EXIT

run_with_heartbeat() {
    local label="$1"
    shift
    "$@" &
    local pid=$!
    local start=$(date +%s)
    while kill -0 "$pid" 2>/dev/null; do
        sleep 60
        if kill -0 "$pid" 2>/dev/null; then
            local elapsed=$(( $(date +%s) - start ))
            printf "[%s] still running: %02d:%02d:%02d\n" \
                "$label" $((elapsed/3600)) $((elapsed%3600/60)) $((elapsed%60))
        fi
    done
    wait "$pid"
}

echo "=== $(date) : start ($N_CORES cores) ==="

# Guards against the futex_wait_queue_me deadlock we hit early on: native
# BLAS/numba libraries inside vdjtools/tcrdist3 otherwise spawn their own
# thread pools on top of --cpus, and the pools fight each other.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# --- system packages ---
SYS_PKG_MARKER="$WORKDIR/.system_packages_installed"
if [ ! -f "$SYS_PKG_MARKER" ]; then
    echo "=== installing system packages ==="
    sudo apt-get update -y
    sudo apt-get install -y build-essential python3-venv unzip curl tmux
    touch "$SYS_PKG_MARKER"
fi

# --- venv ---
if [ ! -d "$VENV_DIR" ]; then
    echo "=== creating venv ==="
    python3 -m venv "$VENV_DIR"
fi
"$VENV_DIR/bin/pip" install -q --upgrade pip
"$VENV_DIR/bin/pip" install -q --prefer-binary \
    vdjtools tcrdist3 polars pyarrow pandas numpy "scipy==1.18.1" statsmodels requests tqdm

# --- accept Anaconda-style Terms of Service where relevant is not needed
#     here (no conda in this venv-only pipeline) ---

# --- step 1: prepare data (run once; skipped on re-run if already done) ---
PREPARE_MARKER="$PREPARED_DIR/.done"
if [ ! -f "$PREPARE_MARKER" ]; then
    echo "=== step 1: prepare data (download + batch-correct + candidates) ==="
    run_with_heartbeat "prepare-data" \
        "$VENV_DIR/bin/python3" "$SCRIPT_DIR/01_prepare_data.py" \
            --data-dir "$DATA_DIR" \
            --output-dir "$PREPARED_DIR" \
            --min-donors 3
    touch "$PREPARE_MARKER"
else
    echo "=== step 1 already done, skipping ==="
fi

# --- step 2-5: run each method (each is idempotent via its own output file) ---
run_method() {
    local name="$1" script="$2" outfile="$3"
    shift 3
    if [ -s "$RESULTS_DIR/$outfile" ]; then
        echo "=== $name already done, skipping ==="
        return
    fi
    echo "=== running $name ==="
    t0=$(date +%s)
    run_with_heartbeat "$name" \
        "$VENV_DIR/bin/python3" "$SCRIPT_DIR/$script" \
            --input-dir "$PREPARED_DIR" \
            --output "$RESULTS_DIR/$outfile" \
            "$@"
    echo "$name: $(( $(date +%s) - t0 )) s"
}

run_method "Fisher"  "02_run_fisher.py"  "fisher_result.tsv"  --cpus "$N_CORES"
run_method "ALICE"   "03_run_alice.py"   "alice_result.tsv"
run_method "TCRnet"  "04_run_tcrnet.py"  "tcrnet_result.tsv"  --cpus "$N_CORES"
run_method "tcrdist3" "05_run_tcrdist3.py" "tcrdist3_result.tsv" --cpus "$N_CORES"

touch "$RESULTS_DIR/DONE"
echo "=== $(date) : all done ==="
echo "Prepared (reusable) data: $PREPARED_DIR"
echo "Method results:           $RESULTS_DIR"
