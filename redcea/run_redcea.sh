#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# REDCEA (official antigenomics/redcea package) on the shared,
# already prepared/batch-corrected data from scripts/01_prepare_data.py.
#
# REDCEA needs conda (a heavier, separate environment from the rest of the
# pipeline), so it is kept on its own here rather than folded into
# scripts/run_all.sh. Before running this, copy the "prepared/" directory
# produced by scripts/01_prepare_data.py onto this machine, e.g.:
#   rsync -avz user@data-host:~/covid_tcr_markers_run/prepared/ ./prepared_data/
# ============================================================

WORKDIR="${WORKDIR:-$HOME/redcea_run}"
PREPARED_DIR="${PREPARED_DIR:-$HOME/covid_tcr_markers_run/prepared}"
DATA_DIR="$WORKDIR/data"
RESULTS_DIR="$WORKDIR/results"
LOG_FILE="$WORKDIR/run.log"
VENV_DIR="$WORKDIR/venv"
CHAIN="TRB"
N_CORES="$(nproc)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$WORKDIR" "$DATA_DIR" "$RESULTS_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

START_TS=$(date +%s)

fail_notify() {
    local exit_code=$?
    echo "=== $(date) : ERROR (code $exit_code) ==="
    echo "See the full log: $LOG_FILE"
    exit "$exit_code"
}
trap fail_notify ERR

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

echo "=== $(date) : start (WORKDIR=$WORKDIR, cores: $N_CORES) ==="

# --- 1. system packages ---
SYS_PKG_MARKER="$WORKDIR/.system_packages_installed"
if [ ! -f "$SYS_PKG_MARKER" ]; then
    echo "=== installing system packages ==="
    sudo apt-get update -y
    sudo apt-get install -y git unzip build-essential curl tmux python3-venv
    touch "$SYS_PKG_MARKER"
fi

# --- 2. Miniconda ---
if [ ! -d "$HOME/miniconda3" ]; then
    echo "=== installing Miniconda ==="
    curl -sL -o "$WORKDIR/miniconda.sh" https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash "$WORKDIR/miniconda.sh" -b -p "$HOME/miniconda3"
    rm "$WORKDIR/miniconda.sh"
fi
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Anaconda requires accepting Terms of Service for the default channels in
# non-interactive mode; without this, `conda create`/`conda install` fails
# with CondaToSNonInteractiveError. Idempotent, safe to re-run.
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# --- 3. venv for the adapter script (isolated from the redcea conda env) ---
if [ ! -d "$VENV_DIR" ]; then
    python3 -m venv "$VENV_DIR"
fi
"$VENV_DIR/bin/pip" install -q --upgrade pip
"$VENV_DIR/bin/pip" install -q pandas

# --- 4. convert the shared prepared data into REDCEA's TSV schema ---
CASE_TSV="$DATA_DIR/redcea_case.tsv"
CONTROL_TSV="$DATA_DIR/redcea_control.tsv"
if [ ! -f "$CASE_TSV" ] || [ ! -f "$CONTROL_TSV" ]; then
    echo "=== converting prepared data for REDCEA ==="
    "$VENV_DIR/bin/python3" "$SCRIPT_DIR/prepare_redcea_input.py" \
        --prepared-dir "$PREPARED_DIR" \
        --out-case "$CASE_TSV" \
        --out-control "$CONTROL_TSV"
else
    echo "=== REDCEA input already prepared, skipping ==="
fi

# --- 5. clone + install redcea ---
REDCEA_DIR="$WORKDIR/redcea"
if [ ! -d "$REDCEA_DIR" ]; then
    git clone https://github.com/antigenomics/redcea.git "$REDCEA_DIR" \
        || git clone https://gitlab.aldan3.itm-rsmu.ru/isagroup/redcea.git "$REDCEA_DIR"
fi
cd "$REDCEA_DIR"

if ! conda env list | grep -qE "^redcea\s"; then
    conda create -n redcea python=3.11 -y
fi
if ! conda run -n redcea python -c "import redcea, tcremp, mir" >/dev/null 2>&1; then
    conda run -n redcea python -m pip install --upgrade pip setuptools wheel
    conda run -n redcea python -m pip install -e .
fi
conda run -n redcea python -c "import redcea, tcremp, mir; print('imports: OK')"

# --- 6. embedding (two-step, checkpointed) ---
CASE_EMB_MARKER="$RESULTS_DIR/.case_embedding_done"
if [ ! -f "$CASE_EMB_MARKER" ]; then
    t0=$(date +%s)
    run_with_heartbeat "case embedding" \
        conda run -n redcea tcremp-run \
            --input "$CASE_TSV" --output "$RESULTS_DIR" --chain "$CHAIN" -np "$N_CORES"
    touch "$CASE_EMB_MARKER"
    echo "case embedding: $(( $(date +%s) - t0 )) s"
fi

CONTROL_EMB_MARKER="$RESULTS_DIR/.control_embedding_done"
if [ ! -f "$CONTROL_EMB_MARKER" ]; then
    t0=$(date +%s)
    run_with_heartbeat "control embedding" \
        conda run -n redcea tcremp-run \
            --input "$CONTROL_TSV" --output "$RESULTS_DIR" --chain "$CHAIN" -np "$N_CORES"
    touch "$CONTROL_EMB_MARKER"
    echo "control embedding: $(( $(date +%s) - t0 )) s"
fi

# --- 7. redcea itself, on the already-computed embeddings ---
SAMPLE_EMB=$(ls "$RESULTS_DIR"/*redcea_case*tcremp.parquet 2>/dev/null | head -n1 || true)
BG_EMB=$(ls "$RESULTS_DIR"/*redcea_control*tcremp.parquet 2>/dev/null | head -n1 || true)

REDCEA_ARGS=(-is "$CASE_TSV" -ib "$CONTROL_TSV" -c "$CHAIN" -o "$RESULTS_DIR" -np "$N_CORES")
[ -n "$SAMPLE_EMB" ] && REDCEA_ARGS+=(-se "$SAMPLE_EMB")
[ -n "$BG_EMB" ] && REDCEA_ARGS+=(-be "$BG_EMB")

t0=$(date +%s)
set +e  # known redcea packaging bug: main() can return a non-empty object,
        # which sys.exit() then treats as an error even on a successful run
        # -- checked properly against actual output files in step 8 below
run_with_heartbeat "redcea" conda run -n redcea redcea "${REDCEA_ARGS[@]}"
REDCEA_EXIT=$?
set -e
echo "redcea: $(( $(date +%s) - t0 )) s (exit code: $REDCEA_EXIT)"
if [ "$REDCEA_EXIT" -ne 0 ]; then
    echo "[warning] redcea returned code $REDCEA_EXIT -- checking success by actual output instead of exit code"
fi

# --- 8. check success by evidence, not by exit code ---
REDCEA_LOG=$(ls -t "$RESULTS_DIR"/*.log 2>/dev/null | head -n1 || true)
if [ -z "$REDCEA_LOG" ] || ! grep -q "TCRempNet pipeline completed." "$REDCEA_LOG"; then
    echo "[error] redcea log does not contain 'TCRempNet pipeline completed.' -- treating the run as failed"
    exit 1
fi
SUMMARY_FILE=$(ls "$RESULTS_DIR"/*_summary_tcrempnet.tsv 2>/dev/null | head -n1 || true)
if [ -z "$SUMMARY_FILE" ]; then
    echo "[error] *_summary_tcrempnet.tsv not found -- treating the run as failed"
    exit 1
fi

ELAPSED=$(( $(date +%s) - START_TS ))
echo "=== $(date) : done, total ${ELAPSED}s ==="
echo "Results: $RESULTS_DIR"

trap - ERR
touch "$RESULTS_DIR/DONE"
