#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="${VENV_DIR:-$ROOT_DIR/.venv}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
DEMO_ROOT="${1:-$ROOT_DIR/demo_workspace}"
DEMO_INPUT="$DEMO_ROOT/demo_input"
DEMO_OUTPUT="$DEMO_ROOT/demo_output"

if [[ ! -d "$VENV_DIR" ]]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"
python -m pip install -e "$ROOT_DIR"

python -m ariadne prepare-demo --output-dir "$DEMO_INPUT"
python -m ariadne run \
  --transcriptomes "$DEMO_INPUT/transcriptomes/demo_sample_transcripts.fasta" \
  --seed-alignment "$DEMO_INPUT/seed_alignment.fasta" \
  --reference-dir "$DEMO_INPUT/references" \
  --output-dir "$DEMO_OUTPUT"

echo "Demo finished: $DEMO_OUTPUT/pipeline_summary.tsv"
