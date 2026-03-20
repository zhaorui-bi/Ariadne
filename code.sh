#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="${VENV_DIR:-$ROOT_DIR/.venv}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
DEMO_ROOT="${1:-$ROOT_DIR/demo_workspace}"
DEMO_INPUT="$DEMO_ROOT/demo_input"
DEMO_OUTPUT_TRANSCRIPT="$DEMO_ROOT/demo_output_transcriptome"
DEMO_OUTPUT_PROTEIN="$DEMO_ROOT/demo_output_protein"

if [[ ! -d "$VENV_DIR" ]]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"
python -m pip install -e "$ROOT_DIR"

ariadne prepare-demo --output-dir "$DEMO_INPUT"
ariadne run \
  --transcriptomes "$DEMO_INPUT/transcriptomes/demo_sample_transcripts.fasta" \
  --reference-dir "$DEMO_INPUT/references" \
  --output-dir "$DEMO_OUTPUT_TRANSCRIPT"

ariadne run \
  --protein-folder "$DEMO_INPUT/proteins" \
  --reference-dir "$DEMO_INPUT/references" \
  --output-dir "$DEMO_OUTPUT_PROTEIN"

echo "Transcriptome demo finished: $DEMO_OUTPUT_TRANSCRIPT/pipeline_summary.tsv"
echo "Protein-folder demo finished: $DEMO_OUTPUT_PROTEIN/pipeline_summary.tsv"
