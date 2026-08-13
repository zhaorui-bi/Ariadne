#!/usr/bin/env bash
#
# Minimal Ariadne README smoke test.
#
# Runs discovery -> filtering -> classification/visualization using a
# user-provided protein input directory and reference FASTA directory.
#
# Set RUN_CEESS=1 to enable the optional ESM2 scoring layer.
#
# Usage:
#   bash tutorial/run_example.sh OUTPUT_DIR PROTEIN_DIR REFERENCE_FASTA_DIR
#
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTPUT_DIR="${1:-}"
INPUT_DIR="${ARIADNE_INPUT_DIR:-${2:-}}"
REFERENCE_DIR="${ARIADNE_REFERENCE_DIR:-${3:-}}"
CEESS_XLSX="${ARIADNE_CEESS_XLSX:-$ROOT_DIR/TPS.xlsx}"

if [[ -z "$OUTPUT_DIR" || -z "$INPUT_DIR" || -z "$REFERENCE_DIR" ]]; then
  echo "Usage: bash tutorial/run_example.sh OUTPUT_DIR PROTEIN_DIR REFERENCE_FASTA_DIR" >&2
  exit 2
fi

ARIADNE="${ARIADNE_BIN:-ariadne}"
if ! command -v "$ARIADNE" >/dev/null 2>&1; then
  if [[ -x "$ROOT_DIR/.venv/bin/ariadne" ]]; then
    ARIADNE="$ROOT_DIR/.venv/bin/ariadne"
  else
    echo "Error: 'ariadne' is not on PATH. Install it first, e.g. 'pip install -e .'." >&2
    exit 1
  fi
fi

run_args=(
  run
  --protein-folder "$INPUT_DIR"
  --reference-dir "$REFERENCE_DIR"
  --output-dir "$OUTPUT_DIR"
)

if [[ "${RUN_CEESS:-0}" == "1" && -f "$CEESS_XLSX" ]]; then
  run_args+=(--ceess-xlsx "$CEESS_XLSX")
else
  run_args+=(--skip-ceess-model)
fi

echo "Running: $ARIADNE ${run_args[*]}"
"$ARIADNE" "${run_args[@]}"

echo
echo "Done. Pipeline summary: $OUTPUT_DIR/pipeline_summary.tsv"
