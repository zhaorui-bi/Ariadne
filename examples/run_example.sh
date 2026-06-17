#!/usr/bin/env bash
#
# Minimal end-to-end Ariadne example.
#
# Runs discovery -> filtering -> classification on the bundled example data
# (input/all_tps.fasta as the protein input, tree/ as the reference backbone).
#
# By default the MAFFT + IQ-TREE phylogeny stage and the optional ESM2 CeeSs
# scoring stage are skipped so the example runs anywhere without extra tools.
# Set RUN_PHYLOGENY=1 (requires mafft + iqtree) or RUN_CEESS=1 (requires the
# `[esm]` extra) to enable them.
#
# Usage:
#   bash examples/run_example.sh [OUTPUT_DIR]
#
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
OUTPUT_DIR="${1:-$ROOT_DIR/results}"

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
  --protein-folder "$ROOT_DIR/input"
  --reference-dir "$ROOT_DIR/tree"
  --output-dir "$OUTPUT_DIR"
)

if [[ "${RUN_PHYLOGENY:-0}" != "1" ]]; then
  run_args+=(--skip-phylogeny)
fi
if [[ "${RUN_CEESS:-0}" != "1" ]]; then
  run_args+=(--skip-ceess-model)
fi

echo "Running: $ARIADNE ${run_args[*]}"
"$ARIADNE" "${run_args[@]}"

echo
echo "Done. Pipeline summary: $OUTPUT_DIR/pipeline_summary.tsv"
