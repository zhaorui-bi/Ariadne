#!/usr/bin/env python3
"""Legacy command-line wrapper for maximum-length FASTA filtering."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parent.parent))

from ariadne.fasta_utils import read_fasta, write_fasta


def main() -> int:
    """Filter a FASTA file by maximum amino-acid length."""
    parser = argparse.ArgumentParser(description="Legacy maximum-length filter wrapper for Ariadne.")
    parser.add_argument("max_length", type=int)
    parser.add_argument("input_fasta", type=Path)
    args = parser.parse_args()
    records = [record for record in read_fasta(args.input_fasta) if len(record.sequence) <= args.max_length]
    output_path = args.input_fasta.with_name(f"{args.input_fasta.stem}.max{args.max_length}{args.input_fasta.suffix}")
    write_fasta(records, output_path)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
