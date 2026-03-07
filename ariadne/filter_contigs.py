#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parent.parent))

from ariadne.fasta_utils import read_fasta, write_fasta
from ariadne.filtering import filter_by_length


def main() -> int:
    parser = argparse.ArgumentParser(description="Legacy minimum-length filter wrapper for Ariadne.")
    parser.add_argument("min_length", type=int)
    parser.add_argument("input_fasta", type=Path)
    args = parser.parse_args()
    records = filter_by_length(read_fasta(args.input_fasta), args.min_length)
    output_path = args.input_fasta.with_name(f"{args.input_fasta.stem}.filter{args.min_length}{args.input_fasta.suffix}")
    write_fasta(records, output_path)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
