#!/usr/bin/env python3
"""Legacy exact-duplicate FASTA remover maintained for compatibility."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parent.parent))

from ariadne.fasta_utils import read_fasta, write_fasta
from ariadne.filtering import deduplicate_exact


def main() -> int:
    """Remove exact duplicate sequences from a FASTA file."""
    parser = argparse.ArgumentParser(description="Exact duplicate remover for Ariadne FASTA files.")
    parser.add_argument("-i", "--input", required=True, type=Path, dest="input_fasta")
    parser.add_argument("-o", "--output", type=Path, default=None)
    parser.add_argument("-v", "--verbose", default="Y")
    args = parser.parse_args()

    input_records = read_fasta(args.input_fasta)
    output_records = deduplicate_exact(input_records)
    output_path = args.output or args.input_fasta.with_name(f"Uniq_{args.input_fasta.name}")
    write_fasta(output_records, output_path)

    if args.verbose.upper() == "Y":
        print(f"[Input seq]\t: {len(input_records)}")
        print(f"[Output seq]\t: {len(output_records)}")
        print(f"[Duplicates]\t: {len(input_records) - len(output_records)}")
        print(f"[Output file]\t: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
