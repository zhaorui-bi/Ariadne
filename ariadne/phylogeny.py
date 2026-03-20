"""Sequence-based phylogeny utilities powered by MAFFT and IQ-TREE."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Union

from ariadne.fasta_utils import FastaRecord, ensure_directory, slugify, write_fasta, write_tsv
from ariadne.references import load_reference_records

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def _resolve_binary(preferred: Optional[str], *fallbacks: str) -> str:
    """Resolve an external executable name from preferred and fallback options."""
    candidates = []
    if preferred:
        candidates.append(preferred)
    candidates.extend(fallbacks)
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    expected = ", ".join(candidates) if candidates else "<none>"
    raise FileNotFoundError(f"Could not find required executable. Tried: {expected}")


def _unique_header(base: str, seen: set[str]) -> str:
    """Generate a unique FASTA header token suitable for MAFFT/IQ-TREE."""
    header = slugify(base)
    if header not in seen:
        seen.add(header)
        return header
    index = 2
    while f"{header}_{index}" in seen:
        index += 1
    unique = f"{header}_{index}"
    seen.add(unique)
    return unique


def prepare_phylogeny_input(
    candidate_fasta: PathLike,
    reference_dir: PathLike,
    output_dir: PathLike,
) -> dict[str, Path]:
    """Combine references and candidates into a clean FASTA for tree building."""
    from ariadne.fasta_utils import read_fasta

    output_root = ensure_directory(output_dir)
    reference_records = load_reference_records(reference_dir)
    candidate_records = read_fasta(candidate_fasta)
    if not reference_records:
        raise ValueError(f"No reference sequences were found in {reference_dir}.")
    if not candidate_records:
        raise ValueError(f"No candidate sequences were found in {candidate_fasta}.")

    combined_records: list[FastaRecord] = []
    mapping_rows: list[dict[str, object]] = []
    seen_headers: set[str] = set()

    for record in reference_records:
        source = record.metadata.get("source", "reference")
        header = _unique_header(f"ref_{source}_{record.id}", seen_headers)
        combined_records.append(FastaRecord(header=header, sequence=record.sequence, metadata=dict(record.metadata)))
        mapping_rows.append(
            {
                "tree_id": header,
                "record_type": "reference",
                "source": source,
                "original_id": record.id,
                "original_header": record.header,
            }
        )

    for record in candidate_records:
        header = _unique_header(f"candidate_{record.id}", seen_headers)
        combined_records.append(FastaRecord(header=header, sequence=record.sequence, metadata=dict(record.metadata)))
        mapping_rows.append(
            {
                "tree_id": header,
                "record_type": "candidate",
                "source": "candidate",
                "original_id": record.id,
                "original_header": record.header,
            }
        )

    combined_fasta = write_fasta(combined_records, output_root / "phylogeny_input.fasta")
    mapping_tsv = write_tsv(mapping_rows, output_root / "phylogeny_sequence_map.tsv")
    return {
        "phylogeny_input": combined_fasta,
        "sequence_map": mapping_tsv,
    }


def run_mafft(
    input_fasta: PathLike,
    output_alignment: PathLike,
    *,
    mafft_bin: Optional[str] = None,
    mode: str = "--auto",
) -> Path:
    """Run MAFFT and write the aligned FASTA output."""
    mafft = _resolve_binary(mafft_bin, "mafft")
    output_path = Path(output_alignment)
    command = [mafft, mode, str(input_fasta)]
    logger.info("Running MAFFT: %s", " ".join(command))
    with output_path.open("w") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)
    return output_path


def run_iqtree(
    alignment_fasta: PathLike,
    output_dir: PathLike,
    *,
    iqtree_bin: Optional[str] = None,
    model: str = "LG",
    threads: str = "AUTO",
    bootstrap: Optional[int] = None,
    fast: bool = True,
) -> dict[str, Path]:
    """Run IQ-TREE on an alignment and collect the main output files."""
    iqtree = _resolve_binary(iqtree_bin, "iqtree2", "iqtree")
    output_root = ensure_directory(output_dir)
    prefix = output_root / "iqtree"
    command = [
        iqtree,
        "-s",
        str(alignment_fasta),
        "--prefix",
        str(prefix),
        "-m",
        model,
        "-T",
        str(threads),
    ]
    if fast:
        command.append("--fast")
    if bootstrap is not None and bootstrap > 0:
        command.extend(["-B", str(bootstrap)])
    logger.info("Running IQ-TREE: %s", " ".join(command))
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    outputs = {
        "iqtree_prefix": prefix,
        "iqtree_treefile": prefix.with_suffix(".treefile"),
        "iqtree_report": prefix.with_suffix(".iqtree"),
        "iqtree_log": prefix.with_suffix(".log"),
    }
    optional = {
        "iqtree_consensus_tree": prefix.with_suffix(".contree"),
        "iqtree_checkpoint": prefix.with_suffix(".ckp.gz"),
    }
    for key, path in optional.items():
        if path.exists():
            outputs[key] = path
    return outputs


def build_phylogeny(
    candidate_fasta: PathLike,
    reference_dir: PathLike,
    output_dir: PathLike,
    *,
    mafft_bin: Optional[str] = None,
    mafft_mode: str = "--auto",
    iqtree_bin: Optional[str] = None,
    iqtree_model: str = "LG",
    iqtree_threads: str = "AUTO",
    iqtree_bootstrap: Optional[int] = None,
    iqtree_fast: bool = True,
) -> dict[str, Path]:
    """Build a sequence alignment and phylogenetic tree from references plus candidates."""
    output_root = ensure_directory(output_dir)
    prepared = prepare_phylogeny_input(candidate_fasta, reference_dir, output_root)
    alignment_path = run_mafft(
        prepared["phylogeny_input"],
        output_root / "phylogeny_alignment.fasta",
        mafft_bin=mafft_bin,
        mode=mafft_mode,
    )
    iqtree_outputs = run_iqtree(
        alignment_path,
        output_root,
        iqtree_bin=iqtree_bin,
        model=iqtree_model,
        threads=iqtree_threads,
        bootstrap=iqtree_bootstrap,
        fast=iqtree_fast,
    )
    outputs = dict(prepared)
    outputs["phylogeny_alignment"] = alignment_path
    outputs.update(iqtree_outputs)
    return outputs
