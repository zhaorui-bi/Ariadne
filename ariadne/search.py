"""Stage 1: build query HMMs and discover candidate TPS proteins.

This module supports two discovery modes:
1. start from transcript FASTA files, predict ORFs with ``pyrodigal``, then
   search translated proteins with a profile HMM;
2. start from precomputed protein FASTA files and search them directly.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import pyhmmer
import pyrodigal

from ariadne.utils import FastaRecord, ensure_directory, parse_coverage, read_fasta, write_fasta, write_tsv

logger = logging.getLogger(__name__)

PathLike = Union[str, Path]


def _resolve_mafft_binary(preferred: Optional[str] = None) -> str:
    """Locate a MAFFT executable for temporary alignment generation."""
    for candidate in [preferred, "mafft"]:
        if not candidate:
            continue
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise FileNotFoundError("Could not find MAFFT executable. Please install `mafft` or pass a prebuilt HMM.")


def _run_mafft_alignment(input_fasta: PathLike, output_fasta: PathLike, *, mafft_bin: Optional[str] = None) -> Path:
    """Align an unaligned protein FASTA with MAFFT."""
    mafft = _resolve_mafft_binary(mafft_bin)
    output_path = Path(output_fasta)
    command = [mafft, "--auto", str(input_fasta)]
    logger.info("Input FASTA is not aligned; running MAFFT before hmmbuild: %s", " ".join(command))
    with output_path.open("w") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)
    return output_path


def build_hmm(alignment_fasta: PathLike, output_hmm: PathLike, *, name: Optional[str] = None) -> Path:
    """Build a profile HMM from an aligned FASTA/MSA file.

    If the input FASTA is not a valid multiple sequence alignment, Ariadne
    temporarily aligns it with MAFFT before running ``hmmbuild``.
    """
    alignment_path = Path(alignment_fasta)
    alphabet = pyhmmer.easel.Alphabet.amino()
    format_name = None
    if alignment_path.suffix.lower() in {".fa", ".faa", ".fasta", ".afa"}:
        format_name = "afa"
    temp_alignment_path: Optional[Path] = None
    try:
        try:
            with pyhmmer.easel.MSAFile(str(alignment_path), digital=True, alphabet=alphabet, format=format_name) as msa_file:
                msa = msa_file.read()
        except ValueError:
            with tempfile.NamedTemporaryFile("w", suffix=".afa", delete=False) as handle:
                temp_alignment_path = Path(handle.name)
            _run_mafft_alignment(alignment_path, temp_alignment_path)
            with pyhmmer.easel.MSAFile(str(temp_alignment_path), digital=True, alphabet=alphabet, format="afa") as msa_file:
                msa = msa_file.read()
    finally:
        if temp_alignment_path is not None and temp_alignment_path.exists():
            temp_alignment_path.unlink(missing_ok=True)
    if not msa.name:
        msa.name = name or Path(output_hmm).stem
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)
    output_path = Path(output_hmm)
    with output_path.open("wb") as handle:
        hmm.write(handle, binary=False)
    return output_path


def _protein_and_nt_records(
    transcriptome_path: PathLike,
    sample_name: Optional[str] = None,
) -> Tuple[List[FastaRecord], Dict[str, FastaRecord]]:
    """Predict ORFs from one transcriptome and return protein/ORF records."""
    transcript_records = read_fasta(transcriptome_path)
    sample = sample_name or Path(transcriptome_path).stem
    finder = pyrodigal.GeneFinder(meta=True)
    proteins: List[FastaRecord] = []
    nucleotide_by_id: Dict[str, FastaRecord] = {}
    for transcript in transcript_records:
        genes = finder.find_genes(transcript.sequence)
        coverage = parse_coverage(transcript.header)
        for index, gene in enumerate(genes, start=1):
            protein_id = f"{sample}|{transcript.id}|orf{index}"
            coverage_text = f"cov_{coverage:g}" if coverage is not None else "cov_na"
            header = (
                f"{protein_id} {coverage_text} begin={gene.begin} end={gene.end} "
                f"strand={gene.strand} transcript={transcript.id}"
            )
            protein_sequence = str(gene.translate()).replace("*", "")
            nucleotide_sequence = str(gene.sequence())
            protein_record = FastaRecord(header=header, sequence=protein_sequence)
            protein_record.metadata["sample"] = sample
            protein_record.metadata["transcript_id"] = transcript.id
            nt_record = FastaRecord(header=header, sequence=nucleotide_sequence)
            nt_record.metadata.update(protein_record.metadata)
            proteins.append(protein_record)
            nucleotide_by_id[protein_record.id] = nt_record
    return proteins, nucleotide_by_id


def search_proteins_with_hmm(
    protein_fasta: PathLike,
    hmm_path: PathLike,
    *,
    min_score: Optional[float] = None,
    max_evalue: Optional[float] = None,
) -> list[dict[str, object]]:
    """Search proteins against one or more HMMs and keep the best hit per ID."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)
    pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
    with pyhmmer.easel.SequenceFile(str(protein_fasta), digital=True, alphabet=alphabet) as sequence_file:
        sequences = sequence_file.read_block()
    best_hits: dict[str, dict[str, object]] = {}
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hmm_file:
        hmms = list(hmm_file)
    for hmm in hmms:
        hits = pipeline.search_hmm(hmm, sequences)
        for hit in hits:
            if min_score is not None and hit.score < min_score:
                continue
            if max_evalue is not None and hit.evalue > max_evalue:
                continue
            domain = hit.best_domain
            row = {
                "sequence_id": hit.name,
                "hmm_name": getattr(hmm, "name", Path(hmm_path).stem),
                "score": round(hit.score, 4),
                "evalue": hit.evalue,
                "bias": round(hit.bias, 4),
                "domain_score": round(domain.score, 4),
                "domain_i_evalue": domain.i_evalue,
                "domain_c_evalue": domain.c_evalue,
                "env_from": domain.env_from,
                "env_to": domain.env_to,
            }
            current = best_hits.get(hit.name)
            if current is None or float(row["score"]) > float(current["score"]):
                best_hits[hit.name] = row
    return sorted(best_hits.values(), key=lambda row: (-float(row["score"]), float(row["evalue"])))


def discover_candidates(
    transcriptome_paths: List[PathLike],
    hmm_path: PathLike,
    output_dir: PathLike,
    *,
    min_score: Optional[float] = None,
    max_evalue: Optional[float] = None,
) -> dict[str, Path]:
    """Run discovery from transcript FASTA files all the way to candidate hits."""
    root = ensure_directory(output_dir)
    combined_proteins: list[FastaRecord] = []
    combined_hits_proteins: list[FastaRecord] = []
    combined_hits_nucleotides: list[FastaRecord] = []
    combined_hit_rows: list[dict[str, object]] = []

    for transcriptome_path in transcriptome_paths:
        sample_name = Path(transcriptome_path).stem
        sample_dir = ensure_directory(root / sample_name)
        proteins, nucleotide_by_id = _protein_and_nt_records(transcriptome_path, sample_name=sample_name)
        proteins_path = sample_dir / f"{sample_name}.proteins.faa"
        nucleotides_path = sample_dir / f"{sample_name}.orfs.fna"
        write_fasta(proteins, proteins_path)
        write_fasta(nucleotide_by_id.values(), nucleotides_path)
        combined_proteins.extend(proteins)

        hit_rows = search_proteins_with_hmm(proteins_path, hmm_path, min_score=min_score, max_evalue=max_evalue)
        protein_by_id = {record.id: record for record in proteins}
        hit_ids = {row["sequence_id"] for row in hit_rows}
        sample_hit_proteins = [protein_by_id[sequence_id] for sequence_id in protein_by_id if sequence_id in hit_ids]
        sample_hit_nucleotides = [nucleotide_by_id[sequence_id] for sequence_id in nucleotide_by_id if sequence_id in hit_ids]
        for row in hit_rows:
            row["sample"] = sample_name
        combined_hit_rows.extend(hit_rows)
        combined_hits_proteins.extend(sample_hit_proteins)
        combined_hits_nucleotides.extend(sample_hit_nucleotides)
        write_fasta(sample_hit_proteins, sample_dir / f"{sample_name}.hits.faa")
        write_fasta(sample_hit_nucleotides, sample_dir / f"{sample_name}.hits.fna")
        write_tsv(hit_rows, sample_dir / f"{sample_name}.hits.tsv")

    combined_proteins_path = write_fasta(combined_proteins, root / "all_predicted_proteins.faa")
    candidate_proteins_path = write_fasta(combined_hits_proteins, root / "candidates.protein.faa")
    candidate_nucleotides_path = write_fasta(combined_hits_nucleotides, root / "candidates.orf.fna")
    hits_tsv_path = write_tsv(combined_hit_rows, root / "candidates.hits.tsv")
    return {
        "all_predicted_proteins": combined_proteins_path,
        "candidate_proteins": candidate_proteins_path,
        "candidate_nucleotides": candidate_nucleotides_path,
        "hits_tsv": hits_tsv_path,
    }


def collect_protein_files(
    protein_dir: PathLike,
    *,
    protein_glob: Sequence[str] = ("*.faa", "*.fa", "*.fasta", "*.pep", "*.prot"),
) -> list[Path]:
    """Recursively collect protein FASTA files from a directory."""
    root = Path(protein_dir)
    if not root.exists():
        raise FileNotFoundError(f"Protein directory does not exist: {root}")
    files: dict[Path, None] = {}
    for pattern in protein_glob:
        for path in root.rglob(pattern):
            if path.is_file():
                files[path] = None
    return sorted(files)


def _protein_records_with_sample_prefix(records: Sequence[FastaRecord], sample_name: str) -> list[FastaRecord]:
    """Prefix record IDs with the sample name to avoid cross-file collisions."""
    normalized: list[FastaRecord] = []
    seen_ids: set[str] = set()
    for index, record in enumerate(records, start=1):
        base_id = f"{sample_name}|{record.id}"
        record_id = base_id
        if record_id in seen_ids:
            record_id = f"{base_id}|dup{index}"
        seen_ids.add(record_id)
        description = record.description
        header = record_id if not description else f"{record_id} {description}"
        normalized_record = FastaRecord(header=header, sequence=record.sequence)
        normalized_record.metadata["sample"] = sample_name
        normalized.append(normalized_record)
    return normalized


def discover_candidates_from_proteins(
    protein_paths: Sequence[PathLike],
    hmm_path: PathLike,
    output_dir: PathLike,
    *,
    min_score: Optional[float] = None,
    max_evalue: Optional[float] = None,
) -> dict[str, Path]:
    """Run discovery directly from protein FASTA files."""
    root = ensure_directory(output_dir)
    combined_proteins: list[FastaRecord] = []
    combined_hits_proteins: list[FastaRecord] = []
    combined_hit_rows: list[dict[str, object]] = []

    for protein_path in protein_paths:
        source_path = Path(protein_path)
        sample_name = source_path.stem
        sample_dir = ensure_directory(root / sample_name)
        records = _protein_records_with_sample_prefix(read_fasta(source_path), sample_name)
        proteins_path = write_fasta(records, sample_dir / f"{sample_name}.proteins.faa")
        combined_proteins.extend(records)

        hit_rows = search_proteins_with_hmm(proteins_path, hmm_path, min_score=min_score, max_evalue=max_evalue)
        protein_by_id = {record.id: record for record in records}
        hit_ids = {row["sequence_id"] for row in hit_rows}
        sample_hit_proteins = [protein_by_id[sequence_id] for sequence_id in protein_by_id if sequence_id in hit_ids]
        for row in hit_rows:
            row["sample"] = sample_name
            row["source_file"] = str(source_path)
        combined_hit_rows.extend(hit_rows)
        combined_hits_proteins.extend(sample_hit_proteins)
        write_fasta(sample_hit_proteins, sample_dir / f"{sample_name}.hits.faa")
        write_tsv(hit_rows, sample_dir / f"{sample_name}.hits.tsv")

    combined_proteins_path = write_fasta(combined_proteins, root / "all_predicted_proteins.faa")
    candidate_proteins_path = write_fasta(combined_hits_proteins, root / "candidates.protein.faa")
    candidate_nucleotides_path = write_fasta([], root / "candidates.orf.fna")
    hits_tsv_path = write_tsv(combined_hit_rows, root / "candidates.hits.tsv")
    return {
        "all_predicted_proteins": combined_proteins_path,
        "candidate_proteins": candidate_proteins_path,
        "candidate_nucleotides": candidate_nucleotides_path,
        "hits_tsv": hits_tsv_path,
    }
