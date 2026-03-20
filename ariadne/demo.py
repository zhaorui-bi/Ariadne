"""Helpers for creating a small reproducible demo workspace."""

from __future__ import annotations

from pathlib import Path
from typing import Union

from ariadne.fasta_utils import FastaRecord, ensure_directory, read_fasta, ungap, write_fasta
from ariadne.references import prepare_insect_reference, write_reference_metadata

PathLike = Union[str, Path]

CODON_TABLE = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}


def reverse_translate(protein_sequence: str) -> str:
    """Generate a simple synthetic CDS for a protein sequence."""
    return "".join(CODON_TABLE[amino_acid] for amino_acid in protein_sequence if amino_acid in CODON_TABLE) + "TAA"


def _pick_records(records: list[FastaRecord], keyword: str, limit: int) -> list[FastaRecord]:
    """Return up to ``limit`` records whose headers contain the keyword."""
    picked = [record for record in records if keyword.lower() in record.header.lower()]
    return picked[:limit]


def prepare_demo_workspace(
    output_dir: PathLike,
    *,
    coral_alignment_fasta: PathLike,
    insect_xlsx: PathLike,
) -> dict[str, Path]:
    """Assemble a compact demo dataset spanning transcripts, proteins, and references."""
    root = ensure_directory(output_dir)
    transcriptome_dir = ensure_directory(root / "transcriptomes")
    protein_dir = ensure_directory(root / "proteins")
    reference_dir = ensure_directory(root / "references")

    aligned_coral = read_fasta(coral_alignment_fasta, keep_gaps=True)
    cembrene_records = [record for record in aligned_coral if "cembrene" in record.header.lower()]
    non_cembrene_records = [record for record in aligned_coral if "cembrene" not in record.header.lower()]
    if len(cembrene_records) < 3:
        raise ValueError("The coral reference file does not contain enough cembrene-labelled sequences for the demo.")

    preferred_candidate_terms = ["DgTC-2-cembreneA", "AbTC-2-cembreneA", "StTC-1-cembreneC"]
    preferred_candidate = None
    for term in preferred_candidate_terms:
        preferred_candidate = next((record for record in aligned_coral if term.lower() in record.header.lower()), None)
        if preferred_candidate is not None:
            break
    if preferred_candidate is None:
        preferred_candidate = cembrene_records[0]
    full_length_protein = ungap(preferred_candidate.sequence).replace("*", "")
    truncated_protein = full_length_protein[:210]
    decoy_protein = "M" + ("GAVLIQNST" * 18)

    transcript_records = [
        FastaRecord(header="demo_full_cov_25", sequence=reverse_translate(full_length_protein)),
        FastaRecord(header="demo_partial_cov_4", sequence=reverse_translate(truncated_protein)),
        FastaRecord(header="demo_decoy_cov_30", sequence=reverse_translate(decoy_protein)),
    ]
    transcriptome_path = write_fasta(transcript_records, transcriptome_dir / "demo_sample_transcripts.fasta")
    protein_records = [
        FastaRecord(header="demo_full_cov_25", sequence=full_length_protein),
        FastaRecord(header="demo_partial_cov_4", sequence=truncated_protein),
        FastaRecord(header="demo_decoy_cov_30", sequence=decoy_protein),
    ]
    protein_path = write_fasta(protein_records, protein_dir / "demo_sample_proteins.faa")

    preferred_cembrene_terms = [
        "DgTC-2-cembreneA",
        "AbTC-2-cembreneA",
        "StTC-1-cembreneC",
        "NA-141_TC4-cembreneA",
        "NA-047_TC3-cembreneA",
        "CdTC-2-cembrene",
    ]
    preferred_non_cembrene_terms = [
        "S_CdTC-2",
        "CdTC-3",
        "VgTC-2",
        "BsTC-1",
        "ErycarTC-1",
        "S_EsTC-2",
    ]

    selected_cembrene: list[FastaRecord] = []
    for term in preferred_cembrene_terms:
        match = next((record for record in cembrene_records if term.lower() in record.header.lower()), None)
        if match is not None and match not in selected_cembrene:
            selected_cembrene.append(match)
    if len(selected_cembrene) < 6:
        for record in cembrene_records:
            if record not in selected_cembrene:
                selected_cembrene.append(record)
            if len(selected_cembrene) >= 6:
                break

    selected_non_cembrene: list[FastaRecord] = []
    for term in preferred_non_cembrene_terms:
        match = next((record for record in non_cembrene_records if term.lower() in record.header.lower()), None)
        if match is not None and match not in selected_non_cembrene:
            selected_non_cembrene.append(match)
    if len(selected_non_cembrene) < 6:
        for record in non_cembrene_records:
            if record not in selected_non_cembrene:
                selected_non_cembrene.append(record)
            if len(selected_non_cembrene) >= 6:
                break

    coral_demo_records: list[FastaRecord] = []
    for record in (selected_cembrene[:6] + selected_non_cembrene[:6]):
        demo_record = FastaRecord(header=record.header, sequence=ungap(record.sequence).replace("*", ""))
        demo_record.metadata["source"] = "coral"
        demo_record.metadata["label"] = "cembrene" if "cembrene" in record.header.lower() else "other"
        demo_record.metadata["is_cembrene"] = "yes" if "cembrene" in record.header.lower() else "no"
        coral_demo_records.append(demo_record)
    coral_reference_path = write_fasta(coral_demo_records, reference_dir / "coral.fasta")

    insect_reference_path, insect_records = prepare_insect_reference(
        insect_xlsx,
        reference_dir,
        filename="insect.fasta",
        limit=12,
    )
    all_reference_records = coral_demo_records + insect_records
    metadata_path = write_reference_metadata(all_reference_records, reference_dir)

    return {
        "transcriptome": transcriptome_path,
        "protein_fasta": protein_path,
        "protein_dir": protein_dir,
        "reference_dir": reference_dir,
        "coral_reference": coral_reference_path,
        "insect_reference": insect_reference_path,
        "metadata": metadata_path,
    }
