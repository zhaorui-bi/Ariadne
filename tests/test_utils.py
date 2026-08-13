from pathlib import Path

from ariadne.embed import _display_method, _source_display_name
from ariadne.utils import (
    FastaRecord,
    as_text,
    first_existing,
    is_fasta_path,
    pairwise_identity,
    parse_coverage,
    read_fasta,
    slugify,
    write_fasta,
    write_tsv,
)


def test_as_text_decodes_bytes_and_preserves_str() -> None:
    assert as_text(b"hit_name") == "hit_name"
    assert as_text("already") == "already"


def test_fasta_path_recognition() -> None:
    assert is_fasta_path("coral.fasta")
    assert is_fasta_path("plant.faa")
    assert not is_fasta_path("notes.tsv")
    assert not is_fasta_path("reads.fastq")


def test_fasta_roundtrip_and_empty_header_id(tmp_path: Path) -> None:
    records = [
        FastaRecord(header="seq1 coral", sequence="MKTIIALSYIFCLVFA"),
        FastaRecord(header="", sequence="ACDE"),
    ]
    path = write_fasta(records, tmp_path / "out.faa")
    loaded = read_fasta(path)
    assert loaded[0].id == "seq1"
    assert loaded[1].id == "unnamed"
    assert loaded[0].sequence == "MKTIIALSYIFCLVFA"


def test_tsv_and_path_helpers(tmp_path: Path) -> None:
    missing = tmp_path / "missing.txt"
    present = tmp_path / "present.txt"
    present.write_text("ok", encoding="utf-8")
    assert first_existing(None, missing, present) == present.resolve()
    assert first_existing(missing) is None
    path = write_tsv([{"a": 1, "b": "x"}, {"b": "y", "c": 3}], tmp_path / "table.tsv")
    text = path.read_text(encoding="utf-8")
    assert text.splitlines()[0] == "a\tb\tc"


def test_coverage_slug_and_identity() -> None:
    assert parse_coverage("NODE_1_length_100_cov_12.5") == 12.5
    assert parse_coverage("no coverage here") is None
    assert slugify("CeeSs group!") == "CeeSs_group"
    assert pairwise_identity("ACDE", "ACDE") == 1.0
    assert pairwise_identity("AAAA", "AAAT") == 0.75


def test_display_helpers_match_figure_s2() -> None:
    assert _display_method("lda_reference_subclade_spread") == "LDA"
    assert _display_method("pca_reference_subclade_spread") == "PCA"
    assert _source_display_name("fungal") == "fungi"
    assert _source_display_name("coral") == "coral"
