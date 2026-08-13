from ariadne.cli import main


def test_missing_subcommand_exits_cleanly() -> None:
    assert main([]) == 0


def test_expected_errors_return_one() -> None:
    assert main(["filter", "--input-fasta", "does-not-exist.faa", "--output-dir", "tmp_out"]) == 1
