"""Tests for ariadne.cli — argument parsing and the program entry point."""

from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from ariadne.cli import (
    _first_existing_path,
    _parse_extra_reference,
    build_parser,
    cmd_run,
    main,
)


class TestParser:
    def test_build_parser_returns_parser(self):
        assert isinstance(build_parser(), argparse.ArgumentParser)

    def test_run_subcommand_dispatches_to_cmd_run(self):
        parser = build_parser()
        args = parser.parse_args(
            ["run", "--protein-folder", "input", "--reference-dir", "tree", "--output-dir", "out"]
        )
        assert args.func is cmd_run
        assert args.skip_phylogeny is False

    def test_filter_defaults(self):
        parser = build_parser()
        args = parser.parse_args(["filter", "--input-fasta", "x.faa", "--output-dir", "out"])
        assert args.min_length == 300
        assert args.identity_threshold == 0.95


class TestHelpers:
    def test_parse_extra_reference(self):
        assert _parse_extra_reference("plant=path/to.fasta") == ("plant", "path/to.fasta")

    def test_parse_extra_reference_requires_equals(self):
        with pytest.raises(ValueError):
            _parse_extra_reference("noequalshere")

    def test_first_existing_path(self, tmp_path: Path):
        present = tmp_path / "f.txt"
        present.write_text("x")
        assert _first_existing_path(tmp_path / "missing", present) == present.resolve()
        assert _first_existing_path(tmp_path / "a", tmp_path / "b") is None


class TestEntryPoint:
    def test_version_exits_zero(self):
        with pytest.raises(SystemExit) as excinfo:
            main(["--version"])
        assert excinfo.value.code == 0

    def test_no_subcommand_returns_zero(self):
        assert main([]) == 0
