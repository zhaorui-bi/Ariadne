"""Tests for the top-level ``ariadne`` package API (lazy re-exports)."""

from __future__ import annotations

import subprocess
import sys

import pytest

import ariadne


def test_version_and_metadata_exposed():
    assert ariadne.__version__
    assert ariadne.__author__
    assert ariadne.__url__.startswith("http")


def test_all_matches_lazy_export_table():
    public = {name for name in ariadne.__all__ if not name.startswith("__")}
    assert public, "expected a non-empty public API"
    assert public == set(ariadne._LAZY_EXPORTS)


def test_dir_lists_public_api():
    listed = set(ariadne.__dir__())
    assert set(ariadne.__all__) <= listed


def test_public_names_resolve_and_are_callable():
    for name in ariadne._LAZY_EXPORTS:
        obj = getattr(ariadne, name)
        assert callable(obj), f"{name} should resolve to a callable/class"
        # Every export lives in the submodule the table claims.
        assert obj.__module__ == f"ariadne.{ariadne._LAZY_EXPORTS[name]}"


def test_resolved_attribute_is_cached():
    # First access goes through __getattr__; afterwards the name is a real global.
    _ = ariadne.filter_candidates
    assert "filter_candidates" in vars(ariadne)


def test_from_import_forms_work():
    from ariadne import FastaRecord, classify_candidates, filter_candidates

    assert isinstance(FastaRecord, type)
    assert callable(classify_candidates)
    assert callable(filter_candidates)


def test_unknown_attribute_raises_attribute_error():
    with pytest.raises(AttributeError):
        _ = ariadne.definitely_not_a_real_export


def test_bare_import_does_not_eagerly_load_heavy_submodules():
    """``import ariadne`` must stay light: no torch / pyhmmer / pyrodigal pulled in."""
    code = (
        "import sys, ariadne\n"
        "for mod in ('torch', 'pyhmmer', 'pyrodigal',\n"
        "            'ariadne.search', 'ariadne.embed', 'ariadne.model', 'ariadne.tree'):\n"
        "    assert mod not in sys.modules, mod\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "ok"
