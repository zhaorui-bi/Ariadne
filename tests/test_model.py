"""Tests for ariadne.model — ESM helpers and TPS workbook loading.

Tests that require torch/transformers are skipped automatically when those
optional dependencies are not installed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ariadne.model import (
    CEESS_POSITIVE_LABELS,
    _collapse_labels_to_ceess_group,
    _color_map,
    _type_output_stem,
    _type_probability_column,
    _type_rawscore_column,
    esm_model_help_text,
    load_tps_xlsx,
    resolve_esm_model_name,
)


class TestResolveEsmModelName:
    @pytest.mark.parametrize(
        "preset,expected",
        [
            ("650M", "facebook/esm2_t33_650M_UR50D"),
            ("150M", "facebook/esm2_t30_150M_UR50D"),
            ("3B", "facebook/esm2_t36_3B_UR50D"),
            ("504M", "facebook/esm2_t33_650M_UR50D"),
        ],
    )
    def test_known_presets(self, preset, expected):
        assert resolve_esm_model_name(preset) == expected

    def test_passthrough_for_unknown(self):
        assert resolve_esm_model_name("my/custom-model") == "my/custom-model"

    def test_help_text_mentions_default(self):
        assert "650M" in esm_model_help_text()


class TestColumnNaming:
    def test_probability_column(self):
        assert _type_probability_column("cembrene A") == "esm_type_probability_cembrene_a"

    def test_rawscore_column(self):
        assert _type_rawscore_column("Sesquiterpenes") == "esm_type_rawscore_sesquiterpenes"

    def test_output_stem(self):
        assert _type_output_stem("klysimplexin R") == "klysimplexin_r"
        assert _type_output_stem("###") == "unknown"


class TestCollapseLabels:
    def test_collapse(self):
        labels = ["cembrene A", "Sesquiterpenes", "cembrene B"]
        out = _collapse_labels_to_ceess_group(labels, set(CEESS_POSITIVE_LABELS))
        assert out == ["CeeSs", "non-CeeSs", "CeeSs"]


class TestColorMap:
    def test_stable_and_complete(self):
        labels = ["b", "a", "c", "a"]
        mapping = _color_map(labels)
        assert set(mapping) == {"a", "b", "c"}
        # deterministic: sorted labels mapped in palette order
        assert mapping == _color_map(["c", "b", "a"])


class TestLoadTpsXlsx:
    def test_loads_bundled_workbook(self, repo_root: Path):
        xlsx = repo_root / "TPS" / "TPS.xlsx"
        if not xlsx.exists():
            pytest.skip("bundled TPS.xlsx not available")
        records = load_tps_xlsx(xlsx)
        assert len(records) > 0
        # every record carries a CeeSs grouping derived from its type label
        groups = {r.ceess_group for r in records}
        assert groups <= {"CeeSs", "non-CeeSs"}
        # cembrene types must be grouped as CeeSs-positive
        for record in records:
            if record.label in CEESS_POSITIVE_LABELS:
                assert record.ceess_group == "CeeSs"


class TestTorchDependentSmoke:
    def test_mlp_classifier_importable_only_with_torch(self):
        torch = pytest.importorskip("torch")
        import numpy as np

        from ariadne.model import _TorchMLPClassifier

        x = np.random.RandomState(0).rand(12, 8).astype("float32")
        y = ["CeeSs"] * 6 + ["non-CeeSs"] * 6
        clf = _TorchMLPClassifier(epochs=3, hidden_dim=8, batch_size=4, random_state=0)
        clf.fit(x, y)
        proba = clf.predict_proba(x)
        assert proba.shape == (12, 2)
        assert torch.is_tensor(torch.tensor(proba))
