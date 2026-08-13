"""Module entry point for ``python -m ariadne``."""

from __future__ import annotations

from ariadne.cli import main

if __name__ == "__main__":  # pragma: no cover - exercised by the interpreter
    raise SystemExit(main())
