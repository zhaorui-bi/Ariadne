"""Allow ``python -m ariadne`` to behave like the ``ariadne`` CLI entrypoint."""

from ariadne.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
