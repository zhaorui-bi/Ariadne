"""Centralised logging and terminal output for the Ariadne pipeline.

Usage
-----
In every module::

    import logging
    logger = logging.getLogger(__name__)

At program entry point::

    from ariadne import log
    log.setup_logging(verbose=args.verbose)
    log.print_banner(__version__)
"""
from __future__ import annotations

import logging
import sys

# ---------------------------------------------------------------------------
# ANSI colour helpers
# ---------------------------------------------------------------------------


def _supports_colour() -> bool:
    """Return True when *stderr* is a colour-capable TTY."""
    return hasattr(sys.stderr, "isatty") and sys.stderr.isatty()


class _Colours:
    """Container for ANSI escape sequences; empty strings when unsupported."""

    def __init__(self) -> None:
        on = _supports_colour()
        self.RESET   = "\033[0m"   if on else ""
        self.BOLD    = "\033[1m"   if on else ""
        self.DIM     = "\033[2m"   if on else ""
        self.RED     = "\033[91m"  if on else ""
        self.YELLOW  = "\033[93m"  if on else ""
        self.GREEN   = "\033[92m"  if on else ""
        self.CYAN    = "\033[96m"  if on else ""
        self.BLUE    = "\033[94m"  if on else ""
        self.MAGENTA = "\033[95m"  if on else ""


# Singleton evaluated once at import time.
C = _Colours()

# ---------------------------------------------------------------------------
# ASCII banner
# ---------------------------------------------------------------------------

# "ARIADNE" rendered in the "Big" figlet style.
_BANNER_ART = (
    r"    _    ____  ___    _    ____  _   _ _____  ",
    r"   / \  |  _ \|_ _|  / \  |  _ \| \ | | ____|",
    r"  / _ \ | |_) || |  / _ \ | | | |  \| |  _|  ",
    r" / ___ \|  _ < | | / ___ \| |_| | |\  | |___ ",
    r"/_/   \_\_| \_\___/_/   \_\____/|_| \_|_____|",
)

_SUBTITLE = "Terpene Synthase Discovery & Annotation Pipeline"

# Total visible width of the box, including the border characters (╭ and ╮).
_BOX_TOTAL = 62


def _box_rule() -> str:
    """Return the horizontal rule string for the top/bottom box border."""
    return "─" * (_BOX_TOTAL - 2)


def _box_line(content: str = "", visible_len: int | None = None) -> str:
    """Return one │-bounded box row, right-padded to *_BOX_TOTAL* width.

    Parameters
    ----------
    content:
        The (potentially ANSI-coloured) string to embed.
    visible_len:
        Explicit *visible* character count when *content* contains ANSI codes
        that inflate ``len(content)``.  Defaults to ``len(content)``.
    """
    if visible_len is None:
        visible_len = len(content)
    # -4 accounts for "│ " prefix and " │" suffix.
    padding = " " * max(0, _BOX_TOTAL - 4 - visible_len)
    return f"│ {content}{padding} │"


def _render_banner(version: str) -> str:
    """Compose the full banner string."""
    rule = _box_rule()
    rows: list[str] = [
        f"╭{rule}╮",
        _box_line(),
    ]
    for art_line in _BANNER_ART:
        coloured = f"  {C.CYAN}{art_line}{C.RESET}"
        rows.append(_box_line(coloured, visible_len=2 + len(art_line)))

    rows.append(_box_line())

    subtitle_coloured = f"  {C.DIM}{_SUBTITLE}{C.RESET}"
    rows.append(_box_line(subtitle_coloured, visible_len=2 + len(_SUBTITLE)))

    ver_str = f"version {version}"
    ver_coloured = f"  {C.BOLD}{ver_str}{C.RESET}"
    rows.append(_box_line(ver_coloured, visible_len=2 + len(ver_str)))

    rows.append(_box_line())
    rows.append(f"╰{rule}╯")
    return "\n".join(rows)


def print_banner(version: str) -> None:
    """Print the Ariadne banner to *stderr*."""
    print(_render_banner(version), file=sys.stderr)
    print(file=sys.stderr)


# ---------------------------------------------------------------------------
# Coloured log formatter
# ---------------------------------------------------------------------------


class _ColourFormatter(logging.Formatter):
    """Log formatter that adds ANSI colours to level names and module names."""

    _LEVEL_LABEL: dict[int, str] = {}

    def __init__(self) -> None:
        super().__init__()
        self._LEVEL_LABEL = {
            logging.DEBUG:    C.DIM    + "DEBUG" + C.RESET,
            logging.INFO:     C.CYAN   + "INFO " + C.RESET,
            logging.WARNING:  C.YELLOW + "WARN " + C.RESET,
            logging.ERROR:    C.RED    + "ERROR" + C.RESET,
            logging.CRITICAL: C.BOLD + C.RED + "CRIT " + C.RESET,
        }

    def format(self, record: logging.LogRecord) -> str:
        level  = self._LEVEL_LABEL.get(record.levelno, record.levelname)
        module = C.MAGENTA + record.name.split(".")[-1] + C.RESET
        return f"[{level}] [{module}] {record.getMessage()}"


# ---------------------------------------------------------------------------
# Public setup
# ---------------------------------------------------------------------------


def setup_logging(verbose: bool = False) -> None:
    """Configure the *ariadne* logger hierarchy.

    Parameters
    ----------
    verbose:
        When *True*, emit DEBUG messages; otherwise INFO and above only.
    """
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_ColourFormatter())

    root = logging.getLogger("ariadne")
    root.setLevel(level)
    root.handlers.clear()
    root.addHandler(handler)
    root.propagate = False


logger = logging.getLogger(__name__)
