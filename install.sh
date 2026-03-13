#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-}"
VENV_DIR="${VENV_DIR:-$ROOT_DIR/.venv}"
INSTALL_BIN_DIR="${INSTALL_BIN_DIR:-$HOME/.local/bin}"
PIP_INDEX_URL="${PIP_INDEX_URL:-https://pypi.tuna.tsinghua.edu.cn/simple}"
PIP_TRUSTED_HOST="${PIP_TRUSTED_HOST:-pypi.tuna.tsinghua.edu.cn}"

pick_python() {
  local candidate
  if [[ -n "$PYTHON_BIN" ]]; then
    if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
      echo "Error: PYTHON_BIN='$PYTHON_BIN' was not found." >&2
      exit 1
    fi
    echo "$PYTHON_BIN"
    return 0
  fi

  for candidate in python3.12 python3.11 python3.10 python3.9 python python3; do
    if command -v "$candidate" >/dev/null 2>&1 && python_version_ok "$candidate"; then
      echo "$candidate"
      return 0
    fi
  done

  echo "Error: cannot find a usable Python interpreter." >&2
  echo "Please install Python 3.9+ and rerun, for example:" >&2
  echo "  PYTHON_BIN=python bash install.sh" >&2
  echo "  PYTHON_BIN=python3.9 bash install.sh" >&2
  exit 1
}

python_version_ok() {
  "$1" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 9) else 1)'
}

PYTHON_BIN="$(pick_python)"

if ! python_version_ok "$PYTHON_BIN"; then
  detected_version="$("$PYTHON_BIN" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')"
  echo "Error: Ariadne requires Python 3.9+, but '$PYTHON_BIN' is $detected_version." >&2
  echo "Please rerun with a newer interpreter, for example:" >&2
  echo "  PYTHON_BIN=python bash install.sh" >&2
  echo "  PYTHON_BIN=python3.9 bash install.sh" >&2
  echo "  PYTHON_BIN=python3.10 bash install.sh" >&2
  exit 1
fi

detected_version="$("$PYTHON_BIN" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')"
echo "Using Python: $PYTHON_BIN ($detected_version)"
echo "Using pip index: $PIP_INDEX_URL"

if [[ -d "$VENV_DIR" ]]; then
  echo "[0/4] Removing existing virtual environment: $VENV_DIR"
  rm -rf "$VENV_DIR"
fi

echo "[1/4] Creating virtual environment: $VENV_DIR"
"$PYTHON_BIN" -m venv "$VENV_DIR"

echo "[2/4] Installing Ariadne into the virtual environment"
PIP_INDEX_URL="$PIP_INDEX_URL" PIP_TRUSTED_HOST="$PIP_TRUSTED_HOST" \
  "$VENV_DIR/bin/python" -m pip install --upgrade pip >/dev/null
PIP_INDEX_URL="$PIP_INDEX_URL" PIP_TRUSTED_HOST="$PIP_TRUSTED_HOST" \
  "$VENV_DIR/bin/python" -m pip install -e "$ROOT_DIR"

echo "[3/4] Creating launcher: $INSTALL_BIN_DIR/ariadne"
mkdir -p "$INSTALL_BIN_DIR"
cat > "$INSTALL_BIN_DIR/ariadne" <<EOF
#!/usr/bin/env bash
exec "$VENV_DIR/bin/python" -m ariadne "\$@"
EOF
chmod +x "$INSTALL_BIN_DIR/ariadne"

echo "[4/4] Done"
echo
echo "Launcher installed to: $INSTALL_BIN_DIR/ariadne"
if [[ ":$PATH:" != *":$INSTALL_BIN_DIR:"* ]]; then
  echo "Your PATH does not include $INSTALL_BIN_DIR."
  echo "Add this line to ~/.bashrc or ~/.zshrc, then reload your shell:"
  echo "export PATH=\"$INSTALL_BIN_DIR:\$PATH\""
  echo
  echo "For the current shell only, run:"
  echo "export PATH=\"$INSTALL_BIN_DIR:\$PATH\""
fi
echo
echo "Then test with:"
echo "ariadne --help"
