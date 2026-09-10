#!/usr/bin/env bash
# Create a Python virtualenv with pairtools 1.1.3 for differential tests.
#
# pairtools 1.1.3 is not on PyPI; it is built from the GitHub tag. Two small
# patches are applied so that it runs on current Python/Cython:
#   * `dedup_cython.pyx`: a legacy code path iterates over memoryviews,
#     which Cython 3 rejects (wrapped in np.asarray; semantics unchanged);
#   * a `pipes` shim (the stdlib module was removed in Python 3.13) that
#     implements the tiny subset pairtools.lib.fileio uses.
#
# Usage: setup_pairtools_oracle.sh VENV_DIR
set -euo pipefail
VENV="${1:?venv dir}"
PY="${PYTHON:-python3}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
"$PY" -m venv "$VENV"
"$VENV/bin/pip" install -q --upgrade pip setuptools wheel
"$VENV/bin/pip" install -q cython numpy "pandas<3" scipy pysam pyyaml click
curl -sL https://github.com/open2c/pairtools/archive/refs/tags/v1.1.3.tar.gz -o "$WORK/pt.tgz"
tar xzf "$WORK/pt.tgz" -C "$WORK"
SRC="$WORK/pairtools-1.1.3"
"$PY" - "$SRC/pairtools/lib/dedup_cython.pyx" <<'PYEOF'
import sys
p = sys.argv[1]
s = open(p).read()
s = s.replace("for ar in [self.c1, self.c2, self.p1, self.p2, self.s1, self.s2]:",
              "for ar in [np.asarray(self.c1), np.asarray(self.c2), np.asarray(self.p1), np.asarray(self.p2), np.asarray(self.s1), np.asarray(self.s2)]:")
open(p, "w").write(s)
PYEOF
SITE="$("$VENV/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
cat > "$SITE/pipes.py" <<'PYEOF'
"""Minimal shim of the removed `pipes` module for pairtools 1.1.3 on Python 3.13+."""
import shlex
import subprocess


class _ProcFile:
    def __init__(self, proc, stream):
        self._proc, self._stream, self.closed = proc, stream, False

    def __getattr__(self, name):
        return getattr(self._stream, name)

    def __iter__(self):
        return iter(self._stream)

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def close(self):
        if not self.closed:
            self.closed = True
            self._stream.close()
            self._proc.wait()


class Template:
    def __init__(self):
        self.steps = []

    def append(self, cmd, kind):
        self.steps.append((cmd, kind))

    def open(self, file, mode):
        cmd = " | ".join(c for c, _ in self.steps)
        if mode == "w":
            proc = subprocess.Popen(f"{cmd} > {shlex.quote(file)}", shell=True, stdin=subprocess.PIPE, text=True, bufsize=1 << 16)
            return _ProcFile(proc, proc.stdin)
        if mode == "r":
            proc = subprocess.Popen(f"{cmd} < {shlex.quote(file)}", shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 16)
            return _ProcFile(proc, proc.stdout)
        raise ValueError(mode)
PYEOF
"$VENV/bin/pip" install -q --no-build-isolation "$SRC"
"$VENV/bin/pairtools" --version
