"""Writable cache folders for numba and matplotlib.

Runs when the package is imported, before numba or matplotlib are loaded. In a container started
with another user id, or a read-only installation, the default cache folders may not be writable:
numba then refuses to compile functions decorated with ``cache=True``, and matplotlib warns. The
first writable choice is used:

numba:      $NUMBA_CACHE_DIR, the package's __pycache__, <user cache>/numba,
            <tmp>/pleiofdr-numba-<uid>; if none is writable, caching is turned off
matplotlib: $MPLCONFIGDIR, <user config>/matplotlib, <tmp>/pleiofdr-matplotlib-<uid>
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path


def writable(path: str | Path) -> bool:
    """True if path exists or can be created, and a file can be written in it."""
    try:
        os.makedirs(path, exist_ok=True)
        tempfile.TemporaryFile(dir=path).close()
    except OSError:
        return False
    return True


def private_tmp_dir(name: str, tmp: str | Path | None = None) -> Path | None:
    """<tmp>/pleiofdr-<name>-<uid>, owned by and only accessible to the current user."""
    uid = os.getuid() if hasattr(os, "getuid") else os.getpid()
    path = Path(tmp or tempfile.gettempdir()) / f"pleiofdr-{name}-{uid}"
    try:
        path.mkdir(mode=0o700, exist_ok=True)
        if hasattr(os, "getuid") and path.stat().st_uid != uid:
            # someone else created it in a shared /tmp; use a fresh private folder instead
            path = Path(tempfile.mkdtemp(prefix=f"pleiofdr-{name}-", dir=tmp))
    except OSError:
        return None
    return path if writable(path) else None


def _user_dir(xdg_var: str, default: str, home: Path) -> Path:
    base = os.environ.get(xdg_var)
    return Path(base) if base else home / default


def configure_numba_cache(package_dir: Path, home: Path, tmp: str | None = None) -> bool:
    """Point NUMBA_CACHE_DIR at a writable folder; return False if caching must be off."""
    chosen = os.environ.get("NUMBA_CACHE_DIR")
    if not (chosen and writable(chosen)):
        chosen = None
        if not writable(package_dir / "__pycache__"):
            user_cache = _user_dir("XDG_CACHE_HOME", ".cache", home) / "numba"
            fallback = private_tmp_dir("numba", tmp)
            chosen = str(user_cache) if writable(user_cache) else fallback and str(fallback)
            if chosen is None:
                os.environ.pop("NUMBA_CACHE_DIR", None)
                return False
    if chosen:
        os.environ["NUMBA_CACHE_DIR"] = chosen
    else:
        os.environ.pop("NUMBA_CACHE_DIR", None)  # numba uses the package's __pycache__
    if "numba" in sys.modules:  # imported before pleiofdr: numba has already read the variable
        from numba.core import config

        config.CACHE_DIR = chosen or ""
    return True


def configure_matplotlib_dir(home: Path, tmp: str | None = None) -> None:
    """Point MPLCONFIGDIR at a writable folder if the configured or default one is not."""
    current = os.environ.get("MPLCONFIGDIR")
    default = _user_dir("XDG_CONFIG_HOME", ".config", home) / "matplotlib"
    if writable(current or default):
        return
    fallback = private_tmp_dir("matplotlib", tmp)
    if fallback is not None:
        os.environ["MPLCONFIGDIR"] = str(fallback)


def _home() -> Path:
    try:
        return Path.home()
    except (KeyError, RuntimeError):  # no HOME and no passwd entry for the uid
        return Path("/nonexistent")


NUMBA_CACHE = configure_numba_cache(Path(__file__).parent, _home())
configure_matplotlib_dir(_home())
