"""Cache-folder selection for numba and matplotlib (pleiofdr/_caches.py)."""

import os
import stat

import pytest

from pleiofdr import _caches


@pytest.fixture
def readonly(tmp_path):
    folder = tmp_path / "readonly"
    folder.mkdir()
    folder.chmod(stat.S_IRUSR | stat.S_IXUSR)
    yield folder
    folder.chmod(stat.S_IRWXU)


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for var in ("NUMBA_CACHE_DIR", "MPLCONFIGDIR", "XDG_CACHE_HOME", "XDG_CONFIG_HOME"):
        monkeypatch.delenv(var, raising=False)


pytestmark = pytest.mark.skipif(
    hasattr(os, "geteuid") and os.geteuid() == 0, reason="root can write to read-only folders"
)


def test_keeps_writable_numba_cache_dir(tmp_path, monkeypatch, readonly):
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "mine"))
    assert _caches.configure_numba_cache(readonly, readonly, str(tmp_path))
    assert os.environ["NUMBA_CACHE_DIR"] == str(tmp_path / "mine")


def test_writable_package_uses_numba_default(tmp_path, monkeypatch, readonly):
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(readonly / "cache"))  # set but unusable
    assert _caches.configure_numba_cache(tmp_path / "pkg", readonly, str(tmp_path))
    assert "NUMBA_CACHE_DIR" not in os.environ


def test_readonly_install_without_home_falls_back_to_private_tmp(tmp_path, readonly):
    # the situation of `docker run -u UID:GID`: read-only site-packages, HOME=/ not writable
    assert _caches.configure_numba_cache(readonly, readonly, str(tmp_path))
    chosen = os.environ["NUMBA_CACHE_DIR"]
    assert chosen == str(tmp_path / f"pleiofdr-numba-{os.getuid()}")
    assert stat.S_IMODE(os.stat(chosen).st_mode) == 0o700


def test_readonly_install_uses_writable_home_cache(tmp_path, readonly):
    assert _caches.configure_numba_cache(readonly, tmp_path / "home", str(readonly))
    assert os.environ["NUMBA_CACHE_DIR"] == str(tmp_path / "home" / ".cache" / "numba")


def test_nothing_writable_turns_caching_off(readonly):
    assert _caches.configure_numba_cache(readonly, readonly, str(readonly)) is False
    assert "NUMBA_CACHE_DIR" not in os.environ


def test_matplotlib_dir(tmp_path, monkeypatch, readonly):
    monkeypatch.setenv("MPLCONFIGDIR", str(readonly / "mpl"))
    _caches.configure_matplotlib_dir(readonly, str(tmp_path))
    assert os.environ["MPLCONFIGDIR"] == str(tmp_path / f"pleiofdr-matplotlib-{os.getuid()}")

    monkeypatch.delenv("MPLCONFIGDIR")
    _caches.configure_matplotlib_dir(tmp_path / "home", str(tmp_path))
    assert "MPLCONFIGDIR" not in os.environ  # the default ~/.config/matplotlib is writable


def test_private_tmp_dir_not_owned_by_us(tmp_path, monkeypatch):
    # pretend to be uid 4242; the folder below then belongs to "someone else" (the real uid)
    monkeypatch.setattr(os, "getuid", lambda: 4242)
    planted = tmp_path / "pleiofdr-numba-4242"
    planted.mkdir()
    chosen = _caches.private_tmp_dir("numba", str(tmp_path))
    assert chosen is not None
    assert chosen != planted and chosen.parent == tmp_path
    assert chosen.name.startswith("pleiofdr-numba-")
