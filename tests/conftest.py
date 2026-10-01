from __future__ import annotations

import os
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parent.parent
FIXTURES = ROOT / "tests" / "fixtures"
DEMO_FILES = [
    ROOT / "ref_1kgPhase3eur_LDr2p1_DEMO.mat",
    ROOT / "CTG_COG_2018_DEMO.mat",
    ROOT / "SSGAC_EDU_2016_DEMO.mat",
]

DEMO_MISSING = "demo data missing: tar -xzvf pleioFDR_demo_data.tar.gz in the repository root"
# CI sets PLEIOFDR_REQUIRE_DEMO=1 so that missing demo data fails the run instead of skipping
REQUIRE_DEMO = os.environ.get("PLEIOFDR_REQUIRE_DEMO") == "1"
HAVE_DEMO = all(p.exists() for p in DEMO_FILES)

requires_demo = pytest.mark.skipif(not HAVE_DEMO, reason=DEMO_MISSING)


def pytest_sessionstart(session):
    if REQUIRE_DEMO and not HAVE_DEMO:
        raise pytest.UsageError(f"PLEIOFDR_REQUIRE_DEMO=1 but {DEMO_MISSING}")


@pytest.fixture(scope="session")
def demo_reference():
    from pleiofdr.io import load_reference

    if not HAVE_DEMO:
        pytest.skip(DEMO_MISSING)
    return load_reference(DEMO_FILES[0])
