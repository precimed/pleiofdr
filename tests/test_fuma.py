"""pleiofdr.fuma against outputs of the original R scripts and novelty script."""

import math

import numpy as np
import pandas as pd
import pytest

from pleiofdr.fuma.combine import combine, is_true, output_name, r_format_number, write_r_table
from pleiofdr.fuma.novelty import novelty_check

from .conftest import FIXTURES

F = FIXTURES / "fuma"


def _cells_match(ours: str, theirs: str) -> bool:
    if ours == theirs:
        return True
    try:
        a, b = float(ours), float(theirs)
    except ValueError:
        return False
    # data.table::fread is not always correctly rounded; allow a last-digit difference
    return math.isclose(a, b, rel_tol=1e-14)


@pytest.mark.parametrize(
    ("mode", "clump"), [("cond", "lead.csv"), ("conj-lead", "lead.csv"), ("conj-snps", "snps.csv")]
)
def test_combine_matches_r(tmp_path, mode, clump):
    table = combine(
        mode,
        F / clump,
        F / "fuma_snps.txt",
        F / "sumstats1.txt.gz",
        F / "sumstats2.txt",
        "T1",
        "T2",
    )
    name = output_name(mode, "T1", "T2")
    write_r_table(table, tmp_path / name)
    ours = (tmp_path / name).read_text().splitlines()
    theirs = (F / f"expected_{name}").read_text().splitlines()
    assert ours[0] == theirs[0]
    assert len(ours) == len(theirs)
    for a, b in zip(ours[1:], theirs[1:], strict=True):
        ca, cb = a.split("\t"), b.split("\t")
        assert len(ca) == len(cb)
        assert all(_cells_match(x, y) for x, y in zip(ca, cb, strict=True)), (a, b)


def test_novelty_matches_original_script():
    result = novelty_check(
        F / "expected_conjFDR_0.05_T1_vs_T2_lead.csv",
        F / "gwascatalog.txt",
        F / "expected_conjFDR_0.05_T1_vs_T2_snps.csv",
        F / "novelty_db.txt",
        "Depress",
        "T1",
    )
    expected = pd.read_csv(
        F / "expected_conjFDR_0.05_T1_vs_T2_novelty.csv", float_precision="round_trip"
    )
    pd.testing.assert_frame_equal(result, expected, rtol=1e-13, atol=0, check_dtype=False)


def test_is_locus_lead_accepts_strings_and_booleans():
    values = pd.Series(["True", "False", True, False, None, "x"], dtype=object)
    assert is_true(values).tolist() == [True, False, True, False, False, False]


@pytest.mark.parametrize(
    ("x", "text"),
    [
        # expected strings printed by write.table in R 4.6 (aarch64)
        (-0.054773587481354, "-0.054773587481354"),
        (1 / 3, "0.333333333333333"),
        (100000.0, "1e+05"),
        (0.0001, "1e-04"),
        (0.001, "0.001"),
        (2.0**60, "1152921504606846976"),
        (0.007037804558528706, "0.00703780455852871"),
        (99999.99999999999, "1e+05"),
        (-2.5e-120, "-2.5e-120"),
        (np.inf, "Inf"),
        (0.0, "0"),
    ],
)
def test_r_number_format(x, text):
    assert r_format_number(x) == text
