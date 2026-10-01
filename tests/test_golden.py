"""End-to-end parity with MATLAB on the demo data (fixtures from tests/fixtures/make_golden.m).

MATLAB's prune indices are passed through randprune_file, so both implementations see
the same random pruning and results must agree to numerical precision.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import scipy.io

from pleiofdr.cli import run_from_config

from .conftest import FIXTURES, ROOT, requires_demo

RTOL = 1e-6
ATOL = 1e-12  # values that are 0 in one implementation and ~1e-17 in the other

BASE = """\
reffile={root}/ref_1kgPhase3eur_LDr2p1_DEMO.mat
traitfolder={root}
traitfile1=CTG_COG_2018_DEMO.mat
traitname1=COGchr21
traitfiles={{'SSGAC_EDU_2016_DEMO.mat'}}
traitnames={{'EDUchr21'}}
randprune=true
randprune_n=20
exclude_chr_pos=[]
manh_colorlist=[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]
reset_pruneidx=true
randprune_repeats=default
pthresh=1
perform_gc=true
use_standard_gc=false
randprune_gc=true
exclude_from_discovery=false
mafthresh=0.005
exclude_ambiguous_snps=true
onscreen=false
exit_matlab_upon_completion=false
"""

RUNS = {
    "conjfdr": ["stattype=conjfdr", "fdrthresh=0.05"],
    "condfdr": ["stattype=condfdr", "fdrthresh=0.01"],
    "condfdr_excl": [
        "stattype=condfdr",
        "fdrthresh=0.05",
        "exclude_chr_pos=[21 30000000 32000000]",
        "exclude_from_discovery=true",
        "use_standard_gc=true",
        "randprune_gc=false",
        "exclude_ambiguous_snps=false",
        "randprune_repeats=maxout",
    ],
}


def _col(a) -> np.ndarray:
    return np.asarray(a, dtype=float).ravel()


@pytest.fixture(scope="module", params=list(RUNS))
def run(request, tmp_path_factory):
    name = request.param
    golden_path = FIXTURES / f"golden_{name}.mat"
    if not golden_path.exists():
        pytest.skip(f"{golden_path.name} missing; run tests/fixtures/make_golden.m")
    golden = scipy.io.loadmat(golden_path)
    tmp = tmp_path_factory.mktemp(name)
    prune_file = tmp / "pruneidx.mat"
    scipy.io.savemat(prune_file, {"pruneidx": golden["pruneidx"].astype(bool)})
    config = tmp / "config.txt"
    lines = [
        BASE.format(root=ROOT),
        *RUNS[name],
        f"outputdir={tmp / 'out'}",
        f"randprune_file={prune_file}",
    ]
    config.write_text("\n".join(lines) + "\n")
    output = run_from_config(config)
    return name, golden, output, tmp / "out"


@requires_demo
def test_arrays_match_matlab(run):
    name, g, out, _ = run
    r = out.results
    np.testing.assert_array_equal(out.pruneidx, g["pruneidx"].astype(bool))
    np.testing.assert_array_equal(out.excludevec, _col(g["excludevec"]).astype(bool))
    for ours, key in [
        (r.logpvec, "logpvec1"),
        (r.logpmat2, "logpmat2"),
        (r.zvec, "zvec1"),
        (r.zmat2, "zmat2"),
        (r.flp, "flp"),
        (r.fdrmat, "fdrmat"),
        (out.fdrmat12, "fdrmat12"),
        (out.fdrmat21, "fdrmat21"),
        (r.fdrvec0, "fdrvec0"),
        (r.logfdrmat, "logfdrmat"),
    ]:
        np.testing.assert_allclose(
            _col(ours), _col(g[key]), rtol=RTOL, atol=ATOL, err_msg=f"{name}:{key}"
        )
    np.testing.assert_array_equal(r.imat.ravel(), _col(g["imat"]).astype(bool))
    np.testing.assert_array_equal(r.imat2.ravel(), _col(g["imat2"]).astype(bool))
    np.testing.assert_array_equal(r.locusnumvec, _col(g["locusnumvec"]))
    assert out.locusnum == float(g["locusnum"].squeeze())


@requires_demo
def test_lookup_tables_match_matlab(run):
    name, g, out, _ = run
    for ours, key in [(out.lookup12, "lookup12"), (out.lookup21, "lookup21")]:
        cells = g[key]
        if cells.size == 0:
            assert ours == []
            continue
        theirs = cells[0, 0]
        np.testing.assert_allclose(
            ours[0][0], theirs[0, 0], rtol=RTOL, atol=ATOL, err_msg=f"{name}:{key}"
        )
        np.testing.assert_array_equal(ours[0][1], theirs[0, 1])


@requires_demo
def test_qq_and_enrichment_match_matlab(run):
    name, g, out, _ = run
    for key in ["mat_qq", "mat_qq_inv", "mat_tdr", "mat_tdr_inv", "mat_enrich", "mat_enrich_inv"]:
        np.testing.assert_allclose(
            getattr(out, key), g[key], rtol=RTOL, atol=ATOL, equal_nan=True, err_msg=f"{name}:{key}"
        )


def _read_csv(path: Path) -> list[list[str]]:
    with open(path, newline="") as f:
        return list(csv.reader(f))


def _cells_equal(a: str, b: str) -> bool:
    if a == b:
        return True
    try:
        return bool(np.isclose(float(a), float(b), rtol=RTOL, atol=0))
    except ValueError:
        return False


@requires_demo
def test_csv_tables_match_matlab(run):
    name, _, _, outdir = run
    expected = sorted((FIXTURES / f"matlab_{name}").glob("*.csv"))
    assert expected
    for ref in expected:
        ours = outdir / ref.name
        assert ours.exists(), f"missing {ref.name}"
        a, b = _read_csv(ours), _read_csv(ref)
        assert a[0] == b[0], f"{ref.name}: header differs"
        assert len(a) == len(b), f"{ref.name}: {len(a)} rows vs {len(b)}"
        for ra, rb in zip(a[1:], b[1:], strict=True):
            assert len(ra) == len(rb)
            assert all(_cells_equal(x, y) for x, y in zip(ra, rb, strict=True)), (ref.name, ra, rb)


@requires_demo
def test_result_mat_and_figures_written(run):
    name, g, _, outdir = run
    result = scipy.io.loadmat(outdir / "result.mat")
    np.testing.assert_allclose(_col(result["fdrmat"]), _col(g["fdrmat"]), rtol=RTOL)
    assert result["qq_t1breaks"].shape == (1, 3001)
    expected_png = {p.name for p in (FIXTURES / f"matlab_{name}").glob("*.png")}
    written = {p.name for p in outdir.glob("*.png")}
    assert expected_png <= written
