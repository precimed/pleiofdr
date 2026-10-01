"""Per-function tests: hand-made edge cases and MATLAB outputs on synthetic data."""

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

from pleiofdr.ldmatrix import LDMatrix
from pleiofdr.loci import independent_loci, locus_number
from pleiofdr.lookup import cond_fdr, hist3, histc, interp2, lookup_table
from pleiofdr.options import Options, matlab_linspace
from pleiofdr.output import fmt_d, fmt_e, fmt_f, nanmin
from pleiofdr.plotting import filter_points
from pleiofdr.pruning import fast_prune, random_prune_idx
from pleiofdr.smoothing import sparse_smooth_2d
from pleiofdr.stats import binofit, fisher_combined, gc_correct_logp, matlab_prctile

from .conftest import FIXTURES

RTOL = 1e-6


@pytest.fixture(scope="module")
def unit():
    path = FIXTURES / "unit_matlab.mat"
    if not path.exists():
        pytest.skip("unit_matlab.mat missing; run tests/fixtures/make_golden.m")
    data = scipy.io.loadmat(path, squeeze_me=False)
    return {k: v for k, v in data.items() if not k.startswith("__")}


def col(a):
    return np.asarray(a, dtype=float).ravel()


def chain_ld(n: int, blocks: list[list[int]]) -> LDMatrix:
    """Symmetric LD matrix with a full diagonal and complete LD within each block."""
    m = sp.lil_array((n, n), dtype=bool)
    for i in range(n):
        m[i, i] = True
    for block in blocks:
        for a in block:
            for b in block:
                m[a, b] = True
    return LDMatrix.from_sparse(m)


# --- hand-made edge cases -------------------------------------------------------------------


def test_histc_edges():
    edges = np.array([0.0, 1.0, 2.0])
    x = np.array([-0.5, 0, 0.5, 1, 1.5, 2, 2.5, np.nan, np.inf])
    assert histc(x, edges).tolist() == [2, 2, 1]


def test_hist3_matches_matlab_probe():
    # Same probe run in MATLAB R2025b: hist3([y x],'Edges',{[0 1 Inf],[0 1 2 Inf]})
    x = np.array([0, 0.5, 1, 1.5, 2, np.nan, np.inf, -1])
    y = np.array([0, 0, 0, 1, 1, 1, 1, 1.0])
    out = hist3(y, x, np.array([0, 1, np.inf]), np.array([0, 1, 2, np.inf]))
    assert out.tolist() == [[2, 1, 0, 0], [0, 1, 1, 1], [0, 0, 0, 0]]


def test_interp2_matches_matlab_probe():
    # MATLAB: V=[1 2 NaN; 4 5 6]; interp2(V,[1 2 3 1.5 2.5 3],[1 1 1 2 2 2]) -> 1 2 NaN 4.5 5.5 6
    v = np.array([[1, 2, np.nan], [4, 5, 6.0]])
    out = interp2(v, np.array([1, 2, 3, 1.5, 2.5, 3]), np.array([1, 1, 1, 2, 2, 2.0]))
    np.testing.assert_array_equal(out, [1, 2, np.nan, 4.5, 5.5, 6])
    np.testing.assert_array_equal(interp2(v, np.array([1, 2.0]), np.array([2, 2.0])), [4, 5])
    assert np.isnan(interp2(v, np.array([0.5]), np.array([1.0]))[0])


def test_prctile_matches_matlab_probe():
    # MATLAB: prctile([1 2 3 4 NaN]', [0 10 50 90 100]) -> 1 1 2.5 4 4
    out = matlab_prctile(np.array([1, 2, 3, 4, np.nan]), np.array([0, 10, 50, 90, 100]))
    np.testing.assert_allclose(out, [1, 1, 2.5, 4, 4])


def test_binofit_edge_cases():
    phat, pci = binofit(np.array([0, 3, 0]), np.array([10, 3, 0]))
    assert phat[0] == 0 and pci[0, 0] == 0
    assert phat[1] == 1 and pci[1, 1] == 1
    assert np.isnan(phat[2]) and pci[2, 0] == 0 and pci[2, 1] == 1


def test_fast_prune_keeps_strongest_per_block():
    ld = chain_ld(6, [[0, 1, 2], [3, 4]])
    logp = np.array([1.0, 3.0, 2.0, np.nan, 5.0, 0.5])
    out = fast_prune(logp, ld)
    np.testing.assert_array_equal(np.isfinite(out), [False, True, False, False, True, True])
    # NaN inputs stay NaN; ties keep the first index (stable descending sort)
    tie = fast_prune(np.array([2.0, 2.0, 1.0, 0, 0, 0]), ld)
    np.testing.assert_array_equal(np.isfinite(tie)[:3], [True, False, False])


def test_random_prune_respects_defvec_and_repeats():
    ld = chain_ld(8, [[0, 1], [2, 3, 4]])
    defvec = np.array([True] * 7 + [False])
    rng = np.random.default_rng(0)
    for repeats in ("default", "maxout", "none"):
        mask = random_prune_idx(10, ld, defvec, repeats, rng)
        assert mask.shape == (8, 10)
        assert not mask[7].any()
        assert (mask[[0, 1]].sum(axis=0) <= 1).all()
        if repeats == "none":
            assert mask.sum(axis=1).max() <= 1  # no SNP is used in more than one iteration


def test_locus_number_and_loci():
    ld = chain_ld(8, [[0, 1, 2], [2, 3], [5, 6]])
    imat = np.array([True, True, False, True, False, True, True, False])[:, None]
    locus, count = locus_number(imat, ld, np.full(8, 0.3), Options())
    # 0-1 connect, 3 joins through SNP 2 only if 2 passes: it does not, so 3 is separate
    np.testing.assert_array_equal(locus, [1, 1, np.nan, 2, np.nan, 3, 3, np.nan])
    assert count == 3
    empty, n = locus_number(np.zeros((8, 1), bool), ld, np.full(8, 0.3), Options())
    assert np.isnan(n) and np.isnan(empty).all()


def test_independent_loci_with_exclusion_region():
    ld = chain_ld(6, [])
    fdr = np.array([[0.01], [0.02], [0.001], [0.5], [0.03], [np.nan]])
    flp = np.ones((6, 1)) * 5
    opts = Options(fdrthresh=0.05, pthresh=1, exclude_from_fit=[(0, 2)])
    imat, imat2, logfdr = independent_loci(fdr, flp, ld, np.full(6, 0.3), opts)
    np.testing.assert_array_equal(imat[:, 0], [True, True, True, False, True, False])
    # one survivor inside the excluded region (the strongest), others as usual
    np.testing.assert_array_equal(imat2[:, 0], [False, False, True, False, True, False])
    assert np.isnan(logfdr[5, 0])


def test_filter_points():
    x = np.array([0.0, 0.0001, 1.0, np.nan, 0.5])
    y = np.array([0.0, 0.0001, 1.0, 1.0, 0.5])
    fx, fy = filter_points(x, y, (10, 10))
    np.testing.assert_array_equal(fx, [0.0, 1.0, 0.5])
    fx, fy = filter_points(np.array([3.0]), np.array([2.0]), (10, 10))
    assert fx.tolist() == [3.0]
    assert filter_points(x, y, None)[0] is x


def test_matlab_number_formats():
    assert fmt_e(1.5e-5) == "1.500000e-05"
    assert fmt_e(np.nan) == "NaN" and fmt_e(np.inf) == "Inf"
    assert fmt_f(-4.4366314) == "-4.436631"
    assert fmt_d(21) == "21" and fmt_d(True) == "1" and fmt_d(30000.5) == "3.000050e+04"
    assert np.isnan(nanmin(np.array([np.nan, np.nan])))
    assert nanmin(np.array([np.nan, 0.2, 0.1])) == 0.1


def test_matlab_linspace():
    grid = matlab_linspace(0, 30, 3001)
    assert grid[2999] == 29.99  # numpy.linspace gives 29.990000000000002
    assert grid[0] == 0 and grid[-1] == 30 and grid.size == 3001


def test_options_validation():
    with pytest.raises(ValueError):
        Options(stattype="fdr")
    with pytest.raises(ValueError):
        Options(exclude_from_fit=[(5, 5)])
    assert Options(stattype="CONJFDR").stattype == "conjfdr"
    assert Options().hv.size == 301


# --- against MATLAB on synthetic data ---------------------------------------------------------


def test_binofit_vs_matlab(unit):
    phat, pci = binofit(col(unit["binofit_x"]), col(unit["binofit_n"]))
    np.testing.assert_allclose(phat, col(unit["binofit_phat"]), rtol=RTOL)
    np.testing.assert_allclose(pci, unit["binofit_pci"], rtol=RTOL)


def test_histc_vs_matlab(unit):
    counts = histc(col(unit["hist_x"]), matlab_linspace(0, 30, 3001))
    np.testing.assert_array_equal(counts, col(unit["histc_counts"]))


def test_sparse_smooth_vs_matlab(unit):
    out = sparse_smooth_2d(unit["smooth_im"], unit["smooth_wim"], (1e2, 1e2))
    np.testing.assert_allclose(out, unit["smooth_out"], rtol=RTOL, atol=1e-10)


def test_fast_prune_vs_matlab(unit, demo_reference):
    out = fast_prune(col(unit["fastprune_in"]), demo_reference.ld)
    np.testing.assert_array_equal(np.isfinite(out), np.isfinite(col(unit["fastprune_out"])))


def test_lookup_tables_vs_matlab(unit):
    opts = Options(randprune_n=5)
    lp1, lp2 = col(unit["lp1"]), col(unit["lp2"])
    pruneidx = unit["pruneidx"].astype(bool)
    np.testing.assert_allclose(lookup_table(lp1, lp2, opts)[0], unit["look_nopr"], rtol=RTOL)
    look, count, std = lookup_table(lp1, lp2, opts, pruneidx)
    np.testing.assert_allclose(look, unit["look_pr"], rtol=RTOL)
    np.testing.assert_array_equal(count, unit["lookcount_pr"])
    np.testing.assert_allclose(std, np.real(unit["lookstd_pr"]), rtol=1e-4, atol=1e-7)
    fdrvec, _, fdrvec0 = cond_fdr(lp1, lp2, opts, pruneidx, np.zeros(lp1.size, bool))
    np.testing.assert_allclose(fdrvec, col(unit["cond_fdrvec"]), rtol=RTOL)
    np.testing.assert_allclose(fdrvec0, col(unit["cond_fdrvec0"]), rtol=RTOL)


def test_genomic_control_vs_matlab(unit):
    lp1 = col(unit["lp1"])
    ivec0 = col(unit["ivec0"]).astype(bool)
    pruneidx = unit["pruneidx"].astype(bool)
    for key, std, prune in [
        ("gc_inhouse", False, None),
        ("gc_std", True, None),
        ("gc_pruned", False, pruneidx),
    ]:
        logp, sig0 = gc_correct_logp(lp1, ivec0, std, prune)
        np.testing.assert_allclose(logp, col(unit[key]), rtol=RTOL)
        np.testing.assert_allclose(sig0, col(unit[f"{key}_sig0"]), rtol=RTOL)


def test_fisher_vs_matlab(unit):
    out = fisher_combined(col(unit["lp1"]), col(unit["lp2"])[:, None])
    np.testing.assert_allclose(out[:, 0], col(unit["fisher"]), rtol=RTOL)
