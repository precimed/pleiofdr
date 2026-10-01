import numpy as np
import pytest
import scipy.io

from pleiofdr.io import load_gwas, load_mat, load_pruneidx, load_reference, save_result_mat

from .conftest import FIXTURES


@pytest.fixture(scope="module")
def v73_files():
    paths = [FIXTURES / n for n in ("ref_v73.mat", "ref_v5.mat", "trait_v73.mat")]
    if not all(p.exists() for p in paths):
        pytest.skip("v7.3 fixtures missing; run tests/fixtures/make_golden.m")
    return paths


def test_v73_reference_matches_v5(v73_files):
    v73, v5 = load_reference(v73_files[0]), load_reference(v73_files[1])
    assert v73.nsnp == v5.nsnp == 3000
    np.testing.assert_array_equal(v73.ld.indptr, v5.ld.indptr)
    np.testing.assert_array_equal(v73.ld.indices, v5.ld.indices)
    assert v73.ld.indices.dtype == np.int32
    for name in ("chrnumvec", "posvec", "mafvec", "is_intergenic", "is_ambiguous"):
        np.testing.assert_array_equal(getattr(v73, name), getattr(v5, name), err_msg=name)


def test_v73_trait_file(v73_files):
    logp, z = load_gwas([v73_files[2]], 3000)
    assert logp.shape == (3000, 1) and z.shape == (3000, 1)
    raw = load_mat(v73_files[2])
    np.testing.assert_array_equal(logp[:, 0], np.ravel(raw["logpvec"]))
    with pytest.raises(ValueError, match="expected 10 snps"):
        load_gwas(v73_files[2], 10)


def test_load_gwas_dummy_zscore(tmp_path):
    path = tmp_path / "trait.mat"
    scipy.io.savemat(path, {"LOGPVEC": np.array([[1.0], [0.0], [np.nan]])})
    with pytest.raises(ValueError, match="must contain variable z"):
        load_gwas(path, 3)
    logp, z = load_gwas(path, 3, dummy_zscore=True)
    np.testing.assert_allclose(z[:2, 0], [1.6448536269514722, 0.0], atol=1e-12)
    assert np.isnan(z[2, 0])


def test_pruneidx_and_result_roundtrip(tmp_path):
    mask = np.random.default_rng(0).random((50, 4)) > 0.5
    scipy.io.savemat(tmp_path / "prune.mat", {"anything": mask})
    np.testing.assert_array_equal(load_pruneidx(tmp_path / "prune.mat"), mask)
    save_result_mat(tmp_path / "result.mat", fdrmat=np.arange(5.0), excludevec=mask[:, 0])
    back = scipy.io.loadmat(tmp_path / "result.mat")
    assert back["fdrmat"].shape == (5, 1)
    np.testing.assert_array_equal(back["excludevec"].ravel().astype(bool), mask[:, 0])
