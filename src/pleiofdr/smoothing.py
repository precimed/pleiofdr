"""Weighted 2-D smoothing with a second-difference penalty (port of SparseSmooth2d.m)."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


def _second_difference(n: int) -> sp.csr_array:
    """MATLAB diff(speye(n), 2): (n-2) x n rows of [1 -2 1]."""
    return sp.diags_array([1.0, -2.0, 1.0], offsets=[0, 1, 2], shape=(n - 2, n), format="csr")


def sparse_smooth_2d(im: np.ndarray, wim: np.ndarray, lamvec: tuple[float, float]) -> np.ndarray:
    """Minimise sum(w * (y - im)^2) + lam * (|D_rows y|^2 + |D_cols y|^2) over the image y.

    The image is vectorised column-major, as im(:) in MATLAB. SparseSmooth2d.m uses lamvec(1)
    for both directions; that is kept.
    """
    nr, nc = im.shape
    l1 = sp.kron(sp.eye_array(nc), _second_difference(nr))  # along rows, within each column
    l2 = sp.kron(_second_difference(nc), sp.eye_array(nr))  # along columns, within each row
    lam = np.sqrt(lamvec[0])
    ell = sp.vstack([lam * l1, lam * l2]).tocsc()
    ltl = (ell.T @ ell).tocsc()

    wvec = wim.ravel(order="F")
    yvec = im.ravel(order="F")
    h = sp.diags_array(np.fmax(0.001, wvec)) + ltl
    # One Newton step from y = 0: g = w .* (0 - y) + LtL * 0
    g = -wvec * yvec
    ysm = -spsolve(h.tocsc(), g)
    return ysm.reshape((nr, nc), order="F")
