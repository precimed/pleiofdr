"""Binary LD matrix in compressed-sparse-column form.

The full reference has ~2.9e9 nonzeros; scipy.sparse would promote both index arrays to int64
(~26 GB). Row indices always fit in int32 (≤ ~9.5M SNPs), so only the column pointers are int64.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class LDMatrix:
    indptr: np.ndarray  # int64, length n + 1
    indices: np.ndarray  # int32 row indices of nonzero entries, column by column
    n: int

    def __post_init__(self) -> None:
        if self.indptr.shape != (self.n + 1,):
            raise ValueError("LDMatrix: indptr must have length n + 1")

    @classmethod
    def from_sparse(cls, matrix: sp.spmatrix | sp.sparray) -> LDMatrix:
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("LD matrix must be square")
        csc = sp.csc_array(matrix)
        csc.eliminate_zeros()
        csc.sort_indices()
        return cls(
            indptr=csc.indptr.astype(np.int64),
            indices=csc.indices.astype(np.int32),
            n=csc.shape[0],
        )

    @property
    def shape(self) -> tuple[int, int]:
        return (self.n, self.n)

    @property
    def nnz(self) -> int:
        return int(self.indptr[-1])

    def column(self, j: int) -> np.ndarray:
        return self.indices[self.indptr[j] : self.indptr[j + 1]]

    def to_sparse(self) -> sp.csc_array:
        data = np.ones(self.indices.shape, dtype=bool)
        return sp.csc_array((data, self.indices, self.indptr), shape=self.shape)
