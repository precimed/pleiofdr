"""Reading MATLAB .mat inputs (v5 and v7.3) and writing result.mat."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from .ldmatrix import LDMatrix
from .stats import logp_to_absz

_H5_CHUNK = 1 << 26  # entries read per chunk when streaming sparse arrays from HDF5


def _is_hdf5(path: Path) -> bool:
    return h5py.is_hdf5(path)


def _h5_value(obj: h5py.Dataset | h5py.Group) -> np.ndarray | sp.csc_array:
    if isinstance(obj, h5py.Group):
        if "MATLAB_sparse" not in obj.attrs:
            raise ValueError(f"unsupported MATLAB v7.3 variable: {obj.name}")
        nrows = int(obj.attrs["MATLAB_sparse"])
        jc = obj["jc"][()].astype(np.int64)
        ir = obj["ir"][()].astype(np.int64) if "ir" in obj else np.zeros(0, np.int64)
        data = obj["data"][()] if "data" in obj else np.zeros(0)
        return sp.csc_array((data, ir, jc), shape=(nrows, jc.size - 1))
    value = obj[()]
    if value.ndim >= 2:
        value = value.T  # HDF5 stores MATLAB arrays in row-major order
    cls = obj.attrs.get("MATLAB_class", b"")
    if cls == b"logical":
        value = value.astype(bool)
    elif cls == b"char":
        value = "".join(map(chr, np.ravel(value, order="F")))
    return value


def load_mat(path: str | Path, variables: list[str] | None = None) -> dict[str, object]:
    """Load variables of a .mat file in file order; v7.3 (HDF5) files are read with h5py."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Data file {path} does not exist")
    if _is_hdf5(path):
        with h5py.File(path, "r") as f:
            names = [k for k in f if not k.startswith("#")]
            return {k: _h5_value(f[k]) for k in names if variables is None or k in variables}
    data = scipy.io.loadmat(path, variable_names=variables, spmatrix=False)
    return {k: v for k, v in data.items() if not k.startswith("__")}


def _h5_ldmatrix(group: h5py.Group) -> LDMatrix:
    """Stream a MATLAB sparse logical matrix into an LDMatrix without int64 row indices.

    MATLAB sparse matrices never store explicit zeros, so the data array is not read.
    """
    n = int(group.attrs["MATLAB_sparse"])
    jc = group["jc"][()].astype(np.int64)
    nnz = int(jc[-1])
    ir_ds = group["ir"]
    indices = np.empty(nnz, dtype=np.int32)
    for start in range(0, nnz, _H5_CHUNK):
        stop = min(start + _H5_CHUNK, nnz)
        indices[start:stop] = ir_ds[start:stop]
    return LDMatrix(indptr=jc, indices=indices, n=n)


@dataclass
class Reference:
    ld: LDMatrix
    chrnumvec: np.ndarray
    posvec: np.ndarray
    mafvec: np.ndarray
    is_intergenic: np.ndarray
    is_ambiguous: np.ndarray

    @property
    def nsnp(self) -> int:
        return self.chrnumvec.size


def load_reference(path: str | Path) -> Reference:
    """Load LDmat, chrnumvec, posvec, mafvec, is_intergenic, is_ambiguous from the reference."""
    path = Path(path)
    names = ["chrnumvec", "posvec", "mafvec", "is_intergenic", "is_ambiguous"]
    if _is_hdf5(path):
        with h5py.File(path, "r") as f:
            missing = [k for k in ["LDmat", "mafvec", "is_intergenic"] if k not in f]
            if missing:
                raise ValueError(f"error loading {', '.join(missing)} from {path}")
            ld = _h5_ldmatrix(f["LDmat"])
            vectors = {k: _h5_value(f[k]) for k in names if k in f}
    else:
        data = load_mat(path, ["LDmat", *names])
        missing = [k for k in ["LDmat", "mafvec", "is_intergenic"] if k not in data]
        if missing:
            raise ValueError(f"error loading {', '.join(missing)} from {path}")
        ld = LDMatrix.from_sparse(data["LDmat"])
        vectors = {k: data[k] for k in names if k in data}

    def vec(name: str, dtype: type, default: object = None) -> np.ndarray:
        if name not in vectors:
            if default is None:
                raise ValueError(f"error loading {name} from {path}")
            return np.full(ld.n, default, dtype=dtype)
        return np.asarray(vectors[name]).ravel().astype(dtype)

    return Reference(
        ld=ld,
        chrnumvec=vec("chrnumvec", np.int64),
        posvec=vec("posvec", np.int64),
        mafvec=vec("mafvec", float),
        is_intergenic=vec("is_intergenic", bool),
        is_ambiguous=vec("is_ambiguous", bool, False),
    )


def load_gwas(
    traitfiles: str | Path | list[str | Path], nsnps: int, dummy_zscore: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """Read the first logp* and z* variables of each trait file into (nsnps, nfiles) matrices."""
    files = [traitfiles] if isinstance(traitfiles, (str, Path)) else list(traitfiles)
    logpmat = np.empty((nsnps, len(files)))
    zmat = np.empty((nsnps, len(files)))
    for i, file in enumerate(files):
        print(f"\n   {file}... ", end="")
        traits = load_mat(file)
        fields = list(traits)
        lpfield = next((k for k in fields if k.lower().startswith("logp")), None)
        if lpfield is None:
            raise ValueError(f"Data file {file} must contain variable logp*")
        print(f"{lpfield}... ", end="")
        lp = np.asarray(traits[lpfield], dtype=float)
        zfield = next((k for k in fields if k.lower().startswith("z")), None)
        if zfield is None and not dummy_zscore:
            raise ValueError(f"Data file {file} must contain variable z*")
        if zfield is None:
            z = logp_to_absz(lp)
        else:
            print(f"{zfield}... ", end="")
            z = np.asarray(traits[zfield], dtype=float)
        if min(lp.shape, default=1) > 1 or min(z.shape, default=1) > 1:
            raise ValueError("logp and z must be vectors, not matrices")
        lp, z = lp.ravel(), z.ravel()
        if lp.size != nsnps:
            raise ValueError(f"{file} logpvec has {lp.size} snps - expected {nsnps} snps")
        if z.size != nsnps:
            raise ValueError(f"{file} zvec has {z.size} snps - expected {nsnps} snps")
        logpmat[:, i] = lp
        zmat[:, i] = z
    return logpmat, zmat


def load_refinfo(path: str | Path, nsnps: int) -> tuple[list[str], list[str], list[str]]:
    """SNP, A1 and A2 columns of a tab-separated reference info file (e.g. 9545380.ref)."""
    df = pd.read_csv(path, sep="\t", usecols=["SNP", "A1", "A2"], dtype=str, keep_default_na=False)
    if len(df) != nsnps:
        raise ValueError(f"{path} has {len(df)} rows - expected {nsnps}")
    return df["SNP"].tolist(), df["A1"].tolist(), df["A2"].tolist()


def load_pruneidx(path: str | Path) -> np.ndarray:
    """Prune indices from the first variable of a .mat file (randprune_file)."""
    data = load_mat(path)
    if not data:
        raise ValueError(f"{path} contains no variables")
    value = next(iter(data.values()))
    if sp.issparse(value):
        value = value.toarray()
    return np.asarray(value).astype(bool)


def save_result_mat(path: str | Path, **variables: np.ndarray) -> None:
    """Uncompressed MATLAB v5 file (save -v6); 1-D arrays become column vectors."""
    scipy.io.savemat(path, variables, do_compression=False, oned_as="column")
