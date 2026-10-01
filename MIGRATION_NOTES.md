# MATLAB → Python migration notes

pleioFDR was ported from MATLAB/Octave to a Python ≥ 3.13 package. This document maps each
MATLAB function to its replacement, lists the MATLAB behaviours that had to be reproduced
exactly, and records every intentional difference. The original code is in `legacy_matlab/`.

## Package layout

```
pyproject.toml                 requires-python = ">=3.13", console script `pleiofdr`
src/pleiofdr/
    __main__.py                python -m pleiofdr
    cli.py                     runme.m: config → Options, load inputs, run
    config.py                  TextConfig.m + runme.m defaults; safe MATLAB-literal parser
    options.py                 pleioOpt.m as a dataclass; MATLAB-compatible linspace
    io.py                      .mat loading (v5 via scipy.io, v7.3 via h5py), load_gwas.m,
                               refinfo (.ref), randprune_file, result.mat
    ldmatrix.py                LDMatrix: CSC pattern with int64 column pointers, int32 row indices
    pruning.py                 FastPrune.m, random_prune_idx_amd.m, random_prune_idx_amd_fb.m
    stats.py                   binofit, GCcorrect_logpvec.m, fisher_comStats.m, sample overlap
    smoothing.py               SparseSmooth2d.m
    lookup.py                  lookup_table.m, conj_lookup_table.m, cond_FDR_amd.m,
                               MATLAB-compatible interp2 / histc / hist3
    fdr.py                     pleioFDR_amd.m
    loci.py                    ind_loci_idx.m, locusnumber.m
    analysis.py                pleiotropy_analysis.m (log file written like MATLAB's diary)
    output.py                  save_to_csv.m, save_fdr.m, save_zscore.m
    plotting/                  save_figure.m, filter_points_for_plotting.m, plot_qq_amd.m,
                               plot_enrichment_amd.m, plot_lookup.m, plot_Manhattan.m
    fuma/combine.py            fuma/*.R (cond / conj-lead / conj-snps)
    fuma/novelty.py            fuma/conj_fuma_combined_novelty.py
tests/                         unit, golden-parity, CLI, fuma and toolkit tests
tests/fixtures/make_golden.m   produces the MATLAB reference outputs
legacy_matlab/                 the original .m files
```

## Dependencies

| Package    | Why |
|------------|-----|
| numpy      | arrays |
| scipy      | `scipy.io` (v5 .mat), `scipy.sparse` + `spsolve` (smoothing), `scipy.stats` / `scipy.special` |
| pandas     | `.ref` file, fuma tables |
| matplotlib | figures (Agg backend unless `onscreen=true`) |
| h5py       | v7.3 .mat files; the full `ref9545380_1kgPhase3eur_LDr2p1.mat` is v7.3 |
| numba      | FastPrune and locusnumber are greedy, inherently sequential walks over up to ~9.5M SNPs (× `randprune_n` for pruning). They cannot be vectorised, and plain Python loops would be orders of magnitude slower. |

Optional extra `ref`: `intervaltree` for `ref4pleioFDR/toolkit/sLDSC_scripts.py`. Dev group: pytest, ruff.

## Function map

| MATLAB | Python | Notes |
|---|---|---|
| runme.m | `cli.main`, `cli.run_from_config` | same config keys and defaults |
| TextConfig.m | `config.TextConfig` | `get_cell`/`get_mat` parse literals; no `eval` |
| pleioOpt.m | `options.Options` | validation kept; exclusion regions stored 0-based inclusive |
| pleiotropy_analysis.m | `analysis.run_analysis` | |
| load_gwas.m | `io.load_gwas` | first variable starting with `logp` / `z` (case-insensitive), in file order |
| pleioFDR_amd.m | `fdr.pleio_fdr` | |
| cond_FDR_amd.m | `lookup.cond_fdr` | |
| lookup_table.m | `lookup.lookup_table`, `lookup.ind_look` | |
| conj_lookup_table.m | `lookup.conj_lookup_table` | |
| SparseSmooth2d.m | `smoothing.sparse_smooth_2d` | penalty built with Kronecker products |
| binofit (Statistics Toolbox) / binofit_wrap.m | `stats.binofit` | Clopper–Pearson via the F distribution, as `statbinoci` |
| GCcorrect_logpvec.m | `stats.gc_correct_logp` | |
| fisher_comStats.m | `stats.fisher_combined` | |
| check/correct_sample_overlap.m | `stats.check_sample_overlap` / `stats.correct_sample_overlap` | |
| FastPrune.m | `pruning.fast_prune` | numba kernel |
| random_prune_idx_amd(_fb).m | `pruning.random_prune_idx` | `repeats` = default / maxout / none |
| ind_loci_idx.m | `loci.independent_loci` | |
| locusnumber.m | `loci.locus_number` | numba kernel |
| save_to_csv / save_fdr / save_zscore.m | `output.*` | same file names, headers and number formats |
| save_figure.m, filter_points_for_plotting.m | `plotting.save_figure`, `plotting.filter_points` | |
| plot_qq_amd.m, plot_enrichment_amd.m | `plotting.qq.plot_qq`, `plotting.qq.plot_enrichment`, `plotting.qq.qq_matrix` | |
| plot_lookup.m, plot_Manhattan.m | `plotting.lookup.plot_lookup`, `plotting.manhattan.plot_manhattan` | |
| is_octave.m, binofit_dale.m | — | Octave-only; dropped |
| suplabel.m, plot_qq_annot.m | — | never called by the pipeline; not ported |

## MATLAB behaviours reproduced

Each of these was checked against MATLAB (R2025b) with a probe or the golden fixtures.

- **`linspace`** computes `d1 + k*(d2-d1)/(n-1)`, multiplying before dividing. numpy multiplies by a precomputed step, which moves some grid points by one ulp (29.99 vs 29.990000000000002) and with them the histogram bin edges. `options.matlab_linspace` is used for every grid.
- **`histc` / `hist3`**: bin k holds `edges(k) <= x < edges(k+1)`, the last bin holds `x == edges(end)`, and values outside the edges and NaN are ignored.
- **`interp2`** (via griddedInterpolant) skips corners with zero weight, so a NaN neighbour does not spoil a query exactly on a grid line (`interp2([1 2 NaN; 4 5 6], 2, 1)` is 2).
- **NaN in `min`/`max`**: MATLAB ignores NaN. `min(1, max(0, pim./qim))` therefore maps NaN to 0, and in `cond_FDR_amd` an undefined p-value maps to grid index 1 (`fdrvec0` is not NaN-masked, `fdrvec` is). `np.fmin`/`np.fmax` are used throughout. The conjFDR `max(fdrmat12, fdrmat21)` also ignores NaN.
- **`prctile`** uses midpoint interpolation (numpy `method="hazen"`) and ignores NaN.
- **`sort(..., 'descend')`** in FastPrune is stable and puts NaN first (NaN and Inf are then dropped). The port uses a stable argsort over the finite entries.
- **Column-major order** in `SparseSmooth2d` (`im(:)`) and in the `binofit` reshapes.
- **`binofit`** follows `statbinoci`: F-distribution quantiles, lower bound 0 when x = 0, upper bound 1 when x = n.
- **Linear indexing quirk**: in `ind_loci_idx.m`, `logfdrmat(mafmask) = NaN` with an nsnp-long mask only touches the **first** column. The port does the same.
- **Output formats**: file names use `%g` of `fdrthresh`, and CSV cells use C `%e`/`%f`/`%d` with MATLAB's `NaN`/`Inf` spelling (`%d` of a non-integer falls back to `%e`).
- **Config parsing**: `TextConfig` keeps only the text between the first and second `=`, `get_bool` treats any non-zero number as true, and unknown keys are an error.

## Quirks preserved for parity (candidate bugs, not fixed)

1. `fisher_comStats.m` GC-corrects only when `~exist('ivec0','var')`, which never happens, so the Fisher statistic is **never** GC-corrected and `pruneidx_fc` has no effect.
2. In `plot_qq_amd.m` / `plot_enrichment_amd.m`, `qqci_max`/`qqci_min` are built from the already-updated `qqmat`. They only affect the dashed lines drawn when `show_ci=true`, which the port does not draw.
3. `mat_qq*` and `mat_enrich*` in `result.mat` hold the matrices of the **last** conditioning trait only.
4. `SparseSmooth2d` uses `lamvec(1)` for both directions (both entries are `smf` anyway).
5. `plot_Manhattan.m` has `showgenes = false` hard-coded, so gene labels are never drawn. The label code is not ported, and `manh_fontsize_genenames`/`manh_yspace`/`manh_ymargin` have no effect, as before.
6. `conj_lookup_table.m` builds both tables as trait 1 | trait 2. It is only reached when lookup tables are missing, which never happens in the pipeline.

## Intentional differences

- **Figures** are PNG and SVG. MATLAB's `.fig` format is not produced; base names are unchanged.
- **Random numbers**: prune indices are drawn with `numpy.random.Generator` (`--seed` makes runs reproducible), so they differ from MATLAB's `unifrnd` stream. Exact comparisons with MATLAB pass MATLAB's prune indices through `randprune_file`.
- **`randprune=false`**: `pleiotropy_analysis.m` errors because `pruneidx_fc` is undefined; the port runs without pruning.
- **Octave** branches are removed. `exit_matlab_upon_completion` and `mlibrary` are accepted and ignored.
- **Exclusion regions** that match no SNP raise a clear error (MATLAB failed with an indexing error).
- **fuma/combine** (ported R scripts): with current data.table (1.18), `fread` parses `True`/`False` as logical, and the scripts' `is_locus_lead == "True"` then selects **no rows**, so they write empty tables. The port accepts both spellings, as intended. Everything else, including column order after `merge`, R's three-valued `&`/`ifelse` NA logic and `write.table` number formatting, matches R.
- **fuma/novelty** accepts tab- as well as comma-separated tables. The original read comma-separated input only, although the combine step writes tab-separated files. Numbers are parsed exactly (`float_precision="round_trip"`).

## Verification

- **Reference outputs**: MATLAB R2025b with the Statistics Toolbox (`tests/fixtures/make_golden.m`) on the demo data (chr21, 133,222 SNPs, 20 pruning iterations) for three configurations: conjFDR; condFDR; and condFDR with an exclusion region, discovery exclusion, standard GC, `randprune_gc=false` and `randprune_repeats=maxout`. Plus per-function outputs on synthetic data and small v7.3 files.
- **Parity** (`tests/test_golden.py`), with MATLAB's prune indices: GC-corrected p- and z-values, Fisher statistics, cond/conjFDR (both directions), unconditional FDR, lookup tables and counts, all QQ/TDR/enrichment matrices, `imat`/`imat2`/locus numbers and all CSV tables agree to rtol 1e-6 (atol 1e-12 for −log10 FDR values that are exactly 0 in MATLAB). In practice the conjFDR CSV tables are byte-identical.
- **fuma** (`tests/test_fuma.py`): outputs of the original R scripts (R 4.6.1, data.table 1.18.6; `tests/fixtures/fuma/make_expected.sh`) and of the original novelty script, on synthetic inputs.
- **Containers** (`.github/workflows/container.yml`): the Docker and Apptainer images run the synthetic self-test with networking disabled, and the conjFDR demo inside each container reproduces the MATLAB CSV tables byte for byte.
- **Toolkit** (`tests/test_toolkit.py`): the annotation pipeline runs end to end on a tiny synthetic genome.
- **Python**: the full suite passes on Python 3.14.7 (the newest release available) and on Python 3.13.15 (the minimum, installed with `uv sync --python 3.13`), with numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, matplotlib 3.11.2, numba 0.67.0 and h5py 3.16.0.
- **Performance** on the demo conjFDR run (Apple Silicon, 16 GB): MATLAB R2025b `runme` 57 s and 1.2 GB peak memory; `pleiofdr` 8.8 s wall time and 0.6 GB peak memory, including all figures.
- **Not verified**: a run on the full 9.5M-SNP reference. Its LD matrix needs about 12 GB in the port (2.9e9 int32 row indices) versus about 26 GB in MATLAB, more than the 16 GB development machine has for either. The v7.3 reader is tested on a MATLAB-written v7.3 subset.
