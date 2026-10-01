[![tests](https://github.com/precimed/pleiofdr/actions/workflows/tests.yml/badge.svg)](https://github.com/precimed/pleiofdr/actions/workflows/tests.yml)

## Contents

- [Contents](#contents)
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Install pleioFDR](#install-pleiofdr)
- [Containers (incl. airgapped systems)](#containers-incl-airgapped-systems)
- [Data downloads](#data-downloads)
- [Data preparation](#data-preparation)
- [Run pleioFDR](#run-pleiofdr)
  - [Running w. Slurm on HPC (TSD)](#running-w-slurm-on-hpc-tsd)
- [Configuration reference](#configuration-reference)
- [pleioFDR results](#pleiofdr-results)
- [FUMA-defined loci](#fuma-defined-loci)
- [Development and testing](#development-and-testing)
- [MATLAB version](#matlab-version)

## Introduction

Pleiotropy-informed conditional and conjunctional false discovery rate allows to boost loci discovery in low-powered GWAS by levereging pleiotropic enrichment with a larger GWAS on related phenotype, and to identify genetic loci joinly associated with two phenotypes.

If you use pleioFDR software for your research publication, please cite the following paper(s):
* Andreassen, O.A. et al. Improved detection of common variants associated with schizophrenia and bipolar disorder using pleiotropy-informed conditional false discovery rate. PLoS Genet 9, e1003455 (2013). 

For an introduction about pleioFDR, please see
* Olav B Smeland et al, Discovery of shared genomic loci using the conditional false discovery rate approach, Hum Genet. 2020 Jan;139(1)

The pleioFDR software may not be used in medical applications.

pleioFDR is a Python package (Python 3.13 or newer). It was ported from the original MATLAB
implementation, and gives the same results, see [MATLAB version](#matlab-version).

## Quick Start

To install and run pleioFDR on a small example, constrained to chromosome 21:
```
git clone https://github.com/precimed/pleiofdr && cd pleiofdr
pip install .
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/pleioFDR_demo_data.tar.gz
tar -xzvf pleioFDR_demo_data.tar.gz
pleiofdr --config config.txt
```

To install and run pleioFDR using full example:
```
git clone https://github.com/precimed/pleiofdr && cd pleiofdr
pip install .
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/CTG_COG_2018.mat
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/SSGAC_EDU_2016.mat
cp config_default.txt config.txt
pleiofdr --config config.txt
```

For the description of the data, see [here](https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt).
For the results, inspect the ``results`` folder.

## Install pleioFDR

Prerequisites:
 - Python 3.13 or newer
 - for the full reference, a workstation with at least 16GB of RAM
   (the LD matrix alone takes about 12GB in memory)

Get the code, either with the "Code" → "Download ZIP" button on https://github.com/precimed/pleiofdr, or from the command line:
  ```
  git clone https://github.com/precimed/pleiofdr && cd pleiofdr
  ```

Install it into a virtual environment with [uv](https://docs.astral.sh/uv/):
  ```
  uv sync          # creates .venv with pleiofdr and its dependencies
  uv run pleiofdr --help
  ```
or with pip:
  ```
  python3 -m venv .venv && source .venv/bin/activate
  pip install .
  pleiofdr --help
  ```

The dependencies are numpy, scipy, pandas, matplotlib, h5py and numba; they are installed
automatically. The optional extra `pip install .[ref]` adds `intervaltree`, used by some scripts
in `ref4pleioFDR/toolkit`.

## Containers (incl. airgapped systems)

pleioFDR is also distributed as a Docker image and as an Apptainer (Singularity) image for x86_64
Linux. Both contain Python and every dependency, so they run without network access; reference and
trait data stay outside the image and are bind-mounted at run time. Tagged releases publish
``ghcr.io/precimed/pleiofdr:X.Y.Z`` and attach ``pleiofdr-X.Y.Z.sif`` (with a ``.sha256`` checksum) and
``pleiofdr-X.Y.Z-docker.tar.gz`` to the [GitHub release](https://github.com/precimed/pleiofdr/releases).

**Getting the image into an airgapped environment.** On a machine with internet access, download
``pleiofdr-X.Y.Z.sif`` and ``pleiofdr-X.Y.Z.sif.sha256`` from the release page (or build the image
yourself with ``apptainer build pleiofdr.sif Apptainer.def``), check it with
``sha256sum -c pleiofdr-X.Y.Z.sif.sha256``, and copy the ``.sif`` file in through your usual file import
route, together with the reference data (``ref9545380_1kgPhase3eur_LDr2p1.mat``, and ``9545380.ref`` if used).

**Running with Apptainer.** Paths in ``config.txt`` are paths *inside* the container, so bind the data
folders to fixed locations:
```
# a config template to edit
apptainer exec pleiofdr.sif cp /opt/pleiofdr/share/config_default.txt config.txt
#   reffile=/ref/ref9545380_1kgPhase3eur_LDr2p1.mat
#   traitfolder=/data/traits
#   outputdir=/data/results

apptainer run --containall -B /path/to/refdata:/ref -B $PWD:/data --pwd /data \
    pleiofdr.sif --config config.txt

apptainer test pleiofdr.sif                              # self-test on synthetic data, no data needed
apptainer run --app fuma-combine pleiofdr.sif --help     # FUMA helpers, see fuma/readme.md
apptainer run --app fuma-novelty pleiofdr.sif --help
apptainer run-help pleiofdr.sif
```
On a SLURM cluster, put the ``apptainer run`` line in the job script and request enough memory
(about 16GB for the full reference).

**Running with Docker or Podman.** Load the image from the release tarball (``docker load < pleiofdr-X.Y.Z-docker.tar.gz``)
or pull ``ghcr.io/precimed/pleiofdr:X.Y.Z``, then
```
docker run --rm --network none -u $(id -u):$(id -g) \
    -v /path/to/refdata:/ref -v $PWD:/data ghcr.io/precimed/pleiofdr:X.Y.Z --config config.txt
```

**Building.** ``docker build --platform linux/amd64 -t pleiofdr .`` or ``apptainer build pleiofdr.sif Apptainer.def``
(add ``--fakeroot`` to build without root; both need network access while building). Both install exactly the dependency versions in ``uv.lock`` on the
same pinned ``python:3.13-slim-bookworm`` base image.

**Building the Docker container imane on macOS (arm64).** First you need to install Docker Desktop (https://docs.docker.com/desktop/install/mac-install/). Then build the image using the following command:
```{shell}
docker buildx build -t pleiofdr .
# running it
docker run --rm --network none \                 
  -v /path/to/refdata:/ref \
  -v $PWD:/data \
  docker.io/library/pleiofdr --config config.txt
```

**Building the Apptainer image on macOS (arm64).** First you need to install Apptainer (https://apptainer.org/docs/admin/main/installation.html#mac). The following was tested using Lima via Homebrew.
First change permissions to allow writing to local directories:

```{shell}
limactl stop apptainer
limactl edit apptainer
```

The "mounts" section should contain:

```{yaml}
mounts:                                     
- location: "~"
  writable: true  
```

Then start the Lima VM and build the image inside it and run pleioFDR with the reference and trait data bind-mounted:

```{shell}
limactl start apptainer
limactl shell apptainer
# inside the Lima VM, build and run the image:
apptainer build pleiofdr-2.0.0-arm64.sif Apptainer.def
apptainer run --containall -B /path/to/refdata:/ref -B $PWD:/data --pwd /data pleiofdr-2.0.0-arm64.sif --config config.txt
```

## Data downloads

Download reference data from [here](https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr). 
The reference is based on 1000 Genomes phase 3 data (May 2, 2013 release).
Variant calls (vcf files) for 22 autosomes were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502 .
We kept only samples of European ancestry (IBS, TSI, GBR, CEU, FIN populations) 
with missing call rate below 10% and only biallelic variants with non-duplicated ids, 
minor allele frequency above 1%, missing call rate below 10% and Hardy-Weinberg equilibrium
exact test p-values greater than 1.E-20. 
The filtering was performed with PLINK 1.9. Resulted template contained 503 samples and 9,545,380 variants.
Further details are available in [about.txt](https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt).

  ```
  wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt
  wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
  wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/CTG_COG_2018.mat
  wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/SSGAC_EDU_2016.mat
  wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_bfile.tar.gz
  wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/9545380.ref
  ```

Those at NORMENT with access to NIRD can also download these data from ``SUMSTAT/misc/9545380_ref`` and ``SUMSTAT/TMP/mat_9545380``.

The reference and trait files are MATLAB ``.mat`` files; pleioFDR reads both the older format
(v5/v7) and the HDF5-based v7.3 format (used by ``ref9545380_1kgPhase3eur_LDr2p1.mat``).

## Data preparation

Here we explain how to convert raw summary statistics to pleioFDR format.
Feel free to skip this step if you would like to try pleioFDR on ``CTG_COG_2018.mat`` and ``SSGAC_EDU_2016.mat``,
or if you downloaded input data from the internal NORMENT ``SUMSTATS`` inventory.

A trait file is a ``.mat`` file with a ``logpvec`` variable (-log10 p-values, one per reference
variant, in reference order) and a ``zvec`` variable (z-scores). The first variables whose names
start with ``logp`` and ``z`` are used. If z-scores are unavailable, see ``dummy_zscore`` below.

The conversion uses the separate [python_convert](https://github.com/precimed/python_convert) repository;
check its README for the Python version it requires.

Downloads:
 - Download code from https://github.com/precimed/python_convert, either via web browser, or ``git clone https://github.com/precimed/python_convert``.
 - Download educational attainment and subjective well-being summary statistics
from SSGAC consortium to traitfolder:
   ```
   wget http://ssgac.org/documents/EduYears_Main.txt.gz -P traitfolder
   wget http://ssgac.org/documents/SWB_Full.txt.gz -P traitfolder
   ```

Conversion steps:
  - Use sumstats.py script to standartizize downloaded summary statistics (csv)
and prepare input files for cond/conj fdr analysis (mat):
    ```
    python src/converter/sumstats.py csv --auto --sumstats traitfolder/EduYears_Main.txt.gz  --n-val 328917 --out traitfolder/ssgac.edu.csv --force
    python src/converter/sumstats.py csv --auto --sumstats traitfolder/SWB_Full.txt.gz --n-val 298420 --out traitfolder/ssgac.swb.csv --force
    python src/converter/sumstats.py mat --sumstats traitfolder/ssgac.edu.csv --ref 9545380.ref --out traitfolder/ssgac.edu.mat
    python src/converter/sumstats.py mat --sumstats traitfolder/ssgac.swb.csv --ref 9545380.ref --out traitfolder/ssgac.swb.mat
    ```
    In the first and second commands --n-val argument indicates sample size. The number is taken from original papers [Okbay et al. (2016)].
  - For more details on input arguments please check:
    ```
    python src/converter/sumstats.py --help
    python src/converter/sumstats.py csv --help
    python src/converter/sumstats.py mat --help
    ```
 
## Run pleioFDR

  Create a configuration file by copying ``config_default.txt`` file, located in the root of pleioFDR repository.
  ```
  cp config_default.txt config.txt
  ```
  
  Edit ``config.txt`` so that 
  * ``reffile`` points to the ``ref9545380_1kgPhase3eur_LDr2p1.mat`` file
  * ``traitfolder`` points to folder containing ``CTG_COG_2018.mat`` and ``SSGAC_EDU_2016.mat``
  * set ``randprune_n=500`` instead of the default ``randprune_n=20``
  You may also want to change ``traitfile1`` and ``traitfiles`` options.
  
  Then run
  ```
  pleiofdr --config config.txt
  ```
  (``python -m pleiofdr --config config.txt`` works too; without ``--config``, ``config.txt`` in the
  current folder is used.)

  Random pruning draws random numbers; add ``--seed 123`` to make a run reproducible. To reuse a
  fixed set of prune indices (for example exported from an earlier run), save them as a boolean
  ``nsnp × randprune_n`` matrix in a ``.mat`` file and set ``randprune_file`` in the config.

  To run many trait pairs, see ``run_batch.py`` (fills ``config_template.txt`` for each pair) and
  ``MultipleRunUtility.sh``.

### Running w. Slurm on HPC (TSD)

On HPC systems with Slurm (e.g., TSD at University of Oslo), you can run pleioFDR in a batch job. For example, create a file ``pleiofdr.sh`` with the following header:

```bash
#!/bin/bash
#SBATCH --job-name=pleiofdr
#SBATCH --output=pleiofdr.out
#SBATCH --error=pleiofdr.err
#SBATCH --open-mode=truncate
#SBATCH --account=<account_name>
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
```

Then, for invocation with Apptainer (make sure `apptainer` exist in `$PATH`), add the following lines to the script:

```bash
module purge

export SIF=<path_to>/pleiofdr-2.0.0.sif
export MYCONFIG=<path_to>/config.txt
export APPTAINERENV_HDF5_USE_FILE_LOCKING=FALSE  # needed on TSD
apptainer run --cleanenv \
  -B /ess:/ess -B /cluster:/cluster \
  -B $PWD:/data --pwd /data \
  $SIF --config $MYCONFIG
```

Then, the job can be submitted with `sbatch pleiofdr.sh`. The output and error will be written to `pleiofdr.out` and `pleiofdr.err`, respectively.
Output files will be written to the folder specified in the config file.

Similarly, the .venv virtual environment can be used in a Slurm job script, by activating it before running pleioFDR.
First make sure that the virtual environment is created locally and activated in the current shell. 

```bash
python -m venv .venv
source .venv/bin/activate
pip install <path_to>/pleiofdr/.[ref]
```

Then, include in the Slurm job script:

```bash
module purge
module load Python/3.13.5-GCCcore-14.3.0

export MYCONFIG=<path_to>/config.txt
export HDF5_USE_FILE_LOCKING=FALSE  # needed on TSD

source .venv/bin/activate
pleiofdr --config $MYCONFIG
```

## Configuration reference

  Configuration files contain ``key=value`` lines; lines starting with ``#`` are comments. Values
  use MATLAB-style literals: ``true``/``false``, numbers, matrices such as ``[6 25119106 33854733; 8 7200000 12500000]``
  and lists of strings such as ``{'a.mat', 'b.mat'}``. Configuration files written for the MATLAB version work unchanged.

  | Key | Default | Meaning |
  |---|---|---|
  | ``reffile`` | ``ref9545380_1kgPhase3eur_LDr2p1.mat`` | reference with ``LDmat``, ``chrnumvec``, ``posvec``, ``mafvec``, ``is_intergenic``, ``is_ambiguous`` |
  | ``refinfo`` | (empty) | optional tab-separated file with ``SNP``, ``A1``, ``A2`` columns (e.g. ``9545380.ref``), used in the output tables |
  | ``traitfolder`` | ``../example_data_for_pleiotropy`` | folder with the trait files; if empty, trait file names must be full paths |
  | ``traitfile1``, ``traitname1`` | ``PGC2_SCZ.mat``, ``SCZ`` | primary trait |
  | ``traitfiles``, ``traitnames`` | ``{'COG_charge.mat'}``, ``{'COGNITION'}`` | trait(s) to condition on |
  | ``outputdir`` | ``test`` | output folder |
  | ``stattype`` | ``conjfdr`` | ``condfdr`` (conditional) or ``conjfdr`` (conjunctional) |
  | ``fdrthresh`` | ``0.05`` | FDR threshold; 0.01 is recommended for condfdr and 0.05 for conjfdr |
  | ``pthresh`` | ``1`` | threshold on Fisher's combined p-value for the loci table (1 disables it) |
  | ``randprune``, ``randprune_n`` | ``true``, ``20`` | random LD pruning and its number of iterations |
  | ``randprune_repeats`` | ``default`` | resampling of SNPs picked in many iterations: ``default``, ``maxout`` or ``none`` |
  | ``randprune_file`` | (empty) | ``.mat`` file with fixed prune indices |
  | ``reset_pruneidx`` | ``true`` | regenerate prune indices for each run |
  | ``exclude_chr_pos`` | ``[6 25119106 33854733]`` | regions ``[CHR BP_from BP_to]`` (hg19) excluded from the FDR fit, one per row |
  | ``exclude_from_discovery`` | ``false`` | also exclude those regions from discovery |
  | ``exclude_ambiguous_snps`` | ``false`` | exclude A/T and C/G SNPs from fit and discovery |
  | ``mafthresh`` | ``0.005`` | exclude SNPs with MAF at or below this (or undefined); ``nan`` keeps them (not recommended) |
  | ``perform_gc``, ``use_standard_gc``, ``randprune_gc`` | ``true``, ``false``, ``false`` | genomic control; standard (median) or in-house (more conservative) estimator; estimate after pruning |
  | ``dummy_zscore`` | ``false`` | compute (positive) z-scores from p-values when trait files lack them |
  | ``onscreen`` | ``false`` | show figures on screen |
  | ``manh_plot`` | ``true`` | draw the Manhattan plot |
  | ``manh_colorlist``, ``manh_legend`` | see ``config_default.txt`` | Manhattan plot colours (scaled by 0.8) and legend position |
  | ``manh_fontsize_genenames``, ``manh_yspace``, ``manh_ymargin`` | ``12``, ``0.75``, ``0.25`` | accepted for compatibility; gene names are not drawn |
  | ``mlibrary``, ``exit_matlab_upon_completion`` | | accepted for compatibility and ignored |

  ``config_default.txt`` documents each option in more detail.

## pleioFDR results

  Results are placed in an output folder, defined in ``config.txt`` file. By default it is named ``results``.
  
  Results contain:
   * table with LD-independent significant loci (``*_loci.csv``)
   * table with all analyzed variants and their cond/conj FDR values (``*_all.csv``);
     for conjfdr also tables with z-scores (``*_zscore_*.csv``)
   * conditional qq plots, true discovery rate (tdr) plots, enrichment plots and the FDR lookup table (PNG and SVG)
   * Manhattan plot (PNG, and SVG for conjfdr)
   * log file
   * ``result.mat`` file containing condFDR or conjFDR values for all SNPs, readable with
     ``scipy.io.loadmat`` or MATLAB

## FUMA-defined loci

  Loci tables generated in the step above use custom non-standard logic to clump results based on LD structure.
  You may want to re-generate loci using ``sumstats.py clump`` script, which implements the same logic as in FUMA.
  To do so, convert ``result.mat`` into a text file, and then perform ``sumstats.py clump``.
  At this step you may use ``ref9545380_bfile.tar.gz`` as a reference to preform clumping. For example:
  ```python
  import pandas as pd
  import scipy.io

  res = scipy.io.loadmat("results/result.mat")
  ref = pd.read_csv("9545380.ref", sep="\t")
  ref["FDR"] = res["fdrmat"][:, 0]
  ref.dropna(subset=["FDR"]).to_csv("results/fdr.csv", sep="\t", index=False)
  ```

  Clumped results can then be combined with FUMA annotations and summary statistics using
  ``python -m pleiofdr.fuma.combine`` and ``python -m pleiofdr.fuma.novelty``, see [fuma/readme.md](fuma/readme.md).

## Development and testing

  ```
  uv sync --group dev        # or: pip install -e . pytest ruff
  uv run pytest
  uv run ruff check src tests
  ```

  Most tests need the demo data (``tar -xzvf pleioFDR_demo_data.tar.gz`` in the repository root) and are
  skipped without it (set ``PLEIOFDR_REQUIRE_DEMO=1`` to make missing demo data an error instead).
  ``tests/test_golden.py`` runs the whole analysis on the demo data and compares every
  intermediate result, the lookup tables, the QQ/enrichment matrices and the CSV tables with the output of
  the original MATLAB code. The reference outputs in ``tests/fixtures`` were produced by
  ``tests/fixtures/make_golden.m``; see [MIGRATION_NOTES.md](MIGRATION_NOTES.md) for how the port maps to the
  MATLAB code and where it intentionally differs.

  GitHub Actions ([.github/workflows/tests.yml](.github/workflows/tests.yml)) runs ruff and the full test
  suite, with the demo data downloaded, on every push and pull request, on Linux and macOS with
  Python 3.13 and 3.14.
  [.github/workflows/container.yml](.github/workflows/container.yml) builds both container images, runs the
  self-test without network access, and runs the demo inside each container offline, requiring its CSV
  tables to be byte-identical to the MATLAB reference (``tests/container/demo_parity.sh``); on ``v*`` tags it
  publishes the images.

## MATLAB version

  The original MATLAB/Octave implementation is kept in [legacy_matlab/](legacy_matlab/) for reference.
  It is no longer maintained. The Python version reproduces its results; the differences are listed in
  [MIGRATION_NOTES.md](MIGRATION_NOTES.md) (for example, figures are saved as PNG/SVG rather than ``.fig``).
