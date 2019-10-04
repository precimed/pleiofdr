## Contents

* [Introduction](#introduction)
* [Quick Start](#quick-start)
* [Install pleioFDR](#install-pleiofdr)
* [Data downloads](#data-downloads)
* [Data preparation](#data-preparation)
* [Run pleioFDR](#run-pleiofdr)
* [pleioFDR results](#pleiofdr-results)
* [FUMA-defined loci](#fuma-defined-loci)
* [Octave support](#octave-support)

## Introduction

Pleiotropy-informed conditional and conjunctional false discovery rate allows to boost loci discovery in low-powered GWAS by levereging pleiotropic enrichment with a larger GWAS on related phenotype, and to identify genetic loci joinly associated with two phenotypes.

If you use pleioFDR software for your research publication, please cite the following paper(s):
* Andreassen, O.A. et al. Improved detection of common variants associated with schizophrenia and bipolar disorder using pleiotropy-informed conditional false discovery rate. PLoS Genet 9, e1003455 (2013). 

The pleioFDR software may not be used for commercial purpose or in medical applications.
We encourage all users to familiarize themselves with US patent https://www.google.no/patents/US20150356243 "Systems and methods for identifying polymorphisms".

## Quick Start

To install and run pleioFDR on a small example, constrained to chromosome 21:
```
git clone https://github.com/precimed/pleiofdr && cd pleiofdr
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/pleioFDR_demo_data.tar.gz
tar -xzvf pleioFDR_demo_data.tar.gz
matlab -nodisplay -nosplash < runme.m
```

To install and run pleioFDR using full example:
```
git clone https://github.com/precimed/pleiofdr && cd pleiofdr
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/CTG_COG_2018.mat
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/SSGAC_EDU_2016.mat
cp config_default.txt config.txt
matlab -nodisplay -nosplash < runme.m
```

For the description of the data, see [here](https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt).
For the results, inspect the ``results`` folder.

## Install pleioFDR

Prerequisites:
 - matlab (tested with versions >= 2015)
 - workstation with at least 16GB of RAM
 
The following step by step instruction assumes you are using Linux, however the same can be done in Windows or Mac with minimal modifications.

Download pleioFDR software by going to https://github.com/precimed/pleiofdr in your favorite internet browser, use "Clone or download" button , and "Download zip" do get the latest code.
  
Alternatively, you may get the code by cloning git repository from command line:
  ```
  git clone https://github.com/precimed/pleiofdr && cd pleiofdr
  ```

## Data downloads

Download reference data from [here](https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr). 
The reference is based on 1000 Genomes phase 3 data (May 2, 2013 release).
Variant calls (vcf files) for 22 autosomes were downloaded from [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 
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

## Data preparation

Here we explain how to convert raw summary statistics to pleioFDR format.
Feel free to skip this step if you would like to try pleioFDR on ``CTG_COG_2018.mat`` and ``SSGAC_EDU_2016.mat``,
or if you downloaded input data from the internal NORMENT ``SUMSTATS`` inventory.

Prerequisites:
 - python >= 2.7 (but < 3) with numpy, scipy and pandas libraries

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
  You may also want to change ``traitfile1`` and ``traitfiles`` options.
  
  Start matlab.
  
  Change current folder to the root of ``pleiofdr`` repository (i.e. a folder containing ``pleiotropy_analysis.m``).
  
  Execute ``runme`` command, which should trigger pleiofdr analysis.

  To run pleioFDR from console:
    ```
    matlab -nodisplay -nosplash < runme.m
    ```
    
## pleioFDR results

  Results are placed in an output folder, defined in ``config.txt`` file. By default it is named ``results``.
  
  Results contain:
   * table with LD-independent significant loci
   * table with all analyzed variants and their cond/conj FDR values
   * conditional qq plots and enrichment plots
   * Manhattan plot (by default only .fig file, so you need to open it with
     matlab and save in another format separately)
   * log file
   * ``results.mat`` file containing condFDR or conjFDR values for all SNPs

## FUMA-defined loci

  Loci tables generated in the step above use custom non-standard logic to clump results based on LD structure.
  You may want to re-generate loci using ``sumstats.py clump`` script, which implements the same logic as in FUMA.
  To do so,  convert ``results.mat`` into a text file (for example using ``scipy.io.loadmat``
  and ``pandas.DataFrame.to_csv``), and then perform ``sumstats.py clump``.
  At this step you may use ``ref9545380_bfile.tar.gz`` as a reference to preform clumping.

## Octave support

NB! Octave support is experimental and not officially supported.

install additional packages:
  ```
  octave --no-gui <(echo "pkg install -forge io statistics")
  octave --no-gui <(echo "pkg install -forge nan")
  ```

install  gammainc function:
  
  ```
  wget http://savannah.gnu.org/bugs/download.php?file_id=37342 -O gammainc.m
  wget http://savannah.gnu.org/bugs/download.php?file_id=37341 -O __gammainc_lentz.cc
  mkoctfile __gammainc_lentz.cc 
  ```

run:

  ```
  octave --no-gui runme.m
  ```

#### Octave notes:
* last tested with octave version 4.0.2
* *.mat files need to be saved with "-v7" max:
  ```
  save('ref9545380_1kgPhase3eur_LDr2p1.mat', '-v7') 
  ```
