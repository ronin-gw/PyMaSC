PyMaSC
======

[![PyPI version](https://badge.fury.io/py/PyMaSC.svg)](https://badge.fury.io/py/PyMaSC)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/PyMaSC.svg)](https://pypi.python.org/pypi/PyMaSC/)
[![Build Status](https://travis-ci.org/ronin-gw/PyMaSC.svg?branch=master)](https://travis-ci.org/ronin-gw/PyMaSC)

Python implementation to calc mappability-sensitive cross-correlation
for fragment length estimation and quality control for ChIP-Seq.

* * *

Introduction
------------
To estimate mean fragment length for single-ended sequencing data, cross-correlation
between positive- and negative-strand reads is commonly used. One of the problems with
this approach is _phantom_ peak, which is the occasionally observed peak corresponding
to the read length. In the ChIP-Seq guidelines by ENCODE consortia, cross-correlation
at fragment length and read length are used for quality control metrics. Additionally,
library length itself is one of the important parameters for analysises. However,
estimating correct flagment length is not a easy task because of _phantom_ peak occurrence.  
P Ramachandran _et al._ proposed __MaSC__, mappability-sensitive cross-correlation to
remove the bias caused from ununiformity of mappability throughout the genome. This
method provides cross-correlation landscape without _phantom_ peak and much accurate
mean fragment length estimation.  
__PyMaSC__ is a tool implemented by python and cython to visualize (mappability-sensitive)
cross-correlation and estimate ChIP-Seq quality metrics and mean fragment length with
MaSC algorithm.


Install
-------
Python version 2.7 or >=3.4 are recommended.  
C compiler needs to build C sources (recommend GCC).  

### Install using pip
1. `numpy` and `pysam==0.15.1` must be installed ___before___ installing PyMaSC

    $ pip install numpy pysam==0.15.1

2. Install PyMaSC from PyPI

    $ pip install pymasc

    If cython is installed, PyMaSC will be built with Cython native compiled sources instead of pre-compiled C sources.

Usage
-----
After installation, PyMaSC provides `pymasc`, `pymasc-precalc` and `pymasc-plot` command.  
Note that `pymasc-precalc` is not essential to calculate mappability-sensitive
cross-correlation.

### `pymasc` command

    pymasc [-h]
           [-p PROCESS]
           [-v {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
           [--disable-progress]
           [--color {TRUE,FALSE}]
           [--version]
           [--successive]
           [-r READ_LENGTH]
           [--estimation-type {MEAN,MEDIAN,MODE,MIN,MAX}]
           [-m REGION_FILE]
           [--mappability-stats MAPPABILITY_STATS]
           [-l LIBRARY_LENGTH]
           [-d MAX_SHIFT]
           [-q MAPQ]
           [--chi2-pval CHI2_PVAL]
           [-w SMOOTH_WINDOW]
           [--skip-ncc]  
           [-n NAME [NAME ...]]
           [-o OUTDIR]
           [--skip-plots]  
           reads [reads ...]  



#### Usage example

Calculate only naïve cross-correlation for `ENCFF000VPI.bam` with max shift size = 1000.  
Output files are `ENCFF000VPI_stats.tab`, `ENCFF000VPI_cc.tab` and `ENCFF000VPI.pdf`.

    $ pymasc -d 1000 ENCFF000VPI.bam


Calculate both naïve and mappability-sensitive cross-correlation by 4 worker processes.  
Output `ENCFF000VPI_mscc.tab` Additionally.

    $ pymasc -p 4 -d 1000 -m wgEncodeCrgMapabilityAlign36mer.bigWig ENCFF000VPI.bam


#### Main input file
SAM and BAM file format are acceptable.
* Input alignment file must be sorted.  
 * Additionally, for parallel processing, input file must be BAM format and indexed.  
* If multiple files specified, PyMaSC processes each file with common parameters.  
* Unmapped or duplicated reads will be discarded.  
* If input file contains paired-end reads, the last (second) segment read will be discarded.


#### Output files

#### General options

##### -v / --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
Set logging message level. (Default: info)  

##### --disable-progress
Disable progress bars.  
Note that progress bar will be disabled automatically if stderr is not connected to terminal.

##### --color {TRUE,FALSE}
Switch coloring log output. (Default: auto; enable if stderr is connected to terminal)

#### --version
Show program's version number and exit


#### Processing settings

##### -p / --process [int]
Set number of worker process. (Default: 1)  
For indexed BAM file, PyMaSC parallel process each reference (chromosome).

#### --successive
Calc with successive algorithm instead of bitarray implementation (Default: false)
Bitarray implementation is recommended　in most situation. See `Computation details`
for more information.

##### --skip-ncc
Both `-m/--mappability` and `--skip-ncc` specified, PyMaSC skips calculate naïve cross-correlation
and calculates only mappability-sensitive cross-correlation. (Default: False)

##### --skip-plots
Skip output figures. (Default: False)


#### Input alignment file settings

##### -r / --read-length [int]
Specify read length explicitly. (Default: get representative by scanning)  
PyMaSC needs representative value of read length to plot figures and to calc
mappability-sensitive cross-correlation. By default, PyMaSC scans input file
read length to get representative read length. If read length is specified, PyMaSC
skips this step.  
Note that this option must be specified to treat unseekable input (like stdin).

##### --readlen-estimator {MEAN,MEDIAN,MODE,MIN,MAX}
Specify how to get representative value of read length. (Default: median)

##### -q / --mapq [int]
Input reads which mapping quality less than specified score will be discarded. (Default: 1)  
MAPQ >= 1 is recommended because MAPQ=0 contains multiple hit reads.

##### -l / --library-length
Specify expected fragment length. (Default: None)  
PyMaSC supplies additional NSC and RSC values calculated from this value.


#### Input mappability file settings

##### -m / --mappability [BigWig file]
Specify mappability (alignability, uniqueness) track to calculate mappability-sensitive
cross-correlation.  
Input file must be BigWig format and each track's score should indicate mappability
in [0, 1] (1 means uniquely mappable position).  
If BigWig file is not supplied, PyMaSC will calculate only naïve cross-correlation.

##### --mappability-stats [json file]
Read and save path to the json file which contains mappability region statistics.  
(Default: same place, same base name as the mappability BigWig file)  
If there is no statistics file for specified BigWig file, PyMaSC calculate total
length of doubly mappable region for each shift size automatically and save them
to reuse for next calculation and faster computing.  
`pymasc-precalc` performs this calculation for specified BigWig file (this is not
necessary, of course).


#### Analysis Parameters

##### -d / --max-shift [int]
PyMaSC calculate cross-correlation with shift size from 0 to this value. (Default: 1000)

##### --chi2-pval [float]
P-value threshold to check strand specificity. (Default: 0.05)  
PyMaSC performs chi-square test between number of reads mapped to positive- and negative-strand.

##### -w / --smooth-window [int]
Before mean fragment length estimation, PyMaSC applies moving average filter to
mappability-sensitive cross-correlation. This option specify filter's window size.
(Default: 15)


#### Output options

##### -o / --outdir [path]
Specify output directory. (Default: current directory)

##### -n / --name [NAME...]
By default, output files are written to `outdir/input_file_base_name`. This option
overwrite output file base name.


### `pymasc-precalc` command

    pymasc-precalc [-h]
                   [-p PROCESS]
                   [-v {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                   [--disable-progress]
                   [--color {TRUE,FALSE}]
                   [--version]
                   [-m REGION_FILE]
                   [--mappability-stats MAPPABILITY_STATS]
                   [-d MAX_SHIFT]
                   [-r MAX_READLEN]

#### Usage example
Calculate total length of doubly mappable region.  
`wgEncodeCrgMapabilityAlign36mer_mappability.json` will be write.

    $ pymasc -p 4 -r 50 -d 1000 -m wgEncodeCrgMapabilityAlign36mer.bigWig

#### Options
Almost same as `pymasc` command.  
Note that actual max shift size is,
- 0 to `read_length` (if `max_shift` < `read_len` * 2)
- 0 to `max_shift` - `read_len` + 1 (if `max_shift` => `read_len` * 2)


### `pymasc-plot` command

    pymasc-plot [-h]
                [-v {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                [--disable-progress]
                [--color {TRUE,FALSE}]
                [--version]
                [-s STATS]
                [-c CC]
                [-m MASC]
                [--chi2-pval
                CHI2_PVAL]
                [-w SMOOTH_WINDOW]
                [-l LIBRARY_LENGTH]
                [-n NAME]
                [-o OUTDIR]
                [statfile]

#### Usage example
(Re)plot figures from `pymasc` outputs `output/ENCFF000VPI_*` with the specified
smoothing window size and library length.

    $ pymasc-plot -w 50 -l 250 output/ENCFF000VPI

#### Input argument
Specify a prefix to `pymasc` output files. For example, set `output/ENCFF000VPI`
to plot figures from `output/ENCFF000VPI_stats.tab` and `output/ENCFF000VPI_cc.tab`
(and/or `output/ENCFF000VPI_masc.tab`). `*_stats.tab` and either or both of `*_cc.tab`
and `*_masc.tab` must be exist.  
To specify these files individually, use `-s/--stats`, `-c/--cc` and `-m/--masc`
options.


Computation details
-------------------

### BitArray and successive implementation
PyMaSC provides two algorithms for calculation.

#### BitArray
BitArray approach is based on
allocated binary arrays with length of references and computation time mostly
depends on the size of reference genome and the maximum shift size. Typically
(hg19, chr1), PyMaSC consumes about 250MB RAM per worker.  

#### Successive implementation
Successive implementation is based on buffer with length of maximum shift size
and calculate overwrapped bases successively while reading reads and mappable
regions. Computation time mostly depends on number of reads and tracks in input
files and BitArray implementation is faster in most cases.
Successive approach surpasses in processing small input data, quite low memory
efficiency and robustness for shift size.


References
----------
* Ramachandran, Parameswaran, et al. "MaSC: mappability-sensitive cross-correlation
  for estimating mean fragment length of single-end short-read sequencing data."
  Bioinformatics 29.4 (2013): 444-450.
