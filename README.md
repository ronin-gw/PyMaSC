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
Python version 2.7 or >=3.5 are recommended.  
C compiler needs to build C sources (recommend GCC).  

### Install using pip
1. `cython`, `numpy` and `pysam` must be installed ___before___ installing PyMaSC

    $ pip install cython numpy pysam

2. Install PyMaSC from PyPI

    $ pip install pymasc


Usage
-----
After installation, PyMaSC provides `pymasc` and `pymasc-precalc` command.  
Note that `pymasc-precalc` is not essential to calculate mappability-sensitive
cross-correlation.

### `pymasc` command

    pymasc [-h]
           [-p [PROCESS]]
           [-v [{DEBUG,INFO,WARNING,ERROR,CRITICAL}]]  
           [--disable-progress]
           [--color [{TRUE,FALSE}]]
           [--version]
           [-r [READ_LENGTH]]  
           [--estimation-type [{MEAN,MEDIAN,MODE,MIN,MAX}]]  
           [-m REGION_FILE]
           [--mappability-stats [MAPPABILITY_STATS]]  
           [-l [LIBRARY_LENGTH]]
           [-d [MAX_SHIFT]]
           [-q [MAPQ]]  
           [--chi2-pval [CHI2_PVAL]]
           [-w [SMOOTH_WINDOW]]
           [--skip-ncc]  
           [-n [NAME [NAME ...]]]
           [-o [OUTDIR]]
           [--skip-plots]  
           reads [reads ...]  



#### Usage example

Calculate only naïve cross-correlation for `ENCFF000VPI.bam` with max shift size = 1000.  
Output files are `ENCFF000VPI_stats.tab`, `ENCFF000VPI_cc.tab` and `ENCFF000VPI.pdf`.

    $ pymasc -d 1000 ENCFF000VPI.bam


Calculate both naïve and mappability-sensitive cross-correlation by 4 worker processes.  
Output `ENCFF000VPI_mscc.tab` Additionally.

    $ pymasc -p 4 -d 1000 -m wgEncodeCrgMapabilityAlign36mer.bigWig ENCFF000VPI.bam


#### Input file
SAM and BAM file format are acceptable. Input alignment file must be sorted.  
For parallel processing, input file must be BAM format and indexed.  
If multiple files specified, PyMaSC processes each file with common parameters.  
Unmapped or duplicated reads will be discarded.  
If input file contains paired-end reads, the last (second) segment read will be discarded.

#### General options

##### -p / --process [int]
Set number of worker process. (Default: 1)  
For indexed BAM file, PyMaSC parallel process each reference (chromosome).

##### -v / --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
Set logging message level. (Default: info)  

##### --disable-progress
Disable progress bars.  
Note that progress bar will be disabled automatically if stderr is not connected to terminal.

##### --color {TRUE,FALSE}
Switch coloring log output. (Default: auto; enable if stderr is connected to terminal)


#### Read length settings

##### -r / --read-length [int]
Specify read length explicitly. (Default: get representative by scanning)  
PyMaSC needs representative value of read length to plot figures and to calc
mappability-sensitive cross-correlation. By default, PyMaSC scans input file
read length to get representative read length. If read length is specified, PyMaSC
skips this step.  
Note that this option must be specified to treat unseekable input (like stdin).

##### --estimation-type {MEAN,MEDIAN,MODE,MIN,MAX}
Specify how to get representative value of read length. (Default: median)


#### Mappability options

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

##### -q / --mapq [int]
Specify discarding reads which mapping quality less than specified score. (Default: 1)  
MAPQ >= 1 is recommended because MAPQ=0 contains multiple hit reads.

##### -w / --smooth-window [int]
Before mean fragment length estimation, PyMaSC applies moving average filter to
mappability-sensitive cross-correlation. This option specify filter's window size.
(Default: 30)

##### --chi2-pval [float]
P-value threshold to check strand specificity. (Default: 0.05)  
PyMaSC performs chi-square test between number of reads mapped to positive- and negative-strand.

##### --skip-ncc
Both `-m/--mappability` and `--skip-ncc` specified, PyMaSC skips calculate naïve cross-correlation
and calculates only mappability-sensitive cross-correlation. (Default: False)

##### -l / --library-length
Specify expected fragment length. (Default: None)  
If mappability region file is not supplied, NSC and RSC estimation will be perform
with this value instead of estimated mean fragment length.


#### Output options

##### -o / --outdir [path]
Specify output directory. (Default: current directory)

##### -n / --name [NAME...]
By default, output files are written to `outdir/input_file_base_name`. This option
overwrite output file base name.

##### --skip-plots
Skip output figures. (Default: False)


### `pymasc-precalc` command

    pymasc-precalc [-h]
                   [-p [PROCESS]]
                   [-v [{DEBUG,INFO,WARNING,ERROR,CRITICAL}]]
                   [--disable-progress]
                   [--color [{TRUE,FALSE}]]
                   [--version]
                   [-m REGION_FILE]
                   [--mappability-stats [MAPPABILITY_STATS]]
                   [-d [MAX_SHIFT]]
                   [-r [MAX_READLEN]]

#### Usage example
Calculate total length of doubly mappable region.  
`wgEncodeCrgMapabilityAlign36mer_mappability.json` will be write.

    $ pymasc -p 4 -r 50 -d 1000 -m wgEncodeCrgMapabilityAlign36mer.bigWig

#### Options
Almost same as `pymasc` command.  
Note that actual max shift size is,
- 0 to `read_length` (if `max_shift` < `read_len` * 2)
- 0 to `max_shift` - `read_len` + 1 (if `max_shift` => `read_len` * 2)
