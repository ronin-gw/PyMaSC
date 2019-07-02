# Change Log

## v0.3.0 (2019-07-03)
- Drop python3.4 support
**Implemented enhancements:**
- Use Fisher's r-to-z transformation to averaging coefficients
- Add chromosome selecting options

## v0.2.6 (2019-01-09)
**Implemented enhancements:**
- Add `--mask-size` option to avoid calling the phantom peak a library length.
  Warnings messages were also implemented.

**Fixed bugs:**
- Fix read counter for unmapped reads.
- Fix bigwig reader.
- Fix mappability stats loader.

## v0.2.5 (2019-01-05)
**Implemented enhancements:**
- Include pre-built C sources and support both pre-built source installation and
  native Cython compile (cython is no longer mandatory dependency).
- Check compatible python versions. (3.4 and 3.7 additionally)
- Scipy is no longer required.

**Fixed bugs:**
- Omit `functools.partialmethod` to keep python2 compatibility.
- Fix minor argparse issues and a error message.

## v0.2.4 (2018-05-09)
**Fixed bugs:**
- Reduce smooth window default size to 15 (Longer lengths are suitable for
  samples with long fragments (such as broad sources) but it's not robust.)
- Use read lengths calculated from CIGAR alignments instead of sequence lengths
- Fix build issue

## v0.2.3 (2018-05-04)
**Implemented enhancements:**
- Add `pymasc-plot` command
- Optimize the default window size for moving average filter for mappability
  sensitive cross-correlation (change from 15 to 100)
- Add annotations to figures about expected library length

**Fixed bugs:**
- Fix bugs relate to mappability sensitive cross-correlation calculation
  Now `bitarray` and `successive` implementations return exactly same results

## v0.2.2 (2017-12-26)
**Implemented enhancements:**
- Warn reference lenth discordancy between alignment and mappability file

**Fixed bugs:**
- Add forgetted initial rotate for mappability array
- BitArray length shortage if mappability tracks include positions which are out
  of range of alignment file reference length

## v0.2.1 (2017-12-23)
**Fixed bugs:**
- BitArray length shortage for reverse strand reads
- BitArray implementation fails for chromosomes which have no mapped reads

## v0.2.0 (2017-12-05)
**Implemented enhancements:**
- Add BitArray implementation
- Some refactoring

**Fixed bugs:**
- Shift size shorter than read length is not acceptable now
- Default window size for smoothing changes 30 to 15

**Known issues:**
- Progress display for bitarray calculation with multiprocessing is currently
  disabled due to `multiprocessing.Queue` stuffing

## v0.1.1 (2017-12-01)
**Implemented enhancements:**
- matplotlib delay import and now import error is not critical for calculation
- Add `--disable-progress` option
- Check python3.5 support
- Some refactoring

**Fixed bugs:**
- Fix unsorted input check (reference disorder)
- Fix progress bar problems
 - Fix multiple line progress bar erasing
 - Fix read count progress bar is not disabling when input has no index
 - Fix progress bar is not disabled when stderr is not a tty
 - Fix interfering mappability progress bar with multiline proggress bar
- Fix some python2 unicode problems
- Fix json serializing
 - `__whole__` entry becomes null for python2
 - serialize error for python3

## v0.1.0 (2017-10-03)
- First release.
