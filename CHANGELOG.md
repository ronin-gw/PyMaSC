# Change Log

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
