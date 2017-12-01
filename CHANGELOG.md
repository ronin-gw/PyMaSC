# Change Log

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
