"""Constants used throughout PyMaSC for cross-correlation analysis.

This module centralizes important constants that were originally defined
in result.py, making them more accessible and maintainable.
"""

# Quality Assessment Constants
NEAR_READLEN_ERR_CRITERION = 5
"""int: Distance threshold for detecting phantom peak contamination.

When the estimated fragment length is within this many base pairs of the
read length, it suggests potential phantom peak contamination. This triggers
warnings to alert users about potential quality issues in the data.
"""

MERGED_CC_CONFIDENCE_INTERVAL = 0.99
"""float: Confidence level for merged cross-correlation intervals.

This sets the confidence level (99%) for calculating confidence intervals
when merging cross-correlation results from multiple chromosomes using
Fisher z-transformation. Higher values provide wider, more conservative
confidence intervals.
"""

NEAR_ZERO_MIN_CALC_LEN = 10
"""int: Number of positions to examine at the beginning of correlation profile.

Used to detect potential issues with minimum correlation calculation.
If the median of the first NEAR_ZERO_MIN_CALC_LEN positions is smaller
than the calculated minimum, it suggests the shift size may be too small.
"""