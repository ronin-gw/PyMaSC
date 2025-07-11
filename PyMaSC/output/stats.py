"""Statistics file output and loading functions.

This module handles the generation and loading of PyMaSC statistics files,
which contain comprehensive numerical summaries of cross-correlation analysis
results. Statistics files use tab-delimited format for easy parsing and
integration with downstream analysis tools.

The module supports:
- Output of complete analysis statistics to tab-delimited files
- Loading of statistics for plot regeneration and analysis continuation
- Both naive cross-correlation and MSCC statistics
- Quality metrics (NSC, RSC, FWHM, VSN) and fragment length estimates
"""
from __future__ import annotations, print_function

import logging
import os
from pathlib import Path
from typing import Any, Callable, Dict, Tuple, Union

from PyMaSC.utils.output import catch_IOError

STATSFILE_SUFFIX = "_stats.tab"
STAT_ATTR = (
    ("Read length", "read_len"),
    ("Expected library length", "library_len"),
    ("Estimated library length", "est_lib_len")
)
STAT_CC_ATTR = (
    ("Genome length", "genomelen"),
    ("Forward reads", "forward_sum"),
    ("Reverse reads", "reverse_sum"),
    ("Minimum NCC", "cc_min"),
    ("NCC at read length", "ccrl"),
    ("NCC at expected library length", "ccfl"),
    ("NCC at estimated library length", "est_ccfl"),
    ("NSC", "nsc"),
    ("RSC", "rsc"),
    ("Estimated NSC", "est_nsc"),
    ("Estimated RSC", "est_rsc"),
    ("FWHM", "cc_width"),
    ("VSN", "vsn"),
    ("Estimated FWHM", "est_cc_width"),
    ("Estimated VSN", "est_vsn")
)


def get_rl_item_from(attr: str) -> Callable[[Any], Union[str, Any]]:
    """Create accessor function for read-length-indexed arrays.

    Generates a function that extracts values from arrays at the read length
    position, handling None values gracefully.

    Args:
        attr: Attribute name to access from CCStats objects

    Returns:
        Function that extracts read-length-indexed values
    """
    def _func(ccstats: Any) -> Union[str, Any]:
        array = getattr(ccstats, attr, None)
        return "nan" if array is None else array[ccstats.read_len - 1]
    return _func


STAT_MSCC_ATTR = (
    ("DMP length", get_rl_item_from("genomelen")),
    ("Forward reads in DMP", get_rl_item_from("forward_sum")),
    ("Reverse reads in DMP", get_rl_item_from("reverse_sum")),
    ("Minimum MSCC", "cc_min"),
    ("MSCC at read length", "ccrl"),
    ("MSCC at expected library length", "ccfl"),
    ("MSCC at estimated library length", "est_ccfl"),
    ("MSCC NSC", "nsc"),
    ("MSCC RSC", "rsc"),
    ("Estimated MSCC NSC", "est_nsc"),
    ("Estimated MSCC RSC", "est_rsc"),
    ("MSCC FWHM", "cc_width"),
    ("MSCC VSN", get_rl_item_from("vsn")),
    ("Estimated MSCC FWHM", "est_cc_width"),
    ("Estimated MSCC VSN", get_rl_item_from("est_vsn"))
)

logger = logging.getLogger(__name__)


@catch_IOError(logger)
def output_stats(outfile: os.PathLike[str], ccr: Any) -> None:
    """Output comprehensive statistics to tab-delimited file.

    Generates a complete statistics file containing all analysis metrics
    including read counts, cross-correlation values, quality metrics,
    and fragment length estimates for both NCC and MSCC.

    Args:
        outfile: Base output file path (suffix will be added)
        ccr: CCResult object containing analysis results
    """
    outfile_path = Path(outfile)
    basename = outfile_path.name
    outfile_with_suffix = str(outfile_path) + STATSFILE_SUFFIX
    logger.info("Output '{}'".format(outfile_with_suffix))

    with open(outfile_with_suffix, 'w') as f:
        print("{}\t{}".format("Name", basename), file=f)
        for row, attr in STAT_ATTR:
            val = getattr(ccr.whole, attr, None)
            if val is None or (isinstance(val, bool) and not val):
                val = "nan"
            print("{}\t{}".format(row, val), file=f)
        for row, attr in STAT_CC_ATTR:
            val = getattr(ccr.whole.cc, attr, None)
            if val is None or (isinstance(val, bool) and not val):
                val = "nan"
            # Handle array values that may be VSN or other arrays
            elif hasattr(val, '__len__') and not isinstance(val, str):
                # Arrays should not reach here for simple attributes
                val = "nan"  # Fallback to avoid ambiguous array evaluation
            print("{}\t{}".format(row, val), file=f)
        for row, attr in STAT_MSCC_ATTR:  # type: ignore[assignment]
            if callable(attr):
                print("{}\t{}".format(row, attr(ccr.whole.masc)), file=f)
            else:
                print("{}\t{}".format(row, getattr(ccr.whole.masc, attr, False) or "nan"), file=f)


@catch_IOError(logger)
def load_stats(path: os.PathLike[str], names: Tuple[str, ...]) -> Dict[str, Any]:
    """Load statistics from tab-delimited file.

    Reads statistics file and extracts requested attributes for use in
    plot regeneration or analysis continuation.

    Args:
        path: Path to statistics file
        names: List of attribute names to extract

    Returns:
        Dictionary mapping attribute names to values
    """
    logger.info("Load statistics from '{}'.".format(path))
    stat2attr = {k: v for k, v in STAT_ATTR if v in names}
    attrs = {}
    with open(path) as f:
        for l in f:
            row, val = l.split('\t', 1)
            if row == "Name":
                attrs["name"] = val.rstrip()
            if row in stat2attr:
                val_converted: Union[int, None] = int(val) if val.strip().isdigit() else None
                val = val_converted  # type: ignore[assignment]
                attrs[stat2attr[row]] = val
    return attrs
