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
from __future__ import print_function

import os.path
import logging

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


def get_rl_item_from(attr):
    """Create accessor function for read-length-indexed arrays.
    
    Generates a function that extracts values from arrays at the read length
    position, handling None values gracefully.
    
    Args:
        attr: Attribute name to access from CCStats objects
        
    Returns:
        Function that extracts read-length-indexed values
    """
    def _func(ccstats):
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
def output_stats(outfile, ccr):
    """Output comprehensive statistics to tab-delimited file.
    
    Generates a complete statistics file containing all analysis metrics
    including read counts, cross-correlation values, quality metrics,
    and fragment length estimates for both NCC and MSCC.
    
    Args:
        outfile: Base output file path (suffix will be added)
        ccr: CCResult object containing analysis results
    """
    basename = os.path.basename(outfile)
    outfile += STATSFILE_SUFFIX
    logger.info("Output '{}'".format(outfile))

    with open(outfile, 'w') as f:
        print("{}\t{}".format("Name", basename), file=f)
        for row, attr in STAT_ATTR:
            print("{}\t{}".format(row, getattr(ccr.whole, attr, False) or "nan"), file=f)
        for row, attr in STAT_CC_ATTR:
            print("{}\t{}".format(row, getattr(ccr.whole.cc, attr, False) or "nan"), file=f)
        for row, attr in STAT_MSCC_ATTR:
            if callable(attr):
                print("{}\t{}".format(row, attr(ccr.whole.masc)), file=f)
            else:
                print("{}\t{}".format(row, getattr(ccr.whole.masc, attr, False) or "nan"), file=f)


@catch_IOError(logger)
def load_stats(path, names):
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
                val = int(val) if val.strip().isdigit() else None
                attrs[stat2attr[row]] = val
    return attrs
