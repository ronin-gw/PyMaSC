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
from typing import Any, Dict, Tuple, Union, TypeVar, Literal, Optional, TextIO

from PyMaSC.utils.output import catch_IOError
from PyMaSC.core.interfaces.stats import GenomeWideStats, CCStats

logger = logging.getLogger(__name__)

STATSFILE_SUFFIX = "_stats.tab"

STAT_ATTR = (
    "genomelen",
    "forward_reads",
    "reverse_reads",
    "cc_min",
    "ccrl"
)

QC_ATTR = (
    "ccfl",
    "nsc",
    "rsc",
    "fwhm",
    "vsn"
)

NCC_LABELS = (
    "Genome length",
    "Forward reads",
    "Reverse reads",
    "Minimum NCC",
    "NCC at read length",
    "NCC at expected library length",
    "NCC at estimated library length",
    "NSC",
    "RSC",
    "Estimated NSC",
    "Estimated RSC",
    "FWHM",
    "VSN",
    "Estimated FWHM",
    "Estimated VSN"
)


MSCC_LABELS = (
    "DMP length",
    "Forward reads in DMP",
    "Reverse reads in DMP",
    "Minimum MSCC",
    "MSCC at read length",
    "MSCC at expected library length",
    "MSCC at estimated library length",
    "MSCC NSC",
    "MSCC RSC",
    "Estimated MSCC NSC",
    "Estimated MSCC RSC",
    "MSCC FWHM",
    "MSCC VSN",
    "Estimated MSCC FWHM",
    "Estimated MSCC VSN"
)


T = TypeVar('T')


def _none2nan(value: T) -> Union[T, Literal["nan"]]:
    """Convert None to 'nan' string for output."""
    return "nan" if value is None else value


@catch_IOError(logger)
def output_stats(outfile: os.PathLike[str], stats_result: GenomeWideStats) -> None:
    """Output comprehensive statistics to tab-delimited file using new StatisticsResult interface.

    Generates a complete statistics file containing all analysis metrics
    including read counts, cross-correlation values, quality metrics,
    and fragment length estimates for both NCC and MSCC.

    Args:
        outfile: Base output file path (suffix will be added)
        stats_result: StatisticsResult object containing analysis results
    """
    outfile_path = Path(outfile)
    basename = outfile_path.name
    outfile_with_suffix = str(outfile_path) + STATSFILE_SUFFIX
    logger.info("Output '{}'".format(outfile_with_suffix))

    with open(outfile_with_suffix, 'w') as f:
        #
        print(f"Name\t{basename}", file=f)
        print(f"Read length\t{stats_result.read_len}", file=f)
        print(f"Expected library length\t{_none2nan(stats_result.expected_lib_len)}", file=f)
        print(f"Estimated library length\t{_none2nan(stats_result.est_lib_len)}", file=f)

        #
        ncc_stats = None if stats_result.whole_ncc_stats is None else stats_result.whole_ncc_stats.stats
        mscc_stats = None if stats_result.whole_mscc_stats is None else stats_result.whole_mscc_stats.stats
        _print_stats(ncc_stats, NCC_LABELS, f)
        _print_stats(mscc_stats, MSCC_LABELS, f)


def _print_stats(stats: Optional[CCStats], labels: Tuple[str, ...], output: TextIO) -> None:
    """Print statistics for a given CCStats object."""
    #
    if stats is None:
        values = ["nan"] * len(labels)
    else:
        values = [_none2nan(getattr(stats, attr)) for attr in STAT_ATTR]
        ccfl, nsc, rsc, fwhm, vsn = (
            _none2nan(getattr(stats.metrics_at_expected_length, attr)) for attr in QC_ATTR
        )
        ccfl_exp, nsc_exp, rsc_exp, fwhm_exp, vsn_exp = (
            _none2nan(getattr(stats.metrics_at_estimated_length, attr)) for attr in QC_ATTR
        )
        values += [
            ccfl, ccfl_exp,
            nsc, rsc, nsc_exp, rsc_exp,
            fwhm, vsn, fwhm_exp, vsn_exp
        ]

    #
    for label, value in zip(labels, values):
        print(f"{label}\t{value}", file=output)


@catch_IOError(logger)
def load_stats(path: os.PathLike[str], names: Tuple[str, ...]) -> Dict[str, Any]:
    """WIP
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
                val = val_converted
                attrs[stat2attr[row]] = val
    return attrs
