"""Statistics file output and loading functions.

Handles generation and loading of PyMaSC statistics files containing
numerical summaries of cross-correlation analysis results. Statistics
files use tab-delimited format for integration with downstream tools.

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
from typing import overload, Union, TypeVar, Literal

from PyMaSC.interfaces.output import (
    StatLabels,
    SummaryItems,
    CorrelationStats, CorrelationItems, NCC_LABELS, MSCC_LABELS,
    OutputStats
)
from PyMaSC.utils.output import catch_IOError
from PyMaSC.interfaces.stats import GenomeWideStats, CCStats

logger = logging.getLogger(__name__)

STATSFILE_SUFFIX = "_stats.tab"

T = TypeVar('T', str, int, float)


@overload
def _none2nan(value: None) -> Literal["nan"]: ...
@overload
def _none2nan(value: T) -> T: ...


def _none2nan(value: Union[T, None]) -> Union[T, Literal["nan"]]:
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
    outfile_with_suffix = str(outfile_path) + STATSFILE_SUFFIX
    logger.info("Output '{}'".format(outfile_with_suffix))

    output = make_output_from_results(outfile_path.name, stats_result)

    with open(outfile_with_suffix, 'w') as f:
        #
        for attr, value in output.get_with_labels():
            print(attr, value, sep='\t', file=f)

        #
        for attr, value in output.ncc_stats.get_with_labels():
            print(attr, value, sep='\t', file=f)

        #
        for attr, value in output.mscc_stats.get_with_labels():
            print(attr, value, sep='\t', file=f)


def make_output_from_results(name: str, result: GenomeWideStats) -> OutputStats:
    """Create OutputStats from GenomeWideStats for output.

    Args:
        result: GenomeWideStats object containing analysis results

    Returns:
        OutputStats object ready for output
    """
    if result.whole_ncc_stats is not None:
        ncc_stats = _make_corr_output_stats(result.whole_ncc_stats.stats, NCC_LABELS)
    else:
        ncc_stats = _make_nan_output_stats(NCC_LABELS)

    if result.whole_mscc_stats is not None:
        mscc_stats = _make_corr_output_stats(result.whole_mscc_stats.stats, MSCC_LABELS)
    else:
        mscc_stats = _make_nan_output_stats(MSCC_LABELS)

    return OutputStats(
        stats=SummaryItems(
            name,
            result.read_len,
            _none2nan(result.expected_lib_len),
            _none2nan(result.est_lib_len)
        ),
        ncc_stats=ncc_stats,
        mscc_stats=mscc_stats
    )


def _make_corr_output_stats(stats: CCStats, labels: StatLabels) -> CorrelationStats:
    expected_len_qc = stats.metrics_at_expected_length
    estimated_len_qc = stats.metrics_at_estimated_length
    return CorrelationStats(
        stats=CorrelationItems(
            _none2nan(stats.genomelen_repr),
            _none2nan(stats.forward_reads_repr),
            _none2nan(stats.reverse_reads_repr),
            _none2nan(stats.cc_min),
            _none2nan(stats.ccrl),
            _none2nan(expected_len_qc.ccfl),
            _none2nan(estimated_len_qc.ccfl),
            _none2nan(expected_len_qc.nsc),
            _none2nan(expected_len_qc.rsc),
            _none2nan(estimated_len_qc.nsc),
            _none2nan(estimated_len_qc.rsc),
            _none2nan(expected_len_qc.fwhm),
            _none2nan(expected_len_qc.vsn),
            _none2nan(estimated_len_qc.fwhm),
            _none2nan(estimated_len_qc.vsn)
        ),
        labels=labels
    )


def _make_nan_output_stats(labels: StatLabels) -> CorrelationStats:
    return CorrelationStats(
        stats=CorrelationItems(
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan',
            'nan'
        ),
        labels=labels
    )
