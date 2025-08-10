"""Matplotlib-based plotting functions for PyMaSC visualization.

Provides plotting functionality for PyMaSC analysis results including
cross-correlation plots, MSCC plots, and comparison visualizations.
Generates multi-page PDF files for viewing and publication.

Key plot types:
- Naive cross-correlation plots with quality metrics
- Mappability-sensitive cross-correlation (MSCC) plots
- Fragment length estimation visualization
- Comparison plots between NCC and MSCC
- Statistical annotations and quality indicators

The module handles matplotlib import failures gracefully and provides detailed
error logging for troubleshooting visualization issues.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Optional, Tuple, Union

import numpy as np

from PyMaSC.interfaces.stats import (
    GenomeWideStats, WholeGenomeStats,
    ChromosomeStats,
    TStats, NCCStats, MSCCStats
)
from PyMaSC.utils.output import catch_IOError

logger = logging.getLogger(__name__)

try:
    import matplotlib.pyplot as plt  # type: ignore[import-untyped]
    from matplotlib.backends.backend_pdf import PdfPages  # type: ignore[import-untyped]
except Exception:
    logger.error("Failed to import matplotlib.")
    import traceback
    logger.warning("Exception traceback:\n|" +
                   traceback.format_exc().replace('\n', "\n|"))
    raise


def _feed_pdf_page(pp: Any) -> None:
    """Save current plot to PDF and close the figure.

    Args:
        pp: PdfPages object for multi-page PDF output
    """
    pp.savefig()
    plt.close()


@catch_IOError(logger)
def plot_figures(outfile: os.PathLike[str], stats: GenomeWideStats) -> None:
    """Generate all PyMaSC plots and save to PDF file using new StatisticsResult interface.

    Creates a multi-page PDF containing all relevant plots based on the
    analysis results. Includes naive cross-correlation, MSCC, and comparison
    plots depending on what calculations were performed.

    Args:
        outfile: Output PDF file path
        stats_result: StatisticsResult object containing analysis results
    """
    outfile_path = Path(outfile)
    logger.info("Output '{}'".format(outfile_path))
    name = outfile_path.stem

    with PdfPages(os.fspath(outfile_path)) as pp:
        if stats.whole_ncc_stats:
            plot_naive_cc(stats.whole_ncc_stats, name)
            _feed_pdf_page(pp)

        if stats.whole_mscc_stats:
            est_lib_len = stats.whole_mscc_stats.est_lib_len
            if plot_naive_cc_just(stats.whole_ncc_stats, est_lib_len, name):
                _feed_pdf_page(pp)

            plot_masc(stats.whole_mscc_stats, name)
            _feed_pdf_page(pp)

        plot_ncc_vs_masc(pp, stats, name)


def _annotate_point(
    x: Union[int, float],
    color: str,
    axis_y: Union[int, float],
    axis_text: str,
    point_y: Optional[Union[int, float]] = None,
    point_text: Optional[str] = None,
    yoffset: Union[int, float] = 0
) -> None:
    """Add vertical line and annotations to plot at specified x position.

    Args:
        x: X-coordinate for vertical line and annotations
        color: Color for the vertical line
        axis_y: Y-coordinate for axis annotation
        axis_text: Text for axis annotation
        point_y: Y-coordinate for point annotation (optional)
        point_text: Text for point annotation (optional)
        yoffset: Vertical offset for point annotation
    """
    plt.axvline(x, color=color, linestyle="dashed", linewidth=0.5)
    plt.annotate(axis_text, (x, axis_y))
    if point_y and point_text:
        plt.scatter(x, point_y, facecolors="none", edgecolors=color)
        plt.annotate(point_text, (x, point_y + yoffset))


def _annotate_bottom_right_box(text: str) -> None:
    """Add text box annotation to bottom right corner of plot.

    Args:
        text: Text content for the annotation box
    """
    plt.annotate(
        text,
        textcoords="axes fraction", xy=(1, plt.gca().get_ylim()[0]), xytext=(0.95, 0.05),
        bbox=dict(boxstyle="round", fc="w", alpha=0.9), horizontalalignment="right"
    )


def _annotate_params(
    nsc: Optional[Union[int, float]] = None,
    rsc: Optional[Union[int, float]] = None,
    est_nsc: Optional[Union[int, float]] = None,
    est_rsc: Optional[Union[int, float]] = None,
    loc: str = "lower right"
) -> None:
    """Add parameter annotations to plot showing quality metrics.

    Args:
        nsc: Normalized Strand Coefficient value
        rsc: Relative Strand Coefficient value
        est_nsc: Estimated NSC value
        est_rsc: Estimated RSC value
        loc: Location for the annotation box
    """
    anno = []
    for stat, label in zip((nsc, rsc, est_nsc, est_rsc),
                           ("NSC", "RSC", "Est NSC", "Est RSC")):
        if stat:
            anno.append("{} = {:.5f}".format(label, stat))

    if anno:
        _annotate_bottom_right_box('\n'.join(anno))


def _set_ylim() -> Tuple[float, float, float]:
    """Adjust Y-axis limits with appropriate padding.

    Sets Y-axis limits to provide better visualization by adding
    padding around the data range.
    """
    axes = plt.gca()
    lower, upper = axes.get_ylim()
    if upper > 0:
        lower, upper = axes.set_ylim((lower, upper * 1.1))
    else:
        lower, upper = axes.set_ylim((lower, upper * 0.95))
    height = upper - lower
    return lower, upper, height


def plot_naive_cc(
    whole_stat: WholeGenomeStats[TStats],
    name: Optional[str] = None,
    xlim: Optional[Tuple[int, Any]] = None
) -> None:
    """Plot naive cross-correlation results using new StatisticsResult interface.

    Creates a comprehensive plot of naive cross-correlation showing:
    - Cross-correlation curve across all shift distances
    - Read length and fragment length markers
    - Quality metrics (NSC, RSC) annotations
    - Statistical significance indicators

    Args:
        stats_result: StatisticsResult object containing analysis results
        name: Sample name for plot title
        xlim: X-axis limits (optional)
    """
    title = "Cross-Correlation"
    if name:
        title += " for " + name

    if whole_stat is None:
        logger.warning("No CC statistics available for plotting")
        return

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Cross-Correlation")

    # Plot confidence intervals if available
    max_shift = len(whole_stat.cc) - 1
    x_range = range(max_shift + 1)

    # Plot confidence intervals
    plt.fill_between(x_range, whole_stat.cc_lower, whole_stat.cc_upper,
                     color="lightskyblue", alpha=0.5, linewidth=0)

    # Plot main cross-correlation curve
    if whole_stat.cc is not None:
        plt.plot(x_range, whole_stat.cc, color="black", linewidth=0.5)

    axes = plt.gca()
    if xlim:
        axes.set_xlim(xlim)
    lower, upper, height = _set_ylim()

    # Plot minimum correlation line
    cc_stats = whole_stat.stats
    if cc_stats.cc_min is not None:
        plt.axhline(cc_stats.cc_min, linestyle="dashed", linewidth=0.5)
        plt.text(0, cc_stats.cc_min, 'min(cc) = {:.5f}'.format(cc_stats.cc_min))

    # Annotate read length
    read_len = cc_stats.read_len
    if read_len and cc_stats.ccrl is not None:
        _annotate_point(
            read_len - 1, "red",
            upper - height / 25, 'read length: {}'.format(read_len),
            cc_stats.ccrl, " cc(read length) = {:.5f}".format(cc_stats.ccrl), height / 50
        )

    # Annotate estimated library length
    qc_stats = cc_stats.metrics_at_estimated_length
    if qc_stats.fragment_length:
        est_ccfl = qc_stats.ccfl
        if est_ccfl is not None:
            _annotate_point(
                qc_stats.fragment_length - 1, "blue",
                upper - height / 10, 'estimated lib len: {}'.format(qc_stats.fragment_length),
                est_ccfl, " cc(est lib len) = {:.5f}".format(est_ccfl), height / 50
            )

    # Annotate expected library length
    qc_stats = cc_stats.metrics_at_expected_length
    if qc_stats and qc_stats.ccfl is not None:
        library_len = qc_stats.fragment_length
        assert library_len is not None
        _annotate_point(
            library_len - 1, "green",
            upper - height / 6, 'expected lib len: {}'.format(library_len),
            qc_stats.ccfl, " cc(lib length) = {:.5f}".format(qc_stats.ccfl), -height / 25
        )

    # Add quality metrics annotations
    _annotate_params(
        getattr(cc_stats, 'nsc', None),
        getattr(cc_stats, 'rsc', None),
        getattr(cc_stats, 'est_nsc', None),
        getattr(cc_stats, 'est_rsc', None)
    )


def plot_naive_cc_just(
    stats: Optional[WholeGenomeStats[NCCStats]],
    est_lib_len: Optional[int],
    name: Optional[str] = None
) -> bool:
    """Plot zoomed naive cross-correlation around fragment length peak using new StatisticsResult interface.

    Creates a focused plot showing the cross-correlation peak region
    with detailed annotations for fragment length estimation.

    Args:
        stats_result: StatisticsResult object containing analysis results
        name: Sample name for plot title

    Returns:
        True if plot was created, False if insufficient data
    """
    if stats is None:
        return False

    # Check if we can create a zoomed plot
    if stats.cc is None:
        return False

    max_shift = len(stats.cc) - 1
    if est_lib_len is not None and est_lib_len * 2 < max_shift + 1:
        plot_naive_cc(stats, name, (0, est_lib_len * 2))
        return True

    return False


def plot_masc(
    masc_stats: Optional[WholeGenomeStats[MSCCStats]],
    name: Optional[str] = None
) -> None:
    """Plot mappability-sensitive cross-correlation (MSCC) results using new StatisticsResult interface.

    Creates comprehensive MSCC plots showing:
    - MSCC curve with mappability correction
    - Library length estimation with confidence intervals
    - Comparison with naive cross-correlation
    - Quality metrics specific to MSCC analysis

    Args:
        stats_result: StatisticsResult object containing analysis results
        name: Sample name for plot title
    """
    title = "MSCC and Library Length Estimation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Mappability Sensitive Cross-Correlation")

    if masc_stats is None:
        logger.warning("No MSCC statistics available for plotting")
        return

    max_shift = len(masc_stats.cc) - 1
    x_range = range(max_shift + 1)

    # Plot confidence intervals for MSCC
    plt.fill_between(x_range, masc_stats.cc_lower, masc_stats.cc_upper,
                     color="lightskyblue", alpha=0.5, linewidth=0)

    # Plot MSCC curves
    if masc_stats.cc is not None:
        plt.plot(x_range, masc_stats.cc, color="black", linewidth=0.5, label="MSCC")
        plt.plot(x_range, masc_stats.avr_cc, alpha=0.8, label="Smoothed", color="pink")

    lower, upper, height = _set_ylim()

    # Annotate estimated library length
    est_lib_len = masc_stats.est_lib_len
    if est_lib_len and masc_stats.cc is not None and est_lib_len <= len(masc_stats.cc):
        masc_est_ll = masc_stats.cc[est_lib_len - 1]
        _annotate_point(
            est_lib_len - 1, "blue",
            upper - height / 2, 'estimated lib len: {}'.format(est_lib_len),
            masc_est_ll, " cc(est lib len) = {:.5f}".format(masc_est_ll), height / 50
        )

    # Annotate expected library length
    library_len = masc_stats.stats.metrics_at_expected_length.fragment_length
    if library_len and masc_stats.cc is not None and library_len <= len(masc_stats.cc):
        masc_ll = masc_stats.cc[library_len - 1]
        _annotate_point(
            library_len - 1, "green",
            upper - height / 1.75, 'expected lib len: {}'.format(library_len),
            masc_ll, " cc(lib length) = {:.5f}".format(masc_ll), -height / 25
        )

    plt.legend(loc="best")

    # Add moving average window size annotation if available
    mv_avr_filter_len = getattr(masc_stats, 'mv_avr_filter_len', None)
    if mv_avr_filter_len:
        _annotate_bottom_right_box("Mov avr win size = {}".format(mv_avr_filter_len))


def plot_ncc_vs_masc(pp: PdfPages, stats: GenomeWideStats, name: str) -> None:
    """Generate comparison plots between NCC and MSCC using new StatisticsResult interface.

    Creates side-by-side comparison plots for each chromosome showing
    differences between naive cross-correlation and MSCC results.

    Args:
        pp: PdfPages object for multi-page PDF output
        stats_result: StatisticsResult object containing analysis results
        name: Sample name for plot titles
    """
    title = "{} Cross-Correlation"
    if name:
        title += " for " + name

    # Plot genome-wide comparison if MSCC is available
    if stats.has_mscc:
        _plot_ncc_vs_masc(stats.whole_ncc_stats, stats.whole_mscc_stats, "Naive CC vs MSCC")
        _feed_pdf_page(pp)

    # Plot per-chromosome comparisons
    for ref in sorted(stats.references):
        try:
            ncc = None if stats.ncc_stats is None else stats.ncc_stats.get(ref)
            mscc = None if stats.mscc_stats is None else stats.mscc_stats.get(ref)
            _plot_ncc_vs_masc(ncc, mscc, title.format(ref))
            _feed_pdf_page(pp)
        except AssertionError:
            logger.debug("Skip plot for {}, valid reads unable.".format(ref))


def _plot_ncc_vs_masc(
    cc_stats: Optional[ChromosomeStats[NCCStats]],
    masc_stats: Optional[ChromosomeStats[MSCCStats]],
    title: str
) -> None:
    """Create individual NCC vs MSCC comparison plot using new UnifiedStats interface.

    Generates a single comparison plot showing both NCC and MSCC
    curves for a specific chromosome or genome-wide data.

    Args:
        stats: UnifiedStats object containing correlation data
        title: Plot title

    Returns:
        None
    """
    # Check if we have valid data to plot
    has_valid_cc = (
        cc_stats is not None and
        cc_stats.cc is not None and
        not np.all(np.isnan(cc_stats.cc))
    )
    has_valid_masc = (
        masc_stats is not None and
        masc_stats.cc is not None and
        not np.all(np.isnan(masc_stats.cc))
    )

    if not (has_valid_cc or has_valid_masc):
        raise AssertionError("No valid correlation data available for plotting")

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Relative Cross-Correlation from each minimum")

    # Determine max_shift from available data
    max_shift = 300  # Default fallback
    if cc_stats and cc_stats.cc is not None:
        max_shift = len(cc_stats.cc) - 1
    elif masc_stats and masc_stats.cc is not None:
        max_shift = len(masc_stats.cc) - 1

    x_range = range(max_shift + 1)

    # Plot NCC if available
    if cc_stats is not None and cc_stats.stats.cc_min is not None:
        plt.plot(x_range, cc_stats.cc - cc_stats.stats.cc_min,
                 color="black", linewidth=0.5, label="Naive CC")

    # Plot MSCC if available
    if masc_stats is not None and masc_stats.stats.cc_min is not None:
        alpha = 1 if not has_valid_cc else 0.8
        plt.plot(x_range, masc_stats.cc - masc_stats.stats.cc_min,
                 alpha=alpha, linewidth=0.5, label="MSCC")

    lower, upper, height = _set_ylim()

    # Annotate read length
    if cc_stats is not None:
        read_len = cc_stats.stats.read_len
    elif masc_stats is not None:
        read_len = masc_stats.stats.read_len
    else:
        raise AssertionError

    _annotate_point(
        read_len, "red",
        upper - height/25, "read length: {}".format(read_len)
    )

    # Annotate estimated library length from MSCC
    if masc_stats is not None:
        est_lib_len = masc_stats.est_lib_len
        if est_lib_len:
            _annotate_point(
                est_lib_len, "blue",
                upper - height/10, "estimated lib len: {}".format(est_lib_len)
            )
        plt.legend(loc="best")

    # Annotate expected library length
    if cc_stats is not None:
        library_len = cc_stats.stats.metrics_at_expected_length.fragment_length
    elif masc_stats is not None:
        library_len = masc_stats.stats.metrics_at_expected_length.fragment_length
    else:
        raise AssertionError

    if library_len:
        _annotate_point(
            library_len, "green",
            upper - height/6, "expected lib len: {}".format(library_len)
        )
