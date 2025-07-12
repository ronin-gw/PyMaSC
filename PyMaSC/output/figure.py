"""Matplotlib-based plotting functions for PyMaSC visualization.

This module provides comprehensive plotting functionality for PyMaSC analysis results,
including cross-correlation plots, MSCC plots, and comparison visualizations.
All plots are generated as multi-page PDF files for easy viewing and publication.

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

from PyMaSC.utils.output import catch_IOError

logger = logging.getLogger(__name__)

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except:
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
def plot_figures(outfile: os.PathLike[str], stats_result: Any) -> None:
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
        if not stats_result.skip_ncc:
            plot_naive_cc(stats_result, name)
            _feed_pdf_page(pp)

        if stats_result.has_mappability:
            if plot_naive_cc_just(stats_result, name):
                _feed_pdf_page(pp)

            plot_masc(stats_result, name)
            _feed_pdf_page(pp)

        plot_ncc_vs_masc(pp, stats_result, name)


def _annotate_point(x: Union[int, float], color: str, axis_y: Union[int, float], axis_text: str, point_y: Optional[Union[int, float]] = None, point_text: Optional[str] = None, yoffset: Union[int, float] = 0) -> None:
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
    if point_y:
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


def _annotate_params(nsc: Optional[Union[int, float]] = None, rsc: Optional[Union[int, float]] = None, est_nsc: Optional[Union[int, float]] = None, est_rsc: Optional[Union[int, float]] = None, loc: str = "lower right") -> None:
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


def plot_naive_cc(stats_result: Any, name: Optional[str] = None, xlim: Optional[Tuple[int, Any]] = None) -> None:
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

    genome_stats = stats_result.genome_wide_stats
    cc_stats = genome_stats.cc if genome_stats else None
    
    if cc_stats is None:
        logger.warning("No CC statistics available for plotting")
        return

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Cross-Correlation")

    # Plot confidence intervals if available (skip if not available in StatisticsResult)
    max_shift = getattr(genome_stats, 'max_shift', len(cc_stats.cc) - 1 if cc_stats.cc is not None else 300)
    x_range = range(max_shift + 1)
    
    # Plot main cross-correlation curve
    if cc_stats.cc is not None:
        plt.plot(x_range, cc_stats.cc, color="black", linewidth=0.5)
    
    axes = plt.gca()
    if xlim:
        axes.set_xlim(xlim)
    lower, upper, height = _set_ylim()

    # Plot minimum correlation line
    if hasattr(cc_stats, 'cc_min') and cc_stats.cc_min is not None:
        plt.axhline(cc_stats.cc_min, linestyle="dashed", linewidth=0.5)
        plt.text(0, cc_stats.cc_min, 'min(cc) = {:.5f}'.format(cc_stats.cc_min))

    # Annotate read length
    read_len = getattr(genome_stats, 'read_len', None) or getattr(cc_stats, 'read_len', None)
    if read_len and hasattr(cc_stats, 'ccrl') and cc_stats.ccrl is not None:
        _annotate_point(
            read_len - 1, "red",
            upper - height / 25, 'read length: {}'.format(read_len),
            cc_stats.ccrl, " cc(read length) = {:.5f}".format(cc_stats.ccrl), height / 50
        )

    # Annotate estimated library length
    if hasattr(cc_stats, 'est_lib_len') and cc_stats.est_lib_len:
        est_ccfl = getattr(cc_stats, 'est_ccfl', None)
        if est_ccfl is not None:
            _annotate_point(
                cc_stats.est_lib_len - 1, "blue",
                upper - height / 10, 'estimated lib len: {}'.format(cc_stats.est_lib_len),
                est_ccfl, " cc(est lib len) = {:.5f}".format(est_ccfl), height / 50
            )
    
    # Annotate expected library length
    library_len = getattr(genome_stats, 'library_len', None) or getattr(cc_stats, 'library_len', None)
    if library_len and hasattr(cc_stats, 'ccfl') and cc_stats.ccfl is not None:
        _annotate_point(
            library_len - 1, "green",
            upper - height / 6, 'expected lib len: {}'.format(library_len),
            cc_stats.ccfl, " cc(lib length) = {:.5f}".format(cc_stats.ccfl), -height / 25
        )
    
    # Add quality metrics annotations
    _annotate_params(
        getattr(cc_stats, 'nsc', None),
        getattr(cc_stats, 'rsc', None),
        getattr(cc_stats, 'est_nsc', None),
        getattr(cc_stats, 'est_rsc', None)
    )


def plot_naive_cc_just(stats_result: Any, name: Optional[str] = None) -> bool:
    """Plot zoomed naive cross-correlation around fragment length peak using new StatisticsResult interface.

    Creates a focused plot showing the cross-correlation peak region
    with detailed annotations for fragment length estimation.

    Args:
        stats_result: StatisticsResult object containing analysis results
        name: Sample name for plot title

    Returns:
        True if plot was created, False if insufficient data
    """
    genome_stats = stats_result.genome_wide_stats
    if genome_stats is None:
        return False
        
    cc_stats = genome_stats.cc
    if cc_stats is None:
        return False
        
    # Check if we can create a zoomed plot
    calc_ncc = not stats_result.skip_ncc
    max_shift = getattr(genome_stats, 'max_shift', 300)
    est_lib_len = getattr(cc_stats, 'est_lib_len', None)
    
    if calc_ncc and est_lib_len and est_lib_len * 2 < max_shift + 1:
        plot_naive_cc(stats_result, name, (0, est_lib_len * 2))
        return True
    return False


def plot_masc(stats_result: Any, name: Optional[str] = None) -> None:
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

    genome_stats = stats_result.genome_wide_stats
    if genome_stats is None:
        logger.warning("No genome-wide statistics available for MSCC plotting")
        return
        
    masc_stats = genome_stats.masc
    if masc_stats is None:
        logger.warning("No MSCC statistics available for plotting")
        return

    max_shift = getattr(genome_stats, 'max_shift', len(masc_stats.cc) - 1 if masc_stats.cc is not None else 300)
    x_range = range(max_shift + 1)

    # Plot MSCC curves
    if masc_stats.cc is not None:
        plt.plot(x_range, masc_stats.cc, color="black", linewidth=0.5, label="MSCC")
    
    # Plot smoothed MSCC if available
    if hasattr(masc_stats, 'avr_cc') and masc_stats.avr_cc is not None:
        plt.plot(x_range, masc_stats.avr_cc, alpha=0.8, label="Smoothed", color="pink")

    lower, upper, height = _set_ylim()

    # Annotate estimated library length
    est_lib_len = getattr(masc_stats, 'est_lib_len', None)
    if est_lib_len and masc_stats.cc is not None and est_lib_len <= len(masc_stats.cc):
        masc_est_ll = masc_stats.cc[est_lib_len - 1]
        _annotate_point(
            est_lib_len - 1, "blue",
            upper - height / 2, 'estimated lib len: {}'.format(est_lib_len),
            masc_est_ll, " cc(est lib len) = {:.5f}".format(masc_est_ll), height / 50
        )

    # Annotate expected library length
    library_len = getattr(genome_stats, 'library_len', None) or getattr(masc_stats, 'library_len', None)
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


def plot_ncc_vs_masc(pp: Any, stats_result: Any, name: str) -> None:
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
    if stats_result.has_mappability:
        _plot_ncc_vs_masc(stats_result.genome_wide_stats, "Naive CC vs MSCC")
        _feed_pdf_page(pp)

    # Plot per-chromosome comparisons
    for ref in sorted(stats_result.references):
        chrom_stats = stats_result.chromosome_stats.get(ref)
        if chrom_stats is not None:
            try:
                _plot_ncc_vs_masc(chrom_stats, title.format(ref))
                _feed_pdf_page(pp)
            except AssertionError:
                logger.debug("Skip plot for {}, valid reads unable.".format(ref))


def _plot_ncc_vs_masc(stats: Any, title: str) -> None:
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
    cc_stats = getattr(stats, 'cc', None)
    masc_stats = getattr(stats, 'masc', None)
    
    has_valid_cc = (cc_stats is not None and 
                   hasattr(cc_stats, 'cc') and cc_stats.cc is not None and 
                   not np.all(np.isnan(cc_stats.cc)))
    has_valid_masc = (masc_stats is not None and 
                     hasattr(masc_stats, 'cc') and masc_stats.cc is not None and 
                     not np.all(np.isnan(masc_stats.cc)))
    
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
    if has_valid_cc and hasattr(cc_stats, 'cc_min') and cc_stats.cc_min is not None:
        plt.plot(x_range, cc_stats.cc - cc_stats.cc_min,
                 color="black", linewidth=0.5, label="Naive CC")
    
    # Plot MSCC if available
    if has_valid_masc and hasattr(masc_stats, 'cc_min') and masc_stats.cc_min is not None:
        alpha = 1 if not has_valid_cc else 0.8
        plt.plot(x_range, masc_stats.cc - masc_stats.cc_min,
                 alpha=alpha, linewidth=0.5, label="MSCC")

    lower, upper, height = _set_ylim()

    # Annotate read length
    read_len = getattr(stats, 'read_len', None)
    if read_len:
        _annotate_point(
            read_len, "red",
            upper - height/25, "read length: {}".format(read_len)
        )

    # Annotate estimated library length from MSCC
    if has_valid_masc:
        est_lib_len = getattr(masc_stats, 'est_lib_len', None)
        if est_lib_len:
            _annotate_point(
                est_lib_len, "blue",
                upper - height/10, "estimated lib len: {}".format(est_lib_len)
            )
        plt.legend(loc="best")

    # Annotate expected library length
    library_len = getattr(stats, 'library_len', None)
    if library_len:
        _annotate_point(
            library_len, "green",
            upper - height/6, "expected lib len: {}".format(library_len)
        )
