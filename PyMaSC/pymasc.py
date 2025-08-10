"""Main PyMaSC CLI application for cross-correlation analysis.

This module provides the main entry point for PyMaSC.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional, Tuple, cast
from itertools import zip_longest

from . import entrypoint, logging_version
from .utils.logfmt import set_rootlogger
from .utils.parsearg import get_pymasc_parser
from .utils.progress import ProgressBase
from .utils.output import prepare_outdir
from .handler.mappability import MappabilityHandler, BWIOError, JSONIOError
from .handler.calc import CalcHandler
from .interfaces.config import PyMaSCConfig, mappability_configured, stats_configured
from .interfaces.stats import GenomeWideStats
from .stats import make_genome_wide_stat
from .core.exceptions import ReadUnsortedError, ReadsTooFew, InputUnseekable, NothingToCalc
from .output.stats import output_stats, STATSFILE_SUFFIX
from .output.table import (output_cc, output_mscc, output_nreads_table,
                           CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, NREADOUTPUT_SUFFIX)


logger = logging.getLogger(__name__)

PLOTFILE_SUFFIX = ".pdf"
EXPECT_OUTFILE_SUFFIXES: Tuple[str, ...] = (PLOTFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, NREADOUTPUT_SUFFIX, STATSFILE_SUFFIX)


def _parse_args() -> argparse.Namespace:
    """Parse and validate command-line arguments.

    Handles argument parsing, validation, and setup of logging configuration.
    Performs consistency checks between arguments and applies default values.

    Returns:
        Parsed and validated argument namespace

    Raises:
        SystemExit: If argument validation fails
    """
    parser = get_pymasc_parser()
    args = parser.parse_args()

    if args.skip_ncc and args.mappability is None:
        parser.error("argument --skip-ncc: -m/--mappable must be specified.")

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logging_version(logger)

    # check args
    if args.mappability_stats and args.mappability_stats == args.mappability:
        args.mappability_stats = None
    if args.library_length and args.library_length > args.max_shift:
        logger.error("Specified expected library length > max shift. Ignore expected length setting.")
        args.library_length = None

    return args


@entrypoint(logger)
def main() -> None:
    """Main PyMaSC application entry point.

    Orchestrates the complete PyMaSC workflow:
    1. Parse command-line arguments
    2. Prepare output directories and file validation
    3. Initialize calculation handlers
    4. Set up mappability correction if requested
    5. Execute cross-correlation calculations
    6. Generate output files and plots

    Returns:
        None on successful completion
    """
    args = _parse_args()
    config = PyMaSCConfig.from_args(args)

    if sys.stderr.isatty() and not args.disable_progress:
        ProgressBase.global_switch = True

    suffixes: List[str] = list(EXPECT_OUTFILE_SUFFIXES)
    if args.mappability:
        if args.skip_ncc:
            suffixes.remove(CCOUTPUT_SUFFIX)
    else:
        suffixes.remove(MSCCOUTPUT_SUFFIX)
    if args.skip_plots:
        suffixes.remove(PLOTFILE_SUFFIX)
    basenames = prepare_output(args.reads, args.name, args.outdir, tuple(suffixes))

    calc_handlers: List[CalcHandler] = []
    for f in args.reads:
        try:
            # Create handler
            handler = CalcHandler(f, config)
            calc_handlers.append(handler)

        except ValueError:
            logger.error("Failed to open file '{}'".format(f))
        except NothingToCalc:
            logger.error("Check your -i/--include-chrom and/or -e/--exclude-chrom options.")
        except InputUnseekable:
            logger.error("If your input can't reread, specify read length using `-r` option.")

    if not calc_handlers:
        return None

    #
    readlen = set_readlen(args, calc_handlers)
    config.read_length = readlen

    #
    mappability_handler: Optional[MappabilityHandler] = None
    if mappability_configured(config):
        try:
            # Use config-based creation for cleaner interface
            mappability_handler = MappabilityHandler.from_config(config)
        except (BWIOError, JSONIOError):
            # Error logs already generated in MappabilityHandler
            sys.exit(1)

        for handler in calc_handlers:
            handler.set_mappability_handler(mappability_handler)
    config = cast(PyMaSCConfig, config)

    #
    logger.info("Calculate cross-correlation between 0 to {} base shift "
                "with reads MAPQ >= {}".format(args.max_shift, args.mapq))
    for handler, output_basename in zip(calc_handlers, basenames):
        result = run_calculation(config, handler, output_basename)
        output_results(args, output_basename, result)

    #
    if mappability_handler:
        mappability_handler.save_mappability_stats()
        mappability_handler.close()


def prepare_output(reads: List[str], names: List[Optional[str]], outdir: str, suffixes: Tuple[str, ...] = EXPECT_OUTFILE_SUFFIXES) -> List[Path]:
    """Prepare output directory and generate output basenames.

    Creates output directory if needed and generates output basenames for each
    input file. Warns about existing files that will be overwritten.

    Args:
        reads: List of input BAM file paths
        names: List of custom output names (can contain None values)
        outdir: Output directory path
        suffixes: Expected output file suffixes to check for conflicts

    Returns:
        List of output basenames as Path objects (without extensions)

    Raises:
        SystemExit: If output directory cannot be created
    """
    #
    if not prepare_outdir(outdir, logger):
        # Error logs already generated in prepare_outdir
        sys.exit(1)

    #
    basenames: List[Path] = []
    for f, n in zip_longest(reads, names):
        if n is None:
            output_basename = Path(outdir) / Path(f).stem
        else:
            output_basename = Path(outdir) / n

        for suffix in suffixes:
            expect_outfile = Path(str(output_basename) + suffix)
            if expect_outfile.exists():
                logger.warning("Existing file '{}' will be overwritten.".format(expect_outfile))
        basenames.append(output_basename)

    return basenames


def set_readlen(args: argparse.Namespace, calc_handlers: List[CalcHandler]) -> int:
    """Set or estimate read length for all calculation handlers.

    Either uses user-specified read length or estimates it from the data.
    Handles multiple files with different read lengths by using the maximum.

    Args:
        args: Parsed command-line arguments
        calc_handlers: List of calculation handlers to configure

    Returns:
        Maximum read length across all handlers

    Note:
        Handlers that fail read length setting are removed from the list
    """
    if args.read_length is not None:
        for handler in calc_handlers:
            handler.read_len = args.read_length
        return args.read_length

    # Estimate read length from data
    logger.info("Check read length: Get {} from read length distribution".format(args.readlen_estimator.lower()))
    readlens: List[int] = []
    for i, handler in enumerate(calc_handlers[:]):
        try:
            readlens.append(handler.estimate_readlen())
        except ValueError:
            calc_handlers.pop(i)
            continue

    #
    max_readlen = max(readlens)
    if len(set(readlens)) != 1:
        logger.warning("There are multiple read length candidates. Use max length "
                       "({}) for MSCC calculation.".format(max_readlen))

    for handler in calc_handlers:
        handler.read_len = max_readlen

    return max_readlen


def run_calculation(config: PyMaSCConfig, handler: CalcHandler, output_basename: Path) -> Optional[GenomeWideStats]:
    """Execute cross-correlation calculation for a single file using new statistics system.

    Runs the main calculation workflow and creates a result handler
    with the computed cross-correlation data using the new StatisticsResult
    class and UnifiedStats system, bypassing legacy compatibility layers.

    Args:
        args: Parsed command-line arguments
        handler: Calculation handler for the input file
        output_basename: Base path for output files

    Returns:
        StatisticsResult object containing calculation results, or None if failed
    """
    logger.info("Process {}".format(handler.path))

    try:
        result = handler.run_calculation()
    except ReadUnsortedError:
        return None

    try:
        # Use new statistics system with container-based data access
        assert stats_configured(config), "Statistics configuration must be set for new stats system"
        return make_genome_wide_stat(
            result,
            config,
            output_warnings=True
        )
    except ReadsTooFew:
        logger.warning("Failed to process {}. Skip this file.".format(handler.path))
        return None


def output_results(args: argparse.Namespace, output_basename: Path, result: Optional[GenomeWideStats]) -> None:
    """Generate all output files and plots from calculation results using new statistics system.

    Creates statistics files, data tables, and plots based on the
    calculation results using the new StatisticsResult interface directly.
    All output modules have been updated to work with StatisticsResult.

    Args:
        args: Parsed command-line arguments
        output_basename: Base path for output files
        result_handler: StatisticsResult object containing calculation results
    """
    if result is None:
        return

    output_stats(output_basename, result)
    output_nreads_table(output_basename, result)
    if result.whole_ncc_stats is not None:
        output_cc(output_basename, result)
    if result.whole_mscc_stats is not None:
        output_mscc(output_basename, result)
    if not args.skip_plots:
        plotfile_path = Path(str(output_basename) + PLOTFILE_SUFFIX)
        try:
            from PyMaSC.output.figure import plot_figures
        except ImportError:
            logger.error("Skip output plots '{}'".format(plotfile_path))
        else:
            plot_figures(plotfile_path, result)
