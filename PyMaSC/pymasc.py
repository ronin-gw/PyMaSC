"""Main PyMaSC CLI application for cross-correlation analysis.

This module provides the main entry point for PyMaSC.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple, Union
from itertools import zip_longest

from PyMaSC import entrypoint, logging_version
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.parsearg import get_pymasc_parser
from PyMaSC.utils.progress import ProgressBase
from PyMaSC.utils.output import prepare_outdir
from PyMaSC.handler.mappability import MappabilityHandler, BWIOError, JSONIOError
from PyMaSC.handler.unified import InputUnseekable
from PyMaSC.handler.base import NothingToCalc
from PyMaSC.handler.unified import UnifiedCalcHandler
from PyMaSC.core.models import (
    CalculationConfig, ExecutionConfig, CalculationTarget, 
    ImplementationAlgorithm, ExecutionMode, MappabilityConfig
)
from PyMaSC.handler.result import CCResult, ReadsTooFew
from PyMaSC.core.ncc import ReadUnsortedError
from PyMaSC.output.stats import output_stats, STATSFILE_SUFFIX
from PyMaSC.output.table import (output_cc, output_mscc, output_nreads_table,
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

    #
    if sys.stderr.isatty() and not args.disable_progress:
        ProgressBase.global_switch = True

    #
    suffixes: List[str] = list(EXPECT_OUTFILE_SUFFIXES)
    if args.mappability:
        if args.skip_ncc:
            suffixes.remove(CCOUTPUT_SUFFIX)
    else:
        suffixes.remove(MSCCOUTPUT_SUFFIX)
    if args.skip_plots:
        suffixes.remove(PLOTFILE_SUFFIX)
    basenames = prepare_output(args.reads, args.name, args.outdir, tuple(suffixes))

    #
    calc_handlers = make_handlers(args)
    if not calc_handlers:
        return None

    #
    max_readlen = set_readlen(args, calc_handlers)

    #
    mappability_handler: Optional[MappabilityHandler] = None
    if args.mappability:
        try:
            mappability_handler = MappabilityHandler(
                args.mappability, args.max_shift, max_readlen,
                args.mappability_stats, args.process
            )
        except (BWIOError, JSONIOError):
            # Error logs already generated in MappabilityHandler
            sys.exit(1)

        for handler in calc_handlers:
            handler.set_mappability_handler(mappability_handler)

    #
    logger.info("Calculate cross-correlation between 0 to {} base shift"
                "with reads MAOQ >= {}".format(args.max_shift, args.mapq))
    for handler, output_basename in zip(calc_handlers, basenames):
        result_handler = run_calculation(args, handler, output_basename)
        output_results(args, output_basename, result_handler)

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


def make_handlers(args: argparse.Namespace) -> List[UnifiedCalcHandler]:
    """Create calculation handlers for input files.

    Instantiates UnifiedCalcHandler with appropriate configuration objects
    based on algorithm selection and input file specifications.

    Args:
        args: Parsed command-line arguments

    Returns:
        List of successfully created calculation handlers

    Note:
        Failed handlers are logged but not included in returned list
    """
    calc_handlers: List[UnifiedCalcHandler] = []

    # Create configuration objects from arguments
    for f in args.reads:
        try:
            # Determine calculation target and implementation based on existing logic
            # Implementation: Use successive if --successive flag, otherwise BitArray (default)
            implementation = ImplementationAlgorithm.SUCCESSIVE if args.successive else ImplementationAlgorithm.BITARRAY
            
            # Target: Auto-determine based on mappability presence (existing behavior)
            if hasattr(args, 'mappability') and args.mappability:
                # Mappability provided: calculate both NCC and MSCC unless --skip-ncc
                target = CalculationTarget.MSCC if args.skip_ncc else CalculationTarget.BOTH
            else:
                # No mappability: only NCC
                target = CalculationTarget.NCC
            
            # Create calculation configuration using new conceptual approach
            calc_config = CalculationConfig(
                target=target,
                implementation=implementation,
                max_shift=args.max_shift,
                mapq_criteria=args.mapq,
                skip_ncc=args.skip_ncc,
                esttype=args.readlen_estimator,
                chromfilter=args.chromfilter
            )

            # Create execution configuration
            exec_config = ExecutionConfig(
                mode=ExecutionMode.MULTI_PROCESS if args.process > 1 else ExecutionMode.SINGLE_PROCESS,
                worker_count=args.process
            )

            # Create mappability configuration if needed
            # Both SUCCESSIVE and BITARRAY algorithms support mappability
            mappability_config = None
            if hasattr(args, 'mappability') and args.mappability:
                mappability_config = MappabilityConfig(
                    mappability_path=args.mappability,
                    mappability_stats_path=args.mappability_stats if getattr(args, 'mappability_stats', None) else None
                )

            # Create handler
            handler = UnifiedCalcHandler(f, calc_config, exec_config, mappability_config)
            calc_handlers.append(handler)

        except ValueError:
            logger.error("Failed to open file '{}'".format(f))
        except NothingToCalc:
            logger.error("Check your -i/--include-chrom and/or -e/--exclude-chrom options.")
        except InputUnseekable:
            logger.error("If your input can't reread, specify read length using `-r` option.")

    return calc_handlers


def set_readlen(args: argparse.Namespace, calc_handlers: List[UnifiedCalcHandler]) -> int:
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
    #
    if args.read_length is None:
        logger.info("Check read length: Get {} from read length distribution".format(args.readlen_estimator.lower()))

    #
    readlens: List[int] = []
    for i, handler in enumerate(calc_handlers[:]):
        try:
            handler.set_or_estimate_readlen(args.read_length)
        except ValueError:
            calc_handlers.pop(i)
            continue
        if handler.read_len is not None:
            readlens.append(handler.read_len)

    #
    max_readlen = max(readlens)
    if len(set(readlens)) != 1:
        logger.warning("There are multiple read length candidates. Use max length "
                       "({}) for MSCC calculation.".format(max_readlen))
    return max_readlen


def run_calculation(args: argparse.Namespace, handler: UnifiedCalcHandler, output_basename: Path) -> Optional[CCResult]:
    """Execute cross-correlation calculation for a single file.

    Runs the main calculation workflow and creates a result handler
    with the computed cross-correlation data.

    Args:
        args: Parsed command-line arguments
        handler: Calculation handler for the input file
        output_basename: Base path for output files

    Returns:
        CCResult object containing calculation results, or None if failed
    """
    logger.info("Process {}".format(handler.path))

    try:
        handler.run_calcuration()
    except ReadUnsortedError:
        return None

    try:
        return CCResult(
            args.smooth_window, args.chi2_pval, args.mask_size, args.bg_avr_width,
            args.library_length, handler
        )
    except ReadsTooFew:
        logger.warning("Faild to process {}. Skip this file.".format(handler.path))
        return None


def output_results(args: argparse.Namespace, output_basename: Path, result_handler: Optional[CCResult]) -> None:
    """Generate all output files and plots from calculation results.

    Creates statistics files, data tables, and plots based on the
    calculation results and user-specified options.

    Args:
        args: Parsed command-line arguments
        output_basename: Base path for output files
        result_handler: CCResult object containing calculation results

    Note:
        Plot generation may be skipped if matplotlib is not available
    """
    if result_handler is None:
        return

    output_stats(output_basename, result_handler)
    output_nreads_table(output_basename, result_handler)
    if not result_handler.skip_ncc:
        output_cc(output_basename, result_handler)
    if result_handler.calc_masc:
        output_mscc(output_basename, result_handler)
    if not args.skip_plots:
        plotfile_path = Path(str(output_basename) + PLOTFILE_SUFFIX)
        try:
            from PyMaSC.output.figure import plot_figures
        except ImportError:
            logger.error("Skip output plots '{}'".format(plotfile_path))
        else:
            plot_figures(plotfile_path, result_handler)
