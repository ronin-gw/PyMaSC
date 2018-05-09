import logging
import os
import sys

from PyMaSC.utils.compatible import zip_longest

from PyMaSC import entrypoint, logging_version
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.parsearg import get_pymasc_parser
from PyMaSC.utils.progress import ProgressBase
from PyMaSC.utils.output import prepare_outdir
from PyMaSC.handler.mappability import MappabilityHandler, BWIOError, JSONIOError
from PyMaSC.handler.masc import CCCalcHandler, InputUnseekable
from PyMaSC.handler.bamasc import BACalcHandler
from PyMaSC.handler.result import CCResult, ReadsTooFew
from PyMaSC.core.ncc import ReadUnsortedError
from PyMaSC.output.stats import (output_cc, output_mscc, output_stats,
                                 CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, STATSFILE_SUFFIX)

logger = logging.getLogger(__name__)

PLOTFILE_SUFFIX = ".pdf"
EXPECT_OUTFILE_SUFFIXES = (PLOTFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, STATSFILE_SUFFIX)


def _get_output_basename(dirpath, filepath):
    return os.path.join(dirpath, os.path.splitext(os.path.basename(filepath))[0])


def _parse_args():
    parser = get_pymasc_parser()
    args = parser.parse_args()

    if args.skip_ncc and args.mappability is None:
        parser.error("argument --skip-ncc: -m/--mappable must be specified.")

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logging_version(logger)

    # check args
    if args.mappability:
        args.mappability = args.mappability[0]
    if args.mappability_stats and args.mappability_stats == args.mappability:
        args.mappability_stats = None
    if args.library_length and args.library_length > args.max_shift:
        logger.error("Specified expected library length > max shift. Ignore expected length setting.")
        args.library_length = None

    return args


@entrypoint(logger)
def main():
    args = _parse_args()

    #
    if sys.stderr.isatty() and not args.disable_progress:
        ProgressBase.global_switch = True

    #
    suffixes = list(EXPECT_OUTFILE_SUFFIXES)
    if args.mappability:
        if args.skip_ncc:
            suffixes.remove(CCOUTPUT_SUFFIX)
    else:
        suffixes.remove(MSCCOUTPUT_SUFFIX)
    if args.skip_plots:
        suffixes.remove(PLOTFILE_SUFFIX)
    basenames = prepare_output(args.reads, args.name, args.outdir, suffixes)

    #
    calc_handlers = make_handlers(args)
    if not calc_handlers:
        return None

    #
    max_readlen = set_readlen(args, calc_handlers)

    #
    mappability_handler = None
    if args.mappability:
        try:
            mappability_handler = MappabilityHandler(
                args.mappability, args.max_shift, max_readlen,
                args.mappability_stats, args.process
            )
        except (BWIOError, JSONIOError):
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


def prepare_output(reads, names, outdir, suffixes=EXPECT_OUTFILE_SUFFIXES):
    #
    if not prepare_outdir(outdir, logger):
        sys.exit(1)

    #
    basenames = []
    for f, n in zip_longest(reads, names):
        if n is None:
            output_basename = os.path.join(outdir, os.path.splitext(os.path.basename(f))[0])
        else:
            output_basename = os.path.join(outdir, n)

        for suffix in suffixes:
            expect_outfile = output_basename + suffix
            if os.path.exists(expect_outfile):
                logger.warning("Existing file '{}' will be overwritten.".format(expect_outfile))
        basenames.append(output_basename)

    return basenames


def make_handlers(args):
    handler_class = CCCalcHandler if args.successive else BACalcHandler
    calc_handlers = []
    for f in args.reads:
        try:
            calc_handlers.append(
                handler_class(f, args.readlen_estimator, args.max_shift, args.mapq,
                              args.process, args.skip_ncc)
            )
        except ValueError:
            logger.error("Failed to open file '{}'".format(f))
        except InputUnseekable:
            logger.error("If your input can't reread, specify read length using `-r` option.")

    return calc_handlers


def set_readlen(args, calc_handlers):
    #
    if args.read_length is None:
        logger.info("Check read length: Get {} from read length distribution".format(args.readlen_estimator.lower()))

    #
    readlens = []
    for i, handler in enumerate(calc_handlers[:]):
        try:
            handler.set_readlen(args.read_length)
        except ValueError:
            calc_handlers.pop(i)
            continue
        readlens.append(handler.read_len)

    #
    max_readlen = max(readlens)
    if len(set(readlens)) != 1:
        logger.warning("There are multiple read length candidates. Use max length "
                       "() for MSCC calculation.".format(max_readlen))
    return max_readlen


def run_calculation(args, handler, output_basename):
    logger.info("Process {}".format(handler.path))

    try:
        handler.run_calcuration()
    except ReadUnsortedError:
        return

    try:
        return CCResult(
            handler, args.smooth_window, args.chi2_pval, args.library_length
        )
    except ReadsTooFew:
        logger.warning("Faild to process {}. Skip this file.".format(handler.path))


def output_results(args, output_basename, result_handler):
    output_stats(output_basename, result_handler)
    if not result_handler.skip_ncc:
        output_cc(output_basename, result_handler)
    if result_handler.calc_masc:
        output_mscc(output_basename, result_handler)
    if not args.skip_plots:
        plotfile_path = output_basename + PLOTFILE_SUFFIX
        try:
            from PyMaSC.output.figure import plot_figures
        except ImportError:
            logger.error("Skip output plots '{}'".format(plotfile_path))
        else:
            plot_figures(plotfile_path, result_handler)
