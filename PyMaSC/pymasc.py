import logging
import os
import sys

from PyMaSC import VERSION
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.parsearg import get_parser
from PyMaSC.utils.progress import ProgressBar
from PyMaSC.reader.align import InputUnseekable
from PyMaSC.handler.mappability import MappabilityHandler, BWIOError, JSONIOError
from PyMaSC.handler.masc import CCCalcHandler
from PyMaSC.handler.result import CCResult, ReadsTooFew
from PyMaSC.output.stats import output_cc, output_stats
from PyMaSC.output.figure import plot_figures

logger = logging.getLogger(__name__)

PLOTFILE_SUFFIX = ".pdf"
CCOUTPUT_SUFFIX = "_cc.tab"
STATSFILE_SUFFIX = "_stats.tab"
EXPECT_OUTFILE_SUFFIXES = (PLOTFILE_SUFFIX, CCOUTPUT_SUFFIX, STATSFILE_SUFFIX)


def _get_output_basename(dirpath, filepath):
    return os.path.join(dirpath, os.path.splitext(os.path.basename(filepath))[0])


def _main():
    # parse args
    parser = get_parser()
    args = parser.parse_args()

    # set up logging
    if args.color == "TRUE":
        colorize = True
    elif args.color == "FALSE":
        colorize = False
    else:
        colorize = sys.stderr.isatty()

    set_rootlogger(colorize, args.log_level)
    logger.info("PyMaSC version " + VERSION)

    # check args
    if args.mappable:
        args.mappable = args.mappable[0]
    if args.map_path and args.map_path == args.mappable:
        args.map_path = None
    if args.library_length and args.library_length > args.max_shift:
        logger.error("Specified expected library length > max shift. Ignore expected length setting.")
        args.library_length = None
    #
    if sys.stderr.isatty():
        ProgressBar.enable = True

    #
    prepare_output(args)

    #
    calc_handlers = []
    need_logging_unseekable_error = False
    for f in args.reads:
        try:
            calc_handlers.append(
                CCCalcHandler(f, args.format, args.estimation_type, args.max_shift, args.mapq, args.process)
            )
        except InputUnseekable:
            need_logging_unseekable_error = True
    if need_logging_unseekable_error:
        logger.error("If your input can't seek back, specify input file format and "
                     "read length using `-f` and `-r` option.")

    #
    readlens = []
    if args.read_length is None:
        logger.info("Check read length: Get {} from read length distribution".format(args.estimation_type.lower()))
    for handler in calc_handlers:
        handler.set_readlen(args.read_length)
        readlens.append(handler.read_len)
    max_readlen = max(readlens)

    #
    mappability_handler = None
    if args.mappable:
        try:
            mappability_handler = MappabilityHandler(
                args.mappable, args.max_shift, max_readlen, args.map_path, args.process
            )
        except (BWIOError, JSONIOError):
            sys.exit(1)

        for handler in calc_handlers:
            handler.set_mappability_handler(mappability_handler)

    #
    logger.info("Calculate cross-correlation between 1 to {} base shift"
                "with reads MAOQ >= {}".format(args.max_shift, args.mapq))
    for handler in calc_handlers:
        logger.info("Process {}".format(handler.path))

        handler.run_calcuration()
        try:
            result_handler = CCResult(handler, args.smooth_window, args.chi2_pval)
        except ReadsTooFew:
            result_handler = None
        if result_handler is None:
            logger.warning("Faild to process {}. Skip this file.".format(f))
            continue
        output_result(handler.path, result_handler, args.outdir)

    #
    if mappability_handler:
        if mappability_handler.need_save_stats:
            mappability_handler.save_mappability_stats()
        mappability_handler.close()
    logger.info("PyMASC finished.")


def prepare_output(args):
    if os.path.exists(args.outdir):
        if not os.path.isdir(args.outdir):
            logger.critical("Specified path as a output directory is not directory.")
            logger.critical(args.outdir)
            sys.exit(1)
    else:
        logger.info("Make output directory: {}".format(args.outdir))
        try:
            os.mkdir(args.outdir)
        except IOError as e:
            logger.critical("Faild to make output directory: [Errno {}] {}".format(e.errno, e.message))
            sys.exit(1)

    if not os.access(args.outdir, os.W_OK):
        logger.critical("Output directory '{}' is not writable.".format(args.outdir))
        sys.exit(1)

    for f in args.reads:
        for suffix in EXPECT_OUTFILE_SUFFIXES:
            expect_outfile = _get_output_basename(args.outdir, f) + suffix
            if os.path.exists(expect_outfile):
                logger.warning("Existing file '{}' will be overwritten.".format(expect_outfile))


def output_result(sourcepath, ccc, outdir):
    outfile_prefix = _get_output_basename(outdir, sourcepath)
    logger.info("Output results to '{}'".format(outfile_prefix))

    output_cc(outfile_prefix + CCOUTPUT_SUFFIX, ccc)
    output_stats(outfile_prefix + STATSFILE_SUFFIX, ccc)
    plot_figures(outfile_prefix + PLOTFILE_SUFFIX, ccc, os.path.basename(sourcepath))


def exec_entrypoint():
    try:
        _main()
    except KeyboardInterrupt:
        sys.stderr.write("\r\033[K")
        sys.stderr.flush()
        logger.info("Got KeyboardInterrupt. bye")
