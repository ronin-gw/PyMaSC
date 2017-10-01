import logging
import os
import sys

from PyMaSC.utils.compatible import zip_longest

from PyMaSC import VERSION
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.parsearg import get_parser
from PyMaSC.utils.progress import ProgressBar, ReadCountProgressBar, MultiLineProgressManager
from PyMaSC.handler.mappability import MappabilityHandler, BWIOError, JSONIOError
from PyMaSC.handler.masc import CCCalcHandler, InputUnseekable
from PyMaSC.handler.result import CCResult, ReadsTooFew
from PyMaSC.core.ncc import ReadUnsortedError
from PyMaSC.output.stats import output_cc, output_mscc, output_stats
from PyMaSC.output.figure import plot_figures

logger = logging.getLogger(__name__)

PLOTFILE_SUFFIX = ".pdf"
CCOUTPUT_SUFFIX = "_cc.tab"
MSCCOUTPUT_SUFFIX = "_mscc.tab"
STATSFILE_SUFFIX = "_stats.tab"
EXPECT_OUTFILE_SUFFIXES = (PLOTFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, STATSFILE_SUFFIX)


def _get_output_basename(dirpath, filepath):
    return os.path.join(dirpath, os.path.splitext(os.path.basename(filepath))[0])


def _main():
    # parse args
    parser = get_parser()
    args = parser.parse_args()

    if args.skip_ncc and args.mappability is None:
        parser.error("argument --skip-ncc: -m/--mappable must be specified.")

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
    if args.mappability:
        args.mappability = args.mappability[0]
    if args.mappability_stats and args.mappability_stats == args.mappability:
        args.mappability_stats = None
    if args.library_length and args.library_length > args.max_shift:
        logger.error("Specified expected library length > max shift. Ignore expected length setting.")
        args.library_length = None
    #
    if sys.stderr.isatty():
        ProgressBar.enable = True
        ReadCountProgressBar.enable = True
        MultiLineProgressManager.enable = True
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
    calc_handlers = []
    need_logging_unseekable_error = False
    for f in args.reads:
        try:
            calc_handlers.append(
                CCCalcHandler(f, args.estimation_type, args.max_shift, args.mapq,
                              args.process, args.skip_ncc)
            )
        except ValueError:
            logger.error("Failed to open file '{}'".format(f))
        except InputUnseekable:
            need_logging_unseekable_error = True
    if need_logging_unseekable_error:
        logger.error("If your input can't reread, specify read length using `-r` option.")
    if not calc_handlers:
        return None

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
    if args.mappability:
        try:
            mappability_handler = MappabilityHandler(
                args.mappability, args.max_shift, max_readlen, args.mappability_stats, args.process
            )
        except (BWIOError, JSONIOError):
            sys.exit(1)

        for handler in calc_handlers:
            handler.set_mappability_handler(mappability_handler)

    #
    logger.info("Calculate cross-correlation between 1 to {} base shift"
                "with reads MAOQ >= {}".format(args.max_shift, args.mapq))
    for handler, output_basename in zip(calc_handlers, basenames):
        logger.info("Process {}".format(handler.path))

        try:
            handler.run_calcuration()
        except ReadUnsortedError:
            continue

        try:
            result_handler = CCResult(handler, args.smooth_window, args.chi2_pval, args.library_length, args.skip_ncc)
        except ReadsTooFew:
            result_handler = None
        if result_handler is None:
            logger.warning("Faild to process {}. Skip this file.".format(f))
            continue

        #
        output_cc(output_basename + CCOUTPUT_SUFFIX, result_handler)
        output_mscc(output_basename + MSCCOUTPUT_SUFFIX, result_handler)
        output_stats(output_basename + STATSFILE_SUFFIX, result_handler)
        if not args.skip_plots:
            plot_figures(output_basename + PLOTFILE_SUFFIX, result_handler)

    #
    if mappability_handler:
        if mappability_handler.need_save_stats:
            mappability_handler.save_mappability_stats()
        mappability_handler.close()


def prepare_output(reads, names, outdir, suffixes=EXPECT_OUTFILE_SUFFIXES):
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            logger.critical("Specified path as a output directory is not directory.")
            logger.critical(outdir)
            sys.exit(1)
    else:
        logger.info("Make output directory: {}".format(outdir))
        try:
            os.mkdir(outdir)
        except IOError as e:
            logger.critical("Faild to make output directory: [Errno {}] {}".format(e.errno, e.message))
            sys.exit(1)

    if not os.access(outdir, os.W_OK):
        logger.critical("Output directory '{}' is not writable.".format(outdir))
        sys.exit(1)

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


def exec_entrypoint():
    try:
        _main()
        logger.info("PyMASC finished.")
    except KeyboardInterrupt:
        sys.stderr.write("\r\033[K")
        sys.stderr.flush()
        logger.info("Got KeyboardInterrupt. bye")
