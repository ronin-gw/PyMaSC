import logging
import os
import sys

from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.parsearg import get_parser
from PyMaSC.utils.progress import ProgressBar
from PyMaSC.reader.align import InputUnseekable
from PyMaSC.handler.alignability import AlignabilityHandler, BWIOError, JSONIOError
from PyMaSC.handler.masc import CCCalcHandler
from PyMaSC.output import output_cc, output_stats, plot_figures

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

    # check args
    if args.mappable:
        args.mappable = args.mappable[0]
    if args.map_path and args.map_path == args.mappable:
        args.map_path = None

    #
    if sys.stderr.isatty():
        ProgressBar.enable = True

    #
    prepare_output(args)

    #
    if args.mappable:
        try:
            alh = AlignabilityHandler(args.mappable, args.max_shift, args.map_path)
        except (BWIOError, JSONIOError):
            sys.exit(1)
    else:
        alh = None

    #
    logger.info("Calculate cross-correlation between 1 to {} base shift with reads MAOQ >= {}".format(args.max_shift, args.mapq))
    for f in args.reads:
        logger.info("Process {}".format(f))

        try:
            cch = CCCalcHandler(f, args.format, args.max_shift, args.mapq, alh, args.chi2_pval)
        except InputUnseekable:
            logger.error("Specify input file format using `-f` option if your input can't seek back.")
            ccr = None

        ccr = cch.run_calcuration()

        if ccr is None:
            logger.warning("Faild to process {}. Skip this file.".format(f))
            continue

        output_result(f, ccc, args.outdir)

    #
    if alh:
        if alh.need_save_stats:
            alh.save_mappability_stats()
        alh.close()
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

    try:
        output_cc(outfile_prefix + CCOUTPUT_SUFFIX, ccc)
    except IOError as e:
        logger.error("Faild to output cc table: {}\n[Errno {}] {}".format(
                     e.filename, e.errno, e.message))

    try:
        output_stats(outfile_prefix + STATSFILE_SUFFIX, ccc)
    except IOError as e:
        logger.error("Faild to output stats: {}\n[Errno {}] {}".format(
                     e.filename, e.errno, e.message))

    try:
        plot_figures(outfile_prefix + PLOTFILE_SUFFIX, ccc, os.path.basename(sourcepath))
    except IOError as e:
        logger.error("Faild to output figure: {}\n[Errno {}] {}".format(
                     e.filename, e.errno, e.message))


def exec_entrypoint():
    try:
        _main()
    except KeyboardInterrupt:
        sys.stderr.write("\r\033[K")
        sys.stderr.flush()
        logger.info("Got KeyboardInterrupt. bye")
