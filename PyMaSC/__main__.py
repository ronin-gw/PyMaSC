import logging
import os
import sys

from logfmt import ColorfulFormatter
from parsearg import get_parser
from reader import get_read_generator_and_init_target, InputUnseekable
from model import CCCalculator, ReadUnsortedError
from output import EXPECT_OUTFILE_SUFFIXES, get_output_basename, output_result

logger = logging.getLogger(__name__)

rl = logging.getLogger('')
h = logging.StreamHandler()
h.setFormatter(ColorfulFormatter(fmt="[%(asctime)s | %(levelname)s] %(name)10s : %(message)s"))
rl.addHandler(h)


def _main():
    parser = get_parser()
    args = parser.parse_args()

    rl.setLevel(args.log_level)
    if args.progress:
        CCCalculator.need_progress_bar = True

    prepare_output(args)

    logger.info("Calculate cross-correlation between 1 to {} base shift with reads MAOQ >= {}".format(args.max_shift, args.mapq))
    for f in args.reads:
        logger.info("Process {}".format(f))

        ccc = compute_cc(f, args.format, args.max_shift, args.mapq)
        if ccc is None:
            logger.warning("Faild to process {}. Skip this file.".format(f))
        else:
            ccc.finishup_calculation()

        output_result(f, ccc, args.outdir)

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
            expect_outfile = get_output_basename(args.outdir, f) + suffix
            if os.path.exists(expect_outfile):
                logger.warning("Existing file '{}' will be overwritten.".format(expect_outfile))


def compute_cc(path, fmt, max_shift, mapq_criteria):
    try:
        parser, stream = get_read_generator_and_init_target(path, fmt)
    except InputUnseekable:
        logger.error("Specify input file format using `-f` option if your input is a stream.")
        return None

    with parser(stream) as alignfile:
        ccc = CCCalculator(max_shift, alignfile.references, alignfile.lengths)

        for i, (is_reverse, is_duplicate, chrom, pos, mapq, readlen) in enumerate(alignfile):
            if mapq < mapq_criteria or is_duplicate:
                continue

            try:
                if is_reverse:
                    ccc.feed_reverse_read(chrom, pos, readlen)
                else:
                    ccc.feed_forward_read(chrom, pos, readlen)
            except ReadUnsortedError:
                logger.error("Input read must be sorted.")
                return None

    return ccc


try:
    _main()
except KeyboardInterrupt:
    logger.info("Got KeyboardInterrupt. bye")
