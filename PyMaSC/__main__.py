import logging

from logfmt import ColorfulFormatter
from parsearg import get_parser
from reader import get_read_generator_and_init_target, InputUnseekable
from model import CCCalculator, ReadUnsortedError

logger = logging.getLogger(__name__)

rl = logging.getLogger('')
h = logging.StreamHandler()
h.setFormatter(ColorfulFormatter(fmt="[%(asctime)s | %(levelname)s] %(name)10s : %(message)s"))
rl.addHandler(h)


def _main():
    parser = get_parser()
    args = parser.parse_args()

    rl.setLevel(args.log_level)
    logger.info("Calculate cross-correlation shift 1 to {}bp with reads MAOQ >= {}".format(args.max_shift, args.mapq))

    for f in args.reads:
        logger.info("Process {}".format(f))
        ccc = compute_cc(f, args.format, args.max_shift, args.mapq)
        if ccc is None:
            logger.warning("Faild to process {}. Skip this file.".format(f))
        else:
            ccc.finishup_calculation()


    logger.info("PyMASC finished.")


def compute_cc(path, fmt, max_shift, mapq_criteria):
    try:
        parser, stream = get_read_generator_and_init_target(path, fmt)
    except InputUnseekable:
        logger.error("Specify input file format using `-f` option if your input is a stream.")
        return None

    with parser(stream) as alignfile:
        ccc = CCCalculator(max_shift, sum(alignfile.lengths))

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


_main()
