import logging
import argparse
import sys

from PyMaSC.utils.parsearg import add_common_args, add_mappability_args, add_shift_arg, ForceNaturalNumber
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.progress import ProgressBar
from PyMaSC.handler.alignability import AlignabilityHandler

logger = logging.getLogger(__name__)


def _main():
    # parse aregs
    parser = argparse.ArgumentParser(description="Pre-compute mappability stats for PyMaSC.")
    add_common_args(parser)
    add_mappability_args(parser)
    add_shift_arg(parser)
    parser.add_argument("-r", "--max-readlen", nargs='?', type=int, action=ForceNaturalNumber, default=1000,
                        help="Set max read length to calculate mappable region length.")

    args = parser.parse_args()
    if not args.mappable:
        parser.error("argument -m/--mappable: expected 1 argument(s)")

    # set up logging
    if args.color == "TRUE":
        colorize = True
    elif args.color == "FALSE":
        colorize = False
    else:
        colorize = sys.stderr.isatty()

    set_rootlogger(colorize, args.log_level)

    # check args
    args.mappable = args.mappable[0]
    if args.map_path and args.map_path == args.mappable:
        args.map_path = None

    #
    if sys.stderr.isatty():
        ProgressBar.enable = True

    #
    logger.info("Calcurate mappable length with max shift size {}.".format(args.max_shift))
    alh = AlignabilityHandler(args.mappable, max(args.max_readlen, args.max_shift), args.map_path)
    alh._calc_alignability()
    alh.save_mappability_stats()
    alh.close()

    logger.info("Calc mappable length finished.")


def exec_entrypoint():
    try:
        _main()
    except KeyboardInterrupt:
        sys.stderr.write("\r\033[K")
        sys.stderr.flush()
        logger.info("Got KeyboardInterrupt. bye")
