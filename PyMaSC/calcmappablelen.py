import logging
import argparse
import sys

from PyMaSC import VERSION
from PyMaSC.utils.parsearg import add_common_args, add_mappability_args, add_shift_arg, ForceNaturalNumber
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.progress import ProgressBase
from PyMaSC.handler.mappability import MappabilityHandler

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
    if not args.mappability:
        parser.error("argument -m/--mappable: expected 1 argument(s)")

    # set up logging
    if args.color == "TRUE":
        colorize = True
    elif args.color == "FALSE":
        colorize = False
    else:
        colorize = sys.stderr.isatty()

    set_rootlogger(colorize, args.log_level)
    logger.info("PyMaSC version {} with Python{}.{}.{}".format(
                *[VERSION] + list(sys.version_info[:3])))
    logger.debug(sys.version)

    # check args
    args.mappability = args.mappability[0]
    if args.mappability_stats and args.mappability_stats == args.mappability:
        args.mappability_stats = None

    #
    if sys.stderr.isatty() and not args.disable_progress:
        ProgressBase.global_switch = True

    #
    alh = MappabilityHandler(args.mappability, args.max_shift, args.max_readlen,
                             args.mappability_stats, args.process)
    alh.calc_mappability()
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
