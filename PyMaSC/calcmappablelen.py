import logging
import sys

from PyMaSC import entrypoint, logging_version
from PyMaSC.utils.parsearg import get_precalc_parser
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.progress import ProgressBase
from PyMaSC.handler.mappability import MappabilityHandler

logger = logging.getLogger(__name__)


@entrypoint(logger)
def main():
    # parse aregs
    parser = get_precalc_parser()

    args = parser.parse_args()
    if not args.mappability:
        parser.error("argument -m/--mappable: expected 1 argument(s)")

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logging_version(logger)

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
