"""PyMaSC pre-calculation utility for mappability statistics.

This module provides the pymasc-precalc command-line utility for pre-calculating
mappability statistics from BigWig files. Pre-calculation allows for faster
subsequent PyMaSC runs by avoiding repeated computation of mappability data.

The pre-calculation process:
- Reads mappability values from BigWig files
- Calculates mappable lengths for each chromosome
- Saves statistics to JSON files for later use
- Supports multiprocessing for faster computation
"""
import logging
import sys

from PyMaSC import entrypoint, logging_version
from PyMaSC.utils.parsearg import get_precalc_parser
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.progress import ProgressBase
from PyMaSC.handler.mappability import MappabilityHandler

logger = logging.getLogger(__name__)


@entrypoint(logger)
def main() -> None:
    """Main pre-calculation workflow for mappability statistics.

    Orchestrates the complete pre-calculation workflow:
    1. Parse and validate command-line arguments
    2. Configure logging and progress display
    3. Initialize mappability handler
    4. Calculate mappability statistics
    5. Save results to JSON file

    The pre-calculated statistics can be used in subsequent PyMaSC
    runs to avoid re-computing mappability data.
    """
    # parse aregs
    parser = get_precalc_parser()

    args = parser.parse_args()
    if not args.mappability:
        parser.error("argument -m/--mappable: expected 1 argument(s)")

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logging_version(logger)

    # check args
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
