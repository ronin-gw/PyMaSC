import logging
import argparse
import sys
import os
import json

from PyMaSC.utils.parsearg import add_common_args, add_mappability_args, add_shift_arg
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.progress import ProgressBar
from PyMaSC.core.alignability import AlignableLengthCalculator

logger = logging.getLogger(__name__)


def _main():
    # parse aregs
    parser = argparse.ArgumentParser(description="Pre-compute mappability stats for a region file.")
    add_common_args(parser)
    add_mappability_args(parser)
    add_shift_arg(parser)

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
    if args.map_path:
        if args.map_path[0] == args.mappable:
            args.map_path = None
        else:
            args.map_path = args.map_path[0]

    if not args.map_path:
        args.map_path = os.path.splitext(args.mappable)[0] + "_mappability.json"

    if not os.path.isfile(args.mappable):
        logger.critical("Failed to open '{}': no such file.".format(args.mappable))
        sys.exit(1)
    elif not os.access(args.mappable, os.R_OK):
        logger.critical("Failed to open '{}': file is unreadable.".format(args.mappable))
        sys.exit(1)

    if not os.path.exists(args.map_path):
        dirname = os.path.dirname(args.map_path)
        if dirname and not os.access(dirname, os.W_OK):
            logger.critical("Directory is not writable: '{}'".format(dirname))
            sys.exit(1)
    elif not os.path.isfile(args.map_path):
        logger.critical("Specified path is not file: '{}'".format(args.map_path))
        sys.exit(1)
    elif not os.access(args.map_path, os.W_OK):
        logger.critical("Failed to overwrite '{}'".format(args.map_path))
    else:
        logger.warning("Existing file '{}' will be overwritten.".format(args.map_path))

    #
    if sys.stderr.isatty():
        ProgressBar.enable = True

    #
    logger.info("Calcurate mappable length with max shift size {}.".format(args.max_shift))
    alc = AlignableLengthCalculator(args.mappable, args.max_shift)
    alc._calc_alignability()

    #
    logger.info("Save mappable length to '{}'".format(args.map_path))
    with open(args.map_path, 'w') as f:
        json.dump({
            "max_shift": alc.max_shift,
            "__whole__": alc.alignable_len,
            "references": alc.chrom2alignable_len
        }, f, indent=4, sort_keys=True)

    #
    alc.close()
    logger.info("Calc mappable length finished.")


def exec_entrypoint():
    try:
        _main()
    except KeyboardInterrupt:
        sys.stderr.write("\r\033[K")
        sys.stderr.flush()
        logger.info("Got KeyboardInterrupt. bye")
