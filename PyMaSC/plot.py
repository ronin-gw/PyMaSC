import logging
import sys
import os.path
from collections import namedtuple

from PyMaSC import VERSION
from PyMaSC.utils.parsearg import get_plot_parser
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.pymasc import prepare_output, PLOTFILE_SUFFIX
from PyMaSC.output.stats import (load_stats, load_cc, load_masc, output_stats,
                                 CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, STATSFILE_SUFFIX)
from PyMaSC.handler.result import PyMaSCStats
from PyMaSC.output.figure import plot_figures

logger = logging.getLogger(__name__)

CCResult = namedtuple("CCResult", ("references", "skip_ncc", "calc_masc", "ref2stats", "whole"))


def _main():
    # parse args
    parser = get_plot_parser()
    args = parser.parse_args()

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logger.info("PyMaSC version {} with Python{}.{}.{}".format(
                *[VERSION] + list(sys.version_info[:3])))
    for line in sys.version.split('\n'):
        logger.debug(line)

    #
    stat_path = cctable_path = masctable_path = None
    if args.statfile:
        stat_path = args.statfile + STATSFILE_SUFFIX
        cctable_path = args.statfile + CCOUTPUT_SUFFIX
        masctable_path = args.statfile + MSCCOUTPUT_SUFFIX
    stat_path = args.stats if args.stats else stat_path
    cctable_path = args.cc if args.cc else cctable_path
    masctable_path = args.masc if args.masc else masctable_path

    #
    if not stat_path:
        parser.error("Statistics file path is not specified.")
    elif not os.path.exists(stat_path):
        parser.error("Statistics file path does not exist: '{}'".format(stat_path))
    elif not cctable_path and not masctable_path:
        parser.error("Neither cross-correlation table file path nor mappability "
                     "sensitive cross-correlation table file path is specified.")
    elif all((cctable_path, masctable_path,
              not os.path.exists(cctable_path), not os.path.exists(masctable_path))):
        parser.error("Neither cross-correlation table file path '{}' nor "
                     "mappability sensitive cross-correlation table file path "
                     "'{}' exists.".format(cctable_path, masctable_path))
    elif cctable_path and not os.path.exists(cctable_path):
        parser.error("Cross-correlation table file path does not exists: "
                     "'{}'".format(cctable_path))
    elif masctable_path and not os.path.exists(masctable_path):
        parser.error("Mappability sensitive cross-correlation table file path "
                     "does not exists: '{}'".format(cctable_path))

    #
    statattrs = load_stats(stat_path)
    if "library_len" in statattrs:
        statattrs["expected_library_len"] = statattrs["library_len"]
        del statattrs["library_len"]
    else:
        statattrs["expected_library_len"] = None
    if args.library_length:
        statattrs["expected_library_len"] = args.library_length
    if "read_len" not in statattrs:
        logger.critical("")
        sys.exit(1)
    if args.smooth_window:
        statattrs["filter_len"] = args.smooth_window
    #
    if cctable_path:
        cc_whole, cc_table = load_cc(cctable_path)
        references = cc_table.keys()
    else:
        cc_whole = cc_table = None
    if masctable_path:
        masc_whole, masc_table = load_masc(masctable_path)
        references = masc_table.keys()
    else:
        masc_whole = masc_table = None

    #
    name = None
    if "name" in statattrs:
        name = statattrs.pop("name")
    name = args.name if args.name else name
    if name is None:
        logger.critical("")
        sys.exit(1)
    prepare_output([None], [name], args.outdir, [PLOTFILE_SUFFIX, STATSFILE_SUFFIX])

    #
    ccr = CCResult(
        references=references,
        skip_ncc=not bool(cc_table),
        calc_masc=bool(masc_table),
        whole=PyMaSCStats(cc=cc_whole, masc=masc_whole, **statattrs),
        ref2stats={
            ref: PyMaSCStats(
                cc=None if cc_table is None else cc_table[ref],
                masc=None if masc_table is None else masc_table[ref],
                **statattrs
            ) for ref in references
        }
    )

    #
    output_stats(os.path.join(args.outdir, name), ccr)
    plot_figures(os.path.join(args.outdir, name + PLOTFILE_SUFFIX), ccr)


def exec_entrypoint():
    try:
        _main()
        logger.info("PyMASC plot finished.")
    except KeyboardInterrupt:
        sys.stderr.write("\r\033[K")
        sys.stderr.flush()
        logger.info("Got KeyboardInterrupt. bye")
