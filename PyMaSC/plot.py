import logging
import sys
import os.path
from collections import namedtuple

from PyMaSC import entrypoint, logging_version
from PyMaSC.utils.parsearg import get_plot_parser
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.pymasc import prepare_output, PLOTFILE_SUFFIX
from PyMaSC.output.stats import (load_stats, load_cc, load_masc, output_stats,
                                 CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, STATSFILE_SUFFIX)
from PyMaSC.handler.result import PyMaSCStats, chi2_test
from PyMaSC.output.figure import plot_figures

logger = logging.getLogger(__name__)

CCResult = namedtuple("CCResult", ("references", "skip_ncc", "calc_masc", "ref2stats", "whole"))


def _complete_path_arg(args, attr, val):
    if not getattr(args, attr):
        setattr(args, attr, val)


def _parse_args():
    # parse args
    parser = get_plot_parser()
    args = parser.parse_args()

    #
    if args.statfile:
        for attr, suffix in zip(("stats", "cc", "masc"),
                                (STATSFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX)):
            _complete_path_arg(args, attr, args.statfile + suffix)
    else:
        if not args.stats:
            parser.error("Statistics file path is not specified.")
        elif not args.cc and not args.masc:
            parser.error("Neither cross-correlation table file path nor mappability "
                         "sensitive cross-correlation table file path is specified.")

    #
    if not os.path.exists(args.stats):
        parser.error("Statistics file path does not exist: '{}'".format(args.stats))
    elif all((args.cc, args.masc,
              not os.path.exists(args.cc), not os.path.exists(args.masc))):
        parser.error("Neither cross-correlation table file path '{}' nor "
                     "mappability sensitive cross-correlation table file path "
                     "'{}' exists.".format(args.cc, args.masc))
    elif args.cc and not os.path.exists(args.cc):
        parser.error("Cross-correlation table file path does not exists: "
                     "'{}'".format(args.cc))
    elif args.masc and not os.path.exists(args.masc):
        parser.error("Mappability sensitive cross-correlation table file path "
                     "does not exists: '{}'".format(args.cc))

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logging_version(logger)

    #
    return args


@entrypoint(logger)
def main():
    args = _parse_args()
    statattrs = _prepare_stats(args)
    cc_whole, cc_table, masc_whole, masc_table, references = _load_tables(args)
    checked_suffixes = _prepare_outputs(args)

    #
    if "forward_sum" in statattrs and "reverse_sum" in statattrs:
        chi2_test(statattrs["forward_sum"], statattrs["reverse_sum"],
                  args.chi2_pval, "Whole genome")

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
    if STATSFILE_SUFFIX in checked_suffixes:
        output_stats(os.path.join(args.outdir, args.name), ccr)
    plot_figures(os.path.join(args.outdir, args.name + PLOTFILE_SUFFIX), ccr)


def _prepare_stats(args):
    #
    statattrs = load_stats(args.stats)
    if "library_len" in statattrs:
        statattrs["expected_library_len"] = statattrs["library_len"]
        del statattrs["library_len"]
    else:
        statattrs["expected_library_len"] = None
    if args.library_length:
        statattrs["expected_library_len"] = args.library_length
    if "read_len" not in statattrs:
        logger.critical("Mandatory attribute 'Read length' not found in '{}'.".format(args.stats))
        sys.exit(1)
    if args.smooth_window:
        statattrs["filter_len"] = args.smooth_window

    #
    if "name" in statattrs:
        args.name = statattrs.pop("name")
    if not args.name:
        logger.critical("Mandatory attribute 'Name' not found in '{}'. "
                        "Set name manually with -n/--name option.".format(args.stats))
        sys.exit(1)

    return statattrs


def _load_table(path, load_fun):
    whole, table = load_fun(path)
    references = table.keys()
    return whole, table, references


def _load_tables(args):
    if args.cc:
        cc_whole, cc_table, references = _load_table(args.cc, load_cc)
    else:
        cc_whole = cc_table = None
    if args.masc:
        masc_whole, masc_table, references = _load_table(args.masc, load_masc)
    else:
        masc_whole = masc_table = None

    return cc_whole, cc_table, masc_whole, masc_table, references


def _prepare_outputs(args):
    check_suffixes = [PLOTFILE_SUFFIX]
    out_stats_path = os.path.normpath(os.path.join(args.outdir, args.name) + STATSFILE_SUFFIX)
    if os.path.normpath(args.stats) == out_stats_path:
        logger.error("Prevent to overwrite input stats file '{}', "
                     "output stats file will be skipped.".format(out_stats_path))
    else:
        check_suffixes.append(STATSFILE_SUFFIX)
    prepare_output([None], [args.name], args.outdir, check_suffixes)

    return check_suffixes
