import logging
import sys
import os.path
import json

from pysam import AlignmentFile

from PyMaSC import entrypoint, logging_version
from PyMaSC.utils.parsearg import get_plot_parser
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.pymasc import prepare_output, PLOTFILE_SUFFIX
from PyMaSC.output.stats import (load_stats, load_cc, load_masc, output_cc, output_mscc, output_stats,
                                 CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, STATSFILE_SUFFIX)
from PyMaSC.handler.result import CCResult
from PyMaSC.output.figure import plot_figures

logger = logging.getLogger(__name__)


def _complete_path_arg(args, attr, val):
    if not getattr(args, attr) and os.path.exists(val):
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

    #
    if not args.stats:
        parser.error("Statistics file path is not specified.")
    elif not args.cc and not args.masc:
        parser.error("Neither cross-correlation table file path nor mappability "
                     "sensitive cross-correlation table file path is specified.")

    #
    if not os.path.exists(args.stats):
        parser.error("Statistics file path does not exist: '{}'".format(args.stats))
    elif all((args.cc and not os.path.exists(args.cc),
              args.masc and not os.path.exists(args.masc))):
        parser.error("Neither cross-correlation table file path '{}' nor "
                     "mappability sensitive cross-correlation table file path "
                     "'{}' exists.".format(args.cc, args.masc))
    if args.cc:
        if not os.path.exists(args.cc):
            parser.error("Cross-correlation table file path does not exists: "
                         "'{}'".format(args.cc))
        elif not args.sizes or not os.path.exists(args.sizes):
            parser.error("Please specify a chromosome sizes file using -s/--sizes option.")
    if args.masc:
        if not os.path.exists(args.masc):
            parser.error("Mappability sensitive cross-correlation table file path "
                         "does not exists: '{}'".format(args.cc))
        else:
            if args.mappability_stats and not args.mappability_stats.endswith(".json"):
                setattr(args, "mappability_stats",
                        os.path.splitext(args.mappability_stats)[0] + "_mappability.json")
            if not args.mappability_stats or not os.path.exists(args.mappability_stats):
                parser.error("Please specify a JSON file which generated by PyMaSC"
                             "for a BigWig file using -m/--mappability-stats option.")

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logging_version(logger)

    #
    return args


@entrypoint(logger)
def main():
    args = _parse_args()
    statattrs = _prepare_stats(args)
    ref2cc, ref2genomelen, ref2masc, ref2mappable_len, references = _load_tables(args)
    checked_suffixes = _prepare_outputs(args)

    #
    ccr = CCResult(
        read_len=statattrs["read_len"],
        references=references,
        ref2genomelen=ref2genomelen,
        ref2cc=ref2cc,
        ref2mappable_len=ref2mappable_len,
        ref2masc=ref2masc
    )

    #
    for outputfunc, suffix in zip((output_stats, output_cc, output_mscc),
                                  (STATSFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX)):
        if suffix in checked_suffixes:
            outputfunc(os.path.join(args.outdir, args.name), ccr)

    #
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
    statattrs["filter_mask_len"] = max(args.mask_size, 0)

    #
    if "name" in statattrs:
        args.name = statattrs.pop("name")
    if not args.name:
        logger.critical("Mandatory attribute 'Name' not found in '{}'. "
                        "Set name manually with -n/--name option.".format(args.stats))
        sys.exit(1)

    return statattrs


def _load_table(path, load_fun):
    try:
        whole, table = load_fun(path)
    except KeyError:
        logger.critical("Failed to load tables.")
        sys.exit(1)
    references = table.keys()
    return whole, table, references


def _load_tables(args):
    if args.cc:
        _, cc_table, references = _load_table(args.cc, load_cc)
        ref2genomelen = _load_chrom_sizes(args.sizes)
        for ref in references:
            if ref not in ref2genomelen:
                logger.critical("Reference '{}' not found in '{}'.".format(ref, args.sizes))
                sys.exit(1)
    else:
        cc_table = ref2genomelen = None

    if args.masc:
        _, masc_table, references = _load_table(args.masc, load_masc)
        try:
            ref2mappable_len = _load_mappable_lengths(args.mappability_stats)
        except json.JSONDecoder as e:
            logger.critical("Failed to load '{}':".format(args.mappability_stats))
            logger.error(e.msg)
            sys.exit(1)
        for ref in references:
            if ref not in ref2mappable_len:
                logger.critical("Reference '{}' not found in '{}'.".format(ref, args.mappability_stats))
                sys.exit(1)
    else:
        masc_table = ref2mappable_len = None

    return cc_table, ref2genomelen, masc_table, ref2mappable_len, references


def _prepare_outputs(args):
    check_suffixes = [PLOTFILE_SUFFIX]
    output_base = os.path.join(args.outdir, args.name)

    for source, suffix in zip((args.stats, args.cc, args.masc),
                              (STATSFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX)):
        if source and _check_overwrite(source, os.path.normpath(output_base + suffix), args.force_overwrite):
            check_suffixes.append(suffix)

    prepare_output([None], [args.name], args.outdir, check_suffixes)

    return check_suffixes


def _check_overwrite(source, output, force_overwrite):
    if os.path.normpath(source) == output:
        if force_overwrite:
            logger.warning("Overwrite option was specified. '{}' will be overwritten!".format(output))
            return True
        else:
            logger.warning("Prevent to overwrite input stats file '{}', "
                           "output stats file will be skipped.".format(output))
            return False
    return True


def _load_chrom_sizes(path):
    try:
        with AlignmentFile(path) as f:
            return {r: l for r, l in zip(f.references, f.lengths)}
    except ValueError:
        ref2len = {}
        with open(path) as f:
            for l in f:
                chrom, length, _ = l.split('\t', 2)
                ref2len[chrom] = int(length)
            return ref2len


def _load_mappable_lengths(path):
    with open(path) as f:
        return json.load(f)["references"]
