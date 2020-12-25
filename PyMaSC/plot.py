import logging
import sys
import os.path
import json

from pysam import AlignmentFile

from PyMaSC import entrypoint, logging_version
from PyMaSC.utils.parsearg import get_plot_parser
from PyMaSC.utils.logfmt import set_rootlogger
from PyMaSC.utils.calc import filter_chroms
from PyMaSC.pymasc import prepare_output, PLOTFILE_SUFFIX
from PyMaSC.output.stats import load_stats, output_stats, STATSFILE_SUFFIX
from PyMaSC.output.table import (load_cc, load_masc, load_nreads_table, output_cc, output_mscc,
                                 CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, NREADOUTPUT_SUFFIX)
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
        for attr, suffix in zip(("stats", "cc", "masc", "nreads"),
                                (STATSFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX, NREADOUTPUT_SUFFIX)):
            _complete_path_arg(args, attr, args.statfile + suffix)

    #
    if not args.stats:
        parser.error("Statistics file path is not specified.")
    if not args.nreads:
        parser.error("# of reads table file path is not specified.")
    elif not args.cc and not args.masc:
        parser.error("Neither cross-correlation table file path nor mappability "
                     "sensitive cross-correlation table file path is specified.")

    #
    if not os.path.exists(args.stats):
        parser.error("Statistics file path does not exist: '{}'".format(args.stats))
    if not os.path.exists(args.nreads):
        parser.error("# of reads table file path does not exist: '{}'".format(args.nreads))
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
                parser.error("Please specify a JSON file which generated by PyMaSC "
                             "for a BigWig file using -m/--mappability-stats option.")

    #
    if "all" in args.force_overwrite:
        args.force_overwrite = ("stats", "cc", "mscc")

    # set up logging
    set_rootlogger(args.color, args.log_level)
    logging_version(logger)

    #
    return args


@entrypoint(logger)
def main():
    args = _parse_args()
    read_len = _prepare_stats(args)
    (ref2cc, ref2genomelen, ref2masc, ref2mappable_len, references,
     (forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum)) = _load_tables(args)
    references = filter_chroms(references, args.chromfilter)
    checked_suffixes = _prepare_outputs(args)

    #
    ccr = CCResult(
        args.smooth_window, args.chi2_pval, args.mask_size, args.bg_avr_width,
        args.library_length,
        read_len=read_len,
        references=references,
        ref2genomelen=ref2genomelen,
        ref2forward_sum=forward_sum,
        ref2reverse_sum=reverse_sum,
        ref2cc=ref2cc,
        ref2mappable_len=ref2mappable_len,
        mappable_ref2forward_sum=mappable_forward_sum,
        mappable_ref2reverse_sum=mappable_reverse_sum,
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
    try:
        statattrs = load_stats(args.stats, ("name", "library_len", "read_len"))
    except (IOError, ):
        logger.critical("Failed to load stats.")
        sys.exit(1)

    if "library_len" in statattrs and not args.library_length:
        args.library_length = statattrs["library_len"]

    if "read_len" not in statattrs:
        logger.critical("Mandatory attribute 'Read length' not found in '{}'.".format(args.stats))
        sys.exit(1)
    read_len = statattrs["read_len"]

    #
    if "name" in statattrs:
        args.name = statattrs.pop("name")
    if not args.name:
        logger.critical("Mandatory attribute 'Name' not found in '{}'. "
                        "Set name manually with -n/--name option.".format(args.stats))
        sys.exit(1)

    return read_len


def _load_table(path, load_fun):
    try:
        return load_fun(path)
    except (IOError, KeyError, IndexError, StopIteration):
        logger.critical("Failed to load tables.")
        sys.exit(1)


def _load_tables(args):
    if args.cc:
        cc_table = _load_table(args.cc, load_cc)
        cc_references = set(cc_table.keys())
        ref2genomelen = _load_chrom_sizes(args.sizes)
        for ref in cc_references:
            if ref not in ref2genomelen:
                logger.critical("Reference '{}' not found in '{}'.".format(ref, args.sizes))
                sys.exit(1)
    else:
        cc_table = ref2genomelen = None
        cc_references = set()

    if args.masc:
        masc_table = _load_table(args.masc, load_masc)
        masc_references = set(masc_table.keys())
        try:
            ref2mappable_len = _load_mappable_lengths(args.mappability_stats)
        except json.JSONDecoder as e:
            logger.critical("Failed to load '{}':".format(args.mappability_stats))
            logger.error(e.msg)
            sys.exit(1)
        for ref in masc_references:
            if ref not in ref2mappable_len:
                logger.critical("Reference '{}' not found in '{}'.".format(ref, args.mappability_stats))
                sys.exit(1)
    else:
        masc_table = ref2mappable_len = None
        masc_references = set()

    dicts = forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum = _load_table(args.nreads, load_nreads_table)
    if forward_sum and reverse_sum:
        nreads_refs = set(forward_sum.keys())
        assert nreads_refs == set(reverse_sum.keys())
    else:
        nreads_refs = set(mappable_forward_sum.keys())
        assert nreads_refs == set(mappable_reverse_sum.keys())
    nreads_refs.discard("whole")

    refsets = [s for s in (cc_references, masc_references) if s]
    union = nreads_refs.union(*refsets)
    intersection = nreads_refs.intersection(*refsets)

    if union != intersection:
        logger.warning("Chromosome names in tables are unmatched.")
        logger.warning("Trying use common names anyway: {}".format(intersection))

    return cc_table, ref2genomelen, masc_table, ref2mappable_len, sorted(intersection), dicts


def _prepare_outputs(args):
    check_suffixes = [PLOTFILE_SUFFIX]
    output_base = os.path.join(args.outdir, args.name)

    for source, suffix, force_overwrite in zip(
            (args.stats, args.cc, args.masc),
            (STATSFILE_SUFFIX, CCOUTPUT_SUFFIX, MSCCOUTPUT_SUFFIX),
            tuple(n in args.force_overwrite for n in ("stats", "cc", "masc"))):
        if source and _check_overwrite(source, os.path.normpath(output_base + suffix), force_overwrite):
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
                chrom, length = l.split('\t', 2)[:2]
                ref2len[chrom] = int(length)
            return ref2len


def _load_mappable_lengths(path):
    with open(path) as f:
        return json.load(f)["references"]
