import argparse
import logging

import PyMaSC
from PyMaSC.handler.result import NEAR_READLEN_ERR_CRITERION

READLEN_ESTIMATION_TYPES = ("MEAN", "MEDIAN", "MODE", "MIN", "MAX")
EPILOG = (" \nVisit PyMaSC web site for more information and to get human genome "
          "mappability tracks\n" + PyMaSC.WEBSITE_URL + '\n ')


def _make_upper(s):
    return s.upper()


class StoreLoggingLevel(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, getattr(logging, values))


class ForceNaturalNumber(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values < 1:
            parser.error("argument {} must be > 0.".format('/'.join(self.option_strings)))
        setattr(namespace, self.dest, values)


def make_multistate_append_action(key):
    class _MultistateAppendAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            args = getattr(namespace, self.dest)
            args = [] if args is None else args
            args.append((key, values))
            setattr(namespace, self.dest, args)

    return _MultistateAppendAction


def add_common_args(parser):
    parser.add_argument(
        "-v", "--log-level", type=_make_upper, default=logging.INFO,
        action=StoreLoggingLevel, choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Set verbosity. (Default: INFO)"
    )
    parser.add_argument(
        "--disable-progress", action="store_true",
        help="Disable progress bar"
    )
    parser.add_argument(
        "--color", type=_make_upper, default=None, choices=("TRUE", "FALSE"),
        help="Coloring log. (Default: auto)"
    )
    parser.add_argument(
        "--version", action="version", version="PyMaSC " + PyMaSC.VERSION
    )


def add_multiprocess_args(group):
    group.add_argument(
        "-p", "--process", type=int, default=1, action=ForceNaturalNumber,
        help="Number of worker process. (Default: 1)"
    )


def add_mappability_args(group):
    group.add_argument(
        "-m", "--mappability", metavar="REGION_FILE",
        help="BigWig format mappable region file."
    )
    group.add_argument(
        "--mappability-stats",
        help="Read/Save path for mappability stats. "
             "(Default: [REGION_FILE]_mappability.json)"
    )


def add_shift_arg(group):
    group.add_argument(
        "-d", "--max-shift", type=int, action=ForceNaturalNumber, default=1000,
        help="PyMaSC calculate CC with reverse strand shift from 1 to [MAX_SHIFT] bases. (Default: 1000)"
    )


def add_liblen_arg(group):
    group.add_argument(
        "-l", "--library-length", type=int, action=ForceNaturalNumber,
        help="Your expected library length for input sample(s)."
    )


def add_chrom_filter_args(group):
    group.add_argument(
        "-i", "--include-chrom", nargs='+', dest="chromfilter", metavar="CHROM",
        action=make_multistate_append_action(True),
        help="Include chromosomes to calculate. You can use Unix shell-style "
             "wildcards ('.', '*', '[]' and '[!]'). This option can be declared "
             "multiple times to include chromosomes specified in a just before "
             "-e/--exclude-chrom option. Note that this option is case-sensitive."
    )
    group.add_argument(
        "-e", "--exclude-chrom", nargs='+', dest="chromfilter", metavar="CHROM",
        action=make_multistate_append_action(False),
        help="Exclude chromosomes from calculation. You can use Unix shell-style "
             "wildcards ('.', '*', '[]' and '[!]'). This option can be declared "
             "multiple times to exclude chromosomes specified in a just before "
             "-i/--include-chrom option. Note that this option is case-sensitive."
    )


def add_result_proc_args(group):
    group.add_argument(
        "--chi2-pval", type=float, default=0.05,
        help="p-value threshold for Chi-squared test to check strand specificity "
             "using number of reads mapped to each strands. (Default: 0.05)"
    )
    group.add_argument(
        "-w", "--smooth-window", type=int, default=15, action=ForceNaturalNumber,
        help="Moving average window size for smoothing MSCC "
             "to estimate library length. (Default: 15)"
    )
    group.add_argument(
        "--mask-size", type=int, default=NEAR_READLEN_ERR_CRITERION,
        help="If difference between a read length and the estimated library length "
             "is equal or less than the length specified by this option, "
             "PyMaSC masks correlation coefficients in the read length +/- specified length "
             "and try to estimate mean library length again. (Default: {}, Specify < 1 to disable)".format(NEAR_READLEN_ERR_CRITERION)
    )
    group.add_argument(
        "--bg-avr-width", type=int, action=ForceNaturalNumber, default=50,
        help="The minimum of coefficient will be calcurated as the median of the end of specified bases. (Default: 50bp)"
    )


def get_pymasc_parser():
    parser = argparse.ArgumentParser(
        description="Estimation and visualization tool for library length, "
                    "NSC and RSC metrics with\nmappability sensitive cross-correlation calculation.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    add_common_args(parser)

    proc_args = parser.add_argument_group("Processing behaviors")
    add_multiprocess_args(proc_args)
    proc_args.add_argument(
        "--successive", action="store_true",
        help="Calc with successive algorithm instead of bit array implementation"
    )
    proc_args.add_argument(
        "--skip-ncc", action="store_true",
        help="Skip naive cross-correlation calculation. Mappability file must be specified."
    )
    proc_args.add_argument(
        "--skip-plots", action="store_true",
        help="Skip output figures."
    )

    input_args = parser.add_argument_group("Input alignment file arguments")
    input_args.add_argument(
        "reads", nargs="+",
        help="SAM/BAM format mapped reads. Input must be sorted by positions."
    )
    input_args.add_argument(
        "-r", "--read-length", type=int, action=ForceNaturalNumber,
        help="Set read length manually and disable read length estimation."
    )
    input_args.add_argument(
        "--readlen-estimator", type=_make_upper,
        default="MEDIAN", choices=READLEN_ESTIMATION_TYPES,
        help="Select a representative value used to estimate a read length from "
             "observed read lengths. "
             "Choices: mean, median, mode, min, max (Default: median)"
    )
    add_liblen_arg(input_args)

    map_args = parser.add_argument_group("Input mappability file arguments")
    add_mappability_args(map_args)

    filter = parser.add_argument_group("Input file filtering arguments")
    filter.add_argument(
        "-q", "--mapq", type=int, default=1,
        help="Filter out reads which have less than specified "
             "SAM mapping quality score. (Default: 1)"
    )
    add_chrom_filter_args(filter)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_shift_arg(proc_params)
    add_result_proc_args(proc_params)

    output = parser.add_argument_group("Output file arguments")
    output.add_argument(
        "-n", "--name", nargs='*', default=[],
        help="Output file base name(s). (Default: input file name without extension)"
    )
    output.add_argument(
        "-o", "--outdir", default='.',
        help="Output directory. (Default: current directory)"
    )

    return parser


def get_precalc_parser():
    parser = argparse.ArgumentParser(
        description="Pre-calculate mappability region statistics for "
                    "PyMaSC successive algorithm.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    add_common_args(parser)

    proc_args = parser.add_argument_group("Processing behaviors")
    add_multiprocess_args(proc_args)

    map_args = parser.add_argument_group("Input mappability file arguments")
    add_mappability_args(map_args)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_shift_arg(proc_params)
    proc_params.add_argument(
        "-r", "--max-readlen", type=int, action=ForceNaturalNumber, default=1000,
        help="Set max read length to calculate mappable region length."
    )

    return parser


def get_plot_parser():
    parser = argparse.ArgumentParser(
        description="Plot figures from PyMaSC statistic outputs.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    add_common_args(parser)

    input_args = parser.add_argument_group("Input alignment file arguments")
    input_args.add_argument(
        "statfile", nargs='?',
        help="A base path to the statistic files (*_stats.tab, *_cc.tab and *_masc.tab) "
             "to plot figures."
    )
    input_args.add_argument(
        "--stats",
        help="For specifying path to a statistic file (*_stats.tab) separately."
    )
    input_args.add_argument(
        "--cc",
        help="For specifying path to a cross-correlation table file (*_cc.tab) separately."
    )
    input_args.add_argument(
        "--masc",
        help="For specifying path to a mappability sensitive cross-correlation file (*_stats.tab) separately."
    )
    input_args.add_argument(
        "--nreads",
        help="For specifying path to a # of reads file (*_nreads.tab) separately."
    )
    input_args.add_argument(
        "-s", "--sizes",
        help="A file to obtain length of chromosomes. Tab delimited files,"
             "like *.chrom.sizes and *.fai files, or SAM/BAM format files are acceptable."
    )
    input_args.add_argument(
        "-m", "--mappability-stats",
        help="A JSON file to obtain mappable length of chromosomes generated by PyMaSC"
             "for a BigWig file."
    )

    filter = parser.add_argument_group("Chromosome filtering arguments")
    add_chrom_filter_args(filter)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_result_proc_args(proc_params)
    add_liblen_arg(proc_params)

    output = parser.add_argument_group("Output file arguments")
    output.add_argument(
        "-n", "--name",
        help="Change output file base name. (Default: same as name field in input)"
    )
    output.add_argument(
        "-o", "--outdir", default='.',
        help="Output directory. (Default: current directory)"
    )
    output.add_argument(
        "-f", "--force-overwrite", nargs="*", type=str.lower, choices=('all', 'stats', 'cc', 'mscc'), default=[],
        help="Overwrite specified files even if input and output path are same. "
             "(choices: 'all', 'stats', 'cc', 'mscc' (multiple choices are acceptable))"
    )

    return parser
