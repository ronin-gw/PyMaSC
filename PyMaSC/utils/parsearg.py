import argparse
import logging

import PyMaSC

READLEN_ESTIMATION_TYPES = ("MEAN", "MEDIAN", "MODE", "MIN", "MAX")


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


def add_common_args(parser):
    parser.add_argument("-p", "--process", nargs='?', type=int, action=ForceNaturalNumber, default=1,
                        help="Number of worker process. (Default: 1)")
    parser.add_argument("-v", "--log-level", nargs='?', type=_make_upper, default=logging.INFO, action=StoreLoggingLevel,
                        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
                        help="Set verbosity. (Default: INFO)")
    parser.add_argument("--disable-progress", action="store_true",
                        help="Disable progress bar")
    parser.add_argument("--color", nargs='?', type=_make_upper, default=None, choices=("TRUE", "FALSE"),
                        help="Coloring log. (Default: auto)")
    parser.add_argument("--version", action="version", version="PyMaSC " + PyMaSC.VERSION)


def add_mappability_args(group):
    group.add_argument("-m", "--mappability", nargs=1, metavar="REGION_FILE",
                       help="BigWig format mappable region file.")
    group.add_argument("--mappability-stats", nargs='?',
                       help="Read/Save path for mappability stats. (Default: [REGION_FILE].json)")


def add_shift_arg(group):
    group.add_argument("-d", "--max-shift", nargs='?', type=int, action=ForceNaturalNumber, default=1000,
                       help="PyMaSC calculate CC with reverse strand shift from 1 to [MAX_SHIFT] base. (Default: 1000)")


def get_parser():
    parser = argparse.ArgumentParser(
        description="Estimation and visualization tool for library length, "
                    "NSC and RSC metrics with mappability sensitive cross-correlation calculation."
    )

    add_common_args(parser)

    input_args = parser.add_argument_group("Input file arguments")
    input_args.add_argument("reads", nargs="+",
                            help="SAM/BAM format mapped reads. Input must be sorted.")
    input_args.add_argument("-r", "--read-length", nargs='?', type=int, action=ForceNaturalNumber,
                            help="Set read length manually and disable read length estimation.")
    input_args.add_argument("--estimation-type", nargs='?', type=_make_upper, default="MEDIAN", choices=READLEN_ESTIMATION_TYPES,
                            help="Define the statistic used to estimate a read length from observed read lengths. "
                                 "Choices: mean, median, mode, min, max (Default: median)")
    add_mappability_args(input_args)

    params = parser.add_argument_group("Parameters")
    params.add_argument("-l", "--library-length", nargs='?', type=int, action=ForceNaturalNumber,
                        help="Expected library length to clac NSC and RSC if library length estimation couldn't perform.")
    add_shift_arg(params)
    params.add_argument("-q", "--mapq", nargs='?', type=int, default=1,
                        help="Filter out reads which have less than specified MAPQ. (Default: 1)")
    params.add_argument("--chi2-pval", nargs='?', type=float, default=0.05,
                        help="p-value threshold for Chi-squared test to check strand specificity. (Default: 0.05)")
    params.add_argument("-w", "--smooth-window", nargs='?', type=int, action=ForceNaturalNumber, default=30,
                        help="Moving average window size for smoothing MSCC to estimate library length. (Default: 30)")
    params.add_argument("--skip-ncc", action="store_true",
                        help="Skip naive cross-correlation calculation. Mappability region file must be specified.")

    output = parser.add_argument_group("Output file arguments")
    output.add_argument("-n", "--name", nargs='*', default=[],
                        help="Output file base name(s). (Default: input file name without extension)")
    output.add_argument("-o", "--outdir", nargs='?', default='.',
                        help="Output directory. (Default: current directory)")
    output.add_argument("--skip-plots", action="store_true",
                        help="Skip output figures.")

    return parser
