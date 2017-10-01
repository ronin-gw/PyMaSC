from __future__ import print_function

import logging

from PyMaSC.utils.output import catch_IOError

logger = logging.getLogger(__name__)


@catch_IOError(logger, "CC table")
def output_cc(outfile, ccr):
    if ccr.skip_ncc:
        return None

    logger.info("Output '{}'".format(outfile))
    with open(outfile, 'w') as f:
        keys = sorted(ccr.ref2cc)
        print('\t'.join(["shift", "whole"] + keys), file=f)
        fmt = '{}\t' * (len(keys) + 1) + '{}'
        for i, cc in enumerate(ccr.cc):
            print(fmt.format(i, cc, *[ccr.ref2cc[k][i] for k in keys]), file=f)


@catch_IOError(logger, "MSCC table")
def output_mscc(outfile, ccr):
    if not ccr.calc_masc:
        return None

    logger.info("Output '{}'".format(outfile))
    with open(outfile, 'w') as f:
        keys = sorted(ccr.ref2masc)
        print('\t'.join(["shift", "whole"] + keys), file=f)
        fmt = '{}\t' * (len(keys) + 1) + '{}'
        for i, cc in enumerate(ccr.masc):
            print(fmt.format(i, cc, *[ccr.ref2masc[k][i] for k in keys]), file=f)


@catch_IOError(logger, "stats")
def output_stats(outfile, ccr):
    logger.info("Output '{}'".format(outfile))

    with open(outfile, 'w') as f:
        for row, attr in (
            ("genome length", "genomelen"),
            ("forward reads", "forward_sum"),
            ("reverse reads", "reverse_sum"),
            ("read length", "read_len"),
            ("minimum cc", "cc_min"),
            ("cc at read length", "ccrl")
        ):
            print("{}\t{}".format(row, getattr(ccr, attr, "nan")), file=f)

        if ccr.expected_library_len:
            print("{}\t{}".format("expected library length", ccr.expected_library_len), file=f)

        if ccr.calc_masc:
            for row, attr in (
                ("estimated library length", "estimated_library_len"),
                ("cc at library length", "estimated_ccfl"),
                ("NSC", "estimated_nsc"),
                ("RSC", "estimated_rsc")
            ):
                print("{}\t{}".format(row, getattr(ccr, attr, "nan")), file=f)

        elif ccr.expected_library_len:
            for row, attr in (
                ("cc at library length", "ccfl"),
                ("NSC", "nsc"),
                ("RSC", "rsc")
            ):
                print("{}\t{}".format(row, getattr(ccr, attr, "nan")), file=f)
