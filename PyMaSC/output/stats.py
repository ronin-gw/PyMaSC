from __future__ import print_function

import logging

from PyMaSC.utils.output import catch_IOError

logger = logging.getLogger(__name__)


@catch_IOError(logger, "cc table")
def output_cc(outfile, ccr):
    with open(outfile, 'w') as f:
        print("shift\tcc\tmscc", file=f)
        if ccr.calc_masc:
            if ccr.skip_ncc:
                for i, cc in enumerate(ccr.masc):
                    print("{}\tnan\t{}".format(i, ccr.masc[i]), file=f)
            else:
                for i, cc in enumerate(ccr.cc):
                    print("{}\t{}\t{}".format(i, cc, ccr.masc[i]), file=f)
        else:
            for i, cc in enumerate(ccr.cc):
                print("{}\t{}\tnan".format(i, cc), file=f)


@catch_IOError(logger, "cc stats")
def output_stats(outfile, ccr):
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
                ("mappable genome length", "mappable_len"),
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

        if ccr.calc_masc:
            print("\nshift\tforward reads in mappable region\treverse reads in mappable region", file=f)
            for i, (forward, reverse) in enumerate(zip(ccr.mappable_forward_sum, ccr.mappable_reverse_sum)):
                print("{}\t{}\t{}".format(i, forward, reverse), file=f)
