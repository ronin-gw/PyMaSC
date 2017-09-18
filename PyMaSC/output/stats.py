import logging

from PyMaSC.utils.output import catch_IOError

logger = logging.getLogger(__name__)


@catch_IOError(logger, "cc table")
def output_cc(outfile, ccr):
    with open(outfile, 'w') as f:
        # if ccr.calc_masc:
        #     print >>f, "shift\tcc\tmscc"
        #     for i, cc in enumerate(ccr.cc):
        #         print >>f, "{}\t{}".format(i + 1, cc, ccr.masc[i])
        # else:
            print >>f, "shift\tcc"
            for i, cc in enumerate(ccr.cc):
                print >>f, "{}\t{}".format(i + 1, cc)


@catch_IOError(logger, "cc stats")
def output_stats(outfile, ccr):
    with open(outfile, 'w') as f:
        for row, attr in (
            ("genome length", "genomelen"),
            ("forward reads", "forward_sum"),
            ("reverse reads", "reverse_sum"),
            ("total forward read length", "forward_read_len_sum"),
            ("total reverse read length", "reverse_read_len_sum"),
            ("read length", "read_len"),
            ("minimum cc", "cc_min"),
            ("cc at read length", "ccrl")
        ):
            print >>f, "{}\t{}".format(row, getattr(ccr, attr))

        if ccr.expected_library_len:
            print >>f, "{}\t{}".format("expected library length", ccr.expected_library_len)

        if ccr.calc_masc:
            for row, attr in (
                ("mappable genome length", "alignable_len"),
                ("forward reads in mappable region", "mappable_forward_sum"),
                ("reverse reads in mappable region", "mappable_reverse_sum"),
                ("total forward read length in mappable region", "mappable_forward_len_sum"),
                ("total reverse read length in mappable region", "mappable_reverse_len_sum"),
                ("estimated library length", "estimated_library_len"),
                ("cc at library length", "estimated_ccfl"),
                ("NSC", "estimated_nsc"),
                ("RSC", "estimated_rsc")
            ):
                print >>f, "{}\t{}".format(row, getattr(ccr, attr))

        elif ccr.expected_library_len:
            for row, attr in (
                ("cc at library length", "ccfl"),
                ("NSC", "nsc"),
                ("RSC", "rsc")
            ):
                print >>f, "{}\t{}".format(row, getattr(ccr, attr))
