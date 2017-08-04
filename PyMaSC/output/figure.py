import logging

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from PyMaSC.utils.calc import moving_avr_filter
from PyMaSC.utils.output import catch_IOError

logger = logging.getLogger(__name__)


@catch_IOError(logger, "figure")
def plot_figures(outfile, ccr, name):
    with PdfPages(outfile) as pp:
        plot_naive_cc(ccr, name)
        pp.savefig()
        plt.close()

        if ccr.calc_masc:
            plot_naive_cc_just(ccr, name)
            pp.savefig()
            plt.close()

            plot_masc(ccr, name)
            pp.savefig()
            plt.close()

        plot_ncc_vs_masc(pp, ccr, name)


def _annotate_point(x, y, text, axis_y, axis_text, color, yoffset=0):
    plt.axvline(x, color=color, linestyle="dashed", linewidth=0.5)
    plt.scatter(x, y, facecolors="none", edgecolors=color)
    plt.annotate(text, (x, y + yoffset))
    plt.annotate(axis_text, (x, axis_y))


def _annotate_params(nsc, rsc, loc="lower right"):
    plt.annotate(
        '\n'.join(["NSC = {:.5f}".format(nsc), "RSC = {:.5f}".format(rsc)]),
        textcoords="axes fraction", xy=(1, plt.gca().get_ylim()[0]), xytext=(0.75, 0.05),
        bbox=dict(boxstyle="round", fc="w", alpha=0.9)
    )


def plot_naive_cc(ccr, name=None, xlim=None):
    title = "Cross-Correlation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Cross-Correlation")

    plt.plot(xrange(1, ccr.max_shift + 1), ccr.cc, color="black", linewidth=0.5)
    axes = plt.gca()
    if xlim:
        axes.set_xlim(xlim)
    lower, upper = axes.get_ylim()
    lower, upper = axes.set_ylim((lower, upper * 1.1))
    height = upper - lower

    plt.axhline(ccr.cc_min, linestyle="dashed", linewidth=0.5)
    plt.text(0, ccr.cc_min, 'min(cc) = {:.5f}'.format(ccr.cc_min))

    _annotate_point(
        ccr.estimated_read_len, ccr.ccrl, " cc(read length) = {:.5f}".format(ccr.ccrl),
        upper - height/25, 'mean read length: {:.1f}'.format(ccr.read_mean_len),
        "red", height/50
    )

    if ccr.estimated_library_len:
        _annotate_point(
            ccr.estimated_library_len, ccr.estimated_ccfl, " cc(lib length) = {:.5f}".format(ccr.estimated_ccfl),
            upper - height/10, 'estimated lib len: {}'.format(ccr.estimated_library_len),
            "blue", height/50
        )
        _annotate_params(ccr.estimated_nsc, ccr.estimated_rsc)
    elif ccr.expected_library_len:
        _annotate_point(
            ccr.expected_library_len, ccr.ccfl, " cc(lib length) = {:.5f}".format(ccr.ccfl),
            upper - height/10, 'estimated lib len: {}'.format(ccr.expected_library_len),
            "green", height/50
        )
        _annotate_params(ccr.nsc, ccr.rsc)


def plot_naive_cc_just(ccr, name=None):
    if ccr.estimated_library_len * 2 < ccr.max_shift + 1:
        plot_naive_cc(ccr, name, (0, ccr.estimated_library_len * 2))


def plot_masc(ccr, name=None):
    title = "MSCC and Library Length Estimation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Mappability Sensitive Cross-Correlation")

    plt.plot(xrange(1, ccr.max_shift + 1), ccr.masc, color="black", linewidth=0.5, label="MSCC")
    plt.plot(xrange(1, ccr.max_shift + 1), moving_avr_filter(ccr.masc, 10), alpha=0.8, label="Smoothed")

    axes = plt.gca()
    lower, upper = axes.get_ylim()
    lower, upper = axes.set_ylim((lower, upper * 1.1))
    height = upper - lower

    masc_max = max(ccr.masc)

    _annotate_point(
        ccr.estimated_library_len, masc_max, " cc(estimated lib len) = {:.5f}".format(masc_max),
        lower + height/25, ' estimated lib len: {}'.format(ccr.estimated_library_len),
        "blue", height/50
    )

    if ccr.expected_library_len:
        plt.axvline(ccr.expected_library_len, color="green", linestyle="dashed", linewidth=0.5)
        plt.annotate('expected lib len: {}'.format(ccr.expected_library_len),
                     (ccr.expected_library_len, upper - height/25))

    plt.legend(loc="best")


def plot_ncc_vs_masc(pp, ccr, name):
    title = "{} Cross-Correlation"
    if name:
        title += " for " + name

    if ccr.calc_masc:
        plt.title("Naive CC vs MSCC")
        plt.xlabel("Reverse Strand Shift")
        plt.ylabel("Relative Cross-Correlation from Infimum")

        plt.plot(xrange(1, ccr.max_shift + 1), ccr.cc - ccr.cc_min, color="black", linewidth=0.5, label="Naive CC")
        plt.plot(xrange(1, ccr.max_shift + 1), ccr.masc - ccr.masc_min, alpha=0.8, linewidth=0.5, label="MSCC")

        axes = plt.gca()
        lower, upper = axes.get_ylim()
        lower, upper = axes.set_ylim((lower, upper * 1.1))
        height = upper - lower

        plt.axvline(ccr.estimated_read_len, color="red", linestyle="dashed", linewidth=0.5)
        plt.annotate('mean read length: {:.1f}'.format(ccr.read_mean_len),
                     (ccr.estimated_read_len, upper - height/25))

        plt.axvline(ccr.estimated_library_len, color="blue", linestyle="dashed", linewidth=0.5)
        plt.annotate('estimated lib len: {}'.format(ccr.estimated_library_len),
                     (ccr.estimated_library_len, upper - height/10))

        if ccr.expected_library_len:
            plt.axvline(ccr.expected_library_len, color="green", linestyle="dashed", linewidth=0.5)
            plt.annotate('expected lib len: {}'.format(ccr.expected_library_len),
                         (ccr.expected_library_len, upper - height/5))

        plt.legend(loc="best")

        pp.savefig()
        plt.close()

    for ref in sorted(ccr.ref2cc):
        cc = ccr.ref2cc[ref]
        if cc is None:
            continue
        plot_masc_data = ccr.calc_masc and ccr.ref2masc.get(ref, None) is not None

        plt.title(title.format(ref))
        plt.xlabel("Reverse Strand Shift")
        plt.ylabel("Relative Cross-Correlation from Infimum")

        cc_min = ccr.ref2cc_min[ref]
        # ccrl = ccr.ref2ccrl[ref]
        plt.plot(xrange(1, ccr.max_shift + 1), cc - cc_min, color="black", linewidth=0.5, label="Naive CC")

        if plot_masc_data:
            masc = ccr.ref2masc[ref]
            masc_min = ccr.ref2masc_min[ref]

            plt.plot(xrange(1, ccr.max_shift + 1), masc - masc_min, alpha=0.8, linewidth=0.5, label="MSCC")

        axes = plt.gca()
        lower, upper = axes.get_ylim()
        lower, upper = axes.set_ylim((lower, upper * 1.1))
        height = upper - lower
        # yoffset = height / 50

        plt.axvline(ccr.estimated_read_len, color="red", linestyle="dashed", linewidth=0.5)
        plt.annotate('mean read length: {:.1f}'.format(ccr.read_mean_len),
                     (ccr.estimated_read_len, upper - height/25))

        if plot_masc_data:
            # mascrl = ccr.ref2mascrl[ref]
            # ccfl = ccr.ref2est_ccfl[ref]
            nsc = ccr.ref2est_nsc[ref]
            rsc = ccr.ref2est_rsc[ref]

            plt.axvline(ccr.estimated_library_len, color="blue", linestyle="dashed", linewidth=0.5)
            plt.annotate('estimated lib len: {}'.format(ccr.estimated_library_len),
                         (ccr.estimated_library_len, upper - height/10))

            plt.legend(loc="best")
            _annotate_params(nsc, rsc, "best")

        elif ccr.expected_library_len:
            nsc = ccr.ref2nsc[ref]
            rsc = ccr.ref2rsc[ref]

            plt.axvline(ccr.expected_library_len, color="green", linestyle="dashed", linewidth=0.5)
            plt.annotate('expected lib len: {}'.format(ccr.expected_library_len),
                         (ccr.expected_library_len, upper - height/10))

        pp.savefig()
        plt.close()
