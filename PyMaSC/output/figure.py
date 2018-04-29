import logging
import os.path

import numpy as np

from PyMaSC.utils.calc import moving_avr_filter
from PyMaSC.utils.output import catch_IOError
from PyMaSC.utils.compatible import xrange

logger = logging.getLogger(__name__)

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except:
    logger.error("Failed to import matplotlib.")
    import traceback
    logger.warning("Exception traceback:\n|" +
                   traceback.format_exc().replace('\n', "\n|"))
    raise


@catch_IOError(logger, "figure")
def plot_figures(outfile, ccr):
    logger.info("Output '{}'".format(outfile))
    name = os.path.basename(os.path.splitext(outfile)[0])

    with PdfPages(outfile) as pp:
        if not ccr.skip_ncc:
            plot_naive_cc(ccr.whole, name)
            pp.savefig()
            plt.close()

        if ccr.calc_masc:
            if plot_naive_cc_just(ccr.whole, name):
                pp.savefig()
                plt.close()

            plot_masc(ccr.whole, name)
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


def plot_naive_cc(stats, name=None, xlim=None):
    title = "Cross-Correlation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Cross-Correlation")

    plt.plot(xrange(stats.max_shift + 1), stats.cc, color="black", linewidth=0.5)
    axes = plt.gca()
    if xlim:
        axes.set_xlim(xlim)
    lower, upper = axes.get_ylim()
    lower, upper = axes.set_ylim((lower, upper * 1.1))
    height = upper - lower

    plt.axhline(stats.cc_min, linestyle="dashed", linewidth=0.5)
    plt.text(0, stats.cc_min, 'min(cc) = {:.5f}'.format(stats.cc_min))

    _annotate_point(
        stats.read_len - 1, stats.ccrl, " cc(read length) = {:.5f}".format(stats.ccrl),
        upper - height/25, 'read length: {:.1f}'.format(stats.read_len),
        "red", height/50
    )

    if stats.est_lib_len or stats.library_len:
        if stats.est_lib_len:
            library_len = stats.est_lib_len
            ccfl = stats.est_ccfl
            nsc = stats.est_nsc
            rsc = stats.est_rsc
            color = "blue"
        elif stats.library_len:
            library_len = stats.library_len
            ccfl = stats.ccfl
            nsc = stats.nsc
            rsc = stats.rsc
            color = "green"

        _annotate_point(
            library_len - 1, ccfl, " cc(lib length) = {:.5f}".format(ccfl),
            upper - height/10, 'estimated lib len: {}'.format(library_len),
            color, height/50
        )
        _annotate_params(nsc, rsc)


def plot_naive_cc_just(stats, name=None):
    if stats.calc_ncc and stats.est_lib_len * 2 < stats.max_shift + 1:
        plot_naive_cc(stats, name, (0, stats.est_lib_len * 2))
        return True
    return False


def plot_masc(stats, name=None):
    title = "MSCC and Library Length Estimation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Mappability Sensitive Cross-Correlation")

    plt.plot(xrange(stats.max_shift + 1), stats.masc, color="black", linewidth=0.5, label="MSCC")
    plt.plot(xrange(stats.max_shift + 1), moving_avr_filter(stats.masc, stats.filter_len), alpha=0.8, label="Smoothed")

    axes = plt.gca()
    lower, upper = axes.get_ylim()
    lower, upper = axes.set_ylim((lower, upper * 1.1))
    height = upper - lower

    masc_ll = stats.masc[stats.est_lib_len - 1]

    _annotate_point(
        stats.est_lib_len - 1, masc_ll, " cc(estimated lib len) = {:.5f}".format(masc_ll),
        lower + height/25, ' estimated lib len: {}'.format(stats.est_lib_len),
        "blue", height/50
    )

    if stats.library_len:
        plt.axvline(stats.library_len, color="green", linestyle="dashed", linewidth=0.5)
        plt.annotate('expected lib len: {}'.format(stats.library_len),
                     (stats.library_len, upper - height/25))

    plt.legend(loc="best")


def plot_ncc_vs_masc(pp, ccr, name):
    title = "{} Cross-Correlation"
    if name:
        title += " for " + name

    if ccr.calc_masc:
        _plot_ncc_vs_masc(pp, "Naive CC vs MSCC", ccr.whole)

    for ref in sorted(ccr.references):
        try:
            _plot_ncc_vs_masc(pp, title.format(ref), ccr.ref2stats[ref])
        except AssertionError:
            logger.debug("Skip plot for {}, valid reads unable.".format(ref))


def _plot_ncc_vs_masc(pp, title, stats):
    try:
        max_shift = stats.max_shift
        read_len = stats.read_len
        cc = stats.cc
        cc_min = stats.cc_min
        masc = stats.masc
        masc_min = stats.masc_min
        nsc = stats.nsc
        rsc = stats.rsc
        estimated_library_len = stats.est_lib_len
        expected_library_len = stats.library_len
    except AttributeError:
        return None

    assert (
        (cc is not None and not np.all(np.isnan(cc))) or
        (masc is not None and not np.all(np.isnan(masc)))
    )

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Relative Cross-Correlation from each minimum")

    if cc is not None:
        plt.plot(xrange(max_shift + 1), cc - cc_min, color="black", linewidth=0.5, label="Naive CC")
    if masc is not None:
        plt.plot(xrange(max_shift + 1), masc - masc_min, alpha=1 if cc is None else 0.8, linewidth=0.5, label="MSCC")

    axes = plt.gca()
    lower, upper = axes.get_ylim()
    lower, upper = axes.set_ylim((lower, upper * 1.1))
    height = upper - lower
    # yoffset = height / 50

    plt.axvline(read_len, color="red", linestyle="dashed", linewidth=0.5)
    plt.annotate('read length: {:.1f}'.format(read_len),
                 (read_len, upper - height/25))

    if masc is not None:
        plt.axvline(estimated_library_len, color="blue", linestyle="dashed", linewidth=0.5)
        plt.annotate('estimated lib len: {}'.format(estimated_library_len),
                     (estimated_library_len, upper - height/10))

        plt.legend(loc="best")
        if None not in (nsc, rsc):
            _annotate_params(nsc, rsc, "best")

    elif expected_library_len:
        plt.axvline(expected_library_len, color="green", linestyle="dashed", linewidth=0.5)
        plt.annotate('expected lib len: {}'.format(expected_library_len),
                     (expected_library_len, upper - height/10))

    pp.savefig()
    plt.close()
