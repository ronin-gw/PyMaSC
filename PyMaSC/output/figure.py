import logging
import os.path

import numpy as np

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


def _feed_pdf_page(pp):
    pp.savefig()
    plt.close()


@catch_IOError(logger)
def plot_figures(outfile, ccr):
    logger.info("Output '{}'".format(outfile))
    name = os.path.basename(os.path.splitext(outfile)[0])

    with PdfPages(outfile) as pp:
        if not ccr.skip_ncc:
            plot_naive_cc(ccr, name)
            _feed_pdf_page(pp)

        if ccr.calc_masc:
            if plot_naive_cc_just(ccr, name):
                _feed_pdf_page(pp)

            plot_masc(ccr, name)
            _feed_pdf_page(pp)

        plot_ncc_vs_masc(pp, ccr, name)


def _annotate_point(x, color, axis_y, axis_text, point_y=None, point_text=None, yoffset=0):
    plt.axvline(x, color=color, linestyle="dashed", linewidth=0.5)
    plt.annotate(axis_text, (x, axis_y))
    if point_y:
        plt.scatter(x, point_y, facecolors="none", edgecolors=color)
        plt.annotate(point_text, (x, point_y + yoffset))


def _annotate_bottom_right_box(text):
    plt.annotate(
        text,
        textcoords="axes fraction", xy=(1, plt.gca().get_ylim()[0]), xytext=(0.95, 0.05),
        bbox=dict(boxstyle="round", fc="w", alpha=0.9), horizontalalignment="right"
    )


def _annotate_params(nsc=None, rsc=None, est_nsc=None, est_rsc=None, loc="lower right"):
    anno = []
    for stat, label in zip((nsc, rsc, est_nsc, est_rsc),
                           ("NSC", "RSC", "Est NSC", "Est RSC")):
        if stat:
            anno.append("{} = {:.5f}".format(label, stat))

    if anno:
        _annotate_bottom_right_box('\n'.join(anno))


def _set_ylim():
    axes = plt.gca()
    lower, upper = axes.get_ylim()
    if upper > 0:
        lower, upper = axes.set_ylim((lower, upper * 1.1))
    else:
        lower, upper = axes.set_ylim((lower, upper * 0.95))
    height = upper - lower
    return lower, upper, height


def plot_naive_cc(ccr, name=None, xlim=None):
    title = "Cross-Correlation"
    if name:
        title += " for " + name

    stats = ccr.whole

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Cross-Correlation")

    plt.fill_between(xrange(stats.max_shift + 1), ccr.ncc_lower, ccr.ncc_upper,
                     color="lightskyblue", alpha=0.5, linewidth=0)
    plt.plot(xrange(stats.max_shift + 1), stats.cc.cc, color="black", linewidth=0.5)
    axes = plt.gca()
    if xlim:
        axes.set_xlim(xlim)
    lower, upper, height = _set_ylim()

    stats = stats.cc

    plt.axhline(stats.cc_min, linestyle="dashed", linewidth=0.5)
    plt.text(0, stats.cc_min, 'min(cc) = {:.5f}'.format(stats.cc_min))

    _annotate_point(
        stats.read_len - 1, "red",
        upper - height / 25, 'read length: {}'.format(stats.read_len),
        stats.ccrl, " cc(read length) = {:.5f}".format(stats.ccrl), height / 50
    )

    if stats.est_lib_len:
        _annotate_point(
            stats.est_lib_len - 1, "blue",
            upper - height / 10, 'estimated lib len: {}'.format(stats.est_lib_len),
            stats.est_ccfl, " cc(est lib len) = {:.5f}".format(stats.est_ccfl), height / 50
        )
    if stats.library_len:
        _annotate_point(
            stats.library_len - 1, "green",
            upper - height / 6, 'expected lib len: {}'.format(stats.library_len),
            stats.ccfl, " cc(lib length) = {:.5f}".format(stats.ccfl), -height / 25
        )
    _annotate_params(stats.nsc, stats.rsc, stats.est_nsc, stats.est_rsc)


def plot_naive_cc_just(ccr, name=None):
    stats = ccr.whole
    if stats.calc_ncc and stats.cc.est_lib_len * 2 < stats.max_shift + 1:
        plot_naive_cc(ccr, name, (0, stats.cc.est_lib_len * 2))
        return True
    return False


def plot_masc(ccr, name=None):
    title = "MSCC and Library Length Estimation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Mappability Sensitive Cross-Correlation")

    stats = ccr.whole

    plt.fill_between(xrange(stats.max_shift + 1), ccr.mscc_lower, ccr.mscc_upper,
                     color="lightskyblue", alpha=0.5, linewidth=0)
    plt.plot(xrange(stats.max_shift + 1), stats.masc.cc,
             color="black", linewidth=0.5, label="MSCC")
    plt.plot(xrange(stats.max_shift + 1), stats.masc.avr_cc,
             alpha=0.8, label="Smoothed", color="pink")

    lower, upper, height = _set_ylim()
    stats = stats.masc

    masc_est_ll = stats.cc[stats.est_lib_len - 1]
    _annotate_point(
        stats.est_lib_len - 1, "blue",
        upper - height / 2, 'estimated lib len: {}'.format(stats.est_lib_len),
        masc_est_ll, " cc(est lib len) = {:.5f}".format(masc_est_ll), height / 50
    )

    if stats.library_len:
        masc_ll = stats.cc[stats.library_len - 1]
        _annotate_point(
            stats.library_len - 1, "green",
            upper - height / 1.75, 'expected lib len: {}'.format(stats.library_len),
            masc_ll, " cc(lib length) = {:.5f}".format(masc_ll), -height / 25
        )

    plt.legend(loc="best")
    _annotate_bottom_right_box("Mov avr win size = {}".format(stats.mv_avr_filter_len))


def plot_ncc_vs_masc(pp, ccr, name):
    title = "{} Cross-Correlation"
    if name:
        title += " for " + name

    if ccr.calc_masc:
        _plot_ncc_vs_masc(ccr.whole, "Naive CC vs MSCC")
        _feed_pdf_page(pp)

    for ref in sorted(ccr.references):
        try:
            _plot_ncc_vs_masc(ccr.ref2stats[ref], title.format(ref))
            _feed_pdf_page(pp)
        except AssertionError:
            logger.debug("Skip plot for {}, valid reads unable.".format(ref))


def _plot_ncc_vs_masc(stats, title):
    assert (
        (stats.calc_ncc and not np.all(np.isnan(stats.cc.cc))) or
        (stats.calc_masc and not np.all(np.isnan(stats.masc.cc)))
    )

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Relative Cross-Correlation from each minimum")

    if stats.calc_ncc:
        plt.plot(xrange(stats.max_shift + 1), stats.cc.cc - stats.cc.cc_min,
                 color="black", linewidth=0.5, label="Naive CC")
    if stats.calc_masc:
        plt.plot(xrange(stats.max_shift + 1), stats.masc.cc - stats.masc.cc_min,
                 alpha=1 if not stats.calc_ncc else 0.8, linewidth=0.5, label="MSCC")

    lower, upper, height = _set_ylim()

    _annotate_point(
        stats.read_len, "red",
        upper - height/25, "read length: {}".format(stats.read_len)
    )

    if stats.calc_masc:
        _annotate_point(
            stats.masc.est_lib_len, "blue",
            upper - height/10, "estimated lib len: {}".format(stats.masc.est_lib_len)
        )
        plt.legend(loc="best")

    if stats.library_len:
        _annotate_point(
            stats.library_len, "green",
            upper - height/6, "expected lib len: {}".format(stats.library_len)
        )

    # _annotate_params(stats.nsc, stats.rsc, stats.est_nsc, stats.est_rsc, "best")
