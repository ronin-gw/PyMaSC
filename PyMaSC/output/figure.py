import logging
import os.path

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
    name = os.path.basename(outfile)

    with PdfPages(outfile) as pp:
        if not ccr.skip_ncc:
            plot_naive_cc(ccr, name)
            pp.savefig()
            plt.close()

        if ccr.calc_masc:
            if plot_naive_cc_just(ccr, name):
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

    plt.plot(xrange(ccr.max_shift + 1), ccr.cc, color="black", linewidth=0.5)
    axes = plt.gca()
    if xlim:
        axes.set_xlim(xlim)
    lower, upper = axes.get_ylim()
    lower, upper = axes.set_ylim((lower, upper * 1.1))
    height = upper - lower

    plt.axhline(ccr.cc_min, linestyle="dashed", linewidth=0.5)
    plt.text(0, ccr.cc_min, 'min(cc) = {:.5f}'.format(ccr.cc_min))

    _annotate_point(
        ccr.read_len - 1, ccr.ccrl, " cc(read length) = {:.5f}".format(ccr.ccrl),
        upper - height/25, 'read length: {:.1f}'.format(ccr.read_len),
        "red", height/50
    )

    if ccr.estimated_library_len:
        _annotate_point(
            ccr.estimated_library_len - 1, ccr.estimated_ccfl, " cc(lib length) = {:.5f}".format(ccr.estimated_ccfl),
            upper - height/10, 'estimated lib len: {}'.format(ccr.estimated_library_len),
            "blue", height/50
        )
        _annotate_params(ccr.estimated_nsc, ccr.estimated_rsc)
    elif ccr.expected_library_len:
        _annotate_point(
            ccr.expected_library_len - 1, ccr.ccfl, " cc(lib length) = {:.5f}".format(ccr.ccfl),
            upper - height/10, 'estimated lib len: {}'.format(ccr.expected_library_len),
            "green", height/50
        )
        _annotate_params(ccr.nsc, ccr.rsc)


def plot_naive_cc_just(ccr, name=None):
    if not ccr.skip_ncc and ccr.estimated_library_len * 2 < ccr.max_shift + 1:
        plot_naive_cc(ccr, name, (0, ccr.estimated_library_len * 2))
        return True
    return False


def plot_masc(ccr, name=None):
    title = "MSCC and Library Length Estimation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Mappability Sensitive Cross-Correlation")

    plt.plot(xrange(ccr.max_shift + 1), ccr.masc, color="black", linewidth=0.5, label="MSCC")
    plt.plot(xrange(ccr.max_shift + 1), moving_avr_filter(ccr.masc, ccr.filter_len), alpha=0.8, label="Smoothed")

    axes = plt.gca()
    lower, upper = axes.get_ylim()
    lower, upper = axes.set_ylim((lower, upper * 1.1))
    height = upper - lower

    masc_ll = ccr.masc[ccr.estimated_library_len - 1]

    _annotate_point(
        ccr.estimated_library_len - 1, masc_ll, " cc(estimated lib len) = {:.5f}".format(masc_ll),
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
        _plot_ncc_vs_masc(
            pp, "Naive CC vs MSCC", ccr.max_shift, ccr.read_len,
            ccr.cc, ccr.cc_min, ccr.masc, ccr.masc_min,
            ccr.nsc, ccr.rsc, ccr.estimated_library_len, ccr.expected_library_len
        )

    for ref in sorted(ccr.ref2genomelen):
        try:
            _plot_ncc_vs_masc(
                pp, title.format(ref), ccr.max_shift, ccr.read_len,
                ccr.ref2cc.get(ref), ccr.ref2cc_min.get(ref), ccr.ref2masc.get(ref), ccr.ref2masc_min.get(ref),
                ccr.ref2est_nsc.get(ref), ccr.ref2est_rsc.get(ref), ccr.estimated_library_len, ccr.expected_library_len
            )
        except AssertionError:
            logger.debug("Skip plot for {}, valid reads unable.".format(ref))


def _plot_ncc_vs_masc(pp, title, max_shift, read_len,
                      cc=None, cc_min=None, masc=None, masc_min=None,
                      nsc=None, rsc=None, estimated_library_len=None, expected_library_len=None):
    assert (cc is not None) or (masc is not None)

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Relative Cross-Correlation from Infimum")

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
