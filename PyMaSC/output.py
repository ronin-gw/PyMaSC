import logging

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

logger = logging.getLogger(__name__)


def output_cc(outfile, ccr):
    with open(outfile, 'w') as f:
        print >>f, "shift\tcc"
        for i, cc in enumerate(ccr.cc):
            print >>f, "{}\t{}".format(i + 1, cc)


def output_stats(outfile, ccr):
    pass


def plot_figures(outfile, ccr, name):
    with PdfPages(outfile) as pp:
        plot_naive_cc(ccr, name)
        pp.savefig()
    plt.clf()


def plot_naive_cc(ccr, name=None):
    title = "Whole Genome Cross-Correlation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Cross-Correlation")

    plt.plot(xrange(1, ccr.max_shift + 1), ccr.cc, color="black")

    plt.axvline(ccr.read_mean_len, color="red", linestyle="dashed", linewidth=0.5)
    plt.text(ccr.read_mean_len, 0, 'mean read length: {:.1f}'.format(ccr.read_mean_len))

    plt.axhline(min(ccr.cc), linestyle="dashed", linewidth=0.5)
    plt.text(0, min(ccr.cc), 'min(cc) = {:.3f}'.format(min(ccr.cc)))

    plt.scatter(ccr.estimated_read_len, ccr.ccrl, facecolors="none", edgecolors="red")
    plt.annotate(" cc(read length) = {:.3f}".format(ccr.ccrl), (ccr.estimated_read_len, ccr.ccrl))
