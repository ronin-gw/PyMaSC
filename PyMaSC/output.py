import logging

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

logger = logging.getLogger(__name__)


def output_cc(outfile, ccc):
    with open(outfile, 'w') as f:
        print >>f, "shift\tcc"
        for i, cc in enumerate(ccc.cc):
            print >>f, "{}\t{}".format(i + 1, cc)


def output_stats(outfile, ccc):
    pass


def plot_figures(outfile, ccc, name):
    with PdfPages(outfile) as pp:
        plot_naive_cc(ccc, name)
        pp.savefig()


def plot_naive_cc(ccc, name=None):
    title = "Genome Wide Cross-Correlation"
    if name:
        title += " for " + name

    plt.title(title)
    plt.xlabel("Reverse Strand Shift")
    plt.ylabel("Cross-Correlation")

    plt.plot(xrange(1, ccc.max_shift + 1), ccc.cc, color="black")

    plt.axvline(ccc.read_mean_len, color="red", linestyle="dashed", linewidth=0.5)
    plt.text(ccc.read_mean_len, 0, 'mean read length: {:.1f}'.format(ccc.read_mean_len))

    plt.axhline(min(ccc.cc), linestyle="dashed", linewidth=0.5)
    plt.text(0, min(ccc.cc), 'min(cc) = {:.3f}'.format(min(ccc.cc)))

    plt.scatter(ccc.ccrl_x, ccc.ccrl_y, facecolors="none", edgecolors="red")
    plt.annotate(" cc(read length) = {:.3f}".format(ccc.ccrl_y), (ccc.ccrl_x, ccc.ccrl_y))
