import os.path
import logging

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

PLOTFILE_SUFFIX = ".pdf"
CCOUTPUT_SUFFIX = "_cc.txt"
STATSFILE_SUFFIX = "_stats.tab"
EXPECT_OUTFILE_SUFFIXES = (PLOTFILE_SUFFIX, CCOUTPUT_SUFFIX, STATSFILE_SUFFIX)

logger = logging.getLogger(__name__)


def get_output_basename(dirpath, filepath):
    return os.path.join(dirpath, os.path.splitext(os.path.basename(filepath))[0])


def output_result(sourcepath, ccc, outdir):
    outfile_prefix = get_output_basename(outdir, sourcepath)
    logger.info("Output results to '{}'".format(outfile_prefix))

    output_cc(outfile_prefix + CCOUTPUT_SUFFIX, ccc)
    output_stats(outfile_prefix + STATSFILE_SUFFIX, ccc)
    plot_figures(outfile_prefix + PLOTFILE_SUFFIX, ccc, os.path.basename(sourcepath))


def output_cc(outfile, ccc):
    with open(outfile, 'w') as f:
        for cc in ccc.cc:
            print >>f, cc


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
    plt.text(0, min(ccc.cc), 'min(cc): {:.3f}'.format(min(ccc.cc)))

    ccrl_x = int(round(ccc.read_mean_len))
    ccrl_y = ccc.cc[ccrl_x - 1]
    plt.scatter(ccrl_x, ccrl_y, facecolors="none", edgecolors="red")
    plt.annotate(" cc(read length): {:.3f}".format(ccrl_y), (ccrl_x, ccrl_y))
