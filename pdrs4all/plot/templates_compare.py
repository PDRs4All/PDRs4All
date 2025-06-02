from astropy.table import Table
from matplotlib import pyplot as plt
from argparse import ArgumentParser
import numpy as np
from matplotlib import ticker
from pdrs4all.plot.templates_overview import nice_ticks


def main():
    ap = ArgumentParser()
    ap.add_argument("templates1")
    ap.add_argument("templates2")
    ap.add_argument(
        "--keys",
        nargs="+",
        help='Which template to plot. E.g. "Atomic" or "DF3"',
        default=["HII", "Atomic", "DF1", "DF2", "DF3"],
    )
    ap.add_argument("--label1", default="templates1")
    ap.add_argument("--label2", default="templates2")
    args = ap.parse_args()
    # t1 = Table.read('/home/dvandepu/BigData/orion/mirifu/v3.2b_test/templates/default_wcscorr_nostitch.ecsv')
    # t2 = Table.read('/home/dvandepu/BigData/orion/mirifu/vX.X_auto_jwst_1322.pmap_1.17.1/templates/default_wcscorr_nostitch.ecsv')
    t1 = Table.read(args.templates1)
    t2 = Table.read(args.templates2)

    for k in args.keys:
        compare_one_template(t1, t2, args.label1, args.label2, k)

    plt.show()


def compare_one_template(t1, t2, label1, label2, k):
    """
    Parameters
    ----------

    t1: Table
        Table loaded from templates.ecsv

    t2: Table
        Table loaded from templates.ecsv

    k: str
        Template key

    """
    fig, axs = plt.subplots(2, 1, height_ratios=[2, 1], sharex=True)
    for ax in axs:
        nice_ticks(ax)

    w1 = t1["wavelength"]
    w2 = t2["wavelength"]
    f1 = t1[f"flux_{k}"]
    f2 = t2[f"flux_{k}"]

    axs[0].plot(w1, f1, label=label1, drawstyle='steps-mid')
    axs[0].plot(w2, f2, label=label2, drawstyle='steps-mid')

    if len(f2) == len(f1):
        residual = (f2 - f1) / (0.5 * (f1 + f2))
        axs[1].plot(w1, residual * 100, color="k")
        axs[1].set_ylabel("(orange - blue) / (0.5 * (orange + blue)) [percent]", ha="left")
        axs[1].set_yscale("asinh")
        axs[1].tick_params(which="minor", left=True)
    else:
        print("Cannot directly compute residuals because wavelength grids are different")

    axs[0].set_title(k)
    axs[0].legend()

    axs[0].set_xscale("log")
    axs[0].set_yscale("log")
    axs[0].set_xticks(list(range(1, 16)) + [16, 18, 20, 24, 28])
    axs[0].set_xticks([i + 0.5 for i in range(1, 15)], minor=True)
    axs[0].set_xlim(np.amin(w1), np.amax(w2))
    axs[0].xaxis.set_minor_formatter(ticker.NullFormatter())
    axs[0].xaxis.set_major_formatter(ticker.ScalarFormatter())



if __name__ == "__main__":
    main()
