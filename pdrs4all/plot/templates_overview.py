"""Simple way to inspect the templates.ecsv file what was created.

Should support saving the plot, so we can store these as a preview,
together with the data products .

"""

from astropy.table import Table
from argparse import ArgumentParser
from matplotlib import pyplot as plt

COLORS = ["b", "orange", "g", "r", "m"]


def flux_and_snr(ax_f, ax_u, ax_snr, w, f, u, **plot_kwargs):
    ax_f.plot(w, f, **plot_kwargs)
    ax_f.set_ylabel("flux (MJy sr-1)")
    ax_u.plot(w, u, **plot_kwargs)
    ax_u.set_ylabel("unc (MJy sr-1)")
    ax_snr.plot(w, f / u, **plot_kwargs)
    ax_snr.set_ylabel("S/N")

    for ax in (ax_f, ax_u, ax_snr):
        ax.set_xlabel("wavelength (μm)")


if __name__ == "__main__":
    ap = ArgumentParser()
    ap.add_argument("templates_ecsv")
    ap.add_argument("--interactive", action="store_true")
    ap.add_argument("-o", help="output file for plot", default="templates.pdf")
    ap.add_argument(
        "--compare",
        help="""
    Another templates.ecsv file to compare with (not implemented
    yet)""",
    )
    ap.add_argument(
        "--segments_ecsv",
        nargs="+",
        help="""
    Files containing extractions for individual cubes. These segments
    will be plotted to inspect how well the stitching works.""",
    )
    args = ap.parse_args()

    t = Table.read(args.templates_ecsv)
    suffixes = [s.split("_", 1)[1] for s in t.colnames if "flux" in s]

    print("Plotting templates ", suffixes)

    fig, axs = plt.subplots(3, 1, sharex=True)
    wavelength = t["wavelength"]
    for i, k in enumerate(suffixes):
        flux = t[f"flux_{k}"]
        unc = t[f"unc_{k}"]
        flux_and_snr(
            axs[0], axs[1], axs[2], wavelength, flux, unc, label=k, color=COLORS[i]
        )

    if args.segments_ecsv is not None:
        for t_segment in (Table.read(fn) for fn in args.segments_ecsv):
            wavelength = t_segment["wavelength"]
            for i, k in enumerate(suffixes):
                flux = t_segment[f"flux_{k}"]
                unc = t_segment[f"unc_{k}"]
                flux_and_snr(
                    axs[0],
                    axs[1],
                    axs[2],
                    wavelength,
                    flux,
                    unc,
                    label=k,
                    color=COLORS[i],
                    alpha=0.33,
                )

    axs[0].legend()
    # for ax in axs:
    #     ax.set_yscale('log')

    fig.set_size_inches(8, 8)
    fig.tight_layout()

    if args.interactive:
        plt.show()
    else:
        fig.savefig(args.o)
