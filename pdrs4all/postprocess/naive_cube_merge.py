from specutils import Spectrum
from astropy.wcs import WCS
from argparse import ArgumentParser
from pdrs4all.postprocess.spectral_segments import merge_nd, merge_nd_memfriendly
from pdrs4all.postprocess.custom_io import write_cube_s1d_wavetab_jwst_s3d_format

if __name__ == "__main__":
    ap = ArgumentParser(
        prog="Naive _s3d.fits merger",
        description="""Merge cubes that already have the same spatial WCS.

                        Works for default NIRSpec cubes, but to apply this to MIRI cubes, you
                        have to deliberately build them on the same spatial grid.

                        Input: _s3d.fits files (need same xy shape and celestial WCS!)
                        
                        Output: _s3d.fits file (header is rather minimal though)

                        Some header information (pipeline version,
                        dates) will be copied from the first cube given.
                        It cannot not guaranteed however, that the
                        NIRSpec and MIRI data are from the same pipeline
                        run.

                        """,
    )

    ap.add_argument("input_s3d_fits", nargs="+")
    ap.add_argument("-o", "--output_s3d_fits", default="naive_merged_s3d.fits")
    ap.add_argument("--memory_friendly", action="store_true")
    args = ap.parse_args()

    s3ds = [Spectrum.read(fn) for fn in args.input_s3d_fits]
    s3ds.sort(key=lambda x: x.spectral_axis.value[0])

    # Make sure that all cubes are the same shape. From this point on,
    # it will be assumed that all cubes have the same celestial WCS.
    xy_shapes = [s3d.shape[:2] for s3d in s3ds]
    if any(shape != xy_shapes[0] for shape in xy_shapes):
        print("All cubes should have the same shape in X and Y")
        print("Shapes are:", xy_shapes)

    print(
        f"Opened cubes with shape {xy_shapes}. Starting merge. This may take a lot of memory"
    )
    cwcs = WCS(s3ds[0].meta["header"]).celestial
    if args.memory_friendly:
        s3dm = merge_nd_memfriendly(s3ds)
    else:
        s3dm = merge_nd(s3ds)

    print(f"Merge finished. Writing result to {args.output_s3d_fits}.")
    write_cube_s1d_wavetab_jwst_s3d_format(args.output_s3d_fits, s3dm, cwcs)
