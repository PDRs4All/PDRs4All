"""Correct NIRSpec WCS and calibration based on NIRCam data.

This implements the worflow by Ryan to improve the NIRSpec data. By
comparing wavelength-integrated NIRSpec cubes to reprojected NIRCam
images, the WCS and flux calibration are improved. At the end, the
improved products are sitched. Here are main steps, with some
interdependency.

0. Naive stitch.

   To make the synthetic photometry in the next step work, we will need
   to make an initial stitched cube first, so that the integration can
   happen without problems.

1. Create comparable maps

   a. NIRSpec cubes x NIRCam filter curves -> synthetic photometry (synphot)
   b. NIRCam images x NIRSpec WCS -> NIRCam maps on same grid (ncrpj)

2. WCS correction

   a. Measure centroid of proplyd in synphot -> RA, Dec offset
   b. Apply to cubes and synphot -> wcscorr cubes

3. Calibration factor

   a. ncrpj x synphot -> factors
   b. wcscorr cubes x factors -> calcubes (See Peeters 2024 to see which ones to use)

4. Stitch end results

Ryans has a lot of workarounds (e.g. masking certain pixels in the cube)
and variations (single pointing vs mosaic vs single segment vs stitched)
in his code that make it hard to understand what works with what. I
think it will be better to reimplement some things, make sure that all
the output is in a format that I understand, and copy only the parts
that I need.

"""

from specutils import Spectrum1D
from argparse import ArgumentParser
from pdrs4all.postprocess import spectral_segments, custom_io, synth, wcscorr
from astropy.wcs import WCS
from astropy.io import fits
from pathlib import Path
from jwst import datamodels


def parse_args(argv=None):
    ap = ArgumentParser()
    ap.add_argument(
        "nirspec_cubes",
        nargs=3,
        help="nirspec _s3d.fits of the same spatial dimensions, in wavelength order, i.e. 100, 170, 290",
    )
    ap.add_argument("--output_dir", default="./")
    args = ap.parse_args(argv)
    return args


def main(args):
    output_path = Path(args.output_dir)
    output_path.mkdir(exist_ok=True)

    # names for the various outputs (and sometimes re-read as inputs)
    naive_stitch_fits = output_path / "nirspec_naive_stitch_s3d.fits"
    synth_fits = output_path / "nirspec_synth.fits"
    synth_unc_fits = output_path / "nirspec_synth_unc.fits"

    print("Step 0: load cubes and make naive stitched cube")
    s3ds = [Spectrum1D.read(fn) for fn in args.nirspec_cubes]
    # merge algorithm needs sorted list of Spectrum1D cubes
    s3ds.sort(key=lambda x: x.spectral_axis.value[0])
    nirspec_cwcs = WCS(s3ds[0].meta["header"]).celestial

    s3dm = spectral_segments.merge_nd_memfriendly(s3ds)
    custom_io.write_cube_s1d_wavetab_jwst_s3d_format(
        naive_stitch_fits, s3dm, nirspec_cwcs
    )

    print("Step 1a: nirspec synthetic photometry to nircam fiters")
    synth_image_dict, cc_dict = synth.synthesize_nircam_images(s3dm)
    # write them out to see if its working
    synth.write_synth_to_fits(
        synth_image_dict,
        nirspec_cwcs,
        synth_fits,
        synth_unc_fits,
    )

    print("Step 2a: compute offset")

    # open the synth image fits file we just wrote out, get one image to
    # use for nirspec wcscorr. Alternatively, a single slice or
    # collapsed range of the cube could be used.
    with fits.open(synth_fits) as hdul:
        synth_F212N = hdul["F212N"].data

    new_cwcs = wcscorr.nirspec_wcscorr_using_proplyd(synth_F212N, nirspec_cwcs)

    print("Step 2b: apply offset to cube")

    # preserve the following header info
    pipeline_cube_dm = datamodels.open(args.nirspec_cubes[0])
    meta_keys = (
        "bunit_data",
        "bunit_err",
        "cal_step",
        "calibration_software_revision",
        "calibration_software_version",
        "compress",
        "coordinates",
        "data_processing_software_version",
        "date",
        "hga_move",
        "origin",
        "oss_software_version",
        "photometry",
        "prd_software_version",
        "program",
        "target",
        "telescope",
        "visit",
    )

    # open the naive stitched cube and the individual cubes and apply
    # the same offset to all of them
    for fn in [naive_stitch_fits] + args.nirspec_cubes:
        cube_dm = datamodels.open(naive_stitch_fits)
        cube_dm.meta.wcsinfo.crval1 = new_cwcs.wcs.crval[0]
        cube_dm.meta.wcsinfo.crval2 = new_cwcs.wcs.crval[1]
        for k in meta_keys:
            value = getattr(pipeline_cube_dm.meta, k)
            setattr(cube_dm.meta, k, value)

        out_fn = Path(fn).name.replace('_s3d.fits', '_wcscorr_s3d.fits')
        cube_dm.write(output_path / out_fn)


if __name__ == "__main__":
    main(parse_args())
