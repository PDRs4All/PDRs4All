import numpy as np

from astropy.nddata import CCDData
from astropy import units as u
from jwst import datamodels
from jwst.datamodels import IFUCubeModel
from jwst.assign_wcs.pointing import create_fitswcs

def write_cube_s1d_wavetab_jwst_s3d_format(fits_fn, s3d, celestial_wcs):
    """Same as write_cube_wavetab_jwst_s3d_format with Spectrum cube as input.

    Only works for MJy / sr flux and micron wavelengths at the moment.

    """
    # needs (w, y, x), but old spectrum1d has (x, y, w), new one
    # typicaly has (w, y, x). Here I assume that spectral axis 0 always
    # means (w,y,x) and that spectral axis 2 always means (x, y, w)
    flux_array = s3d.flux.value
    unc_array = s3d.uncertainty.array
    if s3d.spectral_axis_index != 0:
        flux_array = np.swapaxes(flux_array, 0, -1)
        unc_array = np.moveaxis(unc_array, 0, -1)

    if "header" in s3d.meta:
        naxis321 = [s3d.meta.header[f"NAXIS{i}"] for i in (3, 2, 1)]
        if s3d.shape != naxis321:
            raise ValueError(
                f"Something is wrong, shape {flux_array.shape}"
                "does not match NAXIS3, NAXIS2, NAXIS1 f{naxis321}"
            )

    write_cube_wavetab_jwst_s3d_format(
        fits_fn,
        flux_array,
        unc_array,
        wav_array=s3d.spectral_axis.value,
        celestial_wcs=celestial_wcs,
    )


def write_cube_wavetab_jwst_s3d_format(
    fits_fn, flux_array, unc_array, wav_array, celestial_wcs
):
    """Write cube in same format as jwst pipeline.

    The resulting fits file will be loadable by
    specutils.Spectrum.read with format="JWST s3d"

    fits_fn: str

    flux_array: array of shape (num_wav, num_y, num_x)
       In units MJy / sr (hardcoded)

    unc_array: array of shape (num_wav, num_y, num_x)
       In units MJy / sr (hardcoded)

    wav_array: array of shape (num_wav,)
       In units micron (hardcoded)

    celestial_wcs: 2D wcs (astropy.wcs.WCS class)

    """
    # first set up all the wcsinfo. This info will be used by a utility
    # function from the jwst package to set up a GWCS object. Once this
    # object and some metadata have set as members of the IFUCubeModel,
    # the latter will be written out in a format that specutils
    # understands.
    # first, the crucial part for the wavelength table

    ifucube_model = IFUCubeModel(
        data=flux_array,
        err=unc_array,
        wavetable=np.array(
            [(wav_array[None].T,)], dtype=[("wavelength", "<f4", (len(wav_array), 1))]
        ),
    )
    ifucube_model.meta.wcsinfo.ctype3 = "WAVE-TAB"
    ifucube_model.meta.wcsinfo.ps3_0 = "WCS-TABLE"
    ifucube_model.meta.wcsinfo.ps3_1 = "wavelength"
    ifucube_model.meta.wcsinfo.crval3 = 1.0
    ifucube_model.meta.wcsinfo.crpix3 = 1.0
    ifucube_model.meta.wcsinfo.cdelt3 = None
    ifucube_model.meta.wcsinfo.cunit3 = "um"
    ifucube_model.meta.wcsinfo.pc3_1 = 0.0
    ifucube_model.meta.wcsinfo.pc1_3 = 0.0
    ifucube_model.meta.wcsinfo.pc3_2 = 0.0
    ifucube_model.meta.wcsinfo.pc2_3 = 0.0
    ifucube_model.meta.wcsinfo.pc3_3 = 1.0
    ifucube_model.meta.wcsinfo.wcsaxes = 3
    ifucube_model.wavedim = "(1,{:d})".format(len(wav_array))

    # then, the celestial wcs info
    ifucube_model.meta.wcsinfo.crval1 = celestial_wcs.wcs.crval[0]
    ifucube_model.meta.wcsinfo.crval2 = celestial_wcs.wcs.crval[1]
    ifucube_model.meta.wcsinfo.crpix1 = celestial_wcs.wcs.crpix[0]
    ifucube_model.meta.wcsinfo.crpix2 = celestial_wcs.wcs.crpix[1]
    ifucube_model.meta.wcsinfo.cdelt1 = celestial_wcs.wcs.cdelt[0]
    ifucube_model.meta.wcsinfo.cdelt2 = celestial_wcs.wcs.cdelt[1]
    ifucube_model.meta.wcsinfo.ctype1 = "RA---TAN"
    ifucube_model.meta.wcsinfo.ctype2 = "DEC--TAN"
    ifucube_model.meta.wcsinfo.cunit1 = "deg"
    ifucube_model.meta.wcsinfo.cunit2 = "deg"

    pc = celestial_wcs.wcs.get_pc()
    ifucube_model.meta.wcsinfo.pc1_1 = pc[0, 0]
    ifucube_model.meta.wcsinfo.pc1_2 = pc[0, 1]
    ifucube_model.meta.wcsinfo.pc2_1 = pc[1, 0]
    ifucube_model.meta.wcsinfo.pc2_2 = pc[1, 1]

    ifucube_model.meta.ifu.flux_extension = "SCI"
    ifucube_model.meta.ifu.error_extension = "ERR"
    ifucube_model.meta.ifu.error_type = "ERR"

    ifucube_model.meta.bunit_data = "MJy / sr"
    ifucube_model.meta.bunit_err = "MJy / sr"

    ifucube_model.meta.wcs = create_fitswcs(ifucube_model)
    # seems like the reader also needs this one vvv. In the jwst package
    # code, it uses NAXIS1, NAXIS2, NAXIS3.
    ifucube_model.meta.wcs.bounding_box = (
        (0, flux_array.shape[2] - 1),
        (0, flux_array.shape[1] - 1),
        (0, flux_array.shape[0] - 1),
    )

    ifucube_model.write(fits_fn)


def write_s3d_with_new_crval(output_fn, original_fn, crval):
    """Load cube, edit WCS, and write again.

    Parameters
    ----------

    output_fn : path
        Output file, should end in s3d.fits

    original_fn : path
        Original cube of which the CRVAL will be edited

    crval : pair of float
        Usually new_wcs.wcs.crval

    """
    cube_dm = datamodels.open(original_fn)
    cube_dm.meta.wcsinfo.crval1 = crval[0]
    cube_dm.meta.wcsinfo.crval2 = crval[1]
    cube_dm.write(output_fn)


def write_i2d(fits_fn, array, wcs, unit=u.MJy / u.sr):
    """Write numpy array and to FITS file with WCS

    Parameters
    ----------

    fits_fn : path
        FITS file name for output

    array : 2D array-like
        Flux data

    wcs : WCS
        WCS to be added
    """
    ccd = CCDData(array, wcs=wcs, unit=unit)
    ccd.write(fits_fn, overwrite=True)
