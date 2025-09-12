from astropy.modeling.models import PowerLaw1D
import numpy as np
from astropy import units as u
from astropy.io import fits
from pdrs4all.postprocess import bandpasses
from itertools import product
from tqdm import tqdm

DEFAULT_MIN_THROUGHPUT = 1.0e-3  # 1e-4  # 5e-2
DEFAULT_MIN_COVERAGE = 0.95  # 1.


def compute_colorcor(wave, bandpass, flux_ref, wave_ref, flux_source):
    """
    Compute the color correction K given the bandpass, reference spectrum,
    and source spectrum.  To use this color correction, divide the flux density
    for a band by K.  Such color corrections are needed to compute the correct
    flux density at the reference wavelength for a source with the flux_source
    spectral shape in the photometric convention that provides the flux density
    at a reference wavelength (convention B, see Gordon et al. 2022 for details).
    Parameters
    ----------
    wave : nd float array
       the wavelengths of the bandpass, flux_ref, and flux_source
    bandpass : nd float array
        end-to-end, total throughput bandpass of filter in fractional units
    flux_ref : nd float array
        reference flux density F(lambda) as a function of wave
    wave_ref : float
        reference wavelength
    flux_source : nd float array
        source flux density F(lambda) as a function of wave
    """
    # get the flux densities at the reference waveength
    flux_source_lambda_ref = np.interp(wave_ref, wave, flux_source)
    flux_ref_lambda_ref = np.interp(wave_ref, wave, flux_ref)

    # compute the top and bottom integrals
    inttop = np.trapz(wave * bandpass * flux_source / flux_source_lambda_ref, wave)
    intbot = np.trapz(wave * bandpass * flux_ref / flux_ref_lambda_ref, wave)

    return inttop / intbot


def trapz_uncertainty(y_err, x):
    """
    Uncertainty in np.trapz(y, x) using error propagation with y_err (uncertainty in y)
    """
    dx = np.diff(x)
    res = 0.5 * np.sqrt(np.dot(dx**2, (y_err[1:] ** 2 + y_err[:-1] ** 2)))
    return res


def synthetic_photometry_on_spectrum(
    waves,
    spectrum,
    bandpass,
    spectrum_unc=None,
    min_throughput=DEFAULT_MIN_THROUGHPUT,
    min_coverage=1.0,
):
    """
    Generic function to calculate JWST synthetic images.
    Will calculate uncertainties if spectrum_unc is provided.

    Parameters
    ----------

    waves: array
        wavelengths in micron

    spectrum: array
        flux

    bandpass: 3-tuple
        (reference wavelength, wavelength array, efficiency array)

    spectrum_unc: array
        flux uncertainty

    """
    flux_ref = PowerLaw1D(amplitude=1.0, x_0=1.0, alpha=0.0)(waves)

    # 1. Effective wavelength
    # 2. Wavelength grid
    # 3. Total throughput of instrument over wavelength grid
    lambda_eff, lambdas_filt, tp_filt = bandpass

    lambdas_filt = lambdas_filt.to(u.micron).value

    # Number of wavelength bins where throughput is > threshold
    good_filt = (tp_filt / tp_filt.max()) >= min_throughput
    lambdas_filt = lambdas_filt[good_filt]
    tp_filt = tp_filt[good_filt]

    lambda_min = lambdas_filt.min()
    lambda_max = lambdas_filt.max()

    # Calculate wavelength coverage (we don't want to consider spaxels without full coverage)
    good_unc = np.full(waves.size, fill_value=True, dtype=bool)
    if spectrum_unc is not None:
        good_unc = np.isfinite(spectrum_unc) & (spectrum_unc != 0.0)

    all_spec = (waves >= lambda_min) & (waves <= lambda_max)
    good_spec = all_spec & good_unc
    waves_i = waves[good_spec]
    flux_source_i = spectrum[good_spec]
    flux_ref_i = flux_ref[good_spec]

    if spectrum_unc is not None:
        unc_flux_source_i = spectrum_unc[good_spec]

    # Total number of wavelengths over the bandpass
    n_waves = np.sum(all_spec)
    coverage = np.sum(spectrum[good_spec] != 0)
    good_coverage = (coverage >= n_waves * min_coverage) and (n_waves >= 2)

    if not good_coverage:
        return (None, None, None, None, None)

    # print(coverage, n_waves)
    # print(f"{lambda_min} um to {lambda_max} um")
    # Interpolate throughput on wavelength grid of stitched cube
    tp_interp = np.interp(waves_i, lambdas_filt, tp_filt)

    colour_corr = compute_colorcor(
        waves_i, tp_interp, flux_ref_i, lambda_eff, flux_source_i
    )

    flux_a_numer = np.trapz(waves_i * flux_source_i * tp_interp, waves_i)
    flux_a_denom = np.trapz(waves_i * tp_interp, waves_i)

    flux_convention_a = flux_a_numer / flux_a_denom
    flux_convention_b = np.interp(lambda_eff, waves_i, flux_source_i)

    # Deal with uncertainties if needed
    unc_flux_convention_a = None
    unc_flux_convention_b = None
    if spectrum_unc is not None:
        ### Convention A
        # Calculate uncertainty in numerator and denominator
        unc_numer = trapz_uncertainty(waves_i * unc_flux_source_i * tp_interp, waves_i)
        unc_denom = 0.0
        # Calculate uncertainty in numerator/denominator
        unc_flux_convention_a = np.abs(flux_convention_a) * np.sqrt(
            (unc_numer / flux_a_numer) ** 2 + (unc_denom / flux_a_denom) ** 2
        )

        ### Convention B
        unc_flux_convention_b = np.interp(lambda_eff, waves_i, unc_flux_source_i)

    return (
        flux_convention_a,
        flux_convention_b,
        colour_corr,
        unc_flux_convention_a,
        unc_flux_convention_b,
    )


def make_synthetic_images_from_cube(
    comb_cube,
    waves,
    filter_names=None,
    unc_comb_cube=None,
    min_throughput=DEFAULT_MIN_THROUGHPUT,
    min_coverage=DEFAULT_MIN_COVERAGE,
):
    """
    Returns
    -------

    synth_image_dict: dict with structure d['convention_a'][<filter name>] = 2D array

    cc_dict: dict with structure d[<filter name>] = color correction
    """

    nircam_bandpasses = bandpasses.read_nircam()
    if filter_names is None:
        # Filters used in ERS 1288
        filter_names = [
            "F140M",
            "F182M",
            "F187N",
            "F210M",
            "F212N",
            "F277W",
            "F300M",
            "F335M",
            "F480M",
            "F162M",
            "F164N",
            "F323N",
            "F405N",
            "F470N",
        ]

    print("Making synthetic NIRCam images from NIRSpec cube.")
    print(f"Filters: {filter_names}")

    cc_dict = {}
    synth_image_dict = {}
    synth_image_dict["convention_a"] = {}
    synth_image_dict["convention_b"] = {}
    if unc_comb_cube is not None:
        synth_image_dict["unc_convention_a"] = {}
        synth_image_dict["unc_convention_b"] = {}

    # initialize arrays
    nw, ny, nx = comb_cube.shape
    for filter_key in filter_names:
        cc_dict[filter_key] = np.full((ny, nx), fill_value=np.nan)
        synth_image_dict["convention_a"][filter_key] = np.full(
            (ny, nx), fill_value=np.nan
        )
        synth_image_dict["convention_b"][filter_key] = np.full(
            (ny, nx), fill_value=np.nan
        )
        if unc_comb_cube is not None:
            synth_image_dict["unc_convention_a"][filter_key] = np.full(
                (ny, nx), fill_value=np.nan
            )
            synth_image_dict["unc_convention_b"][filter_key] = np.full(
                (ny, nx), fill_value=np.nan
            )

    # perform photometry for each spaxel
    for ix, iy in tqdm(product(range(nx), range(ny)), total=nx * ny):
        flux_source = comb_cube[:, iy, ix]
        if unc_comb_cube is not None:
            unc_flux_source = unc_comb_cube[:, iy, ix]

        for filter_key in filter_names:
            min_throughput_temp = min_throughput
            min_coverage_temp = min_coverage

            if filter_key in ["F140M", "F277M"]:
                min_throughput_temp = 1e-3
                min_coverage_temp = 0.8

            (
                flux_convention_a,
                flux_convention_b,
                colour_corr,
                unc_flux_convention_a,
                unc_flux_convention_b,
            ) = synthetic_photometry_on_spectrum(
                waves,
                flux_source,
                nircam_bandpasses[filter_key],
                spectrum_unc=unc_flux_source,
                min_throughput=min_throughput_temp,
                min_coverage=min_coverage_temp,
            )

            cc_dict[filter_key][iy, ix] = colour_corr
            synth_image_dict["convention_a"][filter_key][iy, ix] = flux_convention_a
            synth_image_dict["convention_b"][filter_key][iy, ix] = flux_convention_b

            if unc_comb_cube is not None:
                synth_image_dict["unc_convention_a"][filter_key][
                    iy, ix
                ] = unc_flux_convention_a
                synth_image_dict["unc_convention_b"][filter_key][
                    iy, ix
                ] = unc_flux_convention_b

    return synth_image_dict, cc_dict


def write_synth_to_fits(synth_image_dict, w, fname_out_synth, fname_out_synth_unc):
    # New HDUList for synthetic images
    # Every fits file needs a PrimaryHDU. We'll make a blank one
    hdu0 = fits.PrimaryHDU()
    # Start an HDUList
    hdu_list = [hdu0]
    # Same for uncertainties
    hdu0_unc = fits.PrimaryHDU()
    # Start an HDUList
    hdu_list_unc = [hdu0_unc]

    filter_list = list(synth_image_dict["convention_a"].keys())
    hdr0 = w.to_header()

    for filt_tmp in filter_list:
        arr = synth_image_dict["convention_a"][filt_tmp]
        unc_arr = synth_image_dict["unc_convention_a"][filt_tmp]
        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr["EXTNAME"] = (filt_tmp, "Filter")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=arr, header=hdr)
        hdu_list = hdu_list + [hdu_i]

        # One ImageHDU per image
        hdr = hdr0.copy()
        hdr["EXTNAME"] = (filt_tmp, "Filter")
        q = 1 * u.MJy / u.sr
        hdr["UNIT"] = f"{q.unit:FITS}"
        hdu_i = fits.ImageHDU(data=unc_arr, header=hdr)
        hdu_list_unc = hdu_list_unc + [hdu_i]

    hdul = fits.HDUList(hdus=hdu_list)
    hdul.writeto(fname_out_synth, output_verify="fix", overwrite=True)

    hdul_unc = fits.HDUList(hdus=hdu_list_unc)
    hdul_unc.writeto(fname_out_synth_unc, output_verify="fix", overwrite=True)


def synthesize_nircam_images(nirspec_s3d_merged):
    """I removed the path_to_throughputs argument. It should be defined
    globally, probably.

    Parameters
    ----------

    nirspec_s3d_merged : Spectrum
        Merged nirspec cube

    Returns
    -------

    synth_image_dict: dict where d[<filter name>] = 2D numpy array

    cc_dict: dict where d[<filter name>] = color correction (not sure what format)

    """
    # Specutils uses x, y, w. Change it to w, y, x
    comb_cube = np.swapaxes(nirspec_s3d_merged.flux.value, -1, 0)
    unc_comb_cube = np.swapaxes(nirspec_s3d_merged.uncertainty.array, -1, 0)
    synth_image_dict, cc_dict = make_synthetic_images_from_cube(
        comb_cube,
        nirspec_s3d_merged.spectral_axis.value,
        unc_comb_cube=unc_comb_cube,
        min_throughput=DEFAULT_MIN_THROUGHPUT,
        min_coverage=DEFAULT_MIN_COVERAGE,
    )
    return synth_image_dict, cc_dict
