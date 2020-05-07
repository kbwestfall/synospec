#!/usr/bin/env python3

import os
import time
import warnings
import argparse

from IPython import embed

import numpy
from scipy import interpolate

from matplotlib import pyplot, ticker

from astropy import units

from synospec.etc import source, efficiency, telescopes, spectrum, extract, aperture, detector
from synospec.etc.spectrographs import TMTWFOSBlue, TMTWFOSRed, WFOSGrating
from synospec.etc.observe import Observation

def parse_args(options=None):

    parser = argparse.ArgumentParser(description='WFOS Exposure Time Calculator (v0.1)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--spec_file', default=None, type=str,
                        help='A fits or ascii file with the object spectrum to use')
    parser.add_argument('--spec_wave', default='WAVE',
                        help='Extension or column number with the wavelengths.')
    parser.add_argument('--spec_wave_units', default='angstrom',
                        help='Wavelength units')
    parser.add_argument('--spec_flux', default='FLUX',
                        help='Extension or column number with the flux.')
    parser.add_argument('--spec_flux_units', default=None,
                        help='Input units of the flux density. Must be interpretable by '
                             'astropy.units.Unit.  Code assumes 1e-17 erg / (cm2 s angstrom) '
                             'if units are not provided.')

    res_group = parser.add_mutually_exclusive_group()
    res_group.add_argument('--spec_res_indx', default=None,
                           help='Extension or column number with the flux.')
    res_group.add_argument('--spec_res_value', default=None, type=float,
                           help='Single value for the spectral resolution (R = lambda/dlambda) '
                                'for the full spectrum.')
    parser.add_argument('--spec_table', default=None,
                        help='Extension in the fits file with the binary table data.')

    parser.add_argument('-l', '--emline', default=None, type=str,
                        help='File with emission lines to add to the spectrum.')

    parser.add_argument('-m', '--mag', default=24., type=float,
                        help='Total apparent magnitude of the source')
    parser.add_argument('--mag_band', default='g', type=str,
                        help='Broad-band used for the provided magnitude.  Must be '
                             'u, g, r, i, or z.')
    parser.add_argument('--mag_system', default='AB', type=str,
                        help='Magnitude system.  Must be either AB or Vega.')

    parser.add_argument('--sky_mag', default=None, type=float,
                        help='Surface brightness of the sky in mag/arcsec^2 in the defined '
                             'broadband.  If not provided, default dark-sky spectrum is used.')
    parser.add_argument('--sky_mag_band', default='g', type=str,
                        help='Broad-band used for the provided sky surface brightness.  Must be '
                             'u, g, r, i, or z.')
    parser.add_argument('--sky_mag_system', default='AB', type=str,
                        help='Magnitude system.  Must be either AB or Vega.')

    parser.add_argument('-z', '--redshift', default=0.0, type=float,
                        help='Redshift of the object, z')

    dist_group = parser.add_mutually_exclusive_group()
    dist_group.add_argument('-s', '--sersic', default=None, nargs=4, type=float,
                            help='Use a Sersic profile to describe the object surface-brightness '
                                 'distribution; order must be effective radius, Sersic index, '
                                 'ellipticity (1-b/a), position angle (deg).')
    dist_group.add_argument('-u', '--uniform', default=False, action='store_true',
                            help='Instead of a point source or Sersic profile, assume the surface '
                                 'brightness distribution is uniform over the fiber face.  If '
                                 'set, the provided magnitude is assumed to be a surface '
                                 'brightness.  See the MAG option.')

    parser.add_argument('--refl', default='req', type=str,
                        help='Select the reflectivity curve for TMT.  Must be either \'req\' or '
                             '\'goal\' for the required or goal reflectivity performance.')

    parser.add_argument('--blue_grat', default=TMTWFOSBlue.default_grating,
                        type=str, help='Grating to use in the blue arm.  Options are: {0}'.format(
                        ', '.join([g for g in WFOSGrating.available_gratings.keys() if 'B' in g])))
    bgrat_grp = parser.add_mutually_exclusive_group()
    bgrat_grp.add_argument('--blue_wave', default=None, type=float,
                           help='Central wavelength for the blue arm.  If None, will use the '
                                'peak-efficiency wavelength.')
    bgrat_grp.add_argument('--blue_angle', default=None, type=float,
                           help='Grating angle for blue grating.  If None, will use then angle '
                                'the provides the best efficiency for the on-axis spectrum.')
    parser.add_argument('--blue_binning', default=[1,1], nargs=2, type=int,
                        help='On-chip binning for the blue grating.  Order is spectral then '
                             'spatial.  I.e., to bin 2 pixels spectrally and no binning spatial, '
                             'set --blue_binning 2 1')

    parser.add_argument('--red_grat', default=TMTWFOSRed.default_grating,
                        type=str, help='Grating to use in the red arm.  Options are: {0}'.format(
                        ', '.join([g for g in WFOSGrating.available_gratings.keys() if 'R' in g])))
    rgrat_grp = parser.add_mutually_exclusive_group()
    rgrat_grp.add_argument('--red_wave', default=None, type=float,
                           help='Central wavelength for the red arm.  If None, will use the '
                                'peak-efficiency wavelength.')
    rgrat_grp.add_argument('--red_angle', default=None, type=float,
                           help='Grating angle for red grating.  If None, will use then angle '
                                'the provides the best efficiency for the on-axis spectrum.')
    parser.add_argument('--red_binning', default=[1,1], nargs=2, type=int,
                        help='On-chip binning for the red grating.  Order is spectral then '
                             'spatial.  I.e., to bin 2 pixels spectrally and no binning spatial, '
                             'set --red_binning 2 1')

    parser.add_argument('--slit', default=[0., 0., 0.75, 5.0, 0.0], nargs=5, type=float,
                        help='Slit properties: x field center, y field center, width, length, '
                             'rotation.  The rotation is in degrees, everything else is in '
                             'on-sky arcsec.  The slit width is in the *unrotated* frame, '
                             'meaning the effective slit width for a rotated slit is '
                             'slit_width/cos(rotation).  For the field center, x is along the '
                             'dispersion direction with a valid range of +/- 90 arcsec, and y is '
                             'in the cross-dispersion direction with a valid range of +/- 249 '
                             'arcsec.  Coordinate (0,0) is on axis.')
    parser.add_argument('--extract', default=None, type=float,
                        help='Extraction aperture in arcsec *along the slit* centered on the '
                             'source.  At the detector, the extraction aperture is narrower by '
                             'cos(slit rotation).  If not provided, set to the FWHM of the '
                             'seeing disk; see --fwhm')

    parser.add_argument('-t', '--time', default=3600., type=float, help='Exposure time (s)')
    parser.add_argument('-f', '--fwhm', default=0.65, type=float,
                        help='On-sky PSF FWHM (arcsec)')
    parser.add_argument('-a', '--airmass', default=1.0, type=float, help='Airmass')
    parser.add_argument('-i', '--ipython', default=False, action='store_true',
                        help='After completing the setup, embed in an IPython session.')
    parser.add_argument('-p', '--plot', default=True, action='store_false',
                        help='Do not provide a plot of the components of the calculation.')
    parser.add_argument('--snr_units', type=str, default='pixel',
                        help='The units for the S/N.  Options are pixel, angstrom, resolution.')

    parser.add_argument('--sky_err', type=float, default=0.1,
                        help='The fraction of the Poisson error in the sky incurred when '
                             'subtracting the sky from the observation. Set to 0 for a sky '
                             'subtraction that adds no error to the sky-subtracted spectrum; set '
                             'to 1 for a sky-subtraction error that is the same as the Poisson '
                             'error in the sky spectrum acquired during the observation.')

    return parser.parse_args() if options is None else parser.parse_args(options)


def _emission_line_database_dtype():
    return [('name', '<U20'),
            ('flux', float),
            ('restwave', float),
            ('frame', '<U10'),
            ('fwhm', float),
            ('fwhmu', '<U10')]


def read_emission_line_database(dbfile):
    return numpy.genfromtxt(dbfile, dtype=_emission_line_database_dtype())


def get_wavelength_vector(start, end, step, linear=False):
    """
    Get the wavelength vector
    """
    if linear:
        return numpy.arange(start, end+step, step)
    nwave = int((numpy.log10(end)-numpy.log10(start))/step + 1)
    return numpy.power(10., numpy.arange(nwave)*step + numpy.log10(start))


def set_resolution(wave, resolution):
    if resolution is None:
        return None
    _resolution = numpy.atleast_1d(resolution)
    if _resolution.size == 1:
        _resolution = numpy.full_like(wave, _resolution[0], dtype=float)
    if _resolution.size != wave.size:
        raise ValueError('Provided resolution must be a single number or a vector with the same '
                         'length as the wavelength vector.')
    return _resolution


def read_spectrum(spec_file, spec_wave, spec_wave_units, spec_flux, spec_flux_units, spec_res_indx,
                  spec_res_value, spec_table):
    """
    wave and resolution are for the output spectrum, not what's read
    on input.

    If both are vectors, length of wave and resolution must be the
    same, and the resolution is expected to be at the provided wave.
    """
    if not os.path.isfile(spec_file):
        raise FileNotFoundError('{0} does not exist.'.format(spec_file))

    print('Reading spectrum.')
    # First assume it's a fits file
    try:
        spec = spectrum.Spectrum.from_fits(spec_file, waveext=spec_wave, waveunits=spec_wave_units,
                                           fluxext=spec_flux, fluxunits=spec_flux_units,
                                           resext=spec_res, tblext=spec_table,
                                           resolution=spec_res_value,
                                           use_sampling_assessments=True)
    except:
        spec = None

    if spec is not None:
        return spec

    # Then assume it's an ascii file
    try:
        rescol = None if spec_res_indx is None else int(spec_res_indx)
        spec = spectrum.Spectrum.from_ascii(spec_file, wavecol=int(spec_wave),
                                            waveunits=spec_wave_units, fluxcol=int(spec_flux),
                                            fluxunits=spec_flux_units, rescol=rescol,
                                            resolution=spec_res_value,
                                            use_sampling_assessments=True)
    except:
        raise IOError('Could not read provided file.')

    return spec


def observed_spectrum(spec, wave, resolution, mag=None, mag_band='g', mag_system='AB', 
                      redshift=None, emline_db=None):
    """
    Convert input spectrum to expected (noise-free) observation.

    Operations in order:
        - Redshifts spectrum, if redshift provided
        - Rescales magnitude, if mag provided
        - Matches the spectral resolution, if the input spectrum as a defined resolution vector
        - Matches the sampling to the input wavelength vector
        - Add emission lines, if provided
    """
    if redshift is not None:
        spec.redshift(redshift)
    if mag is not None:
        broadband = efficiency.FilterResponse(band=mag_band)
        # Check that the magnitude can be calculated
        indx = (broadband.wave < spec.wave[0]) & (broadband.eta > 0)
        if numpy.any(indx):
            raise ValueError('Input spectrum (after applying any redshift) does not fully '
                             'overlap selected broad-band filter.  Use a different band or '
                             'different spectrum to rescale to the given magnitude.')
        spec.rescale_magnitude(mag, band=broadband, system=mag_system)

    _resolution = set_resolution(wave, resolution)

    # Force the spectrum to be regularly sampled on a logarithmic grid
    if not spec.regular:
        print('Imposing uniform logarithmic sampling')
        # TODO: I'm not sure that this has to be logarithmic.
        spec = spec.resample(log=True)

    # If the resolution is available match it to the resolution
    # expected for the instrument
    if spec.sres is not None:
        new_sres = interpolate.interp1d(wave, _resolution,
                                        fill_value=(_resolution[0],_resolution[-1]),
                                        bounds_error=False)(spec.wave)
        indx = spec.sres < new_sres
        # TODO: set some tolerance (i.e., some fraction of the spectral range)
        if numpy.any(indx):
            warnings.warn('Spectral resolution of input spectrum is lower than what will be '
                          'provided by the instrument over {0:.1f}% of the spectral range.'.format(
                              numpy.sum(indx)/indx.size))
        print('Convolving to delivered spectral resolution (as best as possible).')
        spec = spec.match_resolution(_resolution, wave=wave)
    else:
        spec.sres = interpolate.interp1d(wave, _resolution,
                                         fill_value=(_resolution[0],_resolution[-1]),
                                         bounds_error=False)(spec.wave)

    # Down-sample to the provided wavelength vector
    print('Resampling to the provided wavelength sampling.')
    spec = spec.resample(wave=wave, log=True)

    if emline_db is None:
        return spec

    warnings.warn('Adding emission lines will change the magnitude of the object.')
    return spectrum.EmissionLineSpectrum(spec.wave, emline_db['flux'], emline_db['restwave'],
                                         emline_db['fwhm'], units=emline_db['fwhmu'],
                                         redshift=redshift, resolution=spec.sres, log=True,
                                         continuum=spec.flux)


def get_sky_spectrum(mag=None, mag_band='g', mag_system='AB'):
    sky_spec = spectrum.MaunakeaSkySpectrum()
    if mag is None:
        return sky_spec

    broadband = efficiency.FilterResponse(band=mag_band)
    sky_spec.rescale_magnitude(mag, band=broadband, system=mag_system)
    return sky_spec


def get_source_distribution(fwhm, uniform, sersic):
    if uniform:
        return None

    # Build the source surface brightness distribution with unity
    # integral; intrinsic is set to 1 for a point source
    intrinsic = 1. if sersic is None \
                    else source.OnSkySersic(1.0, sersic[0], sersic[1], ellipticity=sersic[2],
                                            position_angle=sersic[3], unity_integral=True)

    # Set the size of the map used to render the source
    size = fwhm*5 if sersic is None else max(fwhm*5, sersic[0]*5)
    sampling = fwhm/10 if sersic is None else min(fwhm/10, sersic[0]/10/sersic[1])

    # Check the rendering of the Sersic profile
    if sersic is not None:
        intrinsic.make_map(sampling=sampling, size=size)
        r, theta = intrinsic.semi.polar(intrinsic.X, intrinsic.Y)
        flux_ratio = 2 * numpy.sum(intrinsic.data[r < intrinsic.r_eff.value]) \
                        * numpy.square(sampling) / intrinsic.get_integral()
        if numpy.absolute(numpy.log10(flux_ratio)) > numpy.log10(1.05):
            warnings.warn('Difference in expected vs. map-rendered flux is larger than '
                          '5%: {0}%.'.format((flux_ratio-1)*100))

    # Construct the on-sky source distribution
    return source.OnSkySource(fwhm, intrinsic, sampling=sampling, size=size)


def main(args):

    if args.sky_err < 0 or args.sky_err > 1:
        raise ValueError('--sky_err option must provide a value between 0 and 1.')

    t = time.perf_counter()

    # Extract the slit properties for clarity
    slit_x, slit_y, slit_width, slit_length, slit_rotation = args.slit
    effective_slit_width = slit_width/numpy.cos(numpy.radians(slit_rotation))
    _extract_length = args.fwhm if args.extract is None else args.extract

    # TODO: The slit length is currently not used. Instead, the length
    # of the slit is set to the extraction length. This is mostly just
    # because of the current way the Observation class works.

    # Slit aperture. This representation of the slit is *always*
    # centered at (0,0). Set the aperture based on the extraction
    # length for now.
    slit = aperture.SlitAperture(0., 0., slit_width, _extract_length, rotation=slit_rotation)

    # Get the source distribution.  If the source is uniform, onsky is None.
    onsky = get_source_distribution(args.fwhm, args.uniform, args.sersic)

    # Sky spectrum and atmospheric throughput
    sky_spectrum = spectrum.MaunakeaSkySpectrum()
    atmospheric_throughput = efficiency.AtmosphericThroughput(airmass=args.airmass)

    # Emission lines to add
    emline_db = None if args.emline is None else read_emission_line_database(args.emline)

    # Setup the raw object spectrum
    if args.spec_file is None:
        wavelengths = [3100,10000,1e-5]
        wave = get_wavelength_vector(*wavelengths)
        obj_spectrum = spectrum.ABReferenceSpectrum(wave, log=True)
    else:
        obj_spectrum = read_spectrum(args.spec_file, args.spec_wave, args.spec_wave_units,
                                     args.spec_flux, args.spec_flux_units, args.spec_res_indx,
                                     args.spec_res_value, args.spec_table)

    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    # Setup the instrument arms

    #-------------------------------------------------------------------
    # Blue Arm
    blue_arm = TMTWFOSBlue(reflectivity=args.refl, grating=args.blue_grat, cen_wave=args.blue_wave,
                           grating_angle=args.blue_angle)

    # Pixels per resolution element
    blue_res_pix = blue_arm.resolution_element(slit_width=effective_slit_width, units='pixels') \
                        / args.blue_binning[0]

    # Get the wavelength range for each arm
    blue_wave_lim = blue_arm.wavelength_limits(slit_x, slit_y, add_grating_limits=True)

    # Setup dummy wavelength vectors to get something appropriate for sampling
    max_resolution = blue_arm.resolution(blue_wave_lim[1], x=slit_x,
                                         slit_width=effective_slit_width)

    # Set the wavelength vector to allow for a regular, logarithmic binning
    dw = 1/blue_res_pix/max_resolution/numpy.log(10)
    blue_wave = get_wavelength_vector(blue_wave_lim[0], blue_wave_lim[1], dw)
    resolution = blue_arm.resolution(blue_wave, x=slit_x, slit_width=effective_slit_width)
    blue_spec = observed_spectrum(obj_spectrum, blue_wave, resolution, mag=args.mag,
                                  mag_band=args.mag_band, mag_system=args.mag_system,
                                  redshift=args.redshift, emline_db=emline_db)

    # Resample to linear to better match what's expected for the detector
    blue_ang_per_pix = blue_arm.resolution_element(wave=blue_wave_lim,
                                                   slit_width=effective_slit_width,
                                                   units='angstrom') / blue_res_pix
    blue_wave = get_wavelength_vector(blue_wave_lim[0], blue_wave_lim[1],
                                      numpy.mean(blue_ang_per_pix), linear=True)
    blue_spec = blue_spec.resample(wave=blue_wave, log=False)

    # Spectrograph arm efficiency (this doesn't include the telescope)
    blue_arm_eff = blue_arm.efficiency(blue_spec.wave, x=slit_x, y=slit_y, same_type=False)

    # System efficiency combines the spectrograph and the telescope
    blue_thru = efficiency.SystemThroughput(wave=blue_spec.wave,
                                            spectrograph=blue_arm_eff,
                                            telescope=blue_arm.telescope.throughput)

    # Extraction: makes simple assumptions about the monochromatic
    # image and extracts the flux within the aperture, assuming the
    # flux from both the object and sky is uniformly distributed across
    # all detector pixels (incorrect!).

    # Extraction width in pixels
    spatial_width = slit.length * numpy.cos(numpy.radians(slit.rotation)) / blue_arm.pixelscale \
                        / args.blue_binning[1]
    blue_ext = extract.Extraction(blue_arm.det, spatial_width=spatial_width, profile='uniform')

    # Perform the observation
    blue_obs = Observation(blue_arm.telescope, sky_spectrum, slit, args.time, blue_arm.det,
                           system_throughput=blue_thru,
                           atmospheric_throughput=atmospheric_throughput, airmass=args.airmass,
                           onsky_source_distribution=onsky, source_spectrum=blue_spec,
                           extraction=blue_ext, snr_units=args.snr_units)

    # Construct the S/N spectrum
    blue_snr = blue_obs.snr(sky_sub=True, sky_err=args.sky_err)


    #-------------------------------------------------------------------
    # Red Arm
    red_arm = TMTWFOSRed(reflectivity=args.refl, grating=args.red_grat, cen_wave=args.red_wave,
                         grating_angle=args.red_angle)

    # Pixels per resolution element
    red_res_pix = red_arm.resolution_element(slit_width=effective_slit_width, units='pixels') \
                    / args.red_binning[0]

    # Get the wavelength range for each arm
    red_wave_lim = red_arm.wavelength_limits(slit_x, slit_y, add_grating_limits=True)

    # Setup dummy wavelength vectors to get something appropriate for sampling
    max_resolution = red_arm.resolution(red_wave_lim[1], x=slit_x, slit_width=effective_slit_width)

    # Set the wavelength vector to allow for a regular, logarithmic binning
    dw = 1/red_res_pix/max_resolution/numpy.log(10)
    red_wave = get_wavelength_vector(red_wave_lim[0], red_wave_lim[1], dw)
    resolution = red_arm.resolution(red_wave, x=slit_x, slit_width=effective_slit_width)
    red_spec = observed_spectrum(obj_spectrum, red_wave, resolution, mag=args.mag,
                                 mag_band=args.mag_band, mag_system=args.mag_system,
                                 redshift=args.redshift, emline_db=emline_db)

    # Resample to linear to better match what's expected for the detector
    red_ang_per_pix = red_arm.resolution_element(wave=red_wave_lim,
                                                 slit_width=effective_slit_width,
                                                 units='angstrom') / red_res_pix
    red_wave = get_wavelength_vector(red_wave_lim[0], red_wave_lim[1], numpy.mean(red_ang_per_pix),
                                     linear=True)
    ree_spec = red_spec.resample(wave=red_wave, log=False)

    # Spectrograph arm efficiency (this doesn't include the telescope)
    red_arm_eff = red_arm.efficiency(red_spec.wave, x=slit_x, y=slit_y, same_type=False)

    # System efficiency combines the spectrograph and the telescope
    red_thru = efficiency.SystemThroughput(wave=red_spec.wave, spectrograph=red_arm_eff,
                                           telescope=red_arm.telescope.throughput)

    # Extraction: makes simple assumptions about the monochromatic
    # image and extracts the flux within the aperture, assuming the
    # flux from both the object and sky is uniformly distributed across
    # all detector pixels (incorrect!).

    # Extraction width in pixels
    spatial_width = slit.length * numpy.cos(numpy.radians(slit.rotation)) / red_arm.pixelscale \
                        / args.red_binning[1]
    red_ext = extract.Extraction(red_arm.det, spatial_width=spatial_width, profile='uniform')

    # Perform the observation
    red_obs = Observation(red_arm.telescope, sky_spectrum, slit, args.time, red_arm.det,
                          system_throughput=red_thru,
                          atmospheric_throughput=atmospheric_throughput, airmass=args.airmass,
                          onsky_source_distribution=onsky, source_spectrum=red_spec,
                          extraction=red_ext, snr_units=args.snr_units)

    # Construct the S/N spectrum
    red_snr = red_obs.snr(sky_sub=True, sky_err=args.sky_err)

    # Set the wavelength vector
    dw = 1/(5 if red_res_pix > 5 else red_res_pix)/max_resolution/numpy.log(10)
    red_wave = get_wavelength_vector(red_wave_lim[0], red_wave_lim[1], dw)
    resolution = red_arm.resolution(red_wave, x=slit_x, slit_width=effective_slit_width)
    red_spec = observed_spectrum(obj_spectrum, red_wave, resolution, mag=args.mag,
                                  mag_band=args.mag_band, mag_system=args.mag_system,
                                  redshift=args.redshift, emline_db=emline_db)
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------

    snr_label = 'S/N per {0}'.format('R element' if args.snr_units == 'resolution'
                                     else args.snr_units)

    # Report
    g = efficiency.FilterResponse(band='g')
    r = efficiency.FilterResponse(band='r')
#    iband = efficiency.FilterResponse(band='i')
    print('-'*70)
    print('{0:^70}'.format('WFOS S/N Calculation (v0.1)'))
    print('-'*70)
    print('Compute time: {0} seconds'.format(time.perf_counter() - t))
    print('Object g- and r-band AB magnitude: {0:.1f} {1:.1f}'.format(
                    obj_spectrum.magnitude(band=g), obj_spectrum.magnitude(band=r)))
    print('Sky g- and r-band AB surface brightness: {0:.1f} {1:.1f}'.format(
                    sky_spectrum.magnitude(band=g), sky_spectrum.magnitude(band=r)))
    print('Exposure time: {0:.1f} (s)'.format(args.time))
    if not args.uniform:
        print('Aperture Loss: {0:.1f}%'.format((1-red_obs.aperture_factor)*100))
#    print('Extraction Loss: {0:.1f}%'.format((1-obs.extraction.spatial_efficiency)*100))
#    print('Median {0}: {1:.1f}'.format(snr_label, numpy.median(snr.flux)))
#    print('g-band weighted mean {0} {1:.1f}'.format(snr_label,
#                numpy.sum(g(snr.wave)*snr.flux)/numpy.sum(g(snr.wave))))
#    print('r-band weighted mean {0} {1:.1f}'.format(snr_label,
#                numpy.sum(r(snr.wave)*snr.flux)/numpy.sum(r(snr.wave))))
#    print('i-band weighted mean {0} {1:.1f}'.format(snr_label,
#                numpy.sum(iband(snr.wave)*snr.flux)/numpy.sum(iband(snr.wave))))
 
    if args.plot:
        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        ax = fig.add_axes([0.1, 0.5, 0.8, 0.4])
        ax.set_xlim([obj_spectrum.wave[0], obj_spectrum.wave[-1]])
        ax.minorticks_on()
        ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
        ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.set_yscale('log')

        ax = obj_spectrum.plot(ax=ax, label='Object', color='k', lw=1)
        ax = sky_spectrum.plot(ax=ax, label='Sky', color='0.5', lw=0.5)
        ax.legend()
        ax.text(-0.1, 0.5, r'Flux [10$^{-17}$ erg/s/cm$^2$/${\rm \AA}$]', ha='center', va='center',
                transform=ax.transAxes, rotation='vertical')
        
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.4])
        ax.set_xlim([obj_spectrum.wave[0], obj_spectrum.wave[-1]])
        ax.minorticks_on()
        ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
        ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')

        ax = blue_snr.plot(ax=ax, color='C0', label='Blue Arm')
        ax = red_snr.plot(ax=ax, color='C3', label='Red Arm')

        ax.text(0.5, -0.1, r'Wavelength [${\rm \AA}$]', ha='center', va='center',
                transform=ax.transAxes)
        ax.text(-0.1, 0.5, snr_label, ha='center', va='center',
                transform=ax.transAxes, rotation='vertical')

        pyplot.show()

    if args.ipython:
        embed()
