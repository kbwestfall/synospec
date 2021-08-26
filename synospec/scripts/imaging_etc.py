#!/usr/bin/env python3

import os
import time
import warnings
import argparse

from IPython import embed

import numpy

from matplotlib import pyplot, ticker

from astropy import units

from synospec.etc import source, efficiency, telescopes, spectrum, extract, aperture, detector
from synospec.etc.observe import Observation


def parse_args(options=None):

    parser = argparse.ArgumentParser(description='FOBOS Exposure Time Calculator (v0.2)',
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

    parser.add_argument('--spot_fwhm', default=5.8, type=float,
                        help='FHWM of the monochromatic spot size on the detector in pixels.')

    res_group = parser.add_mutually_exclusive_group()
    res_group.add_argument('--spec_res_indx', default=None,
                           help='Extension or column number with the flux.')
    res_group.add_argument('--spec_res_value', default=None,
                           help='Single value for the spectral resolution (R = lambda/dlambda) '
                                'for the full spectrum.')
    parser.add_argument('--spec_table', default=None,
                        help='Extension in the fits file with the binary table data.')

#    parser.add_argument('-w', '--wavelengths', default=[3100,10000,4e-5], nargs=3, type=float,
#                        help='Wavelength grid: start wave, approx end wave, logarithmic step')
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
    parser.add_argument('-l', '--emline', default=None, type=str,
                        help='File with emission lines to add to the spectrum.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--sersic', default=None, nargs=4, type=float,
                       help='Use a Sersic profile to describe the object surface-brightness '
                            'distribution; order must be effective radius, Sersic index, '
                            'ellipticity (1-b/a), position angle (deg).')
    group.add_argument('-u', '--uniform', default=False, action='store_true',
                       help='Instead of a point source or Sersic profile, assume the surface '
                            'brightness distribution is uniform over the fiber face.  If set, '
                            'the provided magnitude is assumed to be a surface brightness.  See '
                            'the MAG option.')

    parser.add_argument('-t', '--time', default=3600., type=float, help='Exposure time (s)')
    parser.add_argument('-f', '--fwhm', default=0.65, type=float,
                        help='On-sky PSF FWHM (arcsec)')
    parser.add_argument('-a', '--airmass', default=1.0, type=float, help='Airmass')
    parser.add_argument('-i', '--ipython', default=False, action='store_true',
                        help='After completing the setup, embed in an IPython session.')
    parser.add_argument('-p', '--plot', default=False, action='store_true',
                        help='Provide a plot of the components of the calculation.')
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


def get_wavelength_vector(start, end, logstep):
    """
    Get the wavelength vector
    """
    nwave = int((numpy.log10(end)-numpy.log10(start))/logstep + 1)
    return numpy.power(10., numpy.arange(nwave)*logstep + numpy.log10(start))


def read_spectrum(spec_file, spec_wave, spec_wave_units, spec_flux, spec_flux_units, spec_res_indx,
                  spec_res_value, spec_table, wave, resolution):
    """
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

    if spec is None:
        # Then assume it's an ascii file
        try:
            rescol = None if spec_res_indx is None else int(spec_res_indx)
            spec = spectrum.Spectrum.from_ascii(spec_file, wavecol=int(spec_wave),
                                                waveunits=spec_wave_units, fluxcol=int(spec_flux),
                                                fluxunits=spec_flux_units, rescol=rescol,
                                                resolution=spec_res_value,
                                                use_sampling_assessments=True)
        except:
            spec = None

    if spec is None:
        raise IOError('Could not read provided file.')

    # Force the spectrum to be regularly sampled on a logarithmic grid
    if not spec.regular:
        print('Imposing uniform sampling')
        spec = spec.resample(log=True)

    # If the resolution is available match it to the resultion set for
    # the FOBOS spectrum
    if spec.sres is not None:
        if numpy.all(spec.sres > resolution):
            print('Convolving to FOBOS spectral resolution.')
            spec = spec.match_resolution(resolution)
        else:
            warnings.warn('Spectral resolution of the spectrum is below what FOBOS will provide.')

    # Down-sample to the provided wavelength vector
    print('Resampling to requested wavelength range.')
    return spec.resample(wave=wave, log=True)


def get_spectrum(wave, mag, mag_band='g', mag_system='AB', spec_file=None, spec_wave=None,
                 spec_wave_units=None, spec_flux=None, spec_flux_units=None, spec_res_indx=None,
                 spec_res_value=None, spec_table=None, emline_db=None, redshift=0.0,
                 resolution=3500):
    """
    """
    spec = spectrum.ABReferenceSpectrum(wave, resolution=resolution, log=True) \
                if spec_file is None \
                else read_spectrum(spec_file, spec_wave, spec_wave_units, spec_flux,
                                   spec_flux_units, spec_res_indx, spec_res_value, spec_table,
                                   wave, resolution)
    broadband = efficiency.FilterResponse(band=mag_band)
    spec.rescale_magnitude(mag, band=broadband, system=mag_system)
    if emline_db is None:
        return spec
    spec = spectrum.EmissionLineSpectrum(wave, emline_db['flux'], emline_db['restwave'],
                                         emline_db['fwhm'], units=emline_db['fwhmu'],
                                         redshift=redshift, resolution=resolution, log=True,
                                         continuum=spec.flux)
    warnings.warn('Including emission lines, spectrum broadband magnitude changed '
                  'from {0} to {1}.'.format(mag, spec.magnitude(band=broadband)))
    return spec


def get_sky_spectrum(mag=None, mag_band='g', mag_system='AB'):
    sky_spec = spectrum.MaunakeaSkySpectrum()
    if mag is None:
        return sky_spec

    broadband = efficiency.FilterResponse(band=mag_band)
    sky_spec.rescale_magnitude(mag, band=broadband, system=mag_system)
    return sky_spec


def get_source_distribution(fwhm, uniform, sersic, size=None, sampling=None):
    if uniform:
        return None

    # Build the source surface brightness distribution with unity
    # integral; intrinsic is set to 1 for a point source
    intrinsic = 1. if sersic is None \
                    else source.OnSkySersic(1.0, sersic[0], sersic[1], ellipticity=sersic[2],
                                            position_angle=sersic[3], unity_integral=True)

    # Set the size of the map used to render the source
    _size = fwhm*5 if sersic is None else max(fwhm*5, sersic[0]*5)
    if size is not None:
        if size < _size:
            warnings.warn(f'Map size ({size}) smaller than nominal calculation ({_size}).')
        _size = size

    _sampling = fwhm/10 if sersic is None else min(fwhm/10, sersic[0]/10/sersic[1])
    if sampling is not None:
        if sampling > _sampling:
            warnings.warn(f'Sampling ({sampling}) larger than nominal calculation ({_sampling}).')
        _sampling = sampling

    # Check the rendering of the Sersic profile
    if sersic is not None:
        intrinsic.make_map(sampling=_sampling, size=_size)
        r, theta = intrinsic.semi.polar(intrinsic.X, intrinsic.Y)
        flux_ratio = 2 * numpy.sum(intrinsic.data[r < intrinsic.r_eff.value]) \
                        * numpy.square(_sampling) / intrinsic.get_integral()
        if numpy.absolute(numpy.log10(flux_ratio)) > numpy.log10(1.05):
            warnings.warn('Difference in expected vs. map-rendered flux is larger than '
                          '5%: {0}%.'.format((flux_ratio-1)*100))

    # Construct the on-sky source distribution
    return source.OnSkySource(fwhm, intrinsic, sampling=_sampling, size=_size)


def main(args):

    if args.sky_err < 0 or args.sky_err > 1:
        raise ValueError('--sky_err option must provide a value between 0 and 1.')

    t = time.perf_counter()

    broadband = efficiency.FilterResponse(band=args.mag_band)
    dw = numpy.mean(numpy.diff(broadband.wave))
    wave = broadband.wave[0]-10*dw

    embed()
    exit()

    # Constants:
    resolution = 3500.      # lambda/dlambda
    fiber_diameter = 0.8    # Arcsec
    rn = 2.                             # Detector readnoise (e-)
    dark = 0.0                          # Detector dark-current (e-/s)

    # Temporary numbers that assume a given spectrograph PSF and LSF.
    # Assume 3 pixels per spectral and spatial FWHM.
    spatial_fwhm = args.spot_fwhm
    spectral_fwhm = args.spot_fwhm

    # Get source spectrum in 1e-17 erg/s/cm^2/angstrom. Currently, the
    # source spectrum is assumed to be
    #   - normalized by the total integral of the source flux 
    #   - independent of position within the source
    dw = 1/spectral_fwhm/resolution/numpy.log(10)
    wavelengths = [3100,10000,dw]
    wave = get_wavelength_vector(wavelengths[0], wavelengths[1], wavelengths[2])
    emline_db = None if args.emline is None else read_emission_line_database(args.emline)
    spec = get_spectrum(wave, args.mag, mag_band=args.mag_band, mag_system=args.mag_system,
                        spec_file=args.spec_file, spec_wave=args.spec_wave,
                        spec_wave_units=args.spec_wave_units, spec_flux=args.spec_flux,
                        spec_flux_units=args.spec_flux_units, spec_res_indx=args.spec_res_indx,
                        spec_res_value=args.spec_res_value, spec_table=args.spec_table,
                        emline_db=emline_db, redshift=args.redshift, resolution=resolution)

    # Get the source distribution.  If the source is uniform, onsky is None.
    onsky = get_source_distribution(args.fwhm, args.uniform, args.sersic)

    # Show the rendered source
#    if onsky is not None:
#        pyplot.imshow(onsky.data, origin='lower', interpolation='nearest')
#        pyplot.show()

    # Get the sky spectrum
    sky_spectrum = get_sky_spectrum(args.sky_mag, mag_band=args.sky_mag_band,
                                    mag_system=args.sky_mag_system)
    
    # Overplot the source and sky spectrum
#    ax = spec.plot()
#    ax = sky_spectrum.plot(ax=ax, show=True)

    # Get the atmospheric throughput
    atmospheric_throughput = efficiency.AtmosphericThroughput(airmass=args.airmass)

#    spec.rescale_magnitude(mag, band=broadband, system=mag_system)

    embed()
    exit()

    # Set the telescope. Defines the aperture area and throughput
    # (nominally 3 aluminum reflections for Keck)
    telescope = telescopes.KeckTelescope()

    # Define the observing aperture; fiber diameter is in arcseconds,
    # center is 0,0 to put the fiber on the target center. "resolution"
    # sets the resolution of the fiber rendering; it has nothing to do
    # with spatial or spectral resolution of the instrument
    fiber = aperture.FiberAperture(0, 0, fiber_diameter, resolution=100)

    # Get the spectrograph throughput (circa June 2018; TODO: needs to
    # be updated). Includes fibers + foreoptics + FRD + spectrograph +
    # detector QE (not sure about ADC). Because this is the total
    # throughput, define a generic efficiency object.
    thru_db = numpy.genfromtxt(os.path.join(os.environ['SYNOSPEC_DIR'], 'data/efficiency',
                               'fobos_throughput.db'))
    spectrograph_throughput = efficiency.Efficiency(thru_db[:,1], wave=thru_db[:,0])

    # System efficiency combines the spectrograph and the telescope
    system_throughput = efficiency.SystemThroughput(wave=spec.wave,
                                                    spectrograph=spectrograph_throughput,
                                                    telescope=telescope.throughput)

    # Instantiate the detector; really just a container for the rn and
    # dark current for now. QE is included in fobos_throughput.db file,
    # so I set it to 1 here.
    det = detector.Detector(rn=rn, dark=dark, qe=1.0)

    # Extraction: makes simple assumptions about the detector PSF for
    # each fiber spectrum and mimics a "perfect" extraction, including
    # an assumption of no cross-talk between fibers. Ignore the
    # "spectral extraction".
    extraction = extract.Extraction(det, spatial_fwhm=spatial_fwhm, spatial_width=1.5*spatial_fwhm,
                                    spectral_fwhm=spectral_fwhm, spectral_width=spectral_fwhm)

    # Perform the observation
    obs = Observation(telescope, sky_spectrum, fiber, args.time, det,
                      system_throughput=system_throughput,
                      atmospheric_throughput=atmospheric_throughput, airmass=args.airmass,
                      onsky_source_distribution=onsky, source_spectrum=spec, extraction=extraction,
                      snr_units=args.snr_units)

    # Construct the S/N spectrum
    snr = obs.snr(sky_sub=True, sky_err=args.sky_err)
    if args.ipython:
        embed()

    snr_label = 'S/N per {0}'.format('R element' if args.snr_units == 'resolution'
                                     else args.snr_units)

    if args.plot:
        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        ax = fig.add_axes([0.1, 0.5, 0.8, 0.4])
        ax.set_xlim([wave[0], wave[-1]])
        ax.minorticks_on()
        ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
        ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.set_yscale('log')

        ax = spec.plot(ax=ax, label='Object')
        ax = sky_spectrum.plot(ax=ax, label='Sky')
        ax.legend()
        ax.text(-0.1, 0.5, r'Flux [10$^{-17}$ erg/s/cm$^2$/${\rm \AA}$]', ha='center', va='center',
                transform=ax.transAxes, rotation='vertical')
        
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.4])
        ax.set_xlim([wave[0], wave[-1]])
        ax.minorticks_on()
        ax.tick_params(which='major', length=8, direction='in', top=True, right=True)
        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True)
        ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')

        ax = snr.plot(ax=ax)

        ax.text(0.5, -0.1, r'Wavelength [${\rm \AA}$]', ha='center', va='center',
                transform=ax.transAxes)
        ax.text(-0.1, 0.5, snr_label, ha='center', va='center',
                transform=ax.transAxes, rotation='vertical')

        pyplot.show()

    # Report
    g = efficiency.FilterResponse(band='g')
    r = efficiency.FilterResponse(band='r')
    iband = efficiency.FilterResponse(band='i')
    print('-'*70)
    print('{0:^70}'.format('FOBOS S/N Calculation (v0.2)'))
    print('-'*70)
    print('Compute time: {0} seconds'.format(time.perf_counter() - t))
    print('Object g- and r-band AB magnitude: {0:.1f} {1:.1f}'.format(
                    spec.magnitude(band=g), spec.magnitude(band=r)))
    print('Sky g- and r-band AB surface brightness: {0:.1f} {1:.1f}'.format(
                    sky_spectrum.magnitude(band=g), sky_spectrum.magnitude(band=r)))
    print('Exposure time: {0:.1f} (s)'.format(args.time))
    if not args.uniform:
        print('Aperture Loss: {0:.1f}%'.format((1-obs.aperture_factor)*100))
    print('Extraction Loss: {0:.1f}%'.format((1-obs.extraction.spatial_efficiency)*100))
    print('Median {0}: {1:.1f}'.format(snr_label, numpy.median(snr.flux)))
    print('g-band weighted mean {0} {1:.1f}'.format(snr_label,
                numpy.sum(g(snr.wave)*snr.flux)/numpy.sum(g(snr.wave))))
    print('r-band weighted mean {0} {1:.1f}'.format(snr_label,
                numpy.sum(r(snr.wave)*snr.flux)/numpy.sum(r(snr.wave))))
    print('i-band weighted mean {0} {1:.1f}'.format(snr_label,
                numpy.sum(iband(snr.wave)*snr.flux)/numpy.sum(iband(snr.wave))))

