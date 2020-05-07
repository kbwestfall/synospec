#!/usr/bin/env python3

import os
import time
import warnings
import argparse

from IPython import embed

import numpy

from matplotlib import pyplot, ticker, colors

from synospec.etc import source, efficiency, telescopes, spectrum, extract, aperture, detector
from synospec.etc.observe import Observation


def parse_args(options=None):

    parser = argparse.ArgumentParser()

    parser.add_argument('-w', '--wavelengths', default=[3100,10000,4e-5], nargs=3, type=float,
                        help='Wavelength grid: start wave, approx end wave, logarithmic step.')
    parser.add_argument('-m', '--mag', default=[21.,25.,0.5], nargs=3, type=float,
                        help='Object apparent AB magnitude; assumed spectrum has a constant AB '
                             'otherwise, it is assumed to be the total magnitude of the source.')
    parser.add_argument('-t', '--time', default=[300.,3600.,300.], nargs=3, type=float,
                        help='Exposure time (s)')
    parser.add_argument('-f', '--fwhm', default=0.65, type=float,
                        help='On-sky PSF FWHM (arcsec)')
    parser.add_argument('-a', '--airmass', default=1.0, type=float, help='Airmass')
    parser.add_argument('--snr_units', type=str, default='pixel',
                        help='The units for the S/N.  Options are pixel, angstrom, resolution.')

    return parser.parse_args() if options is None else parser.parse_args(options)


def get_wavelength_vector(start, end, logstep):
    """
    Get the wavelength vector
    """
    nwave = int((numpy.log10(end)-numpy.log10(start))/logstep + 1)
    return numpy.power(10., numpy.arange(nwave)*logstep + numpy.log10(start))


def get_spectrum(wave, mag, emline_db=None, redshift=0.0, resolution=3500):
    """
    """
    spec = spectrum.ABReferenceSpectrum(wave, resolution=resolution, log=True)
    g = efficiency.FilterResponse()
    spec.rescale_magnitude(mag, band=g)
    if emline_db is None:
        return spec
    spec = spectrum.EmissionLineSpectrum(wave, emline_db['flux'], emline_db['restwave'],
                                         emline_db['fwhm'], units=emline_db['fwhmu'],
                                         redshift=redshift, resolution=resolution, log=True,
                                         continuum=spec.flux)
    warnings.warn('Including emission lines, spectrum g-band magnitude changed '
                  'from {0} to {1}.'.format(mag, spec.magnitude(band=g)))
    return spec


def main(args):

    t = time.perf_counter()

    # Constants:
    resolution = 3500.      # lambda/dlambda
    fiber_diameter = 0.8    # Arcsec
    rn = 2.                             # Detector readnoise (e-)
    dark = 0.0                          # Detector dark-current (e-/s)

    # Temporary numbers that assume a given spectrograph PSF and LSF.
    # Assume 3 pixels per spectral and spatial FWHM.
    spatial_fwhm = 3.0
    spectral_fwhm = 3.0

    mags = numpy.arange(args.mag[0], args.mag[1]+args.mag[2], args.mag[2])
    times = numpy.arange(args.time[0], args.time[1]+args.time[2], args.time[2])

    # Get source spectrum in 1e-17 erg/s/cm^2/angstrom. Currently, the
    # source spectrum is assumed to be
    #   - normalized by the total integral of the source flux 
    #   - independent of position within the source
    wave = get_wavelength_vector(args.wavelengths[0], args.wavelengths[1], args.wavelengths[2])
    spec = get_spectrum(wave, mags[0], resolution=resolution)

    # Get the source distribution.  If the source is uniform, onsky is None.
    onsky = None

    # Get the sky spectrum
    sky_spectrum = spectrum.MaunakeaSkySpectrum()

    # Overplot the source and sky spectrum
#    ax = spec.plot()
#    ax = sky_spectrum.plot(ax=ax, show=True)

    # Get the atmospheric throughput
    atmospheric_throughput = efficiency.AtmosphericThroughput(airmass=args.airmass)

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

    snr_label = 'S/N per {0}'.format('R element' if args.snr_units == 'resolution'
                                     else args.snr_units)

    # SNR
    g = efficiency.FilterResponse()
    r = efficiency.FilterResponse(band='r')
    snr_g = numpy.empty((mags.size, times.size), dtype=float)
    snr_r = numpy.empty((mags.size, times.size), dtype=float)
    for i in range(mags.size):
        spec.rescale_magnitude(mags[i], band=g)
        for j in range(times.size):
            print('{0}/{1} ; {2}/{3}'.format(i+1,mags.size,j+1,times.size), end='\r')
            # Perform the observation
            obs = Observation(telescope, sky_spectrum, fiber, times[j], det,
                              system_throughput=system_throughput,
                              atmospheric_throughput=atmospheric_throughput, airmass=args.airmass,
                              onsky_source_distribution=onsky, source_spectrum=spec,
                              extraction=extraction, snr_units=args.snr_units)

            # Construct the S/N spectrum
            snr_spec = obs.snr(sky_sub=True)
            snr_g[i,j] = numpy.sum(g(snr_spec.wave)*snr_spec.flux)/numpy.sum(g(snr_spec.wave))
            snr_r[i,j] = numpy.sum(r(snr_spec.wave)*snr_spec.flux)/numpy.sum(r(snr_spec.wave))
    print('{0}/{1} ; {2}/{3}'.format(i+1,mags.size,j+1,times.size))

    extent = [args.time[0] - args.time[2]/2, args.time[1] + args.time[2]/2,
              args.mag[0] - args.mag[2]/2, args.mag[1] + args.mag[2]/2]

    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(1.5*w,1.5*h))
    ax = fig.add_axes([0.15, 0.2, 0.7, 0.7])
    img = ax.imshow(snr_g, origin='lower', interpolation='nearest', extent=extent, aspect='auto',
                    norm=colors.LogNorm(vmin=snr_g.min(), vmax=snr_g.max()))
    cax = fig.add_axes([0.86, 0.2, 0.02, 0.7])
    pyplot.colorbar(img, cax=cax)
    cax.text(4, 0.5, snr_label, ha='center', va='center', transform=cax.transAxes,
             rotation='vertical')
    ax.text(0.5, -0.08, 'Exposure Time [s]', ha='center', va='center', transform=ax.transAxes)
    ax.text(-0.12, 0.5, r'Surface Brightness [AB mag/arcsec$^2$]', ha='center', va='center',
            transform=ax.transAxes, rotation='vertical')
    ax.text(0.5, 1.03, r'$g$-band S/N', ha='center', va='center', transform=ax.transAxes,
            fontsize=12)
    pyplot.show()

    #embed()

