import os
import time
import warnings
import argparse

from IPython import embed

import numpy
from scipy import interpolate
from matplotlib import pyplot, ticker

from astropy import units

from ..etc import source, efficiency, telescopes, spectrum, extract, aperture, detector
from ..etc.observe import Observation
from . import scriptbase
from .. import data_file


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


def get_source_distribution(fwhm, beta=None, sersic=None, size=None, sampling=None):

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
    return source.OnSkySource(fwhm, intrinsic, beta=beta, sampling=_sampling, size=_size)


class ImagingETC(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        import argparse

        parser = super().get_parser(description='Imaging Exposure-Time Calculator', width=width)

#        parser.add_argument('--spec_file', default=None, type=str,
#                            help='A fits or ascii file with the object spectrum to use')
#        parser.add_argument('--spec_wave', default='WAVE',
#                            help='Extension or column number with the wavelengths.')
#        parser.add_argument('--spec_wave_units', default='angstrom',
#                            help='Wavelength units')
#        parser.add_argument('--spec_flux', default='FLUX',
#                            help='Extension or column number with the flux.')
#        parser.add_argument('--spec_flux_units', default=None,
#                            help='Input units of the flux density. Must be interpretable by '
#                                 'astropy.units.Unit.  Code assumes 1e-17 erg / (cm2 s angstrom) '
#                                 'if units are not provided.')
#
#        res_group = parser.add_mutually_exclusive_group()
#        res_group.add_argument('--spec_res_indx', default=None,
#                               help='Extension or column number with the flux.')
#        res_group.add_argument('--spec_res_value', default=None,
#                               help='Single value for the spectral resolution (R = lambda/dlambda)'
#                                    ' for the full spectrum.')
#        parser.add_argument('--spec_table', default=None,
#                            help='Extension in the fits file with the binary table data.')

        parser.add_argument('-m', '--mag', default=24., type=float,
                            help='Total apparent magnitude of the source')
        parser.add_argument('--mag_band', default='g', type=str,
                            help='Broad-band used for the provided magnitude.  Must be '
                                 'u, g, r, i, or z.')
        parser.add_argument('--mag_system', default='AB', type=str,
                            help='Magnitude system.  Must be either AB or Vega.')

#        parser.add_argument('--sky_mag', default=None, type=float,
#                            help='Surface brightness of the sky in mag/arcsec^2 in the defined '
#                                 'broadband.  If not provided, default dark-sky spectrum is used.')
#        parser.add_argument('--sky_mag_band', default='g', type=str,
#                            help='Broad-band used for the provided sky surface brightness.  Must '
#                                 'be u, g, r, i, or z.')
#        parser.add_argument('--sky_mag_system', default='AB', type=str,
#                            help='Magnitude system.  Must be either AB or Vega.')

#        parser.add_argument('-z', '--redshift', default=0.0, type=float,
#                            help='Redshift of the object, z')
#        parser.add_argument('-l', '--emline', default=None, type=str,
#                            help='File with emission lines to add to the spectrum.')

        parser.add_argument('-s', '--sersic', default=None, nargs=4, type=float,
                            help='Use a Sersic profile to describe the object surface-brightness '
                                 'distribution; order must be effective radius, Sersic index, '
                                 'ellipticity (1-b/a), position angle (deg).  If None, assume a '
                                 'point source.')

        parser.add_argument('--binning', default=8, type=int,
                            help='Number of binned pixels along one dimension; e.g., '
                                 'setting "-binning 4" means 4x4 binning.')
        parser.add_argument('-t', '--time', default=30., type=float,
                            help='Total exposure time for all combined exposures (s)')
        parser.add_argument('-n', '--nexp', default=2, type=int, help='Number of exposures')
        parser.add_argument('-f', '--fwhm', default=0.65, type=float,
                            help='On-sky PSF FWHM (arcsec)')
        parser.add_argument('--beta', default=None, type=float,
                            help='Beta for Moffat seeing function.  If None, use a Gaussian.')
        parser.add_argument('-a', '--airmass', default=1.0, type=float, help='Airmass')

        return parser

    @staticmethod
    def main(args):

        warnings.warn('THIS SCRIPT IS IN DEVELOPMENT!')

        # Get the sky spectrum
        sky_spectrum = get_sky_spectrum()
    
        # Get the atmospheric throughput
        atmospheric_throughput = efficiency.AtmosphericThroughput(airmass=args.airmass)

        # Set the telescope. Defines the aperture area and throughput
        # (nominally 3 aluminum reflections for Keck)
        telescope = telescopes.KeckTelescope()

        # Set the detector
        qe_file = data_file(filename='efficiency') / 'detectors' / 'Andor-BEX2-DD_qe.db'
        qe = efficiency.Efficiency.from_file(str(qe_file)) 
        # Assume it's square
        det = detector.Detector(shape=(4096,4096), pixelsize=15., rn=9.8, dark=1.3e-4,
                                fullwell=3.5e5, qe=qe)
        pixelsize = args.binning * det.pixelsize*1e-3 / telescope.platescale
        image_size = det.shape[0] * det.pixelsize*1e-3 / telescope.platescale

        # Get the source distribution.  If the source is uniform, onsky is None.
        # TODO: Enable Moffat
        onsky = get_source_distribution(args.fwhm, beta=args.beta, sersic=args.sersic,
                                        size=min(image_size, 10*args.fwhm), sampling=pixelsize)

        # "Extraction" aperture
        phot_aperture = aperture.CircularAperture(0, 0, 1.5*args.fwhm, resolution=100)
        phot_ap_img = phot_aperture.response(onsky.x, onsky.y)

        n = 30
        mag = numpy.linspace(22, 35, n)
        band = ['u', 'g', 'r', 'i', 'z']
        rn = [2.1, 9.8]

        pix_snr = numpy.zeros((n,5,2), dtype=float)
        phot_snr = numpy.zeros((n,5,2), dtype=float)

        for i in range(n):
            for j in range(5):
                for k in range(2):
                    print(f'{i+1}/{n} ; {j+1}/5 ; {k+1}/2', end='\r')

                    det.rn = rn[k]

                    # Setup the filter
                    broadband = efficiency.FilterResponse(band=band[j])

                    # Setup the wavelength grid for the target spectrum
                    # TODO: Is this correct for sampling a true input spectrum?
                    dw = numpy.mean(numpy.diff(broadband.wave))
                    dlogw = numpy.log10(broadband.wave[0]+dw) - numpy.log10(broadband.wave[0])
                    wave = get_wavelength_vector(broadband.wave[0]-10*dw, broadband.wave[-1]+10*dw,
                                                 dlogw)
                    resolution = numpy.mean(broadband.wave)/dw

                    # Get source spectrum in 1e-17 erg/s/cm^2/angstrom. Currently, the
                    # source spectrum is assumed to be
                    #   - normalized by the total integral of the source flux 
                    #   - independent of position within the source
                    spec = get_spectrum(wave, mag[i], mag_band=band[j], mag_system=args.mag_system,
                                        resolution=resolution)

                    # Total number of object counts per exposure
                    _object_flux = numpy.sum(spec.photon_flux(inplace=False) * spec.wavelength_step()
                                             * atmospheric_throughput(spec.wave)
                                             * telescope.area * telescope.throughput(spec.wave)
                                             * broadband(spec.wave) * det(spec.wave) 
                                             * args.time / args.nexp)

                    # Total number of sky counts detected per exposure per square arcsec
                    _sky_flux = numpy.sum(sky_spectrum.photon_flux(inplace=False) 
                                          * sky_spectrum.wavelength_step()
                                          * telescope.area * telescope.throughput(sky_spectrum.wave)
                                          * broadband(sky_spectrum.wave) * det(sky_spectrum.wave) 
                                          * args.time / args.nexp)

                    # Construct the sum of the signal for all images
                    signal = args.nexp * _object_flux * onsky.data * onsky.sampling**2

                    # Construct the variance in the image sum
                    variance = (signal/args.nexp + _sky_flux * onsky.sampling**2 + det.rn**2 \
                                + args.binning**2 * det.dark * args.time / args.nexp) * args.nexp

                    # S/N image
                    snr = signal / numpy.sqrt(variance)
                    pix_snr[i,j,k] = numpy.amax(snr)

                    ext_signal = numpy.sum(phot_ap_img * signal)
                    ext_variance = numpy.sum(phot_ap_img**2 * variance)

                    phot_snr[i,j,k] = ext_signal/numpy.sqrt(ext_variance)

        print(f'{n}/{n} ; 5/5 ; 2/2')

        p2 = numpy.zeros((5,2), dtype=float)
        m5 = numpy.zeros((5,2), dtype=float)
        for j in range(5):
            for k in range(2):
                p2[j,k] = interpolate.interp1d(pix_snr[:,j,k], mag)(2.0)
                m5[j,k] = interpolate.interp1d(phot_snr[:,j,k], mag)(5.0)

        embed()
        exit()




