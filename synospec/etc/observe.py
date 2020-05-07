"""
Module used to mimic observations.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

#Extraction:
#    - fiber-extraction aperture (pixels, optimal?)
#    - fiber PSF FWHM on detector (pixels)
#
#    - spectral resolution (R)
#    - dispersion on detector (A per pixel)
#
#Detector:
#    - arcsec / pixel (spatial)
#    - angstrom / pixel (spectral)
#    - detector readnoise (e-)
#    - detector darkcurrent (e-/hour)
#
#Throughput
#    - spectrograph throughput
#        - detector QE
#        - camera efficiency
#        - grating efficiency
#        - focal-plane to grating tansmission efficiency
#            - foreoptics, lenslet, coupling, fiber transmission, FRD,
#              collimator, dichroic(s)
#
#    - top-of-telescope throughput
#        - telescope efficiency
#        - spectrograph efficiency
#
#Point-source Aperture losses
#    - fiber diameter (arcsec, mm)
#    - focal-plane plate scale (arcsec/mm)
#    - seeing (arcsec)
#    - pointing error (arcsec)
#
#*Source:
#    - source spectrum (1e-17 erg/s/cm^2/angstrom)
#
#    - point source
#        - source mag (at wavelength/in band above)
#
#    - extended source
#        - source surface brightness (at wavelength/in band above)
#        - source size (kpc/arcsec)
#        - cosmology
#
#    - source velocity/redshift
#
#*Sky:
#    - sky spectrum (1e-17 erg/s/cm^2/angstrom)
#
#Observation:
#    - obs date
#    - lunar phase/illumination (days, fraction)
#    - sky coordinates
#    - wavelength for calculation (angstroms)
#    - band for calculation
#    - sky surface brightness (mag/arcsec^2, at wavelength/in band above)
#    - airmass
#    - exposure time (s)
#
#Telescope:
#    - location
#    - telescope diameter (or effective) (m^2)
#    - central obstruction

import os
import warnings

from IPython import embed

import numpy

from scipy import signal

from matplotlib import pyplot

from . import source, spectrum, util

class Observation:
    """
    Observation of the sky with or without a source.

    Args:
        telescope (:class:`~synospec.etc.telescopes.Telescope`):
            Telescope used for the observation.
        sky_spectrum (:class:`~synospec.etc.spectrum.Spectrum`):
            Sky spectrum
        spec_aperture (:class:`~synospec.etc.aperture.Aperture`):
            On-sky entrance aperture.
        exposure_time (:obj:`float`):
            Exposure time in seconds.
        detector (:class:`~synospec.etc.detector.Detector`):
            Instrument detector object.
        system_throughput (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            System throughput; i.e., ratio of the number of photons
            detected to the total number incident on the telescope
            primary.  If None, system throughput is unity.
        atmospheric_throughput (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Atmospheric throughput; i.e., the ratio of the number of
            the photons incident on the telescope primary to the
            number of photons hitting the top of the atmosphere. If
            None, assume there is no atmosphere.
        airmass (:obj:`float`, optional):
            Airmass of the observation. If None, assume
            ``atmospheric_throughput`` is at the correct airmass or
            there is no atmosphere.
        onsky_source_distribution (:class:`~synospec.etc.source.Source`, optional):
            On-sky source surface-brightness distribution, including
            any seeing effects. If None, assume observation is only
            of the sky. Normalization is irrelevant; used to
            calculate aperture losses.
        source_spectrum (:class:`~synospec.etc.spectrum.Spectrum`, optional):
            Spectrum of the source with the flux normalized to the
            *total* source flux.
        extraction (:class:`~synospec.etc.extract.Extraction`, optional):
            Object used to calculate the extraction losses and tally
            the number of pixels included in the extraction to get
            the read-noise hit. TODO: Currently *cannot* be None;
            make a required argument.
        snr_units (:obj:`str`, optional):
            Units of the S/N calculation. Must be ``'pixel'``,
            ``'angstrom'``, or ``'resolution'`` for S/N per pixel,
            per angstrom, or per resolution element.

    Raises:
        ValueError:
            Raised if the requested S/N units cannot be calculated
            (i.e., the provided spectrum does not define the
            resolution meaning that the S/N per resolution element is
            undefined) or if the unit string for the S/N is not
            recognized.

    """
    def __init__(self, telescope, sky_spectrum, spec_aperture, exposure_time, detector,
                 system_throughput=None, atmospheric_throughput=None, airmass=None,
                 onsky_source_distribution=None, source_spectrum=None, extraction=None,
                 snr_units='pixel'):

        # Save the input or use defaults
        self.source = onsky_source_distribution

        # Sky spectrum is expected to be independent of position within
        # the aperture and be the sky flux density per unit area, where
        # the unit area is defined by the aperture object (arcsec^2)
        # i.e., the units are, e.g., erg/s/cm^2/angstrom/arcsec^2
        self.source_spectrum = source_spectrum

        # Match the sampling of the sky and source spectrum, if possible
        self.sky_spectrum = sky_spectrum if self.source_spectrum is None else \
                                spectrum.Spectrum(self.source_spectrum.wave,
                                                  sky_spectrum.interp(self.source_spectrum.wave),
                                                  log=self.source_spectrum.log)

        # In the current implementation, the source spectrum is expected
        # to be:
        #   - independent of position within the source
        #   - the source flux density integrated over the full source
        #     distribution; i.e., units are, e.g., erg/s/cm^2/angstrom
        self.wave = self.sky_spectrum.wave.copy()
        self.sres = self.sky_spectrum.sres.copy() if sky_spectrum.sres is not None \
                        else (None if self.source_spectrum is None 
                                else self.source_spectrum.sres.copy())

        self.atmospheric_throughput = atmospheric_throughput
        self.telescope = telescope
        self.aperture = spec_aperture
        self.system_throughput = system_throughput
        self.detector = detector
        self.exptime = exposure_time
        self.extraction = extraction

        # Get the "aperture factor". If the source distribution is not
        # provided, the source surface brightness is assumed to be
        # uniform within the aperture (like the sky). In this case, the
        # aperture factor is the area of the aperture itself so that
        # the object flux is the integral of the surface brightness
        # over the aperture size (like the sky).
        self.aperture_factor = self.aperture.area if self.source is None \
                                     else self.aperture.integrate_over_source(self.source) \
                                                / self.source.integral

        # Get the total object flux incident on the focal plane in
        # electrons per second per angstrom
        _object_flux = numpy.zeros(self.wave.size, dtype=float) if self.source_spectrum is None \
                            else self.source_spectrum.photon_flux(inplace=False) \
                                    * self.telescope.area * self.aperture_factor \
                                    * self.detector(self.wave)
        if self.atmospheric_throughput is not None:
            _object_flux *= self.atmospheric_throughput(self.wave)
        if self.system_throughput is not None:
            _object_flux *= self.system_throughput(self.wave)

        # Total sky flux in electrons per second per angstrom; the
        # provided sky spectrum is always assumed to be uniform over the
        # aperture
        _sky_flux = self.sky_spectrum.photon_flux(inplace=False) \
                        * self.telescope.area * self.aperture.area \
                        * self.detector(self.wave)
        if self.system_throughput is not None:
            _sky_flux *= self.system_throughput(self.wave)

        # Set the units for the output:
        dw = self.sky_spectrum.wavelength_step()
        if snr_units == 'pixel':
            _object_flux *= dw
            _sky_flux *= dw
            spectral_width = 1.
        elif snr_units == 'angstrom':
            spectral_width = 1./dw
        elif snr_units == 'resolution':
            if self.sres is None:
                raise ValueError('Cannot compute S/N per resolution element without resolution '
                                 'vector.')
            _object_flux *= self.wave/self.sres
            _sky_flux *= self.wave/self.sres
            spectral_width = self.wave/self.sres/dw
        else:
            raise ValueError('Unknown S/N units requested.')

        # Observe and extract the source
        # TODO: Set spectral_pixels...
        self.object_flux, self.obj_shot_var, self.sky_flux, self.sky_shot_var, self.read_var \
                = self.extraction.sum_signal_and_noise(_object_flux, _sky_flux, self.exptime,
                                                       spectral_width=spectral_width)

    def simulate(self, sky_only=False, sky_sub=False, sky_err=0.1):
        """
        Return a simulated spectrum.

        Args:
            sky_only (:obj:`bool`, optional):
                Only include the sky flux in the simulated spectrum
                (no object flux)
            sky_sub (:obj:`bool`, optional):
                Provide the object spectrum only, but include the sky
                flux shot noise and additional error from the sky
                subtraction.
            sky_err (:obj:`float`, optional):
                The fraction of the total sky error incurred due to
                the sky subtraction. Should be between 0 and 1; 0
                means no additional error is incurred, 1 means that
                the sky noise from the sky subtration is the same as
                the sky noise from the observation itself.

        Returns:
            :class:`~synospec.etc.spectrum.Spectrum`: The simulated
            spectrum.
        """
        if sky_only:
            shot_var = self.sky_shot_var
            flux = self.sky_flux
        elif sky_sub:
            shot_var = self.obj_shot_var + (1 + numpy.square(sky_err))*self.sky_shot_var
            flux = self.object_flux
        else:
            shot_var = self.obj_shot_var + self.sky_shot_var
            flux = self.object_flux + self.sky_flux

#        error=numpy.sqrt(shot_var + self.read_var)
#        draw = numpy.random.normal(scale=error)
#        return spectrum.Spectrum(self.wave, flux + draw, error=error,
#                                 log=self.sky_spectrum.log if self.source_spectrum is None
#                                         else self.source_spectrum.log)

        # Draw from a Poisson distribution for the shot noise,
        # subtracted the expectation value of the distribution so that
        # only the noise is added
        shot_draw = numpy.random.poisson(lam=shot_var)-shot_var
        # Draw from a Gaussian distribution for the read noise
        read_draw = numpy.random.normal(scale=numpy.sqrt(self.read_var))
        return spectrum.Spectrum(self.wave, flux + shot_draw + read_draw,
                                 error=numpy.sqrt(shot_var + self.read_var),
                                 log=self.sky_spectrum.log if self.source_spectrum is None
                                         else self.source_spectrum.log)

    def snr(self, sky_sub=False, sky_err=0.1):
        """
        Calculate the S/N.

        Args:
            sky_sub (:obj:`bool`, optional):
                Provide the object spectrum only, but include the sky
                flux shot noise and additional error from the sky
                subtraction.
            sky_err (:obj:`float`, optional):
                The fraction of the total sky error incurred due to
                the sky subtraction. Should be between 0 and 1; 0
                means no additional error is incurred, 1 means that
                the sky noise from the sky subtration is the same as
                the sky noise from the observation itself.

        Returns:
            :class:`~synospec.etc.spectrum.Spectrum`: The S/N spectrum.
            WARNING: Here, :class:`~synospec.etc.spectrum.Spectrum` is
            used as a container class for the S/N vector; some
            functionality of the class will not be valid!
        """
        flux = self.object_flux + self.sky_flux
        var = self.obj_shot_var + self.sky_shot_var + self.read_var
        if sky_sub:
            flux -= self.sky_flux
            var += numpy.square(sky_err)*self.sky_shot_var
        # TODO: add additional noise from sky subtraction
        return spectrum.Spectrum(self.wave, flux / numpy.sqrt(var),
                                 log=self.sky_spectrum.log if self.source_spectrum is None
                                         else self.source_spectrum.log)


def monochromatic_image(sky, spec_aperture, spec_kernel, platescale, pixelsize, onsky_source=None,
                        scramble=False):
    """
    Generate a monochromatic image of the sky, with or with out a
    source, taken by a spectrograph through an aperture.

    .. warning::

        - May resample the `source` and `spec_kernel` maps.

    .. todo::

        - Add effect of differential atmospheric refraction
        - Allow a force map size and pixel sampling

    Args:
        sky (:class:`synospec.etc.source.OnSkySource`):
            Sky flux distribution.
        spec_aperture (:class:`synospec.etc.aperture.Aperture`):
            Spectrograph aperture. Aperture is expected to be
            oriented with the dispersion along the first axis (e.g.,
            the slit width is along the abcissa).
        spec_kernel (:class:`synospec.etc.kernel.SpectrographGaussianKernel`):
            Convolution kernel describing the point-spread function
            of the spectrograph.
        platescale (:obj:`float`):
            Platescale in mm/arcsec at the detector
        pixelsize (:obj:`float`):
            Size of the detector pixels in mm.
        onsky_source (:class:`synospec.etc.source.OnSkySource`, optional):
            On-sky distribution of the source flux. If None, only sky
            is observed through the aperture.
        scramble (:obj:`bool`, optional):
            Fully scramble the source light passing through the
            aperture. This should be False for slit observations. For
            fiber observations, this should be True and makes the
            nominal assumption that the focal plane incident on the
            fiber face is perfectly scrambled.
    """
    # Check input
    if onsky_source is not None:
        if onsky_source.sampling is None or onsky_source.size is None:
            warnings.warn('Source was not provided with an initial map sampling; doing so now '
                          'with default sampling and size.')
            onsky_source.make_map()

    # Detector pixel scale in arcsec/pixel
    pixelscale = pixelsize/platescale

    # Assume the sampling of the source is provided with the maximum
    # allowed pixel size. Determine a pixel size that is no more than
    # this, up to an integer number of detector pixels. Use an integer
    # number specifically so that the result of the kernel convolution
    # can be simply rebinned to match the detector pixel size.
    # `sampling` is in arcsec per pixel
    sampling = pixelscale if onsky_source is None else min(onsky_source.sampling, pixelscale)
    oversample = int(pixelscale/sampling)+1 if sampling < pixelscale else 1
    sampling = pixelscale/oversample

    # Assume the size of the image properly samples the source.
    # Determine a map size that at least encompasses the input source
    # and the input aperture. The factor of 1.5 is ad hoc; could likely
    # be lower. `size` is in arcsec
    dx, dy = 1.5*numpy.diff(numpy.asarray(spec_aperture.bounds).reshape(2,-1), axis=0).ravel()
    size = max(dx,dy) #if onsky_source is None else max(onsky_source.size, dx, dy)

    # TODO: Below alters `source` and `spec_kernel`. Should maybe
    # instead save the old sampling and size and then resample back to
    # the input before returning.

    # Resample the distribution maps. Classes OnSkySource and Aperture
    # use arcsecond units.
    if onsky_source is not None:
        onsky_source.make_map(sampling=sampling, size=size)
    sky.make_map(sampling=sampling, size=size)
    ap_img = spec_aperture.response(sky.x, sky.y)

    # However, SpectrographGaussianKernel uses mm, so we need to
    # convert sampling from arcsec/pixel to mm/pixel using the
    # platescale.
    spec_kernel.resample(pixelscale=platescale*sampling)

    # Construct the image. Note that scrambling simply scales the
    # aperture image by the flux that enters it; otherwise, the source
    # image is just attenuated by the aperture response map.
    source_img = sky.data if onsky_source is None else onsky_source.data + sky.data
    input_img = ap_img*numpy.sum(source_img*ap_img)/numpy.sum(ap_img) if scramble \
                        else source_img*ap_img

    # Convolve it with the spectrograph imaging kernel
    mono_img = signal.fftconvolve(input_img, spec_kernel.array, mode='same')
    # Return the image, downsampling if necessary
    return util.boxcar_average(mono_img, oversample) if oversample > 1 else mono_img


def twod_spectrum(sky_spectrum, spec_aperture, spec_kernel, platescale, linear_dispersion,
                  pixelsize, source_distribution=None, source_spectrum=None, thresh=None,
                  scramble=False, wave_lim=None, field_coo=None, opticalmodel=None):
    """
    Documentation TBW.

    if optical model is not provided:
        - ignore field_coo
        - return rectilinear 2D spectrum

    platescale is in mm/arcsec
    linear_dispersion is in A/mm
    pixelsize is in mm

    """
    # Ensure that the sky spectrum and source spectrum will have the same wavelength limits
    if source_spectrum is not None and wave_lim is None:
        wave_lim = source_spectrum.wave[[0,-1]]

    # Get the sky-only monochromatic image
    sky = source.OnSkyConstant(1.0)
    sky_img = monochromatic_image(sky, spec_aperture, spec_kernel, platescale, pixelsize,
                                  scramble=scramble)

    # Renormalize the sky slit image such that the integral is the area of the aperture.
    sky_img *= spec_aperture.area / numpy.sum(sky_img)/numpy.square(sky.sampling)

    source_img = None
    if source_distribution is not None and source_spectrum is not None:
        # Reset the source distribution map
        source_distribution.reset_map()

        # Get the source-only monochromatic image
        sky = source.OnSkyConstant(0.0)
        source_img = monochromatic_image(sky, spec_aperture, spec_kernel, platescale, pixelsize,
                                         onsky_source=source_distribution, scramble=scramble)
        # Renormalize the source image by the integral of the onsky-source;
        # this maintains the effects of aperture losses
        source_img /= source_distribution.integral

    s = numpy.array([0,0])
    e = numpy.array([*sky_img.shape])
    if thresh is not None:
        indx = sky_img > thresh
        s, e = numpy.append(numpy.where(numpy.any(indx, axis=1))[0][[0,-1]],
                            numpy.where(numpy.any(indx, axis=0))[0][[0,-1]]).reshape(2,2).T
        if source_img is not None:
            indx = source_img > thresh
            _s, _e = numpy.append(numpy.where(numpy.any(indx, axis=1))[0][[0,-1]],
                                  numpy.where(numpy.any(indx, axis=0))[0][[0,-1]]).reshape(2,2).T
            s = numpy.minimum(s, _s)
            e = numpy.maximum(e, _e)

    sky_img = sky_img[s[0]:e[0],s[1]:e[1]]
    if source_img is not None:
        source_img = source_img[s[0]:e[0],s[1]:e[1]]

    # Get the 2D spectrum
    dispscale = linear_dispersion*pixelsize
    wave0, sky_2d_spec = rectilinear_twod_spectrum(sky_spectrum, sky_img, dispscale,
                                                   wave_lim=wave_lim)

    if source_distribution is None or source_spectrum is None:
        return sky_2d_spec if opticalmodel is None \
                    else opticalmodel.project_2d_spectrum(sky_2d_spec, platescale,
                                                          linear_dispersion, pixelsize, wave0,
                                                          field_coo=field_coo)

    # Get the 2D spectrum
    wave0, source_2d_spec = rectilinear_twod_spectrum(source_spectrum, source_img, dispscale,
                                                      wave_lim=wave_lim)

    return sky_2d_spec + source_2d_spec if opticalmodel is None \
                    else opticalmodel.project_2d_spectrum(sky_2d_spec + source_2d_spec, platescale,
                                                          linear_dispersion, pixelsize, wave0,
                                                          field_coo=field_coo)

#def rectilinear_twod_spectrum(spectrum, aperture_image, dispscale, wave_lim=None, oversample=1):
#    """
#    Construct a rectilinear 2D spectrum.
#
#    spectral dimension of aperture_image is along the first axis
#
#    aperture_image has to be odd?
#    """
#
#    # TODO: Let dispersion scale be non-linear?
#
#    # Resample the spectrum to the appropriate dispersion scale up to some constant
#    if wave_lim is None:
#        # TODO: This should be the pixel boundaries, not the pixel centers!
#        wave_lim = spectrum.wave[[0,-1]]
#    resamp_wave = numpy.arange(wave_lim[0], wave_lim[1] + dispscale/oversample,
#                               dispscale/oversample)
#    resampled_spectrum = spectrum.resample(resamp_wave)
##    pyplot.plot(spectrum.wave, spectrum.flux)
##    pyplot.plot(resampled_spectrum.wave, resampled_spectrum.flux)
##    pyplot.show()
#
#    # Oversample the image spectrally
#    _aperture_image = util.block_replicate(aperture_image, (oversample,1)) \
#                            if oversample > 1 else aperture_image
#
#    # Number of spectral and spatial channels in the aperture image
#    nspec, nspat = _aperture_image.shape
#    width = nspec//2
#    # Length of the spectrum
#    npix = len(resampled_spectrum)
#
#    # Pad the spectrum with zeros and create one shifted copy per spectral channel
#    flux = numpy.zeros((nspec,npix), dtype=float)
#    for i in range(nspec):
#        flux[i,width:-width] = resampled_spectrum.flux[i:npix-nspec+i+1]
#
#    import time
#    t = time.perf_counter()
#    flux = numpy.tile(flux, (nspat,1))
#    print('Tile: {0}s'.format(time.perf_counter()-t))
#
#    # Do the convolution
#    t = time.perf_counter()
#    twodspec = flux[:,:] * _aperture_image.ravel()[:,None]
#    indx = numpy.arange(0,nspec*nspat,nspec)
#    twodspec = numpy.add.reduceat(twodspec, indx, axis=0)
#    print('Mult: {0}s'.format(time.perf_counter()-t))
#
#    return util.block_average(twodspec, (oversample,1)) if oversample > 1 else twodspec


def rectilinear_twod_spectrum(spectrum, aperture_image, dispscale, wave_lim=None, oversample=1):
    """
    Construct a rectilinear 2D spectrum.

    Documentation TBW.

    spectral dimension of aperture_image is along the first axis

    aperture_image has to be odd?
    dispscale is A/pixel

    remove rows columns with no pixels above thresh

    """

    # TODO: Let dispersion scale be non-linear?

    # Resample the spectrum to the appropriate dispersion scale up to some constant
    if wave_lim is None:
        # TODO: This should be the pixel boundaries, not the pixel centers!
        wave_lim = spectrum.wave[[0,-1]]
    resamp_wave = numpy.arange(wave_lim[0], wave_lim[1] + dispscale/oversample,
                               dispscale/oversample)
    resampled_spectrum = spectrum.resample(resamp_wave)
    wave0 = numpy.mean(resampled_spectrum.wave[:oversample])
#    pyplot.plot(spectrum.wave, spectrum.flux)
#    pyplot.plot(resampled_spectrum.wave, resampled_spectrum.flux)
#    pyplot.show()

    # Oversample the image spectrally
    _aperture_image = util.block_replicate(aperture_image, (oversample,1)) \
                            if oversample > 1 else aperture_image.copy()

    # Number of spectral and spatial channels in the aperture image;
    # number of convolution kernels is nspat
    nspat, nspec = _aperture_image.shape

    twodspec = numpy.zeros((len(resampled_spectrum),nspat), dtype=float)
    for i in range(nspat):
        twodspec[:,i] = signal.fftconvolve(resampled_spectrum.flux, _aperture_image[i],
                                           mode='same')
    
    return wave0, (util.block_average(twodspec, (oversample,1)) if oversample > 1 else twodspec)

