"""
Simple extraction object

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy
from scipy import special

# TODO: Include an optimal extraction class or method

class Extraction:
    """
    Perform some nominal extraction calculations.

    The primary purpose of this object is to account for how flux is
    distributed at the detector for a given spectral element. In both
    the spatial and spectral dimensions, the ``fwhm`` argument is the
    the FWHM of the intensity distribution and ``width`` is the
    extraction aperture, both in pixels.

    The profile describes how the light is distributed across the
    detector pixels. For now (i.e., until this is improved), use
    ``profile='gaussian'`` for a fiber aperture, and
    ``profile='uniform'`` for a slit aperture. For a uniform profile,
    the FWHM and width should be the same, but the FWHM is never
    used. The spectral profile is constructed, but never used, the
    idea being the loss of one monochromatic element is replaced by
    the adjacent one.

    .. todo::

        Allow the 2D profile to be provided directly.

    Args:
        detector (:class:`~synospec.etc.detector.Detector`):
            Object with the detector data.
        spatial_fwhm (:obj:`float`, optional):
            The FWHM in pixels of the spatial PSF at the spectrograph
            detector.
        spatial_width (:obj:`int`, optional):
            The number of pixels in the extraction.
        spectral_fwhm (:obj:`float`, optional):
            The FWHM in pixels of the spectral LSF at the
            spectrograph detector. **NOT USED**.
        spectral_width (:obj:`int`, optional):
            The number of pixels in the extraction. **NOT USED**.
        profile (:obj:`str`, optional):
            The profile shape. Must be ``'gaussian'`` or
            ``'uniform'``.
    """
    def __init__(self, detector, spatial_fwhm=3, spatial_width=3,
                 spectral_fwhm=3, spectral_width=3, profile='gaussian'):
        self.detector = detector

        # All in pixel units
        self.spatial_fwhm = spatial_fwhm
        self.spatial_width = int(spatial_width)
        self.spectral_fwhm = spectral_fwhm
        self.spectral_width = int(spectral_width)

        self._get_profile(profile)

    @staticmethod
    def _pixelated_gaussian(fwhm, width):
        """
        Construct the profile of the monochromatic element as a
        pixelated Gaussian.

        Args:
            fwhm (:obj:`float`):
                The Gaussican FWHM in pixels.
            width (:obj:`int`, optional):
                The number of pixels to include in the profile. The
                center of the Gaussian is always at ``width/2``
                pixels.

        Returns:
            :obj:`tuple`: Two float vectors are returned: (1) a
            vector with shape ``(width+1,)`` with the coordinates of
            the pixel edges in units of the Gaussian dispersion and
            (2) a vector with shape ``(width,)`` with the profile.
            The sum of the profile is the integral of the Gaussian
            over the window ``[-width/2,width/2]``.
        """
        # Pixel edges in units of the profile dispersion
        edges = width/fwhm * numpy.sqrt(8*numpy.log(2)) * numpy.linspace(-1, 1, width+1)/2
        profile = (special.erf(edges[1:]/numpy.sqrt(2))
                        - special.erf(edges[:-1]/numpy.sqrt(2)))/2.
        return edges, profile

    @staticmethod
    def _pixelated_uniform(width):
        """
        Construct the profile of the monochromatic element as a
        uniform top-hat function.

        Args:
            width (:obj:`int`):
                The number of pixels to include in the profile.

        Returns:
            :obj:`tuple`: Two float vectors are returned: (1) a
            vector with shape ``(width+1,)`` with the coordinates of
            the pixel edges in pixels and (2) a vector with shape
            ``(width,)`` with the profile. The sum of the profile is
            set to unity (i.e., the extraction is perfect).
        """
        # Pixel edges in units of the profile dispersion
        edges = width * numpy.linspace(-1, 1, width+1)/2
        profile = numpy.full(width, 1/width, dtype=float)
        return edges, profile

    def _get_profile(self, profile):
        """
        Construct the spatial profile over the extraction aperture.

        Args:
            profile (:obj:`str`):
                The profile shape. Must be ``'gaussian'`` or
                ``'uniform'``.

        Raises:
            ValueError:
                Raised if the profile string is not recognized.
        """
        if profile == 'gaussian':
            self.spatial_edges, self.spatial_profile \
                    = Extraction._pixelated_gaussian(self.spatial_fwhm, self.spatial_width)
            self.spectral_edges, self.spectral_profile \
                    = Extraction._pixelated_gaussian(self.spectral_fwhm, self.spectral_width)
            return

        if profile == 'uniform':
            self.spatial_fwhm = self.spatial_width
            self.spatial_edges, self.spatial_profile \
                    = Extraction._pixelated_uniform(self.spatial_width)
            self.spectral_fwhm = self.spectral_width
            self.spectral_edges, self.spectral_profile \
                    = Extraction._pixelated_uniform(self.spectral_width)
            return

        raise ValueError('{0} profile not recognized.  Must be \'gaussian\' or \'uniform\'')

        # 2D profile? ...
#        self.profile = self.spectral_profile[:,None]*self.spatial_profile[None,:]

    def sum_signal_and_noise(self, object_flux, sky_flux, exposure_time, spectral_width=1.):
        """
        Get the signal and noise from a summed extraction.

        Fluxes should be in electrons per second per resolution
        element, per angstrom, or per pixel.

        The primary operation is to assume the input spectrum is
        distributed in a spatial and spectral profile that is summed.
        This operation *only* sums over the spatial profile. The
        spectral width provided is used to set the number of pixels
        extracted to construct the desired S/N units (per resolution
        element, angstrom, or pixel). To get the S/N per pixel, use
        the default. Otherwise, ``spectral_width`` is, e.g., the
        number of pixels per angstrom.

        Args:
            object_flux (:obj:`float`, `numpy.ndarray`_):
                Object flux in electrons per second per resolution
                element, per angstrom or per pixel. Type/Shape must
                match ``sky_flux``.
            sky_flux (:obj:`float`, `numpy.ndarray`_):
                Sky flux in electrons per second per resolution
                element, per angstrom or per pixel. Type/Shape must
                match ``sky_flux``.
            exposure_time (:obj:`float`):
                Exposure time in seconds.
            spectral_width (:obj:`float`, optional):
                The number of pixels in the element relevant to the
                flux units. This sets how many read-noise hits to
                include in the total noise calculation.

        Returns:
            :obj:`tuple`: Returns five objects: (1) the extracted
            object flux, (2) the photon shot variance from the object
            (includes any dark current), (3) the extracted sky flux,
            (4) the photon shot variance from the sky (includes any
            dark current), and (5) the read-noise variance. The
            returned fluxes are in units of electrons per resolution
            element, per angstrom, or per pixel, depending on the
            input.
        """
        # Get the variance in each pixel
        extract_box_n = len(self.spatial_profile) * spectral_width
        read_var = numpy.square(self.detector.rn) * extract_box_n
        if isinstance(object_flux, numpy.ndarray):
            ext_obj_flux = numpy.sum(object_flux[None,:]*self.spatial_profile[:,None], axis=0)
            ext_sky_flux = numpy.sum(sky_flux[None,:]*self.spatial_profile[:,None], axis=0)
            if not isinstance(spectral_width, numpy.ndarray):
                read_var = numpy.full_like(ext_obj_flux, read_var, dtype=float)
        else:
            if isinstance(spectral_width, numpy.ndarray):
                raise TypeError('Cannot use vector of spectral widths with individual fluxes.')
            ext_obj_flux = numpy.sum(object_flux*self.spatial_profile)
            ext_sky_flux = numpy.sum(sky_flux*self.spatial_profile)

        ext_obj_flux *= exposure_time
        ext_sky_flux *= exposure_time
        obj_shot_var = ext_obj_flux + self.detector.dark * extract_box_n * exposure_time
        sky_shot_var = ext_sky_flux + self.detector.dark * extract_box_n * exposure_time

        return ext_obj_flux, obj_shot_var, ext_sky_flux, sky_shot_var, read_var

    # TODO: Add optimal_signal_and_noise()

#    def sum_noise_budget(self, object_flux, sky_flux):
#        """
#        Get the noise budget from a summed extraction.
#
#        Fluxes should be in electrons per resolution element
#
#        object_flux and sky_flux can be spectra
#        """
#        # Get the variance in each pixel
#        n_pixels = numpy.nprod(numself.profile.shape)
#        if isinstance(object_flux, numpy.ndarray):
#            object_err = numpy.sqrt(numpy.sum(object_flux[None,None,:] * self.profile[:,:,None],
#                                              axis=0))
#            sky_err = numpy.sqrt(numpy.sum(sky_flux[None,None,:] * self.profile[:,:,None], axis=0))
#
#            dark_err = numpy.full_like(object_err, numpy.sqrt(dark*n_pixels))
#            read_err = numpy.full_like(object_err, detector.rn*numpy.sqrt(n_pixels))
#            return object_err, sky_err, dark_err, read_err
#
#        return numpy.sqrt(object_flux*self.efficiency), numpy.sqrt(sky_flux*self.efficiency), \
#                numpy.sqrt(dark*n_pixels), detector.rn*numpy.sqrt(n_pixels)
    
    @property
    def spatial_efficiency(self):
        """Return the extraction efficiency."""
        return numpy.sum(self.spatial_profile)


