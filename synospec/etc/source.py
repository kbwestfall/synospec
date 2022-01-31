"""
Module defining source surface-brightness distributions.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import warnings
import numpy

from scipy import signal, special, interpolate

from astropy.modeling import functional_models

from ..util.frame import SemiMajorAxisCoo

# TODO: Need to learn how to use abstract classes!  Inherit from
# numpy.ndarray?
class Source:
    """
    This is an abstract class and should not be instantiated on it's
    own!

    Attributes:
        x (vector):
            1D vector with x coordinates
        X (array):
            2D map of x coordinates
        y (vector):
            1D vector with y coordinates
        Y (array):
            2D map of y coordinates
        data (array):
            Map of the surface brightness distribution
        sampling (float):
            Sampling of the square map in arcsec/pixel.
        size (float):
            Size of the square map in arcsec.
    """
    def __init__(self):
        self.x = None
        self.X = None
        self.y = None
        self.Y = None
        self.data = None
        self.sampling = None
        self.size = None

    def __call__(self, x, y):
        pass
    
    def minimum_sampling(self):
        pass

    def minimum_size(self):
        pass

    def make_map(self, sampling=None, size=None):
        """
        Generate a square map of the surface-brightness distribution.

        Object is modified internally with the constructed map saved
        to :attr:`data`.

        Args:
            sampling (:obj:`float`, optional):
                The size of each pixel in the map, typically in
                arcsec per pixel. If None, sampling will be set by
                :func:`minimum_sampling`.
            size (:obj:`float`, optional):
                Length of one axis of the *square* map, typically in
                arcsec. Will be adjusted to ensure an integer number
                of pixels in the map based on :attr:`sampling`. If
                None, set by :func:`minimum_size`.
        """
        if sampling is None and self.sampling is None:
            self.sampling = self.minimum_sampling()
        elif sampling is not None:
            self.sampling = sampling
            
        if size is None and self.size is None:
            self.size = self.minimum_size()
        elif size is not None:
            self.size = size

        # Set the number of pixels that spans the full size requested
        pixsize = numpy.ceil(self.size/self.sampling).astype(int)
        # Force the pixel size to be odd
        if pixsize % 2 == 0:
            pixsize += 1
        # Adjust the size to be an integer number of pixels
        _size = pixsize*self.sampling
        if _size - self.size > 0.1*self.sampling:
            warnings.warn('Size reset to an integer number of pixels: '
                          ' {0} -> {1} arcsec'.format(self.size, _size))
        self.size = _size
        self.y = (pixsize-1)*numpy.linspace(-0.5,0.5,pixsize)*self.sampling
        self.x = self.y.copy()[::-1]

        # Sample it
        self.X, self.Y = numpy.meshgrid(self.x, self.y)
        self.data = self.__call__(self.X, self.Y)

    def __getitem__(self, s):
        """Slice the map."""
        if self.data is None:
            raise ValueError('Distribution data is not defined!')
        return self.data[s]

    @property
    def shape(self):
        """The shape of the current map."""
        return () if self.data is None else self.data.shape

    def reset_map(self):
        """
        Reset mapping attributes for a fresh determination of the
        sampling, size, and surface-brightness map.

        This is mainly used for resetting the internals when using
        the default sampling and size to set the map. To reconstruct
        the map after calling this method, run :func:`make_map`.
        This::

            self.reset_map()
            self.make_map()

        is equivalent to this::

            self.make_map(sampling=self.minimum_sampling(), size=self.minimum_size())

        """
        self.x = None
        self.X = None
        self.y = None
        self.Y = None
        self.data = None
        self.sampling = None
        self.size = None 


class OnSkyConstant(Source):
    """
    An on-sky constant surface brightness.

    Args:
        surfbrightness (scalar-like):
            The constant surface brightness for the source in linear
            flux units per square arcsecond.
        sampling (scalar-like, optional):
            Sampling of a generated map in arcseconds per pixel.
            Default is set by :func:`minimum_sampling`.
        size (scalar-like, optional):
            Size of the image to generate of the distribution in
            *arceconds* along one of the axes.  The map is square.
            Default is defined by :func:`minimum_size`.
    """
    def __init__(self, surfbrightness, sampling=None, size=None):

        # Define internals
        self.surfbrightness = float(surfbrightness)
        # Instantiate the functional_models.Gaussian2D object
        super(OnSkyConstant, self).__init__()

        # Set the map sampling and size
        self.sampling = sampling
        self.size = size

        # Set the map if requested
        if sampling is not None or size is not None:
            self.make_map()

    def minimum_sampling(self):
        r"""
        Return the minimum sampling in arcseconds per pixels.
        Currently just set to 1.
        """
        # TODO: Allow `Source` to understand when this returns None?
        return 1.

    def minimum_size(self):
        r"""
        The minimum size that should be used for the distribution map in
        arcseconds.  Currently just set to 3.
        """
        # TODO: Allow `Source` to understand when this returns None?
        return 3.

    def __call__(self, x, y):
        """
        Return the surface brightness at a given location.
        """
        return numpy.full_like(x, self.surfbrightness, dtype=float) \
                    if isinstance(x, numpy.ndarray) else self.surfbrightness


class OnSkyGaussian(functional_models.Gaussian2D, Source):
    """
    An on-sky Gaussian distribution.

    Args:
        fwhm (scalar-like):
            The FWHM of the Gaussian in *arcseconds*.
        center (array-like, optional):
            The two coordinates of the Gaussian center in *arcseconds*.
        ellipticity (scalar-like, optional):
            The ellipticity (1-b/a) of an elliptical Gaussian
            distribution.
        position_angle (scalar-like, optional):
            The position angle for the elliptical Gaussian distribution,
            defined as the angle from N through E.  The coordinate
            system is defined with positive offsets (in RA) toward the
            east, meaning lower pixel indices.
        sampling (scalar-like, optional):
            Sampling of a generated map in arcseconds per pixel.
            Default is set by :func:`minimum_sampling`.
        size (scalar-like, optional):
            Size of the image to generate of the distribution in
            *arceconds* along one of the axes.  The map is square.
            Default is defined by :func:`minimum_size`.
    """
    def __init__(self, fwhm, center=None, ellipticity=None, position_angle=None, sampling=None,
                 size=None):

        # Define internals
        self.fwhm = float(fwhm)
        self.ellipticity = 0 if ellipticity is None else ellipticity
        self.position_angle = 0 if position_angle is None else position_angle
        sig2fwhm = numpy.sqrt(8*numpy.log(2))
        major_sigma = self.fwhm/sig2fwhm
        minor_sigma = major_sigma * (1-self.ellipticity)

        # Instantiate the functional_models.Gaussian2D object
        super(OnSkyGaussian, self).__init__(amplitude=1/(2*major_sigma*minor_sigma*numpy.pi),
                                            x_mean=0 if center is None else center[0],
                                            y_mean=0 if center is None else center[1],
                                            x_stddev=minor_sigma, y_stddev=major_sigma,
                                            theta=-numpy.radians(self.position_angle))

        # Set the integral to be normalized
        self.integral = 1.0

        # Set the map sampling and size
        self.sampling = sampling
        self.size = size

        # Set the map if requested
        if sampling is not None or size is not None:
            self.make_map()

    def get_integral(self):
        """Return the analytic integral of the source."""
        sig2fwhm = numpy.sqrt(8*numpy.log(2))
        major_sigma = self.fwhm/sig2fwhm
        minor_sigma = major_sigma * (1-self.ellipticity)
        return self.amplitude*2*numpy.pi*major_sigma*minor_sigma

    def minimum_sampling(self):
        r"""
        Return the minimum sampling in arcseconds per pixels.  Currently
        :math:`{\rm FWHM}/2`.
        """
        return self.fwhm/2.

    def minimum_size(self):
        r"""
        The minimum size that should be used for the distribution map in
        arcseconds.  Currently :math:`2\ {\rm FWHM}`.
        """
        return self.fwhm*2.


# TODO: Allow this to have an ellipticity and position angle?
class OnSkyMoffat(functional_models.Moffat2D, Source):
    """
    An on-sky Moffat distribution.

    Args:
        fwhm (scalar-like):
            The FWHM of the Moffat in *arcseconds*.
        beta (scalar-like):
            The power-law index for the Moffat tails (called alpha by astropy)
        center (array-like, optional):
            The two coordinates of the Moffat center in *arcseconds*.
        sampling (scalar-like, optional):
            Sampling of a generated map in arcseconds per pixel.
            Default is set by :func:`minimum_sampling`.
        size (scalar-like, optional):
            Size of the image to generate of the distribution in
            *arceconds* along one of the axes.  The map is square.
            Default is defined by :func:`minimum_size`.
    """
    def __init__(self, fwhm, beta, center=None, sampling=None, size=None):

        # Define internals
        self._fwhm = float(fwhm)
        alpha = float(beta)
        gamma = self._fwhm / 2 / numpy.sqrt(2**(1/alpha)-1)

        # Instantiate the functional_models.Gaussian2D object
        super(OnSkyMoffat, self).__init__(amplitude=(alpha-1)/numpy.pi/gamma**2,
                         x_0=0 if center is None else center[0],
                         y_0=0 if center is None else center[1], gamma=gamma, alpha=alpha)

        # Set the integral to be normalized
        self.integral = 1.0

        # Set the map sampling and size
        self.sampling = sampling
        self.size = size

        # Set the map if requested
        if sampling is not None or size is not None:
            self.make_map()

    def get_integral(self):
        """Return the analytic integral of the source."""
        return self.amplitude * numpy.pi * self.gamma**2 / (self.alpha - 1)

    def minimum_sampling(self):
        r"""
        Return the minimum sampling in arcseconds per pixels.  Currently
        :math:`{\rm FWHM}/2`.
        """
        return self._fwhm/2.

    def minimum_size(self):
        r"""
        The minimum size that should be used for the distribution map in
        arcseconds.  Currently :math:`2\ {\rm FWHM}`.
        """
        return self._fwhm*2.


class OnSkySersic(functional_models.Sersic2D, Source):
    """
    An on-sky Sersic distribution.

    Args:
        sb_eff (scalar-like):
            The surface brightness at 1 effective (half-light) radius.
        r_eff (scalar-like):
            The effective (half-light) radius in *arcseconds*.
        n (scalar-like):
            The Sersic index.
        center (scalar-like, optional):
            The coordinates of the Sersic center in *arcseconds*
            relative to the image center.
        ellipticity (scalar-like, optional):
            The ellipticity (1-b/a) of an elliptical Sersic
            distribution.
        position_angle (scalar-like, optional):
            The position angle for the elliptical Sersic distribution,
            defined as the angle from N through E.  The coordinate
            system is defined with positive offsets (in RA) toward the
            east, meaning lower pixel indices.
        sampling (scalar-like, optional):
            Sampling of a generated map in arcseconds per pixel.
            Default is set by :func:`minimum_sampling`.
        size (scalar-like, optional):
            Size of the image to generate of the distribution in
            *arceconds* along one of the axes.  The map is square.
            Default is defined by :func:`minimum_size`.
        unity_integral (:obj:`bool`, optional):
            Renormalize the distribution so that the integral is unity.
    """
    def __init__(self, sb_eff, r_eff, n, center=None, ellipticity=None, position_angle=None,
                 sampling=None, size=None, unity_integral=False):

        self.position_angle = 0 if position_angle is None else position_angle
        super(OnSkySersic, self).__init__(amplitude=sb_eff, r_eff=r_eff, n=n,
                                          x_0=0 if center is None else center[0],
                                          y_0=0 if center is None else center[1],
                                          ellip=ellipticity,
                                          theta=numpy.radians(90-self.position_angle))

        self.semi = SemiMajorAxisCoo(xc=self.x_0.value, yc=self.y_0.value, pa=self.position_angle,
                                     ell=self.ellip.value)

        self.bn = None
        self.integral = self.get_integral()
        
        if unity_integral:
            self.amplitude /= self.integral
            self.integral = self.get_integral()

        # Set the map sampling and size
        self.sampling = sampling
        self.size = size

        # Set the map if requested
        if sampling is not None or size is not None:
            self.make_map()

    def get_integral(self):
        """
        The analytic integral of the Sersic profile projected on the
        sky.
        """
        # Note the (1-ellipticity) factor.
        self.bn = special.gammaincinv(2. * self.n, 0.5)
        return 2 * numpy.pi * self.n * numpy.exp(self.bn) * self.amplitude \
                            * numpy.square(self.r_eff) * (1-self.ellip) \
                            * special.gamma(2*self.n) * numpy.power(self.bn, -2*self.n)

    def minimum_sampling(self):
        r"""
        Return the minimum sampling in arcseconds per pixels.  Currently
        :math:`R_{\rm eff}/3`.
        """
        return self.r_eff/3.

    def minimum_size(self):
        r"""
        The minimum size that should be used for the distribution map in
        arcseconds.  Currently :math:`3\ R_{\rm eff}`.
        """
        return self.r_eff*3


# TODO: Add an input image distribution
#class OnSkyImage:
#    def __init__(self, fitsfile):


class OnSkySource(Source):
    """
    Container class for an on-sky source convolved with the seeing
    disk.

    Unlike the other :class:`Source` objects, this requires a map to
    work.

    Args:
        seeing (:obj:`float`, :class:`Source`):
            The FWHM of a Gaussian seeing distribution in arcseconds
            or an object used to define the seeing kernel directly.
            If a float is provided, the sampling of the Gaussian
            seeing kernel is set by
            :func:`OnSkyGaussian.minimum_sampling` unless adjusted by
            the intrinsic source object or the ``sampling`` keyword.
            If a :class:`Source` object, the object is used to
            generate a map of the source surface brightness
            distribution. The integral of the seeing kernel should be
            unity!
        intrinsic (:obj:`float`, :class:`Source`):
            The intrinsic surface brightness distribution of the
            source. Can be the total flux of a point source (in,
            e.g., 10^-17 erg/s/cm^2/angstrom) or an object. If a
            :class:`Source` object, the object is used to generate a
            map of the source surface-brightness distribution.
        beta (:obj:`float`, optional):
            The beta for the Moffat function.  If None, use a Gaussian seeing
            distribution.  If seeing is a :class:`Source` object, this is
            ignored.
        sampling (scalar-like, optional):
            Sampling of a generated map in arcseconds per pixel.
            Default is set by :func:`minimum_sampling`.
        size (scalar-like, optional):
            Size of the image to generate of the distribution in
            *arceconds* along one of the axes.  The map is square.
            Default is defined by :func:`minimum_size`.
    """
    def __init__(self, seeing, intrinsic, beta=None, sampling=None, size=None):

        if isinstance(seeing, float):
            self.seeing = OnSkyGaussian(seeing) if beta is None else OnSkyMoffat(seeing, beta)
        else:
            self.seeing = seeing

        # The intrinsic source distribution
        self.intrinsic = intrinsic

        # Get the sampling
        self.sampling = self.minimum_sampling() if sampling is None else sampling
        self.size = self.minimum_size() if size is None else size

        # Make the map
        self.interp = None
        self.make_map()

    def minimum_sampling(self):
        r"""
        Return the minimum sampling in arcseconds per pixels.

        This is determined by the minimum of the seeing disk sampling
        (:attr:`seeing`) and the sampling for the intrinsic
        distribution (if the latter is defined).
        """
        # Sampling in arcsec / pixel
        sampling = self.seeing.minimum_sampling()
        try:
            # Try using `intrinsic` as an object
            sampling = min(self.intrinsic.minimum_sampling(), sampling)
        except AttributeError:
            pass
        return sampling

    def minimum_size(self):
        """
        Return the minimum size of the rendered source map in
        arcseconds.

        This is determined by the maximum of the seeing disk map size
        (:attr:`seeing`) and the intrinsic source map size (if the
        latter is defined).
        """
        # Size in arcsec
        size = self.seeing.minimum_size()
        try:
            # Try using `intrinsic` as an object
            size = max(self.intrinsic.minimum_size(), size)
        except AttributeError:
            pass
        return size

    def make_map(self, sampling=None, size=None):
        """
        Generate a square map of the source surface-brightness
        distribution.

        Object is modified internally with the constructed map saved
        to :attr:`data`.

        Args:
            sampling (:obj:`float`, optional):
                The size of each pixel in the map in arcsec per
                pixel. If None, sampling will be set by
                :func:`minimum_sampling`.
            size (:obj:`float`, optional):
                Length of one axis of the *square* map in arcsec.
                Will be adjusted to ensure an integer number of
                pixels in the map based on :attr:`sampling`. If None,
                set by :func:`minimum_size`.
        """
        if sampling is None and self.sampling is None:
            self.sampling = self.minimum_sampling()
        elif sampling is not None:
            self.sampling = sampling
            
        if size is None and self.size is None:
            self.size = self.minimum_size()
        elif size is not None:
            self.size = size

        # Build the on-sky source distribution
        self.seeing.make_map(sampling=self.sampling, size=self.size)
        self.x = self.seeing.x
        self.X = self.seeing.X
        self.y = self.seeing.y
        self.Y = self.seeing.Y
        try:
            # Construct the intrinsic map of the source
            self.intrinsic.make_map(sampling=self.sampling, size=self.size)
            # Convolve with the seeing distribution, conserving the
            # integral of the intrinsic source
            self.data = signal.fftconvolve(self.intrinsic.data,
                                           self.seeing.data * numpy.square(self.sampling),
                                           mode='same')
        except AttributeError:
            # Renormalize the unity-integral seeing kernal for to
            # represent a point source
            self.data = self.intrinsic*self.seeing.data

        # Get the integral
        try:
            # After convolving with the seeing kernel, the total
            # integral should be the same, up to some tolerance
            self.integral = self.intrinsic.integral
            tolerance = 1e-3
            diff = numpy.absolute(self.integral - numpy.square(self.sampling)*numpy.sum(self.data))
            if diff > tolerance:
                warnings.warn('Map and analytic integrals are discrepant by {0} ({1} %)'.format(
                                    diff, 100*diff/self.integral))
        except AttributeError:
            self.integral = numpy.square(self.sampling) * numpy.sum(self.data)

        # Prep for interpolation
        self.interp = interpolate.interp2d(self.x, self.y, self.data, bounds_error=True)
#        self.interp = interpolate.RectBivariateSpline(self.x, self.y, self.data)

    def __call__(self, x, y):
        """
        Sample the source.

        This interpolates the pre-calculated source at the requested
        coordinate. A `ValueError` will be thrown (see
        `scipy.interpolate.interp2d`_) if the coordinate is outside
        the bounds of the calculated map.

        Args:
            x (:obj:`float`):
                The position in arcsec relative to the field center.
                Positive x is toward the East (positive RA, smaller
                pixel number).
            y (:obj:`float`):
                The position in arcsec relative to the field center.
                Positive x is toward the North (larger pixel number).

        Returns:
            :obj:`float`: The surface brightness at (x,y).
        """
        return self.interp(x,y)


