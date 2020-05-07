"""
Define a series of convolution kernels.

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

from astropy.convolution import Gaussian2DKernel

# TODO: Need to abstract this to any kind of kernel function, as well
# as something generic based on a spot diagram from Zemacs

class SpectrographGaussianKernel(Gaussian2DKernel):
    """
    Define a Gaussian kernel describing the image quality of a
    spectrograph.

    Args:
        spatial_fwhm (scalar-like):
            Spatial FWHM of the Gaussian (point-spread function) in
            *mm*. Should be *independent* of the focal-plane
            entrance aperture. The spatial dimension is assumed to be
            along the ordinate of any map convolved with this kernel.
        spectral_fwhm (scalar-like):
            Spectral FWHM of the Gaussian (line-spread function) in
            *mm*. Should be *independent* of the focal-plane
            entrance aperture. The spectral dimension is assumed to
            be along the abscissa of any map convolved with this
            kernel.
        pixelscale (:obj:`float`, optional):
            The pixel scale for the kernel in mm/pixel.

    Attributes:
        spatial_fwhm (:obj:`float`):
            See class argument.
        spectral_fwhm (:obj:`float`):
            See class argument.
        pixelscale (:obj:`float`):
            See class argument.
    """
    def __init__(self, spatial_fwhm, spectral_fwhm, pixelscale=1.):
        self._sig2fwhm = numpy.sqrt(8*numpy.log(2))
        self.spatial_fwhm = float(spatial_fwhm)
        self.spectral_fwhm = float(spectral_fwhm)
        self.pixelscale = float(pixelscale)

        # Instantiate the functional_models.Gaussian2D object
        super(SpectrographGaussianKernel, self).__init__(
                self.spectral_fwhm/self._sig2fwhm/self.pixelscale,
                y_stddev=self.spatial_fwhm/self._sig2fwhm/self.pixelscale)

    def resample(self, pixelscale):
        """
        Resample the kernel to a new pixel scale.
        """
        self.pixelscale = pixelscale

        # Instantiate the functional_models.Gaussian2D object
        super(SpectrographGaussianKernel, self).__init__(
                self.spectral_fwhm/self._sig2fwhm/self.pixelscale,
                y_stddev=self.spatial_fwhm/self._sig2fwhm/self.pixelscale)
