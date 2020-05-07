# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Implements a set of line profile parameterizations.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import numpy
from scipy import special

class GaussianLSF:
    r"""
    Define a Gaussian line profile, *sampled* over the width of the
    sampling step, parameterized by its integral (:math:`F`), center
    (:math:`\mu`), and standard deviation (:math:`\sigma`). I.e:

    .. math::

        \mathcal{N}(x|f,\mu,\sigma) = \frac{f}{\sqrt{2\pi}\sigma}
        \exp\left(\frac{-\Delta^2}{2\sigma^2}\right)

    where :math:`\Delta = x-\mu`.  The coordinate vector :math:`x` does
    not need to be uniformly sampled.

    Args:
        p (array-like, optional):
            Input parameters ordered as the total integral of the
            profile, the profile center, and the profile standard
            deviation.  Assumed to be (1.0, 0.0, 1.0) by default.

    Attributes:
        p (numpy.ndarray):
            Most recently used parameters
    """
    def __init__(self, p=None):
        self.set_par(p)

    def __call__(self, x, p):
        """Calculate the profile.

        Args:
            x (array-like):
                Independent variable.
            p (array-like):
                LSF parameters.
        """
        self.set_par(p)
        return self.sample(x)

    @staticmethod
    def npar():
        return 3

    def set_par(self, p):
        """
        Set the internal parameters to the provided set.

        Args:
            p (array-like):
                LSF parameters.

        Raises:
            ValueError:
                Raised if the provided parameter vector is not 3
                elements long.
        """
        if p is None:
            self.p = numpy.array([1.0, 0.0, 1.0])
            return
        if len(p) != GaussianLSF.npar():
            raise ValueError('Must provide {0} parameters.'.format(GaussianLSF.npar()))
        self.p = numpy.asarray(p)

    def sample(self, x):
        """
        Sample the profile.

        Args:
            x (array-like):
                Independent variable.
        """
        return self.p[0] * numpy.exp(-numpy.square((x-self.p[1])/self.p[2])/2.) \
                    / numpy.sqrt(2.0*numpy.pi) / self.p[2]

    def parameters_from_moments(self, mom0, mom1, mom2):
        """
        Provided the 0th, 1st, and 2nd moments, produce a set of
        parameters for the profile.
        """
        return numpy.array([mom0, mom1, mom2])


class IntegratedGaussianLSF(GaussianLSF):
    r"""
    Define a Gaussian line profile, *integrated* over the width of the
    sampling step, parameterized by its integral (:math:`F`), center
    (:math:`\mu`), and standard deviation (:math:`\sigma`). I.e:

    .. math::

        \mathcal{N}(x|F,\mu,\sigma) = \frac{F}{2} \left[
        {\rm erf}\left(\frac{\Delta+\delta_x/2}{\sqrt{2}\sigma}\right) - 
        {\rm erf}\left(\frac{\Delta-\delta_x/2}{\sqrt{2}\sigma}\right)\right]

    where :math:`{\rm erf}(x)` is the error function, :math:`\Delta =
    x-\mu`, and :math:`\delta_x` is the sampling step.  The sampling
    *must* be uniform in :math:`x`.

    Args:
        p (array-like, optional): 
            Input parameters ordered as the total integral of the
            profile, the profile center, and the profile standard
            deviation.  Assumed to be (1.0, 0.0, 1.0) by default.
        dx (:obj:`float`, optional):
            Sampling width.  Default is 1.

    Attributes:
        p (numpy.ndarray):
            Most recently used parameters
        dx (:obj:`float`):
            Assumed sampling.
    """
    def __init__(self, p=None, dx=None):
        self.set_par(p)
        self.dx = 1.0 if dx is None else dx

    def sample(self, x):
        """
        Sample the profile.

        .. warning::
            Does **not** check if the provided :math:`x` values are
            sampled at :attr:`dx`.

        Args:
            x (array-like):
                Independent variable.
        """
        n = numpy.sqrt(2.)*self.p[2]
        d = numpy.asarray(x)-self.p[1]
        return self.p[0] * (special.erf((d+self.dx/2.)/n) - special.erf((d-self.dx/2.)/n))/2.

   
class FFTGaussianLSF(GaussianLSF):
    r"""
    Define a Gaussian line profile by first constructing the analytic
    FFT of the profile and then returning the inverse real FFT.  See
    ppxf_util.emline by M. Cappellari.  The sampling *must* be uniform
    in :math:`x`.
    
    Args:
        p (array-like, optional):
            Input parameters ordered as the total integral of the
            profile, the profile center, and the profile standard
            deviation.  Assumed to be (1.0, 0.0, 1.0) by default.
        dx (:obj:`float`, optional):
            Sampling width.  Default is 1.
        pixel (:obj:`bool`, optional):
            Flag to produce profile integrated over the sampling width.
        
    Attributes:
        p (numpy.ndarray):
            Most recently used parameters
        dx (:obj:`float`):
            Assumed sampling.
        pixel (:obj:`bool`):
            Flag to produce profile integrated over the sampling width.
    """
    def __init__(self, p=None, dx=None, pixel=True):
        self.set_par(p)
        self.dx = 1.0 if dx is None else dx
        self.pixel = pixel

    def sample(self, x):
        """
        Sample the profile.

        .. warning::
            Does **not** check if the provided :math:`x` values are
            sampled at :attr:`dx`.

        Args:
            x (array-like):
                Independent variable.
        """
        xsig = self.p[2]/self.dx
        x0 = (self.p[1]-x[0])/self.dx
        npad = 2**int(numpy.ceil(numpy.log2(x.size)))
        w = numpy.linspace(0,numpy.pi,npad//2+1)
        rfft = self.p[0]*numpy.exp(-0.5*numpy.square(w*xsig) - 1j*w*x0)
        if self.pixel:
            rfft *= numpy.sinc(w/(2*numpy.pi))
        lsf = numpy.fft.irfft(rfft, n=npad)[:x.size]
        return lsf if self.pixel else lsf/self.dx


