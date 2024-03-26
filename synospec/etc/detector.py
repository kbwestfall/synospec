"""
Detector class

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

from IPython import embed
import numpy

from .efficiency import Efficiency

class Detector(Efficiency):
    r"""
    Define the detector technical properties.

    Args:
        shape (:obj:`tuple`, optional):
            Dimensions of the **unbinned** detector in number of pixels along
            the spectral axis and number of pixels along the spatial axis. Can
            be None, but limits use if it is.
        pixelsize (:obj:`float`, optional):
            The size of the **unbinned** (square) detector pixels in *micron*.
        binning (:obj:`int`, optional):
            The number of unbinned pixels along one axis included in a square
            binned pixel.  I.e., if the binning is 4x4, this should be 4.
            Currently cannot accommodate binned that is not square.
        rn (:obj:`float`, optional):
            Read-noise in electrons.
        dark (:obj:`float`, optional):
            Dark current in electrons per pixel per second.
        cic (:obj:`float`, optional):
            Clock-induced charge in electrons per pixel.
        gain (:obj:`float`, optional):
            Gain of detector amplifier in electrons per ADU.
        fullwell (:obj:`float`, optional):
            The full well of the pixels in electrons.
        nonlinear (:obj:`float`, optional):
            The fraction of the fullwell above which the detector
            response is nonlinear.
        qe (:obj:`float`, :class:`Efficiency`, optional):
            Detector quantum efficiency.
        em_fac (:obj:`float`, optional):
            For EMCCDs, this is the multiplicative noise factor due to the
            stochasticity in the gain at high gain.  The theoretical value for
            this is :math:`\sqrt{2}`.  Average values for EMCCDs is quoted as
            1.3 by `Hamamatsu
            <https://hamamatsu.magnet.fsu.edu/articles/emccds.html>`__.  This
            should be 1 for a normal CCD.
        em_gain (:obj:`float`, optional):
            For EMCCDs, the readnoise can be provided as either the readnoise
            *without* electron multiplication or *with* multiplication.  If
            provided *with* multiplication, this parameter should be set to 1
            and the inclusion of the readnoise in the detector noise is
            identical to a normal CCD.  If provided *without* electron
            multiplication, ``em_gain`` should provide the electron
            multiplication gain setting, and the readnoise term included in the
            error budget is then ``rn/em_gain``.
        fps (:obj:`float`, optional):
            Readout time expressed in number of full-frame readouts per second.
    """
    # TODO: Allow for multiple amplifiers per detector? Would also need
    # to define amplifier section.
    # TODO: Define overscan and data sections
    def __init__(self, shape=None, pixelsize=15., binning=1, rn=1., dark=0., cic=0., gain=1.,
                 fullwell=1e4, nonlinear=1., qe=0.9, em_fac=1., em_gain=1., fps=1.):
        if not isinstance(qe, (Efficiency, float)):
            raise TypeError('Provided quantum efficiency must be type `float` or `Efficiency`.')
        if isinstance(qe, float):
            super(Detector, self).__init__(qe)
        else:
            super(Detector, self).__init__(qe.eta, wave=qe.wave)

        if shape is not None and len(shape) != 2:
            raise ValueError('Shape must contain two integers.')
        self.binning = binning
        # Shape in binned pixels
        self.shape = None if shape is None else numpy.asarray(shape)/self.binning
        # Size of binned pixels in microns
        self.pixelsize = pixelsize*self.binning
        # Size of the detector in both axis in mm
        self.size = None if self.shape is None else self.shape*self.pixelsize*1e-3
        # Readnoise
        self.rn = rn
        # Dark current per binned pixel per second
        self.dark = dark*self.binning**2
        # Clock-induced charge per binned pixel per *exposure*
        self.cic = cic*self.binning**2
        self.gain = gain
        self.fullwell = fullwell
        self.nonlinear = nonlinear
        self.em_fac = em_fac
        self.em_gain = em_gain
        self.fps = fps

    def count_rate(self, wave, photon_rate):
        """
        Apply the quantum efficiency to convert a photon rate to an electron
        count rate.

        Args:
            wave (:obj:`float`, array-like):
                The wavelength for the incident photons.
            photon_rate (:obj:`float`, array-like):
                The photon flux in number per second at the provided wavelength.

        Returns:
            :obj:`float`, `numpy.ndarray`:
                Count rate.
        """        
        single_value = isinstance(photon_rate, float)
        _rate = numpy.atleast_1d(photon_rate)
        _wave = numpy.atleast_1d(wave)
        if _wave.size == 1 and _rate.size > 1:
            _wave = numpy.full(_rate.shape, wave)
        if _wave.shape != _rate.shape:
            raise ValueError('Wavelength and flux array shapes must match.')
        cnts = self(_wave)*_rate
        return cnts[0] if single_value else cnts

    # TODO: Add digitization noise?
    def variance(self, count_rate, exptime=1.):
        """
        Compute the detector variance for a given count rate.

        Args:
            count_rate (:obj:`float`, `numpy.ndarray`_):
                Number of counts per second.  The quantum efficiency should
                already have been accounted for; see :func:`count_rate`.
            exptime (:obj:`float`, optional):
                The exposure time in seconds.

        Returns:
            :obj:`float`, `numpy.ndarray`_: The expected statistical error in
            the detector counts.  Shape (and type) match ``photon_flux``.
        """
        # NOTE: self.dark and self.cic already account for the binning.
        return self.em_fac**2 * (count_rate*exptime + self.dark*exptime + self.cic) \
                    + (self.rn/self.em_gain)**2

    def t_exp(self, n_exp, t_tot):
        """
        Provided a total observing time and a number of exposures, calculate the
        maximum integration time per exposure given the detector readout time.

        Args:
            n_exp (:obj:`int`):
                Number of exposures.
            t_tot (:obj:`float`):
                Total observing allocation in seconds.

        Returns:
            :obj:`float`: Integration time per exposure to fill to total
            observing allocation.
        """
        return t_tot/n_exp - 1/self.fps

    def n_exp(self, t_exp, t_tot):
        """
        Provided a total observing time and an integration time per exposure,
        calculate the maximum number of exposures that be performed given the
        detector readout time.

        Args:
            t_exp (:obj:`int`):
                Integration time per exposure.
            t_tot (:obj:`float`):
                Total observing allocation in seconds.

        Returns:
            :obj:`float`: Number of exposures possible within the total
            observing allocation.
        """
        return int(numpy.round(t_tot/(t_exp + 1/self.fps)))

    def t_tot(self, n_exp, t_exp):
        """
        Get the total observing time including readout overhead.
        """
        return n_exp*(t_exp + 1/self.fps)




