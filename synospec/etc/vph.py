"""
Various efficiency calculations

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import warnings
import inspect

from IPython import embed

import numpy

from matplotlib import pyplot

from scipy import interpolate

from astropy import units

from . import efficiency

class VPHGrating:
    """
    Object for generic computations for a volume-phase holographic
    grating.
    
    Args:
        line_density (:obj:`float`):
            Line density of the grating in lines per mm.
        n_bulk (:obj:`float`):
            Bulk index of refraction.
        n_mod (:obj:`float`):
            Modulation in the index of refraction.
        thickness (:obj:`float`):
            Grating thickness in micron.
        peak_wave (:obj:`float`, optional):
            Wavelength in angstroms of peak efficiency. This can be
            calculated using the properties of the grating, but a
            reasonable set of wavelength limits (see ``wave_lim``)
            are required, both minimum and maximum. If provided
            directly, this is *not* checked against this direct
            calculation.
        wave_lim (:obj:`tuple`, optional):
            Two-tuple with the lower and upper wavelength limits in
            angstroms to allow for any calculations for this grating.
            If None, no limits are imposed. Elements of the tuple can
            be None, meaning the specific limit is undefined; i.e.,
            to define an upper limit only, set
            ``wave_lim=(None,5400)``.
    """
    def __init__(self, line_density, n_bulk, n_mod, thickness, peak_wave=None, wave_lim=None):
        self.line_density = line_density
        self.n_bulk = n_bulk
        self.n_mod = n_mod
        self.thickness = thickness
        self._peak_wave = peak_wave
        self.wave_lim = (None, None) if wave_lim is None else tuple(wave_lim)
        if not isinstance(self.wave_lim, tuple):
            raise TypeError('Input wavelength limits is not (or could not be converted to) '
                            'a tuple.')
        if len(self.wave_lim) != 2:
            raise ValueError('Wavelength limits must be provided as a two-tuple, the upper and '
                             'lower limit.')

    def _check_wavelengths(self, wave):
        """
        Ensure there are valid wavelengths.

        Args:
            wave (:obj:`float`, `numpy.ndarray`_):
                Wavelength in angstroms.
        """
        _wave = numpy.atleast_1d(wave)
        indx = numpy.ones(_wave.size, dtype=bool)
        if self.wave_lim[0] is not None:
            indx &= (numpy.greater(_wave, self.wave_lim[0]) 
                        | numpy.isclose(_wave, self.wave_lim[0]))
        if self.wave_lim[1] is not None:
            indx &= (numpy.less(_wave, self.wave_lim[1]) 
                        | numpy.isclose(_wave, self.wave_lim[1]))
        if numpy.sum(indx) == 0:
            raise ValueError('Wavelength(s) provided are not in valid range for grating.')
        return _wave, indx

    def _return(self, wave, value, indx):
        """
        For calculations involving wavelength, this is a convenience
        function that does the type jiggering.
        """
        if isinstance(wave, (int, numpy.integer, float, numpy.floating)):
            return value[0]
        if numpy.sum(indx) == indx.size:
            return value
        _value = numpy.ma.MaskedArray(numpy.zeros(indx.size, dtype=float),
                                      mask=numpy.invert(indx))
        _value[indx] = value
        return _value

    @property
    def peak_wave(self):
        """
        The wavelength with the peak super-blaze efficiency.

        Computation is approximate based on sampling the super-blaze
        and picking the sample with the highest efficiency.
        """
        if self._peak_wave is not None:
            return self._peak_wave
        if self.wave_lim is None or self.wave_lim[0] is None or self.wave_lim[1] is None:
            raise ValueError('To calculate wavelength at peak efficiency, wavelength limits '
                             'must be provided.  Use wave_lim attribute.')
        wave = numpy.linspace(self.wave_lim[0], self.wave_lim[1], 1000)
        eff = self.efficiency(wave)
        indx = numpy.argmax(eff)
        self._peak_wave = wave[indx]
        return self._peak_wave

    def littrow_wavelength(self, alpha, m=1):
        """
        Compute the Littrow wavelength at the provided incidence
        angle.

        Args:
            alpha (:obj:`float`):
                Incidence angle in degrees.
            m (:obj:`int`, optional):
                Grating order.

        Returns:
            :obj:`float`: Wavelength in angstroms.
        """
        # Factor of 1e7 converts from mm to angstrom
        return 2e7*numpy.sin(numpy.radians(alpha))/m/self.line_density

    def littrow_angle(self, wave, m=1):
        """
        Calculate the Littrow indicence angle for the provided
        wavelength.

        Args:
            wave (:obj:`float`):
                Wavelength in angstroms.
            m (:obj:`int`, optional):
                Grating order.

        Returns:
            :obj:`float`: Input/output angle in degrees.
        """
        # Factor of 1e7 converts from angstrom to mm
        return numpy.degrees(numpy.arcsin(wave * m * self.line_density / 2e7))

    def bragg_angle_dcg(self, wave, sine_of=False):
        """
        Calculate the grating Bragg angle in degrees in the DCG
        layer.

        Based on IDL code provided Jason Fucik (CalTech), 7 Apr 2020.

        Args:
            wave (:obj:`float`, `numpy.ndarray`_):
                Wavelength at which to calculate the Bragg angle
                (angstroms in air).
            sine_of (:obj:`bool`, optional):
                Return the sine of the angle instead of the angle
                itself.

        Returns:
            :obj:`float`, `numpy.ndarray`, `numpy.ma.MaskedArray`:
            The (sine of) the Bragg angle.
        """
        # NOTE: This will fault if wave was a float and outside the
        # valid wavelength range.
        _wave, indx = self._check_wavelengths(wave)
        # Calculate the (sine of the) Bragg angle, where the angle of
        # incidence = the angle of diffraction, in the DCG layer. The
        # factor of 1e-7 converts the wavelength from angstroms to mm.
        sin_bragg = 1e-7*_wave[indx] * self.line_density / 2 / self.n_bulk
        out_bragg = sin_bragg if sine_of else numpy.degrees(numpy.arcsin(sin_bragg))
        return self._return(wave, out_bragg, indx)

    def bragg_angle(self, wave, radians=False):
        """
        Compute the Bragg incidence angle for the grating.

        Based on IDL code provided Jason Fucik (CalTech), 7 Apr 2020.

        Args:
            wave (:obj:`float`, `numpy.ndarray`_):
                Wavelength at which to calculate the Bragg angle
                (angstroms in air).
            radians (:obj:`bool`, optional):
                Provide the angle in radians, instead of degrees.

        Returns:
            :obj:`float`, `numpy.ndarray`_, `numpy.ma.MaskedArray`_:
            The super-blaze incidence angle.
        """
        # Return the angle of incidence using Snell's Law in air
        alpha = numpy.arcsin(self.n_bulk * self.bragg_angle_dcg(wave, sine_of=True))
        return alpha if radians else numpy.degrees(alpha)

    def efficiency(self, wave, tilt=0.0, alpha=None):
        """
        Calculate the efficiency at the provided wavelengths.

        Based on IDL code provided Jason Fucik (CalTech), 7 Apr 2020.

        Args:
            wave (:obj:`float`, `numpy.ndarray`_):
                Wavelengths in angstroms at which to calculate the
                efficiency.
            tilt (:obj:`float`, optional):
                Grating tilt in degrees.  Only used if 
            alpha (:obj:`float`, optional):
                Angle of incidence onto the grating. If None, return
                the super-blaze function.

        Returns:
            :obj:`float`, `numpy.ndarray`_: The grating efficiency.
            The type returned matches the type provided for ``wave``.
        """
        # TODO: Perform type checking.  Allow list objects for wave.
        if alpha is None:
            # Return the super-blaze function
            bragg = numpy.radians(self.bragg_angle_dcg(wave))
            # Factors of 1e4 are to convert thickness from micron to angstrom
            s = numpy.square(numpy.sin(numpy.pi * self.n_mod * self.thickness * 1e4
                                       / wave / numpy.cos(bragg))) / 2
            p = numpy.square(numpy.sin(numpy.pi * self.n_mod * self.thickness * 1e4
                                       * numpy.cos(2*bragg) / wave / numpy.cos(bragg))) / 2
            return s+p

        # NOTE: This will fault if wave was a float and outside the
        # valid wavelength range.
        _wave, indx = self._check_wavelengths(wave)

        # Return the Bragg envelopes
        _alpha = numpy.arcsin(numpy.sin(numpy.radians(alpha))/self.n_bulk)  # AOI at the DCG layer
        cosa = numpy.cos(_alpha)
        psi = numpy.radians(tilt) + numpy.pi/2.
        tpnbulk = 2 * numpy.pi * self.n_bulk

        # Convert the units of line density from lines/mm to lines/angstrom
        _line_density = self.line_density / 1e7
        # Convert the units of thickness from micron to angstrom
        _thickness = self.thickness * 1e4

        k = 2 * numpy.pi * _line_density
        ct = cosa - _wave[indx] * k * numpy.cos(psi) / tpnbulk
        eta = _thickness * (k * numpy.cos(_alpha - psi) - _wave[indx] * k*k / 2 / tpnbulk) / ct / 2
        ks = numpy.pi * self.n_mod / _wave[indx]
        vs = ks * _thickness / numpy.sqrt(ct * cosa)
        vp = -ks * numpy.cos(2*(_alpha - psi)) * _thickness / numpy.sqrt(ct * cosa)
        s = numpy.square(numpy.sin(numpy.sqrt(vs*vs + eta*eta))) / (1 + numpy.square(eta/vs)) / 2
        p = numpy.square(numpy.sin(numpy.sqrt(vp*vp + eta*eta))) / (1 + numpy.square(eta/vp)) / 2
        return self._return(wave, s+p, indx)

    def diffraction_angle(self, wave, alpha=None, m=1, radians=False):
        """
        Compute the diffraction angle for the provided wavelengths.

        Args:
            wave (:obj:`float`, `numpy.ndarray`_):
                Wavelengths for the computation.
            alpha (:obj:`float`, optional):
                The angle of incidence on the grating. If None, the
                angle is the Bragg angle **for each wavelength**;
                ie., you should basically always be providing
                ``alpha``.
            m (:obj:`int`, optional):
                Grating order.
            radians (:obj:`bool`, optional):
                Provide the angle in radians, instead of degrees.

        Returns:
            :obj:`float`, `numpy.ndarray`_, `numpy.ma.MaskedArray`_:
            The angle of diffraction for each wavelength.
        """
        _wave, indx = self._check_wavelengths(wave)
        _alpha = self.bragg_angle(_wave[indx], radians=True) \
                    if alpha is None else numpy.radians(alpha)
        beta = numpy.arcsin(m * _wave[indx] * self.line_density/1e7 - numpy.sin(_alpha))
        return self._return(wave, beta if radians else numpy.degrees(beta), indx)

    def diffracted_wave(self, beta, alpha=None, m=1):
        """
        Compute the wavelength at the provided angle of diffraction.

        Args:
            beta (:obj:`float`):
                Angle of diffraction.
            alpha (:obj:`float`, optional):
                The angle of incidence on the grating. If None,
                assumed to be the same as the diffraction angle and
                the Littrow wavelength is returned; see
                :func:`littrow_wavelength`.
            m (:obj:`int`, optional):
                Grating order.

        Returns:
            :obj:`float`: Wavelength in angstroms at the provided
            angle of diffraction.
        """
        if alpha is None:
            return self.littrow_wavelength(beta, m=m)
        return (numpy.sin(numpy.radians(alpha)) + numpy.sin(numpy.radians(beta))) * 1e7 \
                    / m / self.line_density

    def angular_dispersion(self, wave, alpha=None, m=1):
        """
        Return the angular dispersion in radians per angstrom.

        Args:
            wave (:obj:`float`, `numpy.ndarray`_):
                Wavelengths for the computation.
            alpha (:obj:`float`, optional):
                The angle of incidence on the grating. If None, the
                angle is the Bragg angle **for each wavelength**;
                ie., you should basically always be providing
                ``alpha``.
            m (:obj:`int`, optional):
                Grating order.

        Returns:
            :obj:`float`, `numpy.ndarray`_, `numpy.ma.MaskedArray`_:
            The angular dispersion in radians per angstrom.
        """
        _wave, indx = self._check_wavelengths(wave)
        beta = self.diffraction_angle(_wave[indx], alpha=alpha, m=m, radians=True)
        return self._return(wave, m * self.line_density / 1e7 / numpy.cos(beta), indx)



