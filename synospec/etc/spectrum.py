"""
Spectrum utilities

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import warnings

from IPython import embed

import numpy
from scipy import interpolate
from astropy.io import fits
from astropy.wcs import WCS
import astropy.constants
import astropy.units

from matplotlib import pyplot

from pydl.goddard.astro import airtovac

from ..util.lineprofiles import IntegratedGaussianLSF
from .sampling import Resample
from .resolution import match_spectral_resolution

from .. import data_file

def spectral_coordinate_step(wave, log=False, base=10.0):
    """
    Return the sampling step for the input wavelength vector.

    If the sampling is logarithmic, return the change in the logarithm
    of the wavelength; otherwise, return the linear step in angstroms.

    Args: 
        wave (array-like):
            Wavelength coordinates of each spectral channel in
            angstroms.
        log (:obj:`bool`, optional):
            Input spectrum has been sampled geometrically.
        base (:obj:`float`, optional):
            If sampled geometrically, the sampling is done using a
            logarithm with this base.  For natural logarithm, use
            numpy.exp(1).

    Returns:
        :obj:`float`: Spectral sampling step in either angstroms
        (`log=False`) or the step in log(angstroms).

    Raises:
        ValueError:
            Raised if the wavelength vector is not uniformly (either
            linearly or log-linearly) sampled to numerical accuracy.
    """
    dw = numpy.diff(numpy.log(wave))/numpy.log(base) if log else numpy.diff(wave)
    if numpy.any(numpy.absolute(numpy.diff(dw)) > 1e-11): #00*numpy.finfo(dw.dtype).eps):
        raise ValueError('Wavelength vector is not uniformly sampled to numerical accuracy.')
    return numpy.mean(dw)


def angstroms_per_pixel(wave, log=False, base=10.0, regular=True):
    """
    Return a vector with the angstroms per pixel at each channel.

    When `regular=True`, the function assumes that the wavelengths are
    either sampled linearly or geometrically.  Otherwise, it calculates
    the size of each pixel as the difference between the wavelength
    coordinates.  The first and last pixels are assumed to have a width
    as determined by assuming the coordinate is at its center.

    .. note::

        If the regular is False and log is True, the code does *not*
        assume the wavelength coordinates are at the geometric center of
        the pixel.

    Args:
        wave (`numpy.ndarray`_):
            (Geometric) centers of the spectrum pixels in angstroms.
        log (`numpy.ndarray`_, optional):
            The vector is geometrically sampled.
        base (:obj:`float`, optional):
            Base of the logarithm used in the geometric sampling.
        regular (:obj:`bool`, optional):
            The vector is regularly sampled.

    Returns:
        numpy.ndarray: The angstroms per pixel.
    """
    if regular:
        ang_per_pix = spectral_coordinate_step(wave, log=log, base=base)
        return ang_per_pix*wave*numpy.log(base) if log else numpy.repeat(ang_per_pix, len(wave))

    return numpy.diff([(3*wave[0]-wave[1])/2] + ((wave[1:] + wave[:-1])/2).tolist()
                      + [(3*wave[-1]-wave[-2])/2])


def convert_flux_density(wave, flux, error=None, density='ang'):
    r"""
    Convert a spectrum with flux per unit wavelength to per unit
    frequency or vice versa.

    For converting from per unit wavelength, this function returns
    
    .. math::
        
        F_{\nu} = F_{\lambda} \frac{d\lambda}{d\nu} = F_{\lambda}
        \frac{\lambda^2}{c}.

    The spectrum independent variable (`wave`) is always expected to
    be the wavelength in angstroms. The input/output units always
    expect :math:`F_{\lambda}` in :math:`10^{-17}\ {\rm erg\ s}^{-1}\
    {\rm cm}^{-2}\ {\rm A}^{-1}` and :math:`F_{\nu}` in microjanskys
    (:math:`10^{-29} {\rm erg\ s}^{-1}\ {\rm cm}^{-2}\ {\rm
    Hz}^{-1}`). Beyond this, the function is ignorant of the
    input/output units. For example, if you provide the function with
    an input spectrum with :math:`F_{\lambda}` in :math:`10^{-11}\
    {\rm erg\ s}^{-1}\ {\rm cm}^{-2}\ {\rm A}^{-1}`, the output will
    be :math:`F_{\nu}` in Janskys.

    Args:
        wave (:obj:`float`, array-like):
            The vector with the wavelengths in angstroms.
        flux (:obj:`float`, array-like):
            The vector with the flux density; cf. `density`.
        error (:obj:`float`, array-like, optional):
            The error in the flux measurements. If None, no errors
            are returned.
        density (:obj:`str`, optional):
            The density unit of the *input* spectrum.  Must be either
            'ang' or 'Hz'.  If the input spectrum is :math:`F_{\lambda}`
            (`density='ang'`), the returned spectrum is :math:`F_{\nu}`.

    Returns:
        :obj:`float`, `numpy.ndarray`_, :obj:`tuple`: The flux with
        the converted units. If the spectrum errors are not provided,
        only the flux value or array is returned.

    Raises:
        ValueError:
            Raised if the `wave` and `flux` arguments do not have the
            same shape.
    """
    # Set to be at least vectors
    _wave = numpy.atleast_1d(wave)
    _flux = numpy.atleast_1d(flux)
    if _wave.shape != _flux.shape:
        raise ValueError('Wavelength and flux arrays must have the same shape.')
    if error is not None:
        _error = numpy.atleast_1d(error)
        if _error.shape != _flux.shape:
            raise ValueError('Error and flux arrays must have the same shape.')
    if density == 'ang':
        # Convert Flambda to Fnu
        factor = numpy.square(_wave)*1e12/astropy.constants.c.to('angstrom/s').value
        fnu = _flux*factor
        if error is not None:
            fnu_err = _error*factor
            return (fnu[0], fnu_err[0]) if isinstance(flux, float) else (fnu, fnu_err)
        return fnu[0] if isinstance(flux, float) else fnu
    if density == 'Hz':
        # Convert Fnu to Flambda
        factor = astropy.constants.c.to('angstrom/s').value/numpy.square(_wave)/1e12
        flambda = _flux*factor
        if error is not None:
            flambda_err = _error*factor
            return (flambda[0], flambda_err[0]) if isinstance(flux, float) \
                        else (flambda, flambda_err)
        return flambda[0] if isinstance(flux, float) else flambda
    raise ValueError('Density units must be either \'ang\' or \'Hz\'.')


def convert_flux_units(wave, flux, current_units, error=None):
    r"""
    Convert flux density units to ``'1e-17 erg / (cm2 s angstrom)'``.

    The flux density can be per unit frequency (:math:`F_{\nu}`) or
    per unit wavelength (:math:`F_{\lambda}`). If a basic unit
    conversion from the current units to ``'1e-17 erg / (cm2 s
    angstrom)'`` fails, it's assumed that the user provided
    :math:`F_{\nu}`. If the conversion from the current units to
    ``'1e-29 erg / (cm2 Hz s)'`` also fails, the method faults.

    Args:
        wave (`numpy.ndarray`_):
            Wavelength in angstroms.
        flux (`numpy.ndarray`_):
            Flux density measurements.
        current_units (:obj:`str`):
            Input units of the flux density. Must be interpretable by
            `astropy.units.Unit`_.
        error (`numpy.ndarray`_):
            Error vector.  Ignored if None.

    Returns:
        `numpy.ndarray`_, :obj:`tuple`: The converted flux vector. If
        an error vector is provided, it is also returned as the
        second element.
    """
    try:
        # Get unit-conversion factor needed to convert to 1e-17
        # erg/s/cm^2/angstrom. First assume the fluxes are in
        # F_lambda
        conversion = astropy.units.Unit(current_units).to('1e-17 erg / (cm2 s angstrom)')
    except astropy.units.core.UnitConversionError as e:
        # Assume this unit conversion failure is because the flux is
        # F_nu, not F_lambda.
        try:
            # Convert F_nu units to those expected by convert_flux_density
            flux *= astropy.units.Unit(current_units).to('1e-29 erg / (cm2 Hz s)')
        except:
            # If that failed, we're really hosed.
            raise ValueError('Unable to handle flux units') from e
        return convert_flux_density(wave, flux, error=error, density='Hz')
    return flux*conversion if error is None else (flux*conversion, error*conversion)


class Spectrum:
    r"""
    Define a spectrum.

    Units are expected to be wavelength in angstroms and flux in 1e-17
    erg/s/cm^2/angstrom.

    If wavelengths are provided in air (instantiated with
    ``airwave=True``), wavelengths are converted to vacuum and the
    sampling is set to be irregular.

    .. todo::
        - allow for different units?

    Args:
        wave (array-like):
            1D wavelength data in angstroms. Can be sampled
            irregularly.
        flux (array-like):
            1D flux data in 1e-17 erg/s/cm^2/angstrom. Use
            :func:`convert_flux_units`, if necessary.
        error (array-like, optional):
            1-sigma error in the flux data
        mask (array-like, optional):
            Boolean mask for flux data (True=bad)
        resolution (:obj:`float`, array-like, optional):
            1D spectral resolution (:math:`$R=\lambda/\Delta\lambda$`)
        regular (:obj:`bool`, optional):
            Spectrum is regularly sampled, either linearly or log-linearly
        log (:obj:`bool`, optional):
            Spectrum is sampled in steps of log base 10. If regular
            is False, this is ignored.
        airwave (:obj:`bool`, optional):
            Spectrum wavelengths are provided in air. If True, the
            wavelengths will be converted to vacuum, and the
            wavelength step will be set to be irregular (if it isn't
            already).
        check (:obj:`bool`, optional):
            Check sampling to confirm input.
        use_sampling_assessments (:obj:`bool`, optional):
            Override any values provided for ``regular`` and ``log``
            and assess the spectral sampling directly from the
            provided wavelength array.
    """
    def __init__(self, wave, flux, error=None, mask=None, resolution=None, regular=True,
                 log=False, airwave=False, check=True, use_sampling_assessments=False):
        # Check the input
        _wave = numpy.atleast_1d(wave)
        if airwave:
            _wave = airtovac(_wave)
        _flux = numpy.atleast_1d(flux)
        if _wave.ndim != 1:
            raise ValueError('Spectrum can only accommodate single vectors for now.')
        if _wave.shape != _flux.shape:
            raise ValueError('Wavelength and flux vectors do not match.')

        # Sampling
        if use_sampling_assessments:
            self.regular, self.log = Spectrum.assess_sampling(_wave)
        else:
            self.regular = False if airwave else regular
            self.log = log
        if check:
            found_regular, found_log = Spectrum.assess_sampling(_wave)
            if self.regular != found_regular:
                raise ValueError('Expected {0} sampling, found {1} sampling.'.format(
                                    'regular' if self.regular else 'irregular',
                                    'regular' if found_regular else 'irregular'))
            if self.log != found_log:
                raise ValueError('Geometric sampling {0} found, but {1} expected.'.format(
                                    'was' if found_log else 'was not',
                                    'was' if self.log else 'was not'))

        self.sres = None if resolution is None else numpy.atleast_1d(resolution)
        if self.sres is not None:
            # Allow resolution to be a single value
            if self.sres.size == 1:
                self.sres = numpy.repeat(self.sres, _flux.size)
            if self.sres.shape != _flux.shape:
                raise ValueError('Resolution vector must be a single number or a vector that '
                                 'matches the length of flux vector.')
            
        self.error = None if error is None else numpy.atleast_1d(error)
        if self.error is not None and self.error.shape != _flux.shape:
            raise ValueError('Error vector must match length of flux vector.')
            
        self.mask = None if mask is None else numpy.atleast_1d(mask)
        if self.mask is not None and self.mask.shape != _flux.shape:
            raise ValueError('Mask vector must match length of flux vector.')

        # TODO: Need interpolators for the error and mask?

        self.interpolator = interpolate.interp1d(_wave, _flux, assume_sorted=True)
        self.size = _wave.size
        # TODO: Check log against the input wavelength vector
        self.nu = self._frequency()
        self.fnu = None

        # Do this brute force Vega magnitude calculations for the time
        # being
        self._vega = None
        
    @staticmethod
    def assess_sampling(wave):
        """
        Assess the wavelength sampling directly from the wavelength vector.

        Args:
            wave (array-like):
                Array with the wavelengths.

        Returns:
            :obj:`tuple`: Two booleans are returned. The first
            indicates if the wavelength vector is regularly sampled
            (either linearly or log-linearly), the second indicates
            if the sampling is log-linear.
        """
        # Use absolute to allow for a monotonically decreasing wavelength vector
        ddw = numpy.absolute(numpy.diff(numpy.diff(wave)))
        regular = numpy.all(ddw < 1e-11) #100*numpy.finfo(ddw.dtype).eps)
        if regular:
            return regular, False
        ddw = numpy.absolute(numpy.diff(numpy.diff(numpy.log(wave))))
        regular = numpy.all(ddw < 1e-11) #100*numpy.finfo(ddw.dtype).eps)
        if regular:
            return regular, True
        return False, False

    @property
    def wave(self):
        """
        The wavelength data vector.
        """
        return self.interpolator.x

    @property
    def flux(self):
        """
        The flux data vector.
        """
        return self.interpolator.y

    def copy(self):
        """
        Return a deep copy of the spectrum.
        """
        return Spectrum(self.wave.copy(), self.flux.copy(),
                        error=None if self.error is None else self.error.copy(),
                        mask=None if self.mask is None else self.mask.copy(),
                        resolution=None if self.sres is None else self.sres.copy(),
                        regular=self.regular, log=self.log, check=False)

    def _arith(self, other, func):
        if isinstance(other, Spectrum):
            return func(other)

        _other = numpy.atleast_1d(other)
        if _other.size == 1:
            _other = numpy.repeat(_other, self.size)
        if _other.size != self.size:
            raise ValueError('Vector to add has incorrect length.')
        return func(Spectrum(self.interpolator.x, _other))

    def __add__(self, rhs):
        return self._arith(rhs, self._add_spectrum)

    def __radd__(self, lhs):
        return self + lhs

    def __sub__(self, rhs):
        return self + rhs*-1

    def __rsub__(self, lhs):
        return self*-1 + lhs

    def __mul__(self, rhs):
        return self._arith(rhs, self._mul_spectrum)

    def __rmul__(self, lhs):
        return self * lhs

    def __truediv__(self, rhs):
        if isinstance(rhs, Spectrum):
            return self * rhs.inverse()
        # NOTE: Parentheses are important below because they force the
        # call to __mul__ for the correct quantity and avoid an
        # infinite loop!
        return self * (1/float(rhs))

    def __rtruediv__(self, lhs):
        return self.inverse() * lhs

    def inverse(self):
        error = None if self.error is None else self.error/numpy.square(self.interpolator.y)
        return Spectrum(self.wave, 1./self.interpolator.y, error=error,
                        mask=None if self.mask is None else self.mask.copy(), resolution=self.sres,
                        regular=self.regular, log=self.log)

    def _mul_spectrum(self, rhs):
        """
        rhs must have type Spectrum
        """
        if not numpy.array_equal(rhs.wave, self.wave):
            raise NotImplementedError('To perform arithmetic on spectra, their wavelength arrays '
                                      'must be identical.')
        flux = self.flux * rhs.flux

        if self.error is None and rhs.error is None:
            error = None
        else:
            error = numpy.zeros(self.size, dtype=float)
            if self.error is not None:
                error += numpy.square(self.error/self.flux)
            if rhs.error is not None:
                error += numpy.square(rhs.error/rhs.flux)
            error = numpy.sqrt(error) * numpy.absolute(flux)

        if self.mask is None and rhs.mask is None:
            mask = None
        else:
            mask = numpy.zeros(self.size, dtype=bool)
            if self.mask is not None:
                mask |= self.mask
            if rhs.mask is not None:
                mask |= rhs.mask
            
        if self.sres is not None and rhs.sres is not None \
                and not numpy.array_equal(self.sres, rhs.sres):
            warnings.warn('Spectral resolution is not correctly propagated.')
        sres = self.sres if self.sres is not None else rhs.sres
        return Spectrum(self.wave, flux, error=error, mask=mask, resolution=sres, regular=regular,
                        log=self.log)

    def _add_spectrum(self, rhs):
        """
        rhs must have type Spectrum
        """
        if not numpy.array_equal(rhs.wave, self.wave):
            raise NotImplementedError('To perform arithmetic on spectra, their wavelength arrays '
                                      'must be identical.')
        flux = self.flux + rhs.flux

        if self.error is None and rhs.error is None:
            error = None
        else:
            error = numpy.zeros(self.size, dtype=float)
            if self.error is not None:
                error += numpy.square(self.error)
            if rhs.error is not None:
                error += numpy.square(rhs.error)
            error = numpy.sqrt(error)

        if self.mask is None and rhs.mask is None:
            mask = None
        else:
            mask = numpy.zeros(self.size, dtype=bool)
            if self.mask is not None:
                mask |= self.mask
            if rhs.mask is not None:
                mask |= rhs.mask
            
        if self.sres is not None and rhs.sres is not None \
                and not numpy.array_equal(self.sres, rhs.sres):
            warnings.warn('Spectral resolution is not correctly propagated.')
        sres = self.sres if self.sres is not None else rhs.sres
        return Spectrum(self.wave, flux, error=error, mask=mask, resolution=sres,
                        regular=self.regular, log=self.log)

    def __len__(self):
        return self.interpolator.x.size

    def __getitem__(self, s):
        """
        Access the flux data directly via slicing.
        """
        return self.interpolator.y[s]

    def _frequency(self):
        """Calculate the frequency in terahertz."""
        return 10*astropy.constants.c.to('km/s').value/self.wave

    def interp(self, w):
        """
        Linearly interpolate the flux at a provided wavelength.

        Args:
            w (:obj:`float`, `numpy.ndarray`_):
                Wavelength in angstroms. Can be a single wavelength
                or a wavelength array.
        """
        if not isinstance(w, numpy.ndarray) and w > self.interpolator.x[0] \
                and w < self.interpolator.x[-1]:
            return self.interpolator(w)

        indx = (w > self.interpolator.x[0]) & (w < self.interpolator.x[-1])
        sampled = numpy.zeros_like(w, dtype=float)
        sampled[indx] = self.interpolator(w[indx])
        return sampled

    @staticmethod
    def wcs_wavelength_vector(wcs, nwave, naxis=1, axis=0):
        """
        Generate the wavelength vector using the provided
        `astropy.wcs.WCS`_ object.
        
        Returns:
            `numpy.ndarray`_: Vector with wavelengths defined by
            :attr:`wcs`.
        """
        coo = numpy.ones((nwave,naxis), dtype=float)
        coo[...,axis] = numpy.arange(nwave)+1
        wave = wcs.all_pix2world(coo, 1)[...,axis]
        dimensionless = wcs.wcs.cunit[axis] is astropy.units.dimensionless_unscaled
        if dimensionless:
            # NOTE: If the axis is dimensionless, assume that the units
            # are angstroms!
            warnings.warn('No units defined for wavelength axis.  Assuming angstroms.')
        return wave if dimensionless else wave*wcs.wcs.cunit[axis].to('angstrom')

    @classmethod
    def from_fits(cls, ifile, waveext='WAVE', waveunits='angstrom', fluxext='FLUX',
                  fluxunits=None, errext=None, ivarext=None, maskext=None, resext=None,
                  tblext=None, **kwargs):
        """
        Construct a spectrum using data from a fits file.

        Fluxes must be in units of flux density (i.e., flux per unit
        angstrom or Hz).

        .. todo::

            - Allow for a "frequency" extension instead of
              wavelengths?

        Args:
            ifile (:obj:`str`):
                Name of the fits file.
            waveext (:obj:`str`, :obj:`int`, optional):
                Extension with the wavelength data. If set to 'WCS',
                the wavelengths are constructed using the WCS in the
                flux extension header. If the data is instead one
                column in a binary table extension, see ``tblext``.
            waveunits (:obj:`str`, optional):
                Units of the wavelengths in the fits file. WARNING:
                If you try to build the wavelength array using the
                WCS in the header, the wavelength array will be
                converted to angstroms based on the WCS units. If the
                WCS units are not defined, the values read from the
                WCS are not altered, but can be fixed using this
                keyword.
            fluxext (:obj:`str`, :obj:`int`, optional):
                Extension with the flux data. If the data is instead
                one column in a binary table extension, see
                ``tblext``. All fluxes are expected to be in units of
                flux density.
            fluxunits (:obj:`str`, optional):
                Units of the flux density in the fits file. Must be
                interpretable by `astropy.units.Unit`_. If None, the
                code will check the header of the flux extension (or
                the binary table extension; see ``tblext``) for the
                ``BUNIT`` keyword. If not None, this overrides any
                units provided by the fits file. If None and no
                ``BUNIT`` keyword is found, the units are assumed to
                be ``'1e-17 erg / (cm2 s angstrom)'``.
            errext (:obj:`str`, :obj:`int`, optional):
                Extension with the flux error. If None, no errors
                will be read. Keyword is mutually exclusive with
                ``ivarext``. If the data is instead one column in a
                binary table extension, see ``tblext``.
            ivarext (:obj:`str`, :obj:`int`, optional):
                Extension with the flux inverse variance. If None, no
                errors will be read. Otherwise, the data are read and
                then converted to 1-sigma errors. Keyword is mutually
                exclusive with ``errext``. If the data is instead one
                column in a binary table extension, see ``tblext``.
            maskext (:obj:`str`, :obj:`int`, optional):
                Extension with the mask data. If not None, the data
                in this extension must be castable to a boolean
                array. If None, all spectral data is assumed to be
                unmasked. If the data is instead one column in a
                binary table extension, see ``tblext``.
            resext (:obj:`str`, :obj:`int`, optional):
                Extension with the spectral resolution. If None, the
                spectral resolution will be undefined. If the data is
                instead one column in a binary table extension, see
                ``tblext``.
            tblext (:obj:`str`, :obj:`int`, optional):
                Instead of being organized in a series of ImageHDU
                extensions, the spectral data are in a BinTableHDU
                with this extension name/index. If provided (not
                None), the "extension" names/indices provided by the
                other keywords are instead interpreted as the columns
                in this binary table.
            **kwargs:
                Keyword arguments passed directly to the class
                instantiation (see :class:`Spectrum`). NOTE: If
                ``resolution`` is provided as one of the kwargs, it
                is given priority over the ``resext`` keyword in this
                method.

        Returns:
            :class:`Spectrum`: Object read from the fits file.

        Raises:
            ValueError:
                Raised if both ``errext`` and ``ivarext`` are set, or
                if the wavelength vector is supposed to be built from
                a WCS, but the data is organized in a binary table.
        """
        # Check the input
        if isinstance(waveext, str) and waveext.lower() == 'wcs' and tblext is not None:
            raise ValueError('WCS cannot be used to define the wavelength coordinate system '
                             'using a binary table data model.')
        if ivarext is not None and errext is not None:
            raise ValueError('An extensions/column can be provided for the inverse variance or '
                             'the error, not both.')
        _resext = resext
        if 'resolution' in kwargs and kwargs['resolution'] is not None and resext is not None:
            warnings.warn('Resolution and extension provided.  Preference given to the former.')
            _resext = None

        hdu = fits.open(ifile)
        # Get the flux
        flux = hdu[fluxext].data if tblext is None else hdu[tblext].data[fluxext]
        # Check the dimensionality
        if flux.ndim != 1:
            raise ValueError('Flux data must be one-dimensional.')

        # Get the mask. This needs to be instantiated before the error
        # vector, in case the mask incorporates ivar-to-error
        # conversion issues.
        mask = None if maskext is None \
                    else (hdu[maskext].data.astype(bool) if tblext is None 
                          else hdu[tblext].data[errext].astype(bool))

        # Get the error
        error = None if errext is None \
                    else (hdu[errext].data if tblext is None else hdu[tblext].data[errext])
        if ivarext is not None:
            error = numpy.ma.power(hdu[ivarext].data if tblext is None 
                                    else hdu[tblext].data[ivarext], -0.5)
            if mask is not None:
                mask |= error.mask
            error = error.data

        # Get the wavelength vector
        if isinstance(waveext, str) and waveext.lower() == 'wcs':
            wave = Spectrum.wcs_wavelength_vector(WCS(header=hdu[fluxext].header, fix=True),
                                                  flux.size)
        else:
            wave = hdu[waveext].data if tblext is None else hdu[tblext].data[waveext]
        # Impose the wavelength units
        wave *= astropy.units.Unit(waveunits).to('angstrom')

        # Fix the flux units
        if fluxunits is None:
            if tblext is not None and 'BUNIT' in hdu[tblext].header:
                fluxunits = hdu[tblext].header['BUNIT']
            elif 'BUNIT' in hdu[fluxext].header:
                fluxunits = hdu[fluxext].header['BUNIT']
        if fluxunits is not None:
#            print('Converting flux units from {0} to 1e-17 erg/s/cm2/angstrom'.format(fluxunits))
            if error is None:
                flux = convert_flux_units(wave, flux, fluxunits)
            else:
                flux, error = convert_flux_units(wave, flux, fluxunits, error=error)

        # Get the spectral resolution
        sres = None if _resext is None \
                    else (hdu[_resext].data if tblext is None else hdu[tblext].data[_resext])
        return cls(wave, flux, error=error, mask=mask, **kwargs) if sres is None \
                else cls(wave, flux, error=error, mask=mask, resolution=sres, **kwargs)

    @classmethod
    def from_ascii(cls, ifile, wavecol=0, waveunits='angstrom', fluxcol=1,
                   fluxunits='1e-17 erg / (cm2 s angstrom)', errcol=None, ivarcol=None,
                   maskcol=None, rescol=None, **kwargs):
        """
        Construct a spectrum using data from an ascii file.

        .. note::

            - All the columns are 0-indexed. I.e., the first column
              is column 0.
            - The core function used to read the data is
              `numpy.genfromtxt`_; i.e., anything that can be parsed
              by that function (e.g., gzipped files), can be parsed
              by this method. Note that all columns are read as
              floats. This means that any string columns are
              converted to nan, and any mask column must have the
              typical cast of 0 for False and anything else for True.

        Args:
            ifile (:obj:`str`):
                Name of the fits file.
            wavecol (:obj:`int`, optional):
                Column index with the wavelength data.
            fluxcol (:obj:`int`, optional):
                Column index with the flux data.
            fluxunits (:obj:`str`, optional):
                Units of the flux density read from the file. Must be
                interpretable by `astropy.units.Unit`_. If None,
                units are assumed to be ``'1e-17 erg / (cm2 s
                angstrom)'``.
            errcol (:obj:`int`, optional):
                Column index with the flux error. If None, no errors
                will be read. Keyword is mutually exclusive with
                ``ivarcol``.
            ivarcol (:obj:`int`, optional):
                Column index with the flux inverse variance. If None,
                no errors will be read. Otherwise, the data are read
                and then converted to 1-sigma errors. Keyword is
                mutually exclusive with ``errcol``.
            maskcol (:obj:`int`, optional):
                Column index with the mask data. If not None, the
                data in this extension must be castable to a boolean
                array. If None, all spectral data is assumed to be
                unmasked.
            rescol (:obj:`str`, :obj:`int`, optional):
                Column index with the spectral resolution. If None,
                the spectral resolution will be undefined.
            **kwargs:
                Keyword arguments passed directly to the class
                instantiation (see :class:`Spectrum`).

        Returns:
            :class:`Spectrum`: Object read from the ascii file.

        Raises:
            ValueError:
                Raised if both ``errcol`` and ``ivarcol`` are set.
        """
        # Check the input
        if ivarcol is not None and errcol is not None:
            raise ValueError('A column can be provided for the inverse variance or '
                             'the error, not both.')
        _rescol = rescol
        if 'resolution' in kwargs and kwargs['resolution'] is not None and rescol is not None:
            warnings.warn('Resolution and column provided.  Preference given to the former.')
            _rescol = None

        db = numpy.genfromtxt(ifile)
        flux = db[:,fluxcol]
        # Check the dimensionality
        if flux.ndim != 1:
            raise ValueError('Flux data must be one-dimensional.')
        wave = db[:,wavecol] * astropy.units.Unit(waveunits).to('angstrom')
        mask = None if maskcol is None else db[:,maskcol].astype(bool)
        error = None if errcol is None else db[:,errcol]
        if ivarcol is not None:
            error = numpy.ma.power(db[:,ivarcol], -0.5)
            if mask is not None:
                mask |= error.mask
            error = error.data
        if fluxunits is not None:
#            print('Converting flux units from {0} to 1e-17 erg/s/cm2/angstrom'.format(fluxunits))
            if error is None:
                flux = convert_flux_units(wave, flux, fluxunits)
            else:
                flux, error = convert_flux_units(wave, flux, fluxunits, error=error)
        sres = None if _rescol is None else db[:,_rescol]
        return cls(wave, flux, error=error, mask=mask, **kwargs) if sres is None \
                else cls(wave, flux, error=error, mask=mask, resolution=sres, **kwargs)

    def wavelength_step(self):
        """
        Return the wavelength step per pixel.

        Returns:
            `numpy.ndarray`_: Change in angstroms per pixel.
        """
        return angstroms_per_pixel(self.wave, log=self.log, regular=self.regular)

    def frequency_step(self):
        """
        Return the frequency step per pixel in THz.
        """
        return 10*astropy.constants.c.to('km/s').value*self.wavelength_step()/self.wave/self.wave

    def magnitude(self, wavelength=None, band=None, system='AB'):
        """
        Calculate the magnitude of the object in the specified band.

        If no arguments are provided, the magnitude is calculated for
        the entire :attr:`flux` vector. Otherwise, the magnitude is
        calculated at a single wavelength (see ``wavelength``) or
        over a filter (see ``band``).

        Args:
            wavelength (:obj:`float`, optional):
                The wavelength at which to calculate the magnitude.
            band (:class:`~synospec.etc.efficiency.FilterResponse`, optional):
                Object with the filter response function
            system (:obj:`str`, optional):
                Photometric system. Currently must be ``AB`` or
                ``Vega``.

        Returns:
            :obj:`float`, `numpy.ndarray`_: The one or more magnitude
            measurements, depending on the input.

        Raises:
            NotImplementedError:
                Raised if the photometric system is not known.
        """
        if wavelength is not None and band is not None:
            warnings.warn('Provided both wavelength and band; wavelength takes precedence.')

        # TODO: Check input.
        if system == 'AB':
            if self.nu is None:
                self._frequency()
            if self.fnu is None:
                # Flux in microJanskys
                self.fnu = convert_flux_density(self.wave, self.flux)
            if wavelength is not None:
                fnu = self.fnu[numpy.argmin(numpy.absolute(self.wave - wavelength))]
            elif band is not None:
                dnu = self.frequency_step()
                fnu = numpy.sum(band(self.wave)*self.fnu*dnu) / numpy.sum(band(self.wave)*dnu)
            else:
                fnu = self.fnu
            # NOTE: This is identically -2.5*numpy.log10(fnu) - 48.6,
            # just accounting for the units of fnu
            return -2.5*numpy.log10(fnu) + 23.9

        if system == 'Vega':
            # NOTE: 
            #   - Vega mags are done brute force, not using a
            #     precomputed set of zero-points.
            #   - This doesn't account for any resolution differences
            #     between this spectrum and the Vega spectrum.
            if self._vega is None:
                # Lazy load the Vega spectrum
                self._vega = VegaSpectrum()
            if wavelength is not None:
                # At one wavelength
                vega_flux = self._vega.flux[numpy.argmin(numpy.absolute(self._vega.wave
                                                                        - wavelength))]
                flux = self.flux[numpy.argmin(numpy.absolute(self.wave - wavelength))]
            elif band is not None:
                # Over a band
                vega_flux = numpy.sum(band(self.wave) * self.wave * self._vega.interp(self.wave)
                                      * self.wavelength_step())
                flux = numpy.sum(band(self.wave) * self.wave * self.flux * self.wavelength_step())
            else:
                # Full spectrum
                flux = self.flux
                vega_flux = self._vega.interp(self.wave)
            return -2.5 * numpy.log10(flux / vega_flux) + 0.03

        raise NotImplementedError('Photometric system {0} not implemented.'.format(system))

    def rescale(self, factor):
        """
        Rescale the spectrum by the provided factor.

        The spectral data are modified *in-place*; nothing is
        returned.

        Args:
            factor (scalar-like or array-like):
                Factor to multiply the fluxes. If it's a vector, it
                must have the same length as the existing spectrum.
        """
        self.interpolator.y *= factor
        if self.error is not None:
            # Adjust the error
            self.error *= numpy.absolute(factor)
        if self.fnu is not None:
            # Adjust the flux density per Hz
            self.fnu *= factor

    def rescale_flux(self, wave, flux):
        """
        Rescale the spectrum to be specifically the flux value at the
        provided wavelength. The input flux should be in units of
        1e-17 erg/s/cm^2/angstrom.

        Args:
            wave (:obj:`float`):
                The wavelength at which to base the rescaling.
            flux (:obj:`float`):
                The target value of the flux at the provided
                wavelength.
        """
        self.rescale(flux/self.interp(wave))

    def rescale_magnitude(self, new_mag, wavelength=None, band=None, system='AB'):
        """
        Rescale the spectrum to an input magnitude.

        Must provide either ``wavelength`` or ``band``. The object is
        edited in-place.

        Args:
            new_mag (:obj:`float`):
                The target magnitude
            wavelength (:obj:`float`, optional):
                The wavelength at which to calculate the magnitude.
            band (:class:`~synospec.etc.efficiency.FilterResponse`, optional):
                Object with the filter response function
            system (:obj:`str`, optional):
                Photometric system. Currently must be ``AB`` or
                ``Vega``.

        Raises:
            ValueError:
                Raised if neither ``wavelength`` nor ``band`` are provided.
            NotImplementedError:
                Raised if the photometric system is not known.
        """
        if wavelength is None and band is None:
            raise ValueError('Must provide either wavelength or the bandpass filter.')
        if system not in ['AB', 'Vega']:
            raise NotImplementedError('Photometric system {0} not implemented.'.format(system))
        dmag = new_mag - self.magnitude(wavelength=wavelength, band=band, system=system)
        self.rescale(numpy.power(10., -dmag/2.5))

    def photon_flux(self, inplace=True):
        r"""
        Convert the spectrum from 1e-17 erg/s/cm^2/angstrom to
        photons/s/cm^2/angstrom.

        Args:
            inplace (:obj:`bool`, optional):
                If ``inplace is True``, the spectrum is modified in
                place and None is returned; otherwise, the converted
                flux vector is returned.

        Returns:
            `numpy.ndarray`_: Flux vector in photons/s/cm^2/angstrom
            if ``inplace is False``; otherwise, nothing is returned.
        """
        ergs_per_photon = astropy.constants.h.to('erg s') * astropy.constants.c.to('angstrom/s') \
                            / (self.wave * astropy.units.angstrom)
        return self.rescale(1e-17 / ergs_per_photon.value) if inplace \
                    else self.interpolator.y * 1e-17 / ergs_per_photon.value

    def plot(self, ax=None, show=False, **kwargs):
        """
        Construct a plot of the spectrum.

        Args:
            ax (`matplotlib.axes.Axes`_):
                Axes for the plot. If None, a new instance is
                generated using ``pyplot.subplot``.
            show (:obj:`bool`, optional):
                Show the plot instead of just returning the modified
                `matplotlib.axes.Axes`_ instance.
            **kwargs:
                Other keywords (like ``color``) passed directly to
                ``pyplot.plot``.

        Returns:
            `matplotlib.axes.Axes`_: Modified or new instance;
            returned only if ``show is False``.
        """
        _ax = pyplot.subplot() if ax is None else ax
        _ax.plot(self.wave, self.flux, **kwargs)
        if show:
            pyplot.show()
        return _ax

    def resample(self, wave=None, dwave=None, log=False, step=True):
        """
        Resample the spectrum to a uniform wavelength grid.

        Args:
            wave (`numpy.ndarray`_, optional):
                New wavelength array. Must be linearly or
                log-linearly sampled. If None, the new wavelength
                grid covers the full extent of the current grid with
                a sampling based in the minimum pixel separation
                (geometrically or linearly) in the current spectrum.
            dwave (:obj:`float`, optional):
                Step for each pixel. If ``log`` is True, this must be
                the geometric step; i.e., the difference in the
                base-10 log of the wavelength of adjacent pixels.
            log (:obj:`bool`, optional):
                Flag that the wavelength array is log-linearly
                sampled.
            step (:obj:`bool`, optional):
                Treat the input function as a step function during
                the resampling integration. If False, use a linear
                interpolation between pixel samples.

        Returns:
            :class:`Spectrum`: Returns a resampled version of itself.
        """
        if wave is None:
            wave = self.wave
        rng = numpy.log10(wave[[0,-1]]) if log else wave[[0,-1]]
        if dwave is None:
            dwave = numpy.amin(numpy.diff(numpy.log10(wave))) \
                        if log else numpy.amin(numpy.diff(wave))
        npix = int(numpy.diff(rng)/dwave)+1
        r = Resample(self.flux, e=self.error, mask=self.mask, x=self.wave,
                     inLog=self.regular and self.log, newRange=wave[[0,-1]], newpix=npix,
                     newLog=log, step=step)
        sres = None if self.sres is None \
                    else interpolate.interp1d(self.wave, self.sres, assume_sorted=True,
                                              bounds_error=False,
                                              fill_value=(self.sres[0],self.sres[-1]))(r.outx)
        return Spectrum(r.outx, r.outy, error=r.oute, mask=r.outf < 0.8, resolution=sres,
                        regular=True, log=log)

    def match_resolution(self, resolution, wave=None):
        r"""
        Match the spectral resolution of the spectrum to the provided
        value.

        The object must be a regularly sampled, either linearly or
        log-linearly (see :attr:`log`). If it isn't, use
        :func:`resample`.

        Args:
            resolution (:obj:`float`, array-like):
                The single value for the spectral resolution
                (:math:`R=\lambda/\Delta\lambda`) or a vector for a
                wavelength-dependent spectral resolution. If a
                vector, the length of the vector must exactly match
                the wavelength vector of the object or the wavelength
                vector must also be provided (see ``wave``).
            wave (`numpy.ndarray`_, optional):
                Wavelength vector for the provided spectral
                resolution. If None, ``resolution`` must either be a
                single number or a vector with the same length as
                :attr:`wave`.

        Returns:
            :class:`Spectrum`: Returns an object with the spectral
            resolution matched (as best as can be done). If
            necessary, the spectrum is masked where the spectral
            resolution could not be matched in detail.
        """
        if not self.regular:
            raise ValueError('The spectrum must be regularly binned; try running resample first.')
        if wave is None:
            wave = self.wave
        new_sres = numpy.atleast_1d(resolution)
        if new_sres.size == 1:
            new_sres = numpy.full(self.wave.size, resolution, dtype=float)
        if new_sres.size != wave.size:
            raise ValueError('Spectral resolution does not match length of wavelength vector.')

        ivar = None if self.error is None else numpy.ma.power(self.error, -2).data
        new_flux, new_sres, _, new_mask, new_ivar \
                = match_spectral_resolution(self.wave, self.flux, self.sres, wave, new_sres,
                                            ivar=ivar, mask=self.mask, log10=self.log,
                                            new_log10=self.log)
        new_mask |= numpy.invert(new_flux > 0)

        if new_ivar is None:
            new_error = None
        else:
            new_error = numpy.ma.power(new_ivar, -0.5)
            new_mask |= new_error.mask
            new_error = new_error.data    
        return Spectrum(self.wave, new_flux, error=new_error, mask=new_mask, resolution=new_sres,
                        use_sampling_assessments=True)

    def redshift(self, z):
        """
        Redshift the spectrum.

        Spectrum is in 1e-17 erg/s/cm^2/angstrom, so this shifts the
        wavelength vector by 1+z and rescales the flux by 1+z to keep
        the flux per *observed* wavelength constant; this is not the
        same as surface-brightness dimming, **which is not accounted
        for**. S/N is kept fixed.

        Args:
            z (:obj:`float`):
                Target redshift
        """
        self.interpolator.x *= (1+z)
        self.rescale(1/(1+z))

    def write(self, ofile, overwrite=False):
        """
        Write the spectrum to a fits file.
        """
        if os.path.isfile(ofile) and not overwrite:
            raise FileExistsError('File exists. Use a different filename or set overwrite=True.')
        
        hdu = fits.HDUList([fits.PrimaryHDU(),
                            fits.ImageHDU(name='WAVE', data=self.wave),
                            fits.ImageHDU(name='FLUX', data=self.flux)])
        if self.error is not None:
            hdu += [fits.ImageHDU(name='ERROR', data=self.error)]
        if self.mask is not None:
            hdu += [fits.ImageHDU(name='MASK', data=self.mask.astype(int))]
        if self.sres is not None:
            hdu += [fits.ImageHDU(name='SPECRES', data=self.sres)]
        hdu.writeto(ofile, overwrite=True)


class EmissionLineSpectrum(Spectrum):
    r"""
    Define an emission-line spectrum.

    .. todo::
        - use MaNGA DAP functions
    
    Flux units are 1e-17 erg/s/cm^2/angstrom.

    Args:
        wave (array-like):
            1D wavelength data in angstroms.  Expected to be sampled
            linearly or geometrically.  These are the *observed*
            wavelengths.
        flux (array-like):
            The total fluxes of one or more emission lines in 1e-17 erg/s/cm^2.
        restwave (array-like):
            The central rest wavelengths of one or emission lines in angstroms.
        fwhm (array-like):
            The FWHM of the Gaussian line profiles.  If the resolution
            vector is provided, these are assumed to be the *intrinsic*
            widths, such that the line included in the spectrum has an
            observed width determined by the quadrature sum of the
            intrinsic and instrumental widths.  If the resolution vector
            is not provided, the line simply has the provided FWHM.
        units (:obj:`str`, array-like, optional):
            The units of the provided FWHM data. Must be either
            'km/s' or 'ang'. Can be a single string or an array with
            the units for each provided value.
        redshift (scalar-like, optional):
            If provided, the emission-line wavelengths are redshifted.
        continuum (array-like, optional):
            If provided, this is the 1D continuum placed below the line,
            which must have the same length as the input wavelength
            vector.  The continuum is 0 if not provided.
        resolution (array-like, optional):
            1D spectral resolution
            (:math:`$R=\lambda/\Delta\lambda$`). If None, the width
            of the lines is set by ``fwhm`` only.
        log (:obj:`bool`, optional):
            Spectrum is sampled in steps of log base 10.
    """
    def __init__(self, wave, flux, restwave, fwhm, units='ang', redshift=0.0, continuum=None,
                 resolution=None, log=False):

        _flux = numpy.atleast_1d(flux).ravel()
        nlines = _flux.size
        _restwave = numpy.atleast_1d(restwave).ravel()
        if _restwave.size != nlines:
            raise ValueError('Number of rest wavelengths does not match the number of fluxes.')
        _fwhm = numpy.atleast_1d(fwhm).ravel()
        if _fwhm.size != nlines:
            raise ValueError('Number of FWHM values does not match the number of fluxes.')
        _units = numpy.atleast_1d(units).ravel()
        # Check the input
        if not numpy.all(numpy.isin(_units, ['km/s', 'ang'])):
            raise ValueError('FWHM units must be \'km/s\' or \'ang\'.')
        if _units.size == 1:
            _units = numpy.repeat(_units, _flux.size)
        if _units.size != nlines:
            raise ValueError('Number of unit values does not match the number of fluxes.')

        if resolution is not None and hasattr(resolution, '__len__') \
                and len(wave) != len(resolution):
            raise ValueError('Resolution vector must match length of wavelength vector.')
        if continuum is not None and len(wave) != len(continuum):
            raise ValueError('Continuum vector must match length of wavelength vector.')

        # Set the line parameters
        _linewave = _restwave*(1+redshift)
        indx = (_linewave > wave[0]) & (_linewave < wave[-1])
        if not numpy.any(indx):
            raise ValueError('Redshifted lines are all outside of the provided wavelength range.')

        sig2fwhm = numpy.sqrt(8.0 * numpy.log(2.0))
        sigma = _fwhm/sig2fwhm
        # Convert the FWHM as needed based on the sampling and units
        in_ang = _units == 'ang'
        in_kms = _units == 'km/s'
        if log and numpy.any(in_ang):
            # Convert to km/s
            sigma[in_ang] = astropy.constants.c.to('km/s').value*sigma[in_ang]/_linewave[in_ang]
        elif not log and numpy.any(in_kms):
            # Convert to angstroms
            sigma[in_kms] = _linewave[in_kms]*sigma[in_kms]/astropy.constants.c.to('km/s').value

        # Add the instrumental resolution in quadrature to the
        # intrinsic width, if the resolution is provided
        if resolution is None:
            _resolution = None
        else:
            _resolution = resolution if hasattr(resolution, '__len__') \
                                else numpy.full_like(wave, resolution, dtype=float)
            sigma_inst = astropy.constants.c.to('km/s').value/_resolution/sig2fwhm if log else \
                            wave/_resolution/sig2fwhm
            interp = interpolate.interp1d(wave, sigma_inst, assume_sorted=True, bounds_error=False,
                                          fill_value=0.)
            sigma = numpy.sqrt(numpy.square(sigma) + numpy.square(interp(_linewave))) 

        # Convert parameters to pixel units
        _dw = spectral_coordinate_step(wave, log=log)
        _linepix = (numpy.log10(_linewave) - numpy.log10(wave[0]) \
                        if log else _linewave - wave[0])/_dw

        # Flux to pixel units so that the spectrum has units of flux
        # density (flux per angstrom); less accurate when spectrum is
        # logarithmically binned
        dl = _linewave*(numpy.power(10.,_dw/2)-numpy.power(10.,-_dw/2)) if log else _dw
        _flux /= dl

        # Convert sigma to pixels
        sigma /= (astropy.constants.c.to('km/s').value*_dw*numpy.log(10.) if log else dl)
        
        # Construct the emission-line spectrum
        pix = numpy.arange(wave.size)
        spectrum = numpy.zeros(wave.size, dtype=float) if continuum is None else continuum.copy()
        profile = IntegratedGaussianLSF()
        for i in range(nlines):
            if not indx[i]:
                continue
            p = profile.parameters_from_moments(_flux[i], _linepix[i], sigma[i])
            spectrum += profile(pix, p)

        # Instantiate
        super().__init__(wave, spectrum, resolution=_resolution, log=log)


# 8329-6104
class BlueGalaxySpectrum(Spectrum):
    """
    An example blue galaxy spectrum pulled from the MaNGA survey.
    """
    def __init__(self, redshift=0.0):
        fitsfile = str(data_file(filename='galaxy') / 'blue_galaxy_8329-6104.fits')
        hdu = fits.open(fitsfile)
        wave = hdu['WAVE'].data * (1+redshift)
        flux = hdu['FLUX'].data
        super().__init__(wave, flux, log=True)

    @classmethod
    def from_file(cls):
        raise NotImplementedError('Spectrum for blue galaxy is fixed.')


# 8131-6102
class RedGalaxySpectrum(Spectrum):
    """
    An example red galaxy spectrum pulled from the MaNGA survey.
    """
    def __init__(self, redshift=0.0):
        fitsfile = str(data_file(filename='galaxy') / 'red_galaxy_8131-6102.fits')
        hdu = fits.open(fitsfile)
        wave = hdu['WAVE'].data * (1+redshift)
        flux = hdu['FLUX'].data
        super().__init__(wave, flux, log=True)

    @classmethod
    def from_file(cls):
        raise NotImplementedError('Spectrum for red galaxy is fixed.')


#class MaunakeaSkySpectrumOLD(Spectrum):
#    def __init__(self):
#        fitsfile = os.path.join(os.environ['SYNOSPEC_DIR'], 'data/sky/manga/apo2maunakeasky.fits')
#        hdu = fits.open(fitsfile)
#        wave = hdu['WAVE'].data
#        flux = hdu['FLUX'].data
#        super(MaunakeaSkySpectrumOLD, self).__init__(wave, flux, log=True)
#
#    @classmethod
#    def from_file(cls):
#        raise NotImplementedError('Maunakea sky spectrum is fixed.')


class MaunakeaSkySpectrum(Spectrum):
    """
    The default, empirical dark night-sky spectrum at Maunakea
    provided by Chuck Steidel.
    """
    def __init__(self):
        fitsfile = str(data_file(filename='sky') / 'lris_esi_skyspec_fnu.fits')
        init = Spectrum.from_fits(fitsfile, waveext='WCS', fluxext=0, airwave=True,
                                  use_sampling_assessments=True)
        init.regular = False
        init.log = False
        wave = numpy.arange(3100., 10500., 2.)
        init = init.resample(wave=wave, dwave=2., log=False)
        wave = init.wave
        flux = init.flux
        indx = numpy.invert(flux > 0) # & (wave < 3210)
        flux[indx] = numpy.median(flux[(wave > 3210) & (wave < 4000)])
        super().__init__(wave, flux)

    @classmethod
    def from_file(cls):
        raise NotImplementedError('Maunakea sky spectrum is fixed.')


class MtHamiltonSkySpectrum(Spectrum):
    """
    An empirical dark night-sky spectrum at Mt. Hamilton provided by X.
    Prochaska/B. Holden.
    """
    def __init__(self):
        fitsfile = str(data_file(filename='sky') / 'lick_sky_d55_2011aug29.fits.gz')
        init = Spectrum.from_fits(fitsfile, waveext=0, fluxext=1, airwave=True,
                                  use_sampling_assessments=True)
        init.rescale(1e17)
        init.regular = False
        init.log = False

        wave = numpy.arange(3160., 8315., 1.)
        init = init.resample(wave=wave, dwave=1., log=False)
        super().__init__(init.wave, init.flux)

#        wave = init.wave
#        flux = init.flux
#        indx = numpy.invert(flux > 0) # & (wave < 3210)
#        flux[indx] = numpy.median(flux[(wave > 3210) & (wave < 4000)])
#        super().__init__(wave, flux)

    @classmethod
    def from_file(cls):
        raise NotImplementedError('Maunakea sky spectrum is fixed.')


class ABReferenceSpectrum(Spectrum):
    """
    Construct a spectrum with a constant flux of 3631 Jy.

    Inherits from :class:`Spectrum`, which we take to mean that the flux
    is always in units of 1e-17 erg/s/cm^2/angstrom.

    Args:
        wave (`numpy.ndarray`_):
            Wavelength vector for the spectrum.
        resolution (:obj:`float`, array-like, optional):
            1D spectral resolution (:math:`$R=\lambda/\Delta\lambda$`)
        regular (:obj:`bool`, optional):
            Spectrum is regularly sampled, either linearly or log-linearly
        log (:obj:`bool`, optional):
            Spectrum is sampled in steps of log base 10. If regular
            is False, this is ignored.
    """
    def __init__(self, wave, resolution=None, regular=True, log=False):
        norm = numpy.power(10., 29 - 48.6/2.5)  # Reference flux in microJanskys
        fnu = numpy.full_like(wave, norm, dtype=float)
        flambda = convert_flux_density(wave, fnu, density='Hz')
        super().__init__(wave, flambda, resolution=resolution, log=log, regular=regular)


class VegaSpectrum(Spectrum):
    """
    Return the spectrum of Vega constructed using data from:

    https://ssb.stsci.edu/cdbs/calspec/alpha_lyr_stis_009.fits

    Downloaded on 3 Apr 2020.  Wavelengths are in vacuum.

    Args:
        waverange (array-like, optional):
            Used to limit the wavelength range of the returned
            spectrum. The maximum wavelength range of the current
            spectrum goes from 900 Angstrom to 300 microns.
    """
    def __init__(self, waverange=None):
        if waverange is not None and numpy.asarray(waverange).size != 2:
            raise ValueError('Wavelength range must be a two-element list, tuple, etc.')
        fitsfile = str(data_file(filename='spectra') / 'alpha_lyr_stis_009.fits')
        hdu = fits.open(fitsfile)
        indx = numpy.ones(hdu[1].data['WAVELENGTH'].size, dtype=bool) if waverange is None \
                    else (hdu[1].data['WAVELENGTH'] > waverange[0]) \
                            & (hdu[1].data['WAVELENGTH'] < waverange[1])
        resolution=hdu[1].data['WAVELENGTH']/hdu[1].data['FWHM']
        super().__init__(hdu[1].data['WAVELENGTH'][indx], hdu[1].data['FLUX'][indx]*1e17,
                         resolution=resolution[indx], regular=False)

