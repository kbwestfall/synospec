"""
Provides Efficiency object, primarily just a linear interpolator for
efficiency data.

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

from .. import data_file

class Efficiency:
    r"""
    Base class for efficiency data.

    Provided a wavelength independent efficiency value (:math:`\eta`;
    ``eta``), or a sampled vector of :math:`\eta` vs. wavelength
    (:math:`\lambda`; ``wave``), this class mainly just allows for
    access and interpolation of those data. When provided as a
    function of wavelength, the efficiency is assumed to be 0 outside
    the bounds of the provided wavelength vector.

    Args:
        eta (:obj:`float`, array-like):
            Constant or 1D efficiency data.
        wave (array-like, optional):
            1D array with wavelengths in angstroms. If not provided,
            ``eta`` must be a constant; if provided, ``eta`` must be
            an array.

    Attributes:
        interpolator (`scipy.interpolate.interp1d`_):
            Linear interpolator used to sample efficiency at any
            wavelength, if vectors are provided.

    Raises:
        TypeError:
            Raised if ``wave`` is not provided and ``eta`` has a
            ``len`` attribute.
        ValueError:
            Raised if ``wave`` is provided and ``eta`` is *not* a
            vector or if ``wave`` and ``eta`` have different shapes.

    """
    def __init__(self, eta, wave=None):
        if wave is None:
            # eta has to be a constant
            if hasattr(eta, '__len__'):
                raise TypeError('If instantiated without a wavelength vector, the efficiency '
                                'must be a wavelength-independent contant.')
            self.interpolator = None
            self._eta = eta
        else:
            _wave = numpy.atleast_1d(wave)
            _eta = numpy.atleast_1d(eta)
            if _eta.ndim > 1 or len(_eta) == 1:
                raise ValueError('When providing wavelengths, efficiency must be a 1D vector '
                                 'with more than one element.')
            if _wave.shape != _eta.shape:
                raise ValueError('Efficiency and wavelengths must have the same shape.')
            self.interpolator = interpolate.interp1d(_wave, _eta, assume_sorted=True,
                                                     bounds_error=False, fill_value=0.0)
            self._eta = None
    
    @classmethod
    def from_file(cls, datafile, wave_units='angstrom'):
        """
        Read from an ascii file.

        The format of the file must be columnated data with the first
        column providing the wavelength and the second column
        providing the efficiency. The data is read using
        `numpy.genfromtxt`_.

        Args:
            datafile (:obj:`str`):
                Ascii file with the data.
            wave_units (:obj:`str`, optional):
                Units of the wavelength in the file. Must be
                interpretable by `astropy.units`_ and convertable to
                angstroms.
        """
        if not os.path.isfile(datafile):
            raise FileNotFoundError('File does not exist: {0}'.format(datafile))
        db = numpy.genfromtxt(datafile)
        return cls(db[:,1], wave=db[:,0]*units.Unit(wave_units).to('angstrom'))

    def __call__(self, wave):
        """
        Return the efficiency interpolated at the provided
        wavelengths.

        Args:
            wave (array-like):
                Wavelengths at which to return the efficiency.
        """
        _wave = numpy.atleast_1d(wave)
        if self.interpolator is None:
            return self._eta if _wave.size == 1 \
                        else numpy.full(_wave.shape, self._eta, dtype=float)
        _eta = self.interpolator(_wave)
        return _eta if hasattr(wave, '__len__') else _eta[0]

    def __getitem__(self, k):
        """
        Return the selected element of the efficiency vector.

        Args:
            k (:obj:`int`):
                Vector element to return.
        """
        if self.interpolator is None:
            # TODO: Handle if k is a slice...
            warnings.warn('Efficiency is not a vector!  Returning constant value.')
            return self._eta
        return self.interpolator.y[k]

    @property
    def wave(self):
        """The wavelength vector with direct efficiency values."""
        if self.interpolator is None:
            warnings.warn('Efficiency is wavelength independent.')
            return None
        return self.interpolator.x

    @property
    def eta(self):
        """The effeciency vector."""
        if self.interpolator is None:
            return self._eta
        return self.interpolator.y

    def rescale(self, scale):
        """
        Rescale the efficiency data.

        The internal attributes are directly modified.

        Args:
            scale (scalar-like, array-like):
                Factor for the efficiency. *No check is performed to
                make sure this doesn't result in a value larger than
                1!* Value must be broadcastable to the shape of the
                efficiency vector.
        """
        if self.interpolator is None:
            self._eta *= scale
        else:
            self.interpolator.y *= scale


class CombinedEfficiency(Efficiency):
    """
    A class that combines multiple efficiencies that can be accessed
    separately or act as a single efficiency.

    Accessing this object as you would an
    :class:`~synospec.etc.efficiency.Efficiency` object always returns
    the total efficiency.

    Args:
        efficiencies (:obj:`list`, :obj:`dict`):
            The set of efficiencies to combine. Nominally this should
            be a dictionary that gives the efficiencies and a keyword
            identifier for each. A list can be entered, meaning that
            the efficiencies can only be access by their index, not a
            keyword. Each list element or dictionary value *must* be
            an :class:`~synospec.etc.efficiency.Efficiency` object.
        wave (array-like, optional):
            Wavelengths of/for efficiency measurements.

    Attributes:
        efficiencies (:obj:`dict`):
            The efficiencies combined. Access to individual
            efficiencies is by keyword; if keywords are not provided,
            access is by single integer index.

    Raises:
        TypeError:
            Raised if ``efficiencies`` is not a list or dictionary,
            or if the elements of either are not
            :class:`~synospec.etc.efficiency.Efficiency` objects.
    """
    def __init__(self, efficiencies, wave=None):
        if isinstance(efficiencies, list):
            self.efficiencies = dict([(i,eff) for i,eff in enumerate(efficiencies)])
        elif isinstance(efficiencies, dict):
            self.efficiencies = efficiencies
        else:
            raise TypeError('Efficiencies to include must provided as a list or dict.')

        # Make sure the components of self are Efficiency objects
        for eff in self.efficiencies.values():
            if not isinstance(eff, Efficiency):
                raise TypeError('Each element of input must be an Efficiency object.')

        if wave is None:
            # Consolidate wavelengths from other efficiencies
            wave = numpy.empty(0, dtype=float)
            for inp in self.efficiencies.values():
                if inp.wave is None:
                    continue
                wave = numpy.append(wave, inp.wave)
            wave = None if len(wave) == 0 else numpy.unique(wave)
            if wave is None:
                warnings.warn('No wavelengths provided for any efficiencies to combine.')

        # Construct the total efficiency
        total = 1. if wave is None else numpy.ones_like(wave, dtype=float)
        for eff in self.efficiencies.values():
            total *= (eff.eta if wave is None else eff(wave))

        super(CombinedEfficiency, self).__init__(total, wave=wave)

    @classmethod
    def from_total(cls, total, wave=None):

        """
        Construct the combined efficiency object from the total
        efficiency only. This is virtually identical to a single
        :class:`~synospec.etc.efficiency.Efficiency` object.

        Args:
            total (:obj:`float`, array-like):
                Constant or 1D total efficiency data.
            wave (array-like, optional):
                1D array with wavelengths in angstroms. If not
                provided, ``eta`` must be a constant; if provided,
                ``eta`` must be an array.
        """
        return cls({'total': Efficiency(total, wave=wave)})

    def keys(self):
        """The iterable with the efficiency keywords."""
        return self.efficiencies.keys()

    def __getitem__(self, key):
        """Return the specified efficiency."""
        return self.efficiencies[key]


class FiberThroughput(Efficiency):
    """
    Fiber efficiency object.

    Currently this only returns efficiencies based on pre-computed
    data. See $SYNOSPEC_DIR/data/efficiency/fibers/.

    Args:
        fiber (:obj:`str`, optional):
            The fiber vendor. Currently this can only be 'polymicro'.
    """
    def __init__(self, fiber='polymicro'):
        datafile = FiberThroughput.select_data_file(fiber)
        if not os.path.isfile(datafile):
            raise FileNotFoundError('No file: {0}'.format(datafile))
        db = numpy.genfromtxt(datafile)
        super(FiberThroughput, self).__init__(db[:,1], wave=db[:,0])

    @staticmethod
    def select_data_file(fiber):
        if fiber == 'polymicro':
            return str(data_file(filename='efficiency') / 'fibers' / 'polymicro.db')
        raise NotImplementedError('Unknown fiber type: {0}'.format(fiber))


class FilterResponse(Efficiency):
    """
    The efficiency of a broad-band filter.

    Args:
        band (:obj:`str`, optional):
            The band to use. Options are for the response functions
            in the $SYNOSPEC_DIR/data/broadband_filters directory.

    Raises:
        FileNotFoundError:
            Raised if the default file for the given band is not
            available.
    """
    def __init__(self, band='g'):
        dfile = data_file(filename='broadband_filters') / f'gunn_2001_{band}_response.db'
        if not dfile.is_file():
            raise FileNotFoundError(f'No file: {str(dfile)}')
        db = numpy.genfromtxt(str(dfile))
        super(FilterResponse, self).__init__(db[:,1], wave=db[:,0])


class AtmosphericThroughput(Efficiency):
    """
    Atmospheric throughput.

    Args:
        airmass (:obj:`float`, optional):
            Airmass of the observation.
        location (:obj:`str`, optional):
            Location of the observations. Currently this can only be
            ``'maunakea'``.
    """
    def __init__(self, airmass=1.0, location='maunakea'):
        if location == 'maunakea':
            dfile = data_file(filename='sky') / 'mauna_kea_extinction.db'
            if not dfile.is_file():
                raise FileNotFoundError(f'No file: {str(dfile)}')
            db = numpy.genfromtxt(str(dfile))
        else:
            raise NotImplementedError('Extinction unknown at {0}.'.format(location))
        self.airmass = airmass
        super(AtmosphericThroughput, self).__init__(numpy.power(10, -0.4*db[:,1]*self.airmass),
                                                    wave=db[:,0])

    def reset_airmass(self, airmass):
        self.rescale(numpy.power(self.interpolator.y, (airmass/self.airmass -1)))


class SpectrographThroughput(CombinedEfficiency):
    """
    Define the spectrograph throughput from the telescope focal plane
    to the detector.

    Args:
        wave (array-like, optional):
            Wavelength vector at which to define the total
            efficiency.
        coupling (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Focal-plain coupling efficiency.
        fibers (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Fiber efficiency.
        grating (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Grating efficiency
        camera (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Camera efficiency
        detector (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Detector efficiency
        other (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Other efficiency terms
    """
    def __init__(self, wave=None, coupling=None, fibers=None, grating=None, camera=None,
                 detector=None, other=None):
        values = inspect.getargvalues(inspect.currentframe())
        keys = numpy.array(values.args[2:])
        objects = numpy.array([values.locals[key] for key in keys])
        indx = numpy.array([o is not None for o in objects])
        efficiencies = dict([(k,o) for k,o in zip(keys[indx], objects[indx])])
        super(SpectrographThroughput, self).__init__(efficiencies, wave=wave)


class SystemThroughput(CombinedEfficiency):
    """
    Define the full system throughput from the top of the telescope
    to the detector; i.e., ratio of detected-to-incident photons.

    Args:
        wave (array-like, optional):
            Wavelength vector at which to define the total
            efficiency.
        spectrograph (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Spectrograph efficiency
        telescope (:class:`~synospec.etc.efficiency.Efficiency`, optional):
            Telescope efficiency
    """
    # TODO: Also allow for 'other' here?
    def __init__(self, wave=None, spectrograph=None, telescope=None):
        values = inspect.getargvalues(inspect.currentframe())
        keys = numpy.array(values.args[2:])
        objects = numpy.array([values.locals[key] for key in keys])
        indx = numpy.array([o is not None for o in objects])
        efficiencies = dict([(k,o) for k,o in zip(keys[indx], objects[indx])])
        super(SystemThroughput, self).__init__(efficiencies, wave=wave)

