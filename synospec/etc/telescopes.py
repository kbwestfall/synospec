#!/bin/env/python3
# -*- encoding utf-8 -*-
"""
Define telesope parameters

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import numpy

from . import efficiency
from .. import data_file

# TODO: Define a `TelescopePort` or `TelescopeFocalPlane` class that
# holds the f-ratio, platescale, and throughput information, and allow
# each `Telescope` definition to have multiple focal planes?

class Telescope:
    """
    Collect useful telescope parameters.

    Args:
        longitude (scalar-like):
            Earth coordinate with the location of the telescope in
            degrees.
        latitude (scalar-like):
            Earth coordinate with the location of the telescope in
            degrees.
        elevation (scalar-like):
            Earth elevation above sea level of the telescope in
            meters.
        fratio (scalar-like):
            F-ratio or beam speed (focal length over diameter) for
            the telescope focal plane.
        platescale (scalar-like):
            Telescope platescale in mm/arcsec.
        throughput (:obj:`float`, :class:`~synospec.etc.efficiency.Efficiency`, :class:`~synospec.etc.efficiency.CombinedEfficiency`, optional):
            The throughput of the telescope from the top of the
            telescope to the focal plane.
        area (scalar-like, optional):
            The true or effective area of the telescope aperture in
            square centimeters. If not provided, calculated using
            ``diameter``. Must be provided if ``diameter`` is not.
        diameter (scalar-like, optional):
            Telescope diameter in meters. If provided, used to set
            the telescope area. Must be provided if `area` is not.
        obstruction (scalar-like, optional):
            The unitless fraction of the total telescope area lost
            due to the central obstruction. If provided, the
            telescope area is multiplied by (1-`obstruction`) to
            obtain its effective area. If not provided, the `area` or
            `diameter` is assumed to account for the central
            obstruction.

    Raises:
        ValueError:
            Raised if both or neither of `diameter` or `area` are
            provoded.
    """
    def __init__(self, longitude, latitude, elevation, fratio, platescale, throughput=1.,
                 area=None, diameter=None, obstruction=None):
        self.longitude = longitude
        self.latitude = latitude
        self.elevation = elevation
        self.fratio = fratio
        self.platescale = platescale
        self._throughput = throughput

        # If area is provided, use it directly:
        if area is not None:
            self.area = area
            if diameter is not None:
                raise ValueError('Cannot provide area and diameter; provide one or the other.')
            self.diameter = numpy.sqrt(area/numpy.pi)*2/100
            if obstruction is not None:
                warnings.warn('Obstruction and area provided, combining to get effective area.')
        elif diameter is not None:
            self.diameter = diameter
            self.area = numpy.pi*numpy.square(self.diameter*100/2)
        else:
            raise ValueError('Must provide area or diameter!')

        # Apply the central obsruction if provided.
        if obstruction is not None:
            self.area *= (1-obstruction)

    @property
    def throughput(self): #, wave=None):
        if self._throughput is None:
            # TODO: Should this really raise an error...
            raise ValueError('Throughput not defined.')
        return self._throughput


class KeckTelescope(Telescope):
    """
    Keck telescopes on MaunaKea
    """
    def __init__(self):
        # This assumes the f/15 secondary
        fratio = 15
        platescale = 0.725
        # Assumes the Nasmyth port for the throughput; three
        # reflections off an aluminum coating.
        eta_file = str(data_file(filename='efficiency') / 'avg_keck_reflectivity.db')
        single_reflection = efficiency.Efficiency.from_file(eta_file, wave_units='nm')
        throughput = efficiency.CombinedEfficiency(dict([(key,single_reflection)
                                                         for key in ['m1', 'm2', 'm3']]))
        super().__init__(155.47833, 19.82833, 4160.0, fratio, platescale,
                         throughput=throughput, area=723674.)


class SDSSTelescope(Telescope):
    """
    SDSS 2.5-meter telescope at Apache Point Observatory.
    """
    def __init__(self):
        super().__init__(105.82028, 32.78028, 2788.0, 5., 0.06048, diameter=2.5, obstruction=0.286)


class LickNickelTelescope(Telescope):
    """
    Nickel Telescope at Lick Observatory
    """
    def __init__(self):
        fratio = 17
        platescale = 0.0815
        eta_file = str(data_file(filename='efficiency') / 'avg_keck_reflectivity.db')
        single_reflection = efficiency.Efficiency.from_file(eta_file, wave_units='nm')
        cassegrain_throughput = efficiency.CombinedEfficiency(dict([(key,single_reflection)
                                                         for key in ['m1', 'm2']]))
        super().__init__(121.64278, 37.34139, 1283.0, fratio, platescale,
                         throughput=cassegrain_throughput, diameter=0.978, obstruction=0.12)


class APFTelescope(Telescope):
    """
    Automated Planet Finder telescope at Lick Observatory.
    """
    def __init__(self):
        super().__init__(121.64278, 37.34139, 1283.0, 15., 0.17452, diameter=2.41, obstruction=0.02)


class TMTTelescope(Telescope):
    """
    Thirty-Meter Telescope, assuming it is constructed on Maunakea.
    """
    def __init__(self, reflectivity='req'):
        """
        TMT telescope parameters.

        Args:
            reflectivity (:obj:`str`, optional):
                The reflectivity curve to use. For the ORD
                requirement or goal, use ``'req'`` or ``'goal'``,
                respectively.
        """
        if reflectivity not in ['req', 'goal']:
            raise ValueError('Must set reflectivity to \'req\' or \'goal\'.')
        # Assumes the Nasymth port for the throughput; three
        # reflections off an aluminum coating.
#        eta_file = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'efficiency',
#                                'uv_enhanced_silver.db')
#        single_reflection = efficiency.Efficiency.from_file(eta_file)
        eta_file = str(data_file(filename='efficiency') / 'tmt_ord.db')
        db = numpy.genfromtxt(eta_file)
        c = 3 if reflectivity == 'req' else 4
        single_reflection = efficiency.Efficiency(db[:,c], wave=db[:,0]) 
        throughput = efficiency.CombinedEfficiency(dict([(key,single_reflection)
                                                         for key in ['m1', 'm2', 'm3']]))
        # TODO: These numbers are all temporary
        # Focal ratio is f/15
        # Area is just taken as 30m diameter with a 10% obstruction
        # Plate scale is 1e3 * 15[f/] * 30 [m] / 206265 [arcsec/rad] = 2.182 [mm/arcsec]"""
        # Location is just taken to be the same as Keck.
        super().__init__(155.47833, 19.82833, 4160.0, 15., 2.182,
                         throughput=throughput, diameter=30, obstruction=0.1)


#class Observation:
#    """Not yet implemented..."""
#    def __init__(self, airmass, sky_brightness, exposure_time, wavelength=None, band=None):
#        if wavelength is None and band is None:
#            raise ValueError('Must provide band or wavelength.')
#
#        self.airmass = airmass
#        self.sky_brightness = sky_brightness
#        self.exposure_time = exposure_time
#        self.wavelength = wavelength
#        self.band = band
#
#    @classmethod
#    def from_date_telescope_target(cls, date, telescope, target, exposure_time, wavelength=None,
#                                   band=None):
#        raise NotImplementedError('Not yet implemented!')

        
#        # Use date to get lunar phase
#        lunar_cycle = 29.53 # days
#        day_since_new_moon = int(lunar_cycle*(0.5 
#                                    - astroplan.moon_phase_angle(Time(date)).value/2/numpy.pi))
#        sky_brightness = interpolate_sky_brightness(day_since_new_moon, wavelength=wavelength,
#                                                    band=band)
#        # Use date, telescope, and target to get airmass
#        targ = SkyCoord(target)
#        tele = EarthLocation(lat=telescope.latitude*units.deg, lon=telescope.longitude*units.deg,
#                             height=telescope.elevation*units.m)
#        targaltaz = targ.transform_to(AltAz(obstime=Time(date), location=tele))
#        airmass = targaltaz.secz
#
#        return cls(airmass, sky_brightness, exposure_time, wavelength=wavelength, band=band)



