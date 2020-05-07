#!/bin/env/python3
# -*- encoding utf-8 -*-
"""
Module with defined spectrographs.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os

from IPython import embed

import numpy

from . import telescopes, observe, kernel, detector, efficiency, vph

# from . import telescopes, observe, kernel, detector, efficiency, optical

#class Spectrograph:
#    """
#    Base class.
#    """
#    def __init__(self):
#        pass
#
#class MultiArmSpectrograph:
#    """
#    Base class.
#    """
#    def __init__(self):
#        pass

# Below is the 2D version in place before Apr 2020.  This needs to be updated.
#class SpectrographArm:
#    """
#    Base class for a spectrograph arm.
#    """
#    def __init__(self, platescale, cwl, dispersion, spectral_range, arm_detector, arm_kernel,
#                 throughput, opticalmodel, scramble=False):
#        self.platescale = platescale            # focal-plane plate scale in mm/arcsec
#        self.cwl = cwl                          # central wavelength in angstroms
#        self.dispersion = dispersion            # linear dispersion in A/mm
#        self.spectral_range = spectral_range    # free spectral range (Delta lambda)
#        self.detector = arm_detector            # Detector instance
#        self.kernel = arm_kernel                # Monochromatic kernel
#        self.throughput = throughput            # Should describe the throughput from the focal
#                                                # plane to the detector, including detector QE;
#                                                # see `SpectrographThroughput`
#        self.opticalmodel = opticalmodel        # Optical model used to propagate rays from the
#                                                # focal plane to the camera detector
#        self.scramble = scramble                # Does the source image get scrambled by the
#                                                # entrance aperture?
#
#    @property
#    def pixelscale(self):
#        return self.detector.pixelsize/self.platescale     # arcsec/pixel
#
#    @property
#    def dispscale(self):
#        return self.detector.pixelsize*self.dispersion      # A/pixel
#
#    def monochromatic_image(self, sky, spec_aperture, onsky_source=None):
#        """
#        Construct a monochromatic image of a source through an
#        aperture as observed by this spectrograph arm.
#        """
#        return observe.monochromatic_image(sky, spec_aperture, self.kernel, self.platescale,
#                                           self.detector.pixelsize, onsky_source=onsky_source,
#                                           scramble=self.scramble)
#
#    def twod_spectrum(self, sky_spectrum, spec_aperture, source_distribution=None,
#                      source_spectrum=None, wave_lim=None, field_coo=None, rectilinear=False):
#        optical_args = {} if rectilinear \
#                            else {'field_coo': field_coo, 'opticalmodel':self.opticalmodel}
#        return observe.twod_spectrum(sky_spectrum, spec_aperture, self.kernel, self.platescale,
#                                     self.dispersion, self.detector.pixelsize,
#                                     source_distribution=source_distribution,
#                                     source_spectrum=source_spectrum, thresh=1e-10,
#                                     wave_lim=wave_lim, **optical_args)
#
#    def observe(sky, sky_spectrum, spec_aperture, exposure_time, airmass, onsky_source=None,
#                source_spectrum=None, extraction=None):
#        """
#        Take an observation through an aperture.
#        """
#        pass
#
#
#class TMTWFOSBlue(SpectrographArm):
#    def __init__(self, setting='lowres'):
#        if setting not in TMTWFOSBlue.valid_settings():
#            raise ValueError('Setting {0} not known.'.format(setting))
#
#        # Setup the telescope
#        telescope = telescopes.TMTTelescope()
#        # Plate-scale in mm/arcsec assuming the camera fratio is 2
#        platescale = telescope.platescale * 2 / telescope.fratio
#        # Assume camera yields 0.2 arcsec FWHM in both dimensions
#        spatial_FWHM, spectral_FWHM = numpy.array([0.2, 0.2])*platescale
#        # Assign the kernel without setting the pixel sampling
#        arm_kernel = kernel.SpectrographGaussianKernel(spatial_FWHM, spectral_FWHM)
#        # The detector
#        qe_file = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'efficiency', 'detectors',
#                               'itl_sta_blue.db')
#        arm_detector = detector.Detector((4*4096, 2*4096), pixelsize=0.015, rn=2.5, dark=0.1,
#                                         qe=efficiency.Efficiency.from_file(qe_file))
#        # Focal-plane up to, but not including, detector throughput
#        throughput_file = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'efficiency', 
#                                       'wfos_throughput.db')
#        pre_detector_eta = efficiency.Efficiency.from_file(throughput_file)
#        # Total throughput
#        throughput = efficiency.SpectrographThroughput(detector=arm_detector,
#                                                       other=pre_detector_eta)
#
#        opticalmodel = TMTWFOSBlueOpticalModel(setting=setting)
#
#        if setting == 'lowres':
#            #   Central wavelength is 4350 angstroms
#            #   Linear dispersion is 13.3 angstroms per mm
#            #   Free spectral range is 2500 angstroms
#            super(TMTWFOSBlue, self).__init__(platescale, 4350., 13.3, 2500., arm_detector,
#                                              arm_kernel, throughput, opticalmodel)
#
#    @staticmethod
#    def valid_settings():
#        return ['lowres']
#
#
#class TMTWFOSRed(SpectrographArm):
#    def __init__(self, setting='lowres'):
#        if setting not in TMTWFOSRed.valid_settings():
#            raise ValueError('Setting {0} not known.'.format(setting))
#
#        # Setup the telescope
#        telescope = telescopes.TMTTelescope()
#        # Plate-scale in mm/arcsec assuming the camera fratio is 2
#        platescale = telescope.platescale * 2 / telescope.fratio
#        # Assume camera yields 0.2 arcsec FWHM in both dimensions
#        spatial_FWHM, spectral_FWHM = numpy.array([0.2, 0.2])*platescale
#        # Assign the kernel without setting the pixel sampling
#        arm_kernel = kernel.SpectrographGaussianKernel(spatial_FWHM, spectral_FWHM)
#        # The detector
#        qe_file = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'efficiency', 'detectors',
#                               'itl_sta_red.db')
#        arm_detector = detector.Detector((4*4096, 2*4096), pixelsize=0.015, rn=2.5, dark=0.1,
#                                         qe=efficiency.Efficiency.from_file(qe_file))
#        # Focal-plane up to, but not including, detector throughput
#        throughput_file = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'efficiency', 
#                                       'wfos_throughput.db')
#        pre_detector_eta = efficiency.Efficiency.from_file(throughput_file)
#        # Total throughput
#        throughput = efficiency.SpectrographThroughput(detector=arm_detector,
#                                                       other=pre_detector_eta)
#
#        opticalmodel = TMTWFOSRedOpticalModel(setting=setting)
#
#        if setting == 'lowres':
#            #   Central wavelength is 7750 angstroms
#            #   Linear dispersion is 23.7 angstroms per mm
#            #   Free spectral range is 4500 angstroms
#            super(TMTWFOSRed, self).__init__(platescale, 7750., 23.7, 4500., arm_detector,
#                                             arm_kernel, throughput, opticalmodel)
#
#    @staticmethod
#    def valid_settings():
#        return ['lowres']
#
#class TMTWFOS:  #(MultiArmSpectrograph)
#    """
#    Instantiate a setting of the WFOS spectrograph on TMT.
#    """
#    # TODO: Sky stuff should default to Maunakea and be defined by target position and moon-phase...
#    # TODO: Allow airmass to be defined by target position and UT start of observation
#    def __init__(self, setting='lowres'):
#        self.telescope = telescopes.TMTTelescope()
#        # TODO: Allow settings to be different for each arm.
#        self.arms = {'blue': TMTWFOSBlue(setting=setting), 'red': TMTWFOSRed(setting=setting)}
#
#    def monochromatic_image(self, sky, spec_aperture, onsky_source=None, arm=None):
#        """
#        Generate monochromatic images of the source through the
#        aperture in one or more of the spectrograph arms.
#        """
#        if arm is not None:
#            return self.arms[arm].monochromatic_image(sky, spec_aperture, onsky_source=onsky_source)
#        return dict([(key, a.monochromatic_image(sky, spec_aperture, onsky_source=onsky_source))
#                         for key,a in self.arms.items()])
#
#    def twod_spectrum(self, sky_spectrum, spec_aperture, source_distribution=None,
#                      source_spectrum=None, wave_lim=None, arm=None, field_coo=None,
#                      rectilinear=False):
#        """
#        Generate a 2D spectrum of the source through the aperture in
#        one or more of the spectrograph arms.
#        """
#        if arm is not None:
#            return self.arms[arm].twod_spectrum(sky_spectrum, spec_aperture,
#                                                source_distribution=source_distribution,
#                                                source_spectrum=source_spectrum, wave_lim=wave_lim,
#                                                field_coo=field_coo, rectilinear=rectilinear)
#        return dict([(key, a.twod_spectrum(sky_spectrum, spec_aperture,
#                                           source_distribution=source_distribution,
#                                            source_spectrum=source_spectrum, wave_lim=wave_lim,
#                                            field_coo=field_coo, rectilinear=rectilinear))
#                         for key,a in self.arms.items()])
#
##    def observe(self, source_distribution, source_spectrum, sky_distribution, sky_spectrum,
##                spec_aperture, airmass, exposure_time, extraction):
##        """
##        Returns the total spectrum and variance, sky spectrum
##        """
##        # TODO: sky_spectrum should default to Maunakea
#
#
#class TMTWFOSBlueOpticalModel(optical.OpticalModelInterpolator):
#    def __init__(self, setting='lowres'):
#        if setting == 'lowres':
#            modelfile = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'instr_models', 'wfos',
#                                     'Blue_Low_Spot_Data_2020_150.txt')
#        else:
#            raise NotImplementedError('Setting {0} not yet recognized.'.format(setting))
#        xf, yf, wave, xc, yc, rays = numpy.genfromtxt(modelfile).T
#        indx = rays > 0
#        # Convert field coordinates from arcmin to arcsec, wavelengths
#        # from micron to angstroms, and percentage of incident rays to
#        # a fraction
#        super(TMTWFOSBlueOpticalModel, self).__init__(60*xf[indx], 60*yf[indx], 10000*wave[indx],
#                                                      xc[indx], yc[indx], vignette=0.01*rays[indx])
#
#
#class TMTWFOSRedOpticalModel(optical.OpticalModelInterpolator):
#    def __init__(self, setting='lowres'):
#        if setting == 'lowres':
#            modelfile = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'instr_models', 'wfos',
#                                     'Red_Low_Spot_Data_2020_150.txt')
#        else:
#            raise NotImplementedError('Setting {0} not yet recognized.'.format(setting))
#        xf, yf, wave, xc, yc, rays = numpy.genfromtxt(modelfile).T
#        indx = rays > 0
#        # Convert field coordinates from arcmin to arcsec, wavelengths
#        # from micron to angstroms, and percentage of incident rays to
#        # a fraction
#        super(TMTWFOSRedOpticalModel, self).__init__(60*xf[indx], 60*yf[indx], 10000*wave[indx],
#                                                     xc[indx], yc[indx], vignette=0.01*rays[indx])

class SlitSpectrographArm:
    """
    Base class for spectrograph arms.

    Cannot be instantiated as a useful class by itself because
    constants are undefined.

    Args:
        telescope (:class:`~synospec.etc.telescopes.Telescope`):
            Spectrograph telescope properties
        grating (:class:`~synospec.etc.vph.VPHGrating`):
            The grating object.
        cen_wave (:obj:`float`, optional):
            Central wavelength in angstroms for an on-axis slit. Note
            this is mutually exclusive with ``grating_angle``. If
            this is provided, the grating angle is set to the Bragg
            angle. If None, ``grating_angle`` must be provided, and
            the central wavelength is set to the Littrow value.
        grating_angle (:obj:`float`, optional):
            The angle of the grating normal with respect to the
            incoming beam in degrees. Note this is mutually exclusive
            with ``cen_wave``. If provided, the central wavelength is
            set to the Littrow value. If None, ``cen_wave`` must be
            provided, and the grating angle is set to the Bragg
            angle.

    Raises:
        NotImplementedError:
            Raised if the derived class (see, e.g.,
            :class:`TMTWFOSBlue`) does not have :attr:`base_file`
            defined.

    """

    det = None
    """
    :class:`~synospec.etc.detector.Detector` instance with the detector
    properties.
    """

    focal_length_coll = None
    """Collimator focal length in mm"""

    focal_length_cam = None
    """Camera focal length in mm"""

#    cam_kernel = None
#    """Kernel instance with the imaging kernel."""

    max_field_angle_cam = None
    """Maximum camera field entrance angle (radius) in deg."""

    focal_plane_limits = None
    """
    Upper and lower limits on focal plane coordinates relative to the
    field center in arcsec.
    """

    base_file = None
    """
    File with the grating-independent efficiencies. Must be defined
    by the derived class.
    """

    def __init__(self, telescope, grating, cen_wave=None, grating_angle=None):
        # Check mutually exclusive input
        if cen_wave is None and grating_angle is None:
            raise ValueError('Must provide cen_wave or grating_angle.')
        if cen_wave is not None and grating_angle is not None:
            raise ValueError('Cannot provide both cen_wave or grating_angle.')

        self.telescope = telescope

        # Setup the grating
        self.grating = grating
        if cen_wave is None:
            self.grating_angle = grating_angle
            self.cen_wave = self.grating.littrow_wavelength(self.grating_angle)
        else:
            self.cen_wave = cen_wave
            self.grating_angle = self.grating.bragg_angle(self.cen_wave)

        # Get the "base-level" efficiency; i.e., everything except the
        # telescope and the grating.
        self.base_efficiency = self._base_efficiency()

        # Check the focal-plane limits
        if self.focal_plane_limits is not None and len(self.focal_plane_limits) != 4:
            raise ValueError('{0} poorly defined; '.format(self.__class__.__name__)
                             + 'focal_plane_limits must be None or a 4-element tuple.')

        # Set the plate scale of the spectrograph arm at the detector
        # in mm/arcsec
        self.platescale = self.telescope.platescale * self.focal_length_cam/self.focal_length_coll

    @property
    def pixelscale(self):
        """Pixelscale of the detector in arcsec/pixel."""
        return self.det.pixelsize/self.platescale/1e3

    def _base_efficiency(self):
        """
        Internal method to construct the baseline efficiency (i.e.,
        the combined efficiency of all spectrograph components for
        this arm, except for the grating and the telescope).
        """
        if self.base_file is None:
            raise NotImplementedError('Must define the base_file for {0}'.format(
                                        self.__class__.__name__))
        db = numpy.genfromtxt(self.base_file)
        # Wavelengths in these files are expected to be in nm.
        # Wavelengths are expected to be in the first column. All other
        # columns are expected to be efficiencies of all other
        # instrument components, the product of which provides the
        # total efficiency of everything but the grating.
        return efficiency.Efficiency(numpy.prod(db[:,1:], axis=1), wave=db[:,0]*10)

# NOT TESTED
#    def monochromatic_image(self, sky, spec_aperture, onsky_source=None):
#        """
#        Construct a monochromatic image of a source through an
#        aperture as observed by this spectrograph arm.
#
#        Args:
#            sky (:class:`synospec.etc.source.OnSkySource`):
#                Sky flux distribution.
#            spec_aperture (:class:`synospec.etc.aperture.Aperture`):
#                Spectrograph aperture. Aperture is expected to be
#                oriented with the dispersion along the first axis
#                (e.g., the slit width is along the abcissa).
#            onsky_source (:class:`synospec.etc.source.OnSkySource`, optional):
#                On-sky distribution of the source flux. If None, only
#                sky is observed through the aperture.
#
#        Returns:
#            `numpy.ndarray`_: 2D array with a representation of the
#            monochromatic image as sampled by the detector pixels.
#        """
#        return observe.monochromatic_image(sky, spec_aperture, self.cam_kernel, self.platescale,
#                                           self.det.pixelsize, onsky_source=onsky_source)

    def wavelength_limits(self, x, y, add_grating_limits=False):
        """
        Calculate the wavelength limits for a spectrum at the
        provided focal-plane position.

        Args:
            x (:obj:`float`):
                Field position along the dispersion direction in
                arcseconds relative to the field center.
            y (:obj:`float`):
                Field position along the cross-dispersed direction in
                arcseconds relative to the field center.
            add_grating_limits (:obj:`bool`, optional):
                The returned limits include any limits set by the
                grating.

        Returns:
            :obj:`tuple`: **Approximate** minimum and maximum
            wavelengths viewable by the camera. There are known
            differences between these results and those returned by a
            Zemax model of the instrument.
        """
        # Pupil Space Field Angles
        theta_x = self.slit_angle(x)
        theta_y = self.slit_angle(y)

        if theta_x > self.max_field_angle_cam or theta_y > self.max_field_angle_cam:
            raise ValueError('Slit not viewable by camera!')

        # Input angle at this X position
        alpha = self.grating_angle - theta_x
        # Assume the camera angle is set by the output angle for the
        # central wavelength for an on-axis slit.
        beta = self.grating.diffraction_angle(self.cen_wave, alpha=self.grating_angle)
        # Maximum camera input angle at this theta_y (relative to the
        # field x center)
        cam_angle_max = numpy.sqrt(numpy.square(self.max_field_angle_cam) - numpy.square(theta_y))

        wave_min = self.grating.diffracted_wave(beta-cam_angle_max, alpha=alpha)
        wave_max = self.grating.diffracted_wave(beta+cam_angle_max, alpha=alpha)

        if not add_grating_limits:
            return wave_min, wave_max

        if self.grating.wave_lim[0] is not None and wave_min < self.grating.wave_lim[0]:
            wave_min = self.grating.wave_lim[0]
        if self.grating.wave_lim[1] is not None and wave_max > self.grating.wave_lim[1]:
            wave_max = self.grating.wave_lim[1]
        return wave_min, wave_max

    def efficiency(self, wave, x=0.0, y=0.0, same_type=True):
        r"""
        Compute the total efficiency for this arm, from the ADC
        through to (and including) the detector.

        .. note::

            - Note that x and y must be scalars for now.

        Args:
            wave (:obj:`float`, `numpy.ndarray`_):
                One or more wavelengths (angstroms in air) at which
                to sample the spectrograph efficiency.
            x (:obj:`float`, optional):
                Field position along the dispersion direction in
                arcseconds relative to the field center.
            y (:obj:`float`, optional):
                Field position along the cross-dispersed direction in
                arcseconds relative to the field center.
            same_type (:obj:`bool`, optional):
                Return the efficiency values relevant with the same
                type as the provided ``wave``. If False, ``wave``
                must be an array and an :class:`Efficiency` object is
                returned.

        Returns:
            :obj:`float`, `numpy.ndarray`_, :class:`Efficiency`: The
            efficiency of the spectrograph arm at one or more
            wavelengths. The returned type matches the input type of
            ``wave``.

        Raises:
            TypeError:
                Raised if ``same_type`` is False and ``wave`` is not
                an array.
            ValueError:
                Raised if ``x`` or ``y`` are outside the field limits
                of the instrument.
        """
        if not same_type and not isinstance(wave, (list, numpy.ndarray)):
            raise TypeError('To return an Efficiency object, must provide a wavelength array.')
        if self.focal_plane_limits is not None:
            if  x < self.focal_plane_limits[0] or x > self.focal_plane_limits[1]:
                raise ValueError('Slit x position must be within {0:.1f} -- {1:.1f}!'.format(
                                   self.focal_plane_limits[0], self.focal_plane_limits[1]))
            if  y < self.focal_plane_limits[2] or y > self.focal_plane_limits[3]:
                raise ValueError('Slit y position must be within {0:.1f} -- {1:.1f}!'.format(
                                   self.focal_plane_limits[2], self.focal_plane_limits[3]))

        # Calculate the efficiencies
        alpha = self.slit_angle(x) + self.grating_angle
        eta = self.base_efficiency(wave) * self.grating.efficiency(wave, alpha=alpha)

        # If a masked array is returned, fill the masked components
        # (outside the defined grating wavelength limits) with 0s
        if isinstance(eta, numpy.ma.MaskedArray):
            eta = eta.filled(0.0)

        # Truncate the spectrum to the limits of the camera field-of-view
        wave_min, wave_max = self.wavelength_limits(x, y)
        indx = wave < wave_min
        eta[indx] = 0.0
        indx = wave > wave_max
        eta[indx] = 0.0
        return eta if same_type else efficiency.Efficiency(eta, wave=wave)

    def slit_angle(self, x):
        """
        Entrance angle relative to the field center for a slit at
        position ``x``.

        Args:
            x (:obj:`float`):
                Field position along the dispersion direction in
                arcseconds relative to the field center.

        Returns:
            :obj:`float`: Slit entrance angle (should be added to the
            grating angle to get the grating incident angle).
        """
        return numpy.degrees(numpy.arctan(x*self.telescope.platescale/self.focal_length_coll))

    def linear_dispersion(self, wave, x=0., m=1):
        """
        Return the linear dispersion in mm per angstrom.

        Args:
            wave (:obj:`float`, array-like):
                Wavelengths at which to calculate the linear
                dispersion.
            x (:obj:`float`, optional):
                Field position along the dispersion direction in
                arcseconds relative to the field center.
            m (:obj:`int`, optional):
                Grating order.

        Returns:
            :obj:`float`, `numpy.ndarray`_: Linear dispersion at each
            wavelength position in mm/angstrom.
        """
        alpha = self.slit_angle(x) + self.grating_angle
        return self.focal_length_cam*self.grating.angular_dispersion(wave, alpha=alpha, m=m)

    def resolution(self, wave, x=0, slit_width=0.75, m=1):
        r"""
        Return the spectral resolution (:math:`R \equiv
        \lambda/\Delta\lambda`) for a given slit width.

        .. warning::

            Assumes perfect image quality and a uniformly illuminated
            aperture.

        Args:
            wave (:obj:`float`, array-like):
                Wavelengths at which to calculate the spectral
                resolution.
            x (:obj:`float`, optional):
                Field position along the dispersion direction in
                arcseconds relative to the field center.
            slit_width (:obj:`float`, optional):
                The width of the slit in arcsec along the dispersion
                direction (i.e., should account for any rotation of
                the slit).
            m (:obj:`int`, optional):
                Grating order.

        Returns:
            :obj:`float`, `numpy.ndarray`_: The spectral resolution,
            :math:`R \equiv \lambda/\delta\lambda`, at each
            wavelength position.
        """
        return wave / self.resolution_element(slit_width=slit_width, units='angstrom', wave=wave,
                                              x=x, m=m)

    def resolution_element(self, slit_width=0.75, units='micron', wave=None, x=0, m=1):
        """
        Return the width of the resolution element (monochromatic
        image) at the camera focal plane.

        Units can be micron, pixels, or angstroms. If in angstroms,
        wave must be provided.

        .. warning::

            Assumes perfect image quality and a uniformly illuminated
            aperture.

        Args:
            slit_width (:obj:`float`, optional):
                The width of the slit in arcsec along the dispersion
                direction (i.e., should account for any rotation of
                the slit).
            units (:obj:`str`, optional):
                Units for the resolution element. Must be
                ``'micron'``, ``'pixels'``, or ``'angstrom'``. Note
                that ``wave`` is ignored (and the returned value is a
                single float) unless angstrom units are requested,
                when ``wave`` is required.
            wave (:obj:`float`, array-like, optional):
                Wavelengths at which to calculate the spectral
                resolution. Ignored unless units are angstroms.
            x (:obj:`float`, optional):
                Field position along the dispersion direction in
                arcseconds relative to the field center. Ignored
                unless units are angstroms.
            m (:obj:`int`, optional):
                Grating order. Ignored unless units are angstroms.

        Returns:
            :obj:`float`, `numpy.ndarray`_: The width of the spectral
            resolution element in the requested units. If ``units ==
            'angstroms``, the returned type matches ``wave``;
            otherwise, a single float is returned.
        
        """
        if units == 'angstrom' and wave is None:
            raise ValueError('For resolution element in angstrom units, must at least provide '
                             'wavelength.')

        # Slit width in micron at the detector
        width = self.platescale * 1e3 * slit_width

        # Return with the correct units
        if units == 'micron':
            return width
        if units == 'pixels':
            return width / self.det.pixelsize
        if units == 'angstrom':
            return width * 1e-3 / self.linear_dispersion(wave, x=x, m=m)

        raise ValueError('Unknown units: {0}'.format(units))

class WFOSGrating(vph.VPHGrating):
    """
    :class:`~synospec.etc.vph.VPHGrating` instances available in the
    current WFOS design.

    To list the available gratings::

        from synospec.etc.spectrographs import WFOSGrating
        print(WFOSGrating.available_gratings.keys())

    Args:
        grating (:obj:`str`):
            Grating name. Must be one of the gratings with parameters
            defined by :attr:`available_gratings`.
    """
    available_gratings \
            = dict(B1210=dict(line_density=1210., n_bulk=1.35, n_mod=0.05, thickness=4.0,
                              wave_lim=(3100., 5600.)),
                   B2479=dict(line_density=2479., n_bulk=1.35, n_mod=0.13, thickness=1.742,
                              wave_lim=(3100., 5600.)),
                   B2700=dict(line_density=2700., n_bulk=1.17, n_mod=0.17, thickness=3.50,
                              wave_lim=(3100., 5600.)),
                   B3600=dict(line_density=3600., n_bulk=1.17, n_mod=0.15, thickness=3.03,
                              wave_lim=(3100., 5600.)),
                   R680=dict(line_density=680., n_bulk=1.35, n_mod=0.07, thickness=5.35,
                              wave_lim=(5400., 10000.)),
                   R1392=dict(line_density=1392., n_bulk=1.35, n_mod=0.14, thickness=2.85,
                              wave_lim=(5400., 10000.)),
                   R1520=dict(line_density=1520., n_bulk=1.17, n_mod=0.23, thickness=4.67,
                              wave_lim=(5400., 10000.)),
                   R2052=dict(line_density=2052., n_bulk=1.17, n_mod=0.20, thickness=4.01,
                              wave_lim=(5400., 10000.)))

    def __init__(self, grating):
        if grating is None:
            raise ValueError('Grating ID must be defined.')
        if grating not in self.available_gratings.keys():
            raise ValueError('Unknown grating ID ({0}).  Options are: {1}'.format(grating,
                                ', '.join(list(self.available_gratings.keys()))))
        super(WFOSGrating, self).__init__(
                    self.available_gratings[grating]['line_density'],
                    self.available_gratings[grating]['n_bulk'],
                    self.available_gratings[grating]['n_mod'],
                    self.available_gratings[grating]['thickness'],
                    wave_lim=self.available_gratings[grating]['wave_lim'])
        self.name = grating


class TMTWFOSArm(SlitSpectrographArm):
    """
    Instance of :class:`SlitSpectrographArm` specific to TMT-WFOS.

    Args:
        reflectivity (:obj:`str`, optional):
            The reflectivity curve to use. For the ORD requirement or
            goal, use ``'req'`` or ``'goal'``, respectively.
        grating (:obj:`str`, optional):
            Grating name. Must be one of the gratings with parameters
            defined by :attr:`available_gratings`. If None, use
            :attr:`default_grating` defined by each arm.
        cen_wave (:obj:`float`, optional):
            Central wavelength in angstroms for an on-axis slit. Note
            this is mutually exclusive with ``grating_angle``. If
            this is provided, the grating angle is set to the Bragg
            angle. If None and ``grating_angle`` is also None, the
            value is set to the wavelength with the peak efficiency
            for the selected grating. If None and ``grating_angle``
            is provided, the central wavelength is set to the Littrow
            value for the provided grating angle.
        grating_angle (:obj:`float`, optional):
            The angle of the grating normal with respect to the
            incoming beam in degrees. Note this is mutually exclusive
            with ``cen_wave``. If provided, the central wavelength is
            set to the Littrow value. If None, the grating angle is
            set to the Bragg angle for the central wavelength.
    """
    focal_length_coll = 4500.
    focal_plane_limits = [-90, 90, -249, 249]
    default_grating = None
    def __init__(self, reflectivity='req', grating=None, cen_wave=None, grating_angle=None):
        telescope = telescopes.TMTTelescope(reflectivity=reflectivity)
        grating = WFOSGrating(self.default_grating if grating is None else grating)
        if cen_wave is None and grating_angle is None:
            cen_wave = grating.peak_wave
        super(TMTWFOSArm, self).__init__(telescope, grating, cen_wave=cen_wave,
                                         grating_angle=grating_angle)


class TMTWFOSBlue(TMTWFOSArm):
    """
    Object used to compute the efficiency of the blue arm of the
    TMT-WFOS spectrograph. See the base class for the description of
    the instantiation arguments and methods.
    """
    det = detector.Detector((4*4096, 2*4096), pixelsize=15., rn=2., dark=0.0)
#    cam_kernel = kernel.SpectrographGaussianKernel(*list(numpy.array([0.2, 0.2])
#                                                            * TMTWFOSArm.telescope.platescale))
    focal_length_cam = 600.
    max_field_angle_cam = 12.

    base_file = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'efficiency',
                             'wfos_blue_efficiency.db')
    """
    The file containing tabulated efficiency data for the elements
    of the spectrograph arm, excluding the grating.
    """

    default_grating = 'B1210'
    """
    The default grating used when one is not specified at
    instantiation.
    """


class TMTWFOSRed(TMTWFOSArm):
    """
    Object used to compute the efficiency of the red arm of the
    TMT-WFOS spectrograph. See the base class for the description of
    the instantiation arguments and methods.
    """

    det = detector.Detector((4*4096, 2*4096), pixelsize=15., rn=2., dark=0.0)
#    cam_kernel = kernel.SpectrographGaussianKernel(*list(numpy.array([0.2, 0.2])
#                                                            * TMTWFOSArm.telescope.platescale))
    focal_length_cam = 600.
    max_field_angle_cam = 12.

    base_file = os.path.join(os.environ['SYNOSPEC_DIR'], 'data', 'efficiency',
                             'wfos_red_efficiency.db')
    """
    The file containing tabulated efficiency data for the elements
    of the spectrograph arm, excluding the grating.
    """

    default_grating = 'R680'
    """
    The default grating used when one is not specified at
    instantiation.
    """
        
#class TMTWFOS:
#    """
#    Specific setup for both arms of TMT WFOS
#    """
#    def __init__(self, blue_cen, red_cen, blue_grating=None, blue_grating_angle=None,
#                 red_grating=None, red_grating_angle=None):
#
#        self.arms = dict(blue=TMTWFOSBlue(blue_cen, grating=blue_grating,
#                                          grating_angle=blue_grating_angle),
#                         red=TMTWFOSRed(red_cen, grating=red_grating,
#                                        grating_angle=red_grating_angle))
#
#    def efficiency(self, wave, x=0.0, y=0.0, same_type=True, combine_arms=True):
#        blue_eta = self.arms['blue'].efficiency(wave, x=x, y=y)
#        red_eta = self.arms['red'].efficiency(wave, x=x, y=y)
#        if combine_arms:
#            eta = blue_eta + red_eta
#            return efficiency.Efficiency(eta, wave=wave) if to_efficiency else eta
#
#        return dict(blue=efficiency.Efficiency(blue_eta, wave=wave) if to_efficiency else blue_eta,
#                    red=efficiency.Efficiency(red_eta, wave=wave) if to_efficiency else red_eta)
#
#    def wavelengths(self):
#        """
#        Construct the wavelength vector 
#
#    def resolution(self, wave):
#
#        
#    def efficiency():

