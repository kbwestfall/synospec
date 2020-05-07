"""
Module with the optical model interpolation class.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy
from scipy import interpolate

from matplotlib import pyplot

from .sampling import Resample

class OpticalModelInterpolator:
    """
    Class used to interpolate between the focal plane and the
    detector positions based on an archived result from a run of an
    optical model (e.g., using Zemax).

    The arguments of the instantiation should be the known values
    from the optical model used to setup the interpolation object.

    Args:
        xf (`numpy.ndarray`_):
            On-sky focal-plane position in arcseconds, along the
            spatial axis, of the incoming rays relative to the field
            center.
        yf (`numpy.ndarray`_):
            On-sky focal-plane position in arcseconds, along the
            dispersive axis, of the incoming rays relative to the
            field center.
        wave (`numpy.ndarray`_):
            Wavelength of the incoming rays in angstroms.
        xc (`numpy.ndarray`_):
            The camera field position in mm along the spatial axis
            relative to the field center where the incoming rays
            land.
        yc (`numpy.ndarray`_):
            The camera field position in mm along the dispersive axis
            relative to the field center where the incoming rays
            land.
        vignette (`numpy.ndarray`_):
            The vignetting of the incoming rays as a fraction of the
            rays that make it to the designated position. The
            vignetting is important because it sets if the
            interpolation should be masked in specific regions of the
            grid space. If the vignetting is 0., the output is
            masked.

    Attributes:
        field (`numpy.ndarray`_):
            Field information of the input grid: the x and y
            positions and wavelengths.
        camera (`numpy.ndarray`_):
            Camera information of the input grid: the x and y
            positions and vignetting.

    """
    def __init__(self, xf, yf, wave, xc, yc, vignette):

        # Check the input
        _xf = numpy.atleast_1d(xf)
        if _xf.ndim != 1:
            raise ValueError('Input field x coordinates must be a 1D vector.')
        _yf = numpy.atleast_1d(yf)
        if _yf.shape != _xf.shape:
            raise ValueError('Input field y coordinates do not match field x coordinates.')
        _wave = numpy.atleast_1d(wave)
        if _wave.shape != _xf.shape:
            raise ValueError('Input wavelengths do not match field x coordinates.')
        _xc = numpy.atleast_1d(xc)
        if _xc.shape != _xf.shape:
            raise ValueError('Input camera x coordinates do not match field x coordinates.')
        _yc = numpy.atleast_1d(yc)
        if _yc.shape != _xf.shape:
            raise ValueError('Input camera y coordinates do not match field x coordinates.')
        _vignette = numpy.atleast_1d(vignette)
        if _vignette.shape != _xf.shape:
            raise ValueError('Input vignetting do not match field x coordinates.')
        if numpy.any(_vignette > 1):
            raise ValueError('Vignetting is an efficiency and cannot be greater than 1.')

        self.field = numpy.hstack((_xf.reshape(-1,1), _yf.reshape(-1,1), _wave.reshape(-1,1)))

        self.gridxf = numpy.unique(_xf)
        self.gridyf = numpy.unique(_yf)
        self.gridwv = numpy.unique(_wave)

        self.shape = (len(self.gridxf), len(self.gridyf), len(self.gridwv))

        self.xc_interp = interpolate.RegularGridInterpolator((self.gridxf, self.gridyf, self.gridl),
                                                             _xc.reshape(self.shape))
        self.yc_interp = interpolate.RegularGridInterpolator((self.gridxf, self.gridyf, self.gridl),
                                                             _yc.reshape(self.shape))
        self.vg_interp = interpolate.RegularGridInterpolator((self.gridxf, self.gridyf, self.gridl),
                                                             _vignette.reshape(self.shape))
        self.mk_interp = interpolate.RegularGridInterpolator((self.gridxf, self.gridyf, self.gridl),
                                                (_vignette > 0).astype(float).reshape(self.shape))

    def field2camera(self, x, y, wave):
        """
        Use the interpolator to propagate field points to the camera
        focal plane.

        Args:
            x (`numpy.ndarray`_):
                On-sky focal-plane position in arcseconds, along the
                spatial axis, of the incoming rays relative to the
                field center.
            y (`numpy.ndarray`_):
                On-sky focal-plane position in arcseconds, along the
                dispersive axis, of the incoming rays relative to the
                field center.
            wave (`numpy.ndarray`_):
                Wavelength of the incoming rays in angstroms.

        Returns:
            tuple: Returns four `numpy.ndarray`_ objects with the
            camera focal-plane coordinates, vignetting, and a
            bad-calculation mask.
        """
        _x = numpy.atleast_1d(x)
        if _x.ndim > 1:
            raise ValueError('Input x coordinates must be a 1D vector.')
        _y = numpy.atleast_1d(y)
        if _y.shape != _x.shape:
            raise ValueError('Input y coordinates must match input x coordinates.')
        _wave = numpy.atleast_1d(wave)
        if _wave.shape != _x.shape:
            raise ValueError('Input wavelengths must match input x coordinates.')
        samp = numpy.hstack((_x.reshape(-1,1), _y.reshape(-1,1), _wave.reshape(-1,1)))
        return self.xc_interp(samp), self.yc_interp(samp), self.vg_interp(samp), \
                self.mk_interp(samp) < 1.

# TODO: This requires some more thought. Can we have Zemax files that
# run in the opposite direction?
#    def camera2field(self, xc, yc, yf):
#        """
#        Use the interpolator to propagate camera focal-plane
#        coordinates to field x and wavelength for a given focal-plane
#        coordinate along the dispersive axis.
#
#        Args:
#            xc (`numpy.ndarray`_):
#                The camera field position in mm along the spatial
#                axis relative to the field center where the incoming
#                rays land.
#            yc (`numpy.ndarray`_):
#                The camera field position in mm along the dispersive
#                axis relative to the field center where the incoming
#                rays land.
#            yf (`numpy.ndarray`_):
#                On-sky focal-plane position in arcseconds, along the
#                dispersive axis, of the incoming rays relative to the
#                field center.
#
#        Returns:
#            tuple: Returns two `numpy.ndarray`_ objects with the
#            on-sky focal-plane position in arcminutes, along the
#            spatial axis, of the incoming rays relative to the field
#            center and their wavelength.
#        """
#        _xc = numpy.atleast_1d(xc)
#        if _xc.ndim > 1:
#            raise ValueError('Input camera x coordinates must be a 1D vector.')
#        _yc = numpy.atleast_1d(yc)
#        if _yc.shape != _xc.shape:
#            raise ValueError('Input camera y coordinates must match camera x coordinates.')
#        _yf = numpy.atleast_1d(yf)
#        if _yf.shape != _x.shape:
#            raise ValueError('Input field y coordinates must match camera x coordinates.')
#
#        if self._camera2field is None:    
#            self._camera2field = interpolate.LinearNDInterpolator(numpy.array([self.camera[:,0],
#                                                                               self.camera[:,1],
#                                                                               self.field[:,1]]).T,
#                                                                  self.field[:,[0,2]],
#                                                                  fill_value=self.fill_value,
#                                                                  rescale=True)
#        result = self._camera2field(numpy.array([_xc, _yc, _yf]).T)
#        return (result[:,0], result[:,1]) if hasattr(x, '__len__') \
#                    else (result[0,0], result[0,1])

    def project_2d_spectrum(self, spec, platescale, linear_dispersion, pixelsize, wave0,
                            field_coo=None):
        """
        Project a 2D spectrum from the field onto the detector.

        shape of spec should be (nspec,nspat)

        TODO: Allow field_coo to be a 2D array for multiple field positions?

        Args:
            pixelsize (:obj:`float`):
                Pixel size in mm.
            field_coo (array-like):
                Field coordinates for the center of the entrance
                aperture in arcseconds relative to the field center.

        Returns:
            tuple: Returns the 2D spectrum projected from the field
            to the detector, the first pixel for the first axis of
            the image (spectral) and the first pixel for the second
            axis of the image (spatial)
        """
        _field_coo = numpy.zeros(2, dtype=float) if field_coo is None \
                        else numpy.atleast_1d(field_coo)
        if _field_coo.ndim > 1:
            raise NotImplementedError('Can only do one coordinate.')
        if _field_coo.size != 2:
            raise ValueError('Field coordinates must be x and y relative to the field center.')

        nspec, nspat = spec.shape

        # Pixels scale in arcsec/pixel
        pixelscale = pixelsize/platescale
        # Dispersion scale in A/pixel
        dispscale = linear_dispersion * pixelsize

        # Generate the coordinate grid for the input.
        #   - Spatial slit coordinates relative to the field center in arcsec
        x = (numpy.arange(nspat) - nspat//2)*pixelscale + _field_coo[0]
        #   - Wavelengths in angstroms
        w = wave0 + numpy.arange(nspec)*dispscale
        x,w = map(lambda x: x.ravel(), numpy.meshgrid(x,w))
        #   - Single coordinate along the dispersing axis relative to the field center in arcsec
        y = numpy.full(x.size, _field_coo[1], dtype=float)

        # Coordinates are returned in mm
        xc, yc, eta, bpm = map(lambda x: x.reshape(nspec,nspat), self.field2camera(x, y, w))
        # Convert mm to pixels
        xc /= pixelsize
        yc /= pixelsize
        # Vignette the spectrum
        _spec = spec*eta
        # Get the pixels with good values
        gpm = numpy.invert(bpm)

        # Re-orient so that pixel coordinates increase with map pixel
        # number
        if numpy.median(numpy.diff(xc[0,:])) < 0:
            vig_spec = vig_spec[:,::-1]
            xc = xc[:,::-1]
            yc = yc[:,::-1]
        if numpy.median(numpy.diff(yc[:,0])) < 0:
            vig_spec = vig_spec[::-1,:]
            xc = xc[::-1,:]
            yc = yc[::-1,:]

        # Resample the spatial dimension at each wavelength
        newx = numpy.arange(numpy.floor(numpy.amin(xc[gpm])), numpy.ceil(numpy.amax(xc[gpm]))+1.,
                            1.)
        xrng = newx[[0,-1]]
        _yc = numpy.zeros((nspec,newx.size), dtype=float)
        _gpm = numpy.zeros((nspec,newx.size), dtype=float)
        _spec = numpy.zeros((nspec,newx.size), dtype=float)
        for i in range(nspec):
            _yc[i,:] = interpolate.interp1d(xc[i,:], yc[i,:], fill_value='extrapolate',
                                            bounds_error=False)(newx)
            _gpm[i,:] = interpolate.interp1d(xc[i,:], gpm[i,:].astype(float),
                                     fill_value='extrapolate', bounds_error=False)(newx) > 0.8
            _spec[i,:] = Resample(vig_spec[i,:], x=xc[i,:], mask=_gpm[i,:], newRange=xrng,
                                  newdx=1., newLog=False).outy

        # Resample the spectral dimension at each spatial coordinate
        newy = numpy.arange(numpy.floor(numpy.amin(_yc)), numpy.ceil(numpy.amax(_yc))+1., 1.)
        yrng = newy[[0,-1]]
        proj_spec = numpy.zeros((newy.size,newx.size), dtype=float)
        proj_spec_frac = numpy.zeros((newy.size,newx.size), dtype=float)
        for i in range(newx.size):
            r = Resample(_spec[:,i], x=_yc[:,i], mask=_gpm[:,i], newRange=yrng, newdx=1.,
                         newLog=False)
            proj_spec[:,i] = r.outy
            proj_spec_frac[:,i] = r.outf

        return proj_spec, proj_spec_frac > 0.5, newy[0], newx[0]
