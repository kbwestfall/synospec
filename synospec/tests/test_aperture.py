
from IPython import embed

import numpy

from synospec.etc.source import OnSkyMoffat
from synospec.etc import aperture


def test_hex():

    ap = aperture.HexagonalAperture(0., 0., 0.8, incircle=True, orientation='vertical')
    
    fwhm = 0.7
    beta = 3.5
    sampling = 0.03
    size = 3.

    mof = OnSkyMoffat(fwhm, beta, sampling=sampling, size=size)
    assert ap.integrate_over_source(mof) < 0.5, 'Integral changed'


def test_fiber():
    fwhm = 0.7
    beta = 3.5
    sampling = 0.03
    size = 3.

    mof = OnSkyMoffat(fwhm, beta, sampling=sampling, size=size)
    dy = numpy.diff(mof.y)
    assert numpy.allclose(dy, dy[0]), 'Sampling should be uniform'
    dy = dy[0]
    assert numpy.allclose(numpy.absolute(numpy.diff(mof.x)), dy), \
            'Sampling should be the same in x and y'
    
    fiber = aperture.FiberAperture(0., 0., 1.)
    r = fiber.response(mof.x, mof.y)
    assert numpy.isclose(fiber.area, numpy.sum(r)*dy*dy), 'Area mismatch'


test_fiber()

