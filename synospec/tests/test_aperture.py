
from IPython import embed

import numpy

from synospec.etc.source import OnSkyMoffat
from synospec.etc.aperture import HexagonalAperture

def test_hex():

    ap = HexagonalAperture(0., 0., 0.8, incircle=True, orientation='vertical')
    
    fwhm = 0.7
    beta = 3.5
    sampling = 0.03
    size = 3.

    mof = OnSkyMoffat(fwhm, beta, sampling=sampling, size=size)
    assert ap.integrate_over_source(mof) < 0.5, 'Integral changed'


