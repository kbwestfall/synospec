
from IPython import embed

import numpy

from synospec.etc.source import OnSkyMoffat, OnSkyGaussian

def test_moffat():

    fwhm = 0.7
    beta = 3.5
    sampling = 0.05
    size = 5.

    gau = OnSkyGaussian(fwhm, sampling=sampling, size=size)
    mean_gau_x = numpy.sum(gau.data*gau.X)/numpy.sum(gau.data)
    mean_gau_y = numpy.sum(gau.data*gau.Y)/numpy.sum(gau.data)
    assert numpy.allclose([mean_gau_x, mean_gau_y], [0,0]), 'Gaussian not centered'

    peak = numpy.amax(gau.data)
    indx = (gau.data > 0.47*peak) & (gau.data < 0.53*peak)
    gau_fwhm = numpy.mean(2*numpy.sqrt(gau.X[indx]**2 + gau.Y[indx]**2))
    assert numpy.absolute(fwhm - gau_fwhm) < sampling, 'Bad FWHM'

    mof = OnSkyMoffat(fwhm, beta, sampling=sampling, size=size)
    mean_mof_x = numpy.sum(mof.data*mof.X)/numpy.sum(mof.data)
    mean_mof_y = numpy.sum(mof.data*mof.Y)/numpy.sum(mof.data)
    assert numpy.allclose([mean_mof_x, mean_mof_y], [0,0]), 'Moffat not centered'

    peak = numpy.amax(mof.data)
    indx = (mof.data > 0.47*peak) & (mof.data < 0.53*peak)
    mof_fwhm = numpy.mean(2*numpy.sqrt(mof.X[indx]**2 + mof.Y[indx]**2))
    assert numpy.absolute(fwhm - mof_fwhm) < sampling, 'Bad FWHM'

    assert numpy.array_equal(gau.x, mof.x), 'Coordinate arrays are different'
    assert numpy.isclose(gau_fwhm, mof_fwhm), 'Gaussian and Moffat FWHM too different'


