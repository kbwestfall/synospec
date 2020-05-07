
import pytest
import warnings
warnings.simplefilter('ignore', UserWarning)
warnings.simplefilter('ignore', RuntimeWarning)

import numpy

from synospec.etc import spectrum 

def test_convert_flux():
    fnu = spectrum.convert_flux_density(5000., 1.)
    assert numpy.isclose(1.0, spectrum.convert_flux_density(5000., fnu, density='Hz')), \
            'Revert flux density failed'

