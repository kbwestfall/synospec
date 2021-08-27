
from IPython import embed

from synospec.etc.source import OnSkyMoffat, OnSkyGaussian

def test_moffat():

    gau = OnSkyGaussian(0.7, sampling=0.05, size=5.)

    mof = OnSkyMoffat(0.7, 3.5, sampling=0.05, size=5.)

    from matplotlib import pyplot

    embed()
    exit()

if __name__ == '__main__':
    test_moffat()
