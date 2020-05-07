"""
Detector class

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import os
import numpy

from .efficiency import Efficiency

class Detector(Efficiency):
    """
    Define the detector statistics.

    Args:
        shape (:obj:`tuple`, optional):
            Dimensions of the detector in number of pixels along the
            spectral axis and number of pixels along the spatial
            axis. Can be None, but limits use if it is.
        pixelsize (:obj:`float`, optional):
            The size of the (square) detector pixels in *micron*.
        rn (:obj:`float`, optional):
            Read-noise in electrons.
        dark (:obj:`float`, optional):
            Dark current in electrons per second.
        gain (:obj:`float`, optional):
            Gain of detector amplifier in e- per ADU.
        fullwell (:obj:`float`, optional):
            The full well of the pixels in e-.
        nonlinear (:obj:`float`, optional):
            The fraction of the fullwell above which the detector
            response is nonlinear.
        qe (:obj:`float`, :class:`Efficiency`, optional):
            Detector quantum efficiency.
    """
    # TODO: Allow for multiple amplifiers per detector? Would also need
    # to define amplifier section.
    # TODO: Define overscan and data sections
    def __init__(self, shape=None, pixelsize=15., rn=1., dark=0., gain=1., fullwell=1e4,
                 nonlinear=1., qe=0.9):
        if shape is not None and len(shape) != 2:
            raise ValueError('Shape must contain two integers.')
        self.shape = shape
        self.pixelsize = pixelsize
        self.rn = rn
        self.dark = dark
        self.gain = gain
        self.fullwell = fullwell
        self.nonlinear = nonlinear
        if not isinstance(qe, (Efficiency, float)):
            raise TypeError('Provided quantum efficiency must be type `float` or `Efficiency`.')
        if isinstance(qe, float):
            super(Detector, self).__init__(qe)
        else:
            super(Detector, self).__init__(qe.eta, wave=qe.wave)
