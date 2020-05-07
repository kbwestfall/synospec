#!/bin/env/python3
# -*- encoding: utf-8 -*-

import os

__version__ = '0.1.0dev'
__license__ = 'BSD3'
__author__ = 'Kyle B. Westfall'
__maintainer__ = 'Kyle B. Westfall'
__email__ = 'westfall@ucolick.org'
__copyright__ = '(c) 2020, Kyle B. Westfall'

def source_dir():
    """Return the root path to the DAP source directory."""
    import pkg_resources
    data_dir = pkg_resources.resource_filename('synospec', 'data')
    return os.path.split(data_dir) [0]

os.environ['SYNOSPEC_DIR'] = source_dir()

def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
    """
    return ' %s: %s\n' % (category.__name__, message)

import warnings
warnings.formatwarning = short_warning
