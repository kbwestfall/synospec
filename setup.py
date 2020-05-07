#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# K. Westfall, 22 May 2018
#   Adapted from Marvin's setup.py file

# Imports
import os
import glob
from setuptools import setup, find_packages

import requests
import warnings

VERSION = '0.1.0dev'
RELEASE = 'dev' not in VERSION
MINIMUM_PYTHON_VERSION = '3.5'

def get_package_data(root='synospec'):
    """Generate the list of package data."""
    return [os.path.relpath(f, root) 
                    for f in glob.glob(os.path.join(root, 'data/*/*'))] \
           + [os.path.relpath(f, root) 
                    for f in glob.glob(os.path.join(root, 'data/README.md'))]

#def get_data_files():
#    """Generate the list of data files."""
#    data_files = []
#    data_roots = [ 'data' ]
#    for root in data_roots:
#        for path, directories, files in os.walk(root):
#            for f in files:
#                data_path = '/'.join(path.split('/')[1:])
#                data_files.append(os.path.join(data_path, f))
#    return data_files


def get_scripts():
    """ Grab all the scripts in the bin directory.  """
    scripts = []
    if os.path.isdir('bin'):
        scripts = [ fname for fname in glob.glob(os.path.join('bin', '*')) ]
    return scripts


def get_requirements():
    """ Get the package requirements from the system file. """
    requirements_file = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']
    return install_requires


def run_setup(package_data, scripts, packages, install_requires):

    setup(name='synospec',
          version=VERSION,
          license='BSD3',
          description='Synthetic Observatory Spectrum',
          long_description=open('README.md').read(),
          author='Kyle B. Westfall',
          author_email='westfall@ucolick.org',
          keywords='astronomy spectroscopy',
          url='https://github.com/kbwestfall/synospec',
          python_requires='>='+MINIMUM_PYTHON_VERSION,
          packages=packages,
          package_dir={'synospec': 'synospec'},
          package_data={'synospec': package_data},
          include_package_data=True,
          install_requires=install_requires,
          scripts=scripts,
          setup_requires=[ 'pytest-runner' ],
          tests_require=[ 'pytest' ],
          classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3.7',
              'Programming Language :: Python :: 3 :: Only',
              'Topic :: Documentation :: Sphinx',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Topic :: Software Development :: Libraries :: Python Modules'
          ])

#def default_paths():
#    return { 'MANGADRP_VER': _MANGADRP_VER,
#             'MANGA_SPECTRO_REDUX': os.path.join(os.environ['HOME'], 'MaNGA', 'redux'),
#             'MANGADAP_VER': _VERSION,
#             'MANGA_SPECTRO_ANALYSIS': os.path.join(os.environ['HOME'], 'MaNGA', 'analysis'),
#             'MANGACORE_VER': _MANGACORE_VER,
#             'MANGACORE_DIR': os.path.join(os.environ['HOME'], 'MaNGA', 'core', _MANGACORE_VER)
#           }
#
#def check_environment():
#    ev = default_paths()
#    for k in ev.keys():
#        if k not in os.environ:
#            warnings.warn('{0} environmental variable undefined.  Default is: {1}'.format(k,ev[k]))
#
#
#def short_warning(message, category, filename, lineno, file=None, line=None):
#    """
#    Return the format for a short warning message.
#    """
#    return ' %s: %s\n' % (category.__name__, message)


#-----------------------------------------------------------------------
if __name__ == '__main__':

#    warnings.formatwarning = short_warning

    # Compile the data files to include
#    data_files = get_data_files()

    # Get the package data (data inside the main product root)
    package_data = get_package_data()

    # Compile the scripts in the bin/ directory
    scripts = get_scripts()

    # Get the packages to include
    packages = find_packages()

    # Collate the dependencies based on the system text file
    install_requires = get_requirements()

    # Run setup from setuptools
#    run_setup(data_files, scripts, packages, install_requires)
    run_setup(package_data, scripts, packages, install_requires)

#    # Check if the environmental variables are found and warn the user
#    # of their defaults
#    check_environment()

