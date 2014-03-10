# -*- coding: utf-8 -*

import ez_setup

ez_setup.use_setuptools()

import os
import imp
import sys
import glob
import tempfile
import unittest
import commands
import platform

from setuptools import setup, Extension, Distribution, find_packages

try:
    import numpy
except ImportError:
    pass

setup(
    name="interfits",
    version="0.0",
    description="Interferometric Data Interchange and Visualization",
    author="Danny Price",
    author_email="",
    long_description="Interferometric Data Interchange and Visualization",
    classifiers=['Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy'],
    packages=find_packages(),
    scripts=glob.glob('scripts/*.py'),
    setup_requires=['numpy>=1.2'],
    install_requires=['pyfits>=3.1', 'numpy>=1.6', 'lxml', 'h5py'],
    dependency_links=['http://www.stsci.edu/resources/software_hardware/pyfits/Download'],
    include_package_data=True,
    ext_package='interfits',
    zip_safe=False,
) 
