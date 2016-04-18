# Needed to build Cython
from distutils.core import setup, Extension
from Cython.Build import cythonize

# Change to newest gcc on Darwin
import os
from sys import platform
if platform.lower() == 'darwin':
    os.environ['CC'] = 'gcc-5'

# Do the build
ext = Extension('fourier', sources=['fourier.pyx', 'tsfourier.c'])
setup(name='fourier', ext_modules=cythonize([ext]))
