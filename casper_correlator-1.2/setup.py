from distutils.core import setup, Extension

import os, glob, numpy, sys

__version__ = '1.4.0'

def indir(dir, files): return [dir+f for f in files]
def globdir(dir, files): 
    rv = []
    for f in files: rv += glob.glob(dir+f)
    return rv

setup(name = 'casper_correlator',
    version = __version__,
    description = 'Interface to CASPER correlators',
    long_description = 'Interface to CASPER correlators.',
    license = 'GPL',
    author = 'Aaron Parsons',
    author_email = 'aparsons at astron.berkeley.edu',
    url = '',
    package_dir = {'casper_correlator':'src'},
    packages = ['casper_correlator'],
    ext_modules = [
        Extension('casper_correlator.rx',
            globdir('src/rx/',
                ['*.cpp','*.c']),
            include_dirs = [numpy.get_include(), 'src/rx/include'],
            libraries=['rt'],
        )
    ],
    scripts=glob.glob('scripts/*'),
)

