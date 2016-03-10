from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
from os import sys, path
import os,shutil,re
from glob import glob
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']


from imp import find_module
try: find_module('numpy')
except: sys.exit('### Error: python module numpy not found')
    
try: find_module('pyfits')
except: sys.exit('### Error: python module pyfits not found')

try: find_module('pyraf')
except: sys.exit('### Error: python module pyraf not found')

try: find_module('matplotlib')
except: sys.exit('### Error: python module matplotlib not found')

try: find_module('scipy')
except: sys.exit('### Error: python module matplotlib not found')


setup(
    name='s3',
    version='1.1.0',
    author='C.Inserra',
    author_email='c.inserra@qub.ac.uk',
    classifiers=[
        # How mature is this project?
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Stable',
        'Intended Audience :: General users',
        'Topic :: Astronomy :: photometric corrections',
        'Programming Language :: Python :: 2.7',
    ],
    scripts=['bin/STHREE','bin/SNAKE','bin/SNAKELOOP','bin/SMS','bin/SNAP'],
    url='https://github.com/cinserra',
    license=open('LICENSE.rst').read(),
    description='S3 is a package for K and P-correction and synthetic mags',
    long_description=open('README.rst').read(),
    keywords='K-correction P-correction magnitudes',
    install_requires = ['numpy','pyfits','pyraf','matplotlib','scipy'],
    packages=['s3'],
    package_dir={'':'src'},
    package_data = {'s3' : ["metadata/*.txt","metadata/NTT/*.txt","metadata/NOT/*.txt",\
                                    "metadata/PS1/*.txt","metadata/ASIAGO/*.txt","metadata/LCOGT/*.txt",\
                                    "metadata/SKYMAPPER/*.txt","metadata/LT/*.txt","metadata/LSQ/*.txt",\
                                    "metadata/WHT/*.txt","metadata/OGLE/*.txt","metadata/VLT/*.txt"]}
    #entry_points={
    #    'console_scripts': [
    #        'sample=sample:main',
    #    ],
    #}
)

