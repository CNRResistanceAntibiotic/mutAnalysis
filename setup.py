#!/usr/bin/python3

# setup.py file install libraries wanted for readmapper
# launch it with pip install -e .

# Install setuptools if not already present.
from setuptools import setup, find_packages
import glob

from scripts.mutAnalysis import version

setup(
    name='mutAnalysis',
    version=version(),
    description='mutAnalysis: count and report specific mutation',
    packages=find_packages(),
    author='Aurélien Birer',
    author_email='abirer@chu-clermontferrand.fr',
    url='https://github.com/CNRResistanceAntibiotic/mutAnalysis',
    scripts=glob.glob('scripts/*'),
    install_requires=[],
    license='GPLv3',
    classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],

)
