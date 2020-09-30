#!/usr/bin/env python3
"""
The version is stored here in a separate file so it can exist in only one place.
https://stackoverflow.com/a/7071358/2438989
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/CNRResistanceAntibiotic/mutAnalysis

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""

# Install setuptools if not already present.
from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = ""
exec(open('mutanalysis/version.py').read())

setup(name='mutanalysis',
      version=__version__,
      description='mutanalysis: mutation analysis tool',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/CNRResistanceAntibiotic/mutAnalysis',
      author='Aurélien Birer',
      author_email='abirer36@gmail.com',
      license='GPLv3',
      packages=["mutanalysis"],
      package_data={'mutanalysis': ['database/mutations.tsv', 'database/sequences.fasta']},
      include_package_data=True,
      install_requires=['biopython', 'pandas'],
      entry_points={"console_scripts": ['mutanalysis = mutanalysis.mutAnalysis:run']},
      zip_safe=False,
      python_requires='>=3.6')
