from setuptools import setup, find_packages
import sys
import os
import shutil

name = 'pyoma'

req_packages = ['numpy>=1.16', 'tables>=3.5.1', 'future', 'fuzzyset>=0.0.17',
                'tqdm', 'pyopa>=0.8', 'pandas>=0.21', 'biopython']
if sys.version_info < (3, 3):
    req_packages.extend(['mock', 'functools32'])

# Create oma2hdf to install
shutil.copyfile('bin/importdata.py', 'bin/oma2hdf')

setup(
    name=name,
    version='0.7.0',
    author='Adrian Altenhoff',
    author_email='adrian.altenhoff@inf.ethz.ch',
    description='todoc',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    scripts=['bin/importdata.py', 'bin/oma2hdf'],
    package_data={'pyoma': ['browser/*.drw']},
    install_requires=req_packages,
    extras_require={
        'create_db': ['PySAIS', 'familyanalyzer', 'matplotlib', 'scikit-learn',
                      'scikit-fuzzy', 'lark-parser'],
    },
)

# Remove local copy of oma2hdf (installed)
os.remove('bin/oma2hdf')
