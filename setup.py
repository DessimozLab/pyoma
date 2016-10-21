from setuptools import setup, find_packages
import os
import shutil

name = 'pyoma'

# Create standalone2hdf
shutil.copyfile('bin/importdata.py', 'bin/oma2hdf')

setup(
    name=name,
    version='0.3.1-dev',
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
    scripts=['bin/oma2hdf'],
    data_files=[('pyoma', ['pyoma/browser/convert.drw'])],
    install_requires=['numpy', 'tables>=3.2', 'future'],
)

# Remove it
os.remove('bin/oma2hdf')
