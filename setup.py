from setuptools import setup, find_packages
import sys

name = 'pyoma'

req_packages = ['numpy', 'tables>=3.2', 'future', 'familyanalyzer']
if sys.version_info < (3, 3):
    req_packages.extend(['mock', 'functools32'])

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
    package_data={'pyoma': ['browser/*.drw']},
    install_requires=req_packages,
)
