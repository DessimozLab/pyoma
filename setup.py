from setuptools import setup, find_packages

name = 'pyoma'

setup(
    name=name,
    version='0.3.0',
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
    data_files=[('pyoma', ['pyoma/browser/convert.drw'])],
    install_requires=['numpy', 'tables>=3.2', 'future'],
)
