from setuptools import setup, find_packages

name = 'pyoma'

setup(
    name=name,
    version='0.0.1',
    author='Adrian Altenhoff',
    author_email='adrian.altenhoff@inf.ethz.ch',
    description='todoc',
    packages=find_packages(),
    data_files=[('pyoma', ['pyoma/convert.drw'])],
    install_requires=['numpy', 'numexpr', 'cython', 'tables>=3.1'],
)
