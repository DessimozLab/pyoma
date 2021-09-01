from setuptools import setup, find_packages
import os
import shutil

name = "pyoma"

req_packages = [
    "numpy >= 1.16",
    "tables >= 3.5.1",
    "future",
    "fuzzyset2 >= 0.1.1",
    "tqdm",
    "pyopa >= 0.8",
    "pandas >= 0.22",
    'biopython >= 1.76 ; python_version >= "3.6"',
    'biopython == 1.76 ; python_version < "3.6"',
    "datasketch",
    "ete3",
    "networkx",
]

cur_dir = os.path.abspath(os.path.dirname(__file__))

# Create oma2hdf to install
shutil.copyfile("bin/importdata.py", "bin/oma2hdf")

__version__ = "Undefined"
for line in open("{}/__init__.py".format(name)):
    if line.startswith("__version__"):
        exec(line.strip())

with open(os.path.join(cur_dir, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name=name,
    version=__version__,
    author="DessimozLab",
    author_email="contact@omabrowser.org",
    description="library to interact and build OMA hdf5 files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    licence="MPL 2.0",
    packages=find_packages(exclude=("tests",)),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    scripts=["bin/importdata.py", "bin/oma2hdf", "bin/map_to_closest_seq"],
    package_data={"pyoma": ["browser/*.drw"]},
    install_requires=req_packages,
    extras_require={
        "create_db": [
            "PySAIS",
            "familyanalyzer>=0.7.3",
            "matplotlib",
            "scikit-learn",
            "scikit-fuzzy",
            "pebble",
            "lark-parser",
            "pyham",
        ],
        "docs": ["sphinx"],
    },
    test_require=["nose"],
    python_requires=">=3.6",
)

# Remove local copy of oma2hdf (installed)
os.remove("bin/oma2hdf")
