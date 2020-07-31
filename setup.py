from setuptools import setup, find_packages
import os
import shutil

name = "pyoma"

req_packages = [
    "numpy >= 1.16",
    "tables >= 3.5.1",
    "future",
    "fuzzyset >= 0.0.17",
    "tqdm",
    "pyopa >= 0.8",
    "pandas >= 0.22",
    'biopython >= 1.76 ; python_version >= "3.6"',
    'biopython == 1.76 ; python_version < "3.6"',
    "datasketch",
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
    author="Adrian Altenhoff",
    author_email="adrian.altenhoff@inf.ethz.ch",
    description="python library to interact and build OMA hdf5 files",
    long_description=long_description,
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    scripts=["bin/importdata.py", "bin/oma2hdf"],
    package_data={"pyoma": ["browser/*.drw"]},
    install_requires=req_packages,
    extras_require={
        "create_db": [
            "PySAIS",
            "familyanalyzer >= 0.6.0",
            "matplotlib",
            "scikit-learn",
            "scikit-fuzzy",
            "lark-parser",
            "pyham",
        ]
    },
    dependency_links=[
        "git+ssh://gitolite@lab.dessimoz.org:2222/family-analyzer@master#egg=familyanalyzer"
    ],
    python_requires=">=3.5",
)

# Remove local copy of oma2hdf (installed)
os.remove("bin/oma2hdf")
