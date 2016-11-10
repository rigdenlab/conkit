"""Python Interface for Residue-Residue Contact Predictions"""

from conkit._version import __version__
from setuptools import setup, find_packages
import glob


setup(
    name='conkit',
    description=__doc__.replace("\n", ""),
    long_description=open('README.rst').read()
    version=__version__,
    author='Felix Simkovic',
    author_email='felixsimkovic@me.com',
    license='GNU General Public License',
    url='https://github.com/fsimkovic/conkit',
    download_url='https://github.com/fsimkovic/conkit/tarball/' + __version__,
    package_dir={'conkit': 'conkit'},
    packages=find_packages(exclude="tests"),
    scripts=[script for script in glob.glob("scripts/*")],
    platforms=['Linux', 'Mac OS-X', 'Unix'],
    install_requires=['biopython >=1.64', 'matplotlib >=1.3.1', 'numpy >=1.8.2', 'scipy >=0.16.0'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    test_suite='nose.collector',
    tests_require=['nose >=1.3.7'],
    include_package_data=True,
    zip_safe=False,
)

