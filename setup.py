"""Python Interface for Residue-Residue Contact Predictions"""

from setuptools import setup, find_packages
from distutils.util import convert_path
import glob

def get_version():
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path('conkit/_version.py')
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns['__version__']

# Obtain the current version of ConKit
__version__ = get_version()

# Do the actual setup below
setup(
    name='conkit',
    description=__doc__.replace("\n", ""),
    long_description=open('README.rst').read(),
    version=__version__,
    author='Felix Simkovic',
    author_email='felixsimkovic@me.com',
    license='BSD License',
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
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    test_suite='nose.collector',
    tests_require=['nose >=1.3.7'],
    include_package_data=True,
    zip_safe=False,
)

