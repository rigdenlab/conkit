"""Python Interface for Residue-Residue Contact Predictions"""

from setuptools import setup, find_packages
from distutils.util import convert_path


def get_version():
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path('conkit/_version.py')
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns['__version__']


def readme():
    with open('README.rst', 'r') as f_in:
        return f_in.read()


# Obtain the current version of ConKit
__version__ = get_version()

# Do the actual setup below
setup(
    name='conkit',
    description=__doc__.replace("\n", ""),
    long_description=readme(),
    version=__version__,
    author='Felix Simkovic',
    author_email='felixsimkovic@me.com',
    license='BSD License',
    url='https://github.com/rigdenlab/conkit',
    package_dir={'conkit': 'conkit'},
    packages=find_packages(exclude="tests"),
    scripts=[
        'bin/conkit-plot', 'bin/conkit-msatool',
        'bin/conkit-predict', 'bin/conkit-precision',
        'bin/conkit-convert',
    ],
    platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
    install_requires=['numpy >=1.8.2', 'scipy >=0.16.0', 'biopython >=1.64', 'matplotlib >=1.3.1'],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    test_suite='nose.collector',
    tests_require=['nose >=1.3.7'],
    include_package_data=True,
    zip_safe=False,
)

