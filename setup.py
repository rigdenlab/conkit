"""Python Interface for Residue-Residue Contact Predictions"""

from setuptools import setup
from distutils.util import convert_path

import os
import sys

# Determine the Python executable
if "--script-python-path" in sys.argv:
    PYTHON_EXE = sys.argv[sys.argv.index("--script-python-path") + 1]
    sys.argv.pop(sys.argv.index("--script-python-path") + 1)
    sys.argv.remove("--script-python-path")
else:
    PYTHON_EXE = sys.executable


def dependencies():
    with open('requirements.txt', 'r') as f_in:
        return [l for l in f_in.read().rsplit(os.linesep) 
                if l and not l.startswith("#")]


def readme():
    with open('README.rst', 'r') as f_in:
        return f_in.read()


def scripts():
    extension = '.bat' if sys.platform.startswith('win') else ''
    header = '' if sys.platform.startswith('win') else '#!/bin/bash'
    bin_dir = 'bin'
    command_dir = convert_path('conkit/command_line')
    scripts = []
    for file in os.listdir(command_dir):
        if not file.startswith('_') and file.endswith('.py'):
            # Make sure we have a workable name
            f_name = os.path.basename(file).rsplit('.', 1)[0]
            for c in ['.', '_']:
                new_f_name = f_name.replace(c, '-')
            # Write the content of the script
            script = os.path.join('bin', new_f_name + extension)
            with open(script, "w") as f_out:
                f_out.write(header + os.linesep)
                # BATCH file
                if sys.platform.startswith('win'):
                    string = "@{0} -m conkit.command_line.{1} %*"
                # BASH file
                else:
                    string = "{0} -m conkit.command_line.{1} \"$@\""
                f_out.write(string.format(PYTHON_EXE, f_name) + os.linesep)
            os.chmod(script, 0o777)
            scripts.append(script)
    return scripts


def version():
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path('conkit/_version.py')
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns['__version__']


AUTHOR = "Felix Simkovic"
AUTHOR_EMAIL = "felixsimkovic@me.com"
DESCRIPTION = __doc__.replace("\n", "")
DEPENDENCIES = dependencies()
LICENSE = "BSD License"
LONG_DESCRIPTION = readme()
PACKAGE_DIR = "conkit"
PACKAGE_NAME = "conkit"
PLATFORMS = ['Linux', 'Mac OS-X', 'Unix', 'Windows']
SCRIPTS = scripts()
URL = "http://www.conkit.org/en/latest/"
VERSION = version()

PACKAGES = [
    'conkit', 
    'conkit/applications',
    'conkit/command_line',
    'conkit/core', 
    'conkit/io',
    'conkit/plot', 
]

CLASSIFIERS = [
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
]


# Do the actual setup below
setup(
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    name=PACKAGE_NAME,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    license=LICENSE,
    version=VERSION,
    url=URL,
    packages=PACKAGES,
    package_dir={PACKAGE_NAME: PACKAGE_DIR},
    install_requires=DEPENDENCIES,
    scripts=SCRIPTS,
    platforms=PLATFORMS,
    classifiers=CLASSIFIERS,
    test_suite='nose.collector',
    tests_require=['nose >=1.3.7'],
    include_package_data=True,
    zip_safe=False,
)

