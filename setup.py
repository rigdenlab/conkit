"""Python Interface for Residue-Residue Contact Predictions"""
import os
import subprocess
import sys

from distutils.command.build import build
from distutils.util import convert_path
from setuptools import setup, Extension

from Cython.Build import cythonize
import numpy as np


# ==============================================================
# Setup.py command extensions
# ==============================================================


# Credits to http://stackoverflow.com/a/33181352
class BuildCommand(build):
    user_options = build.user_options + [
        ("script-python-path=", None, "Path to Python interpreter to be included in the scripts")
    ]

    def initialize_options(self):
        build.initialize_options(self)
        self.script_python_path = None

    def finalize_options(self):
        build.finalize_options(self)

    def run(self):
        global script_python_path
        script_python_path = self.script_python_path
        build.run(self)


# ==============================================================
# Functions, functions, functions ...
# ==============================================================


def dependencies():
    with open("requirements.txt", "r") as f_in:
        deps = f_in.read().splitlines()
    return deps


def extensions():
    return cythonize(
        [
            Extension(name="conkit.misc.ext.c_bandwidth", sources=["conkit/misc/ext/c_bandwidth.c"],
                      include_dirs=["conkit/misc/ext"]),
        ],
        language_level=sys.version_info[0],
    )


def readme():
    with open("README.rst", "r") as f_in:
        return f_in.read()


def scripts():
    extension = ".bat" if sys.platform.startswith("win") else ""
    header = "" if sys.platform.startswith("win") else "#!/bin/sh"
    bin_dir = "bin"
    command_dir = convert_path("conkit/command_line")
    scripts = []
    for file in os.listdir(command_dir):
        if not file.startswith("_") and file.endswith(".py"):
            # Make sure we have a workable name
            f_name = os.path.basename(file).rsplit(".", 1)[0]
            for c in [".", "_"]:
                new_f_name = f_name.replace(c, "-")
            # Write the content of the script
            script = os.path.join(bin_dir, new_f_name + extension)
            with open(script, "w") as f_out:
                f_out.write(header + os.linesep)
                # BATCH file
                if sys.platform.startswith("win"):
                    string = "@{0} -m conkit.command_line.{1} %*"
                # BASH file
                else:
                    string = '{0} -m conkit.command_line.{1} "$@"'
                f_out.write(string.format(PYTHON_EXE, f_name) + os.linesep)
            os.chmod(script, 0o777)
            scripts.append(script)
    return scripts


def version():
    # Credits to http://stackoverflow.com/a/24517154
    main_ns = {}
    ver_path = convert_path("conkit/version.py")
    with open(ver_path) as f_in:
        exec(f_in.read(), main_ns)
    return main_ns["__version__"]


# ==============================================================
# Determine the Python executable
# ==============================================================
PYTHON_EXE = None
for arg in sys.argv:
    if arg[0:20] == "--script-python-path" and len(arg) == 20:
        option, value = arg, sys.argv[sys.argv.index(arg) + 1]
        PYTHON_EXE = value
    elif arg[0:20] == "--script-python-path" and arg[20] == "=":
        option, value = arg[:20], arg[21:]
        PYTHON_EXE = value

if not PYTHON_EXE:
    PYTHON_EXE = sys.executable

# ==============================================================
# Define all the relevant options
# ==============================================================
AUTHOR = "Felix Simkovic"
AUTHOR_EMAIL = "felixsimkovic@me.com"
DESCRIPTION = __doc__.replace("\n", "")
DEPENDENCIES = dependencies()
EXT_MODULES = extensions()
LICENSE = "BSD License"
LONG_DESCRIPTION = readme()
PACKAGE_DIR = "conkit"
PACKAGE_NAME = "conkit"
PLATFORMS = ["POSIX", "Mac OS", "Windows", "Unix"]
SCRIPTS = scripts()
URL = "http://www.conkit.org/en/latest/"
VERSION = version()

PACKAGES = [
    "conkit",
    "conkit/applications",
    "conkit/command_line",
    "conkit/core",
    "conkit/core/ext",
    "conkit/io",
    "conkit/misc",
    "conkit/misc/ext",
    "conkit/plot",
]

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

TEST_REQUIREMENTS = [
    "codecov",
    "coverage",
    "pytest",
    "pytest-cov",
    "pytest-pep8",
    "pytest-helpers-namespace",
]

setup(
    cmdclass={"build": BuildCommand},
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    name=PACKAGE_NAME,
    description=DESCRIPTION,
    ext_modules=EXT_MODULES,
    include_dirs=[np.get_include()],
    long_description=LONG_DESCRIPTION,
    license=LICENSE,
    version=version(),
    url=URL,
    packages=PACKAGES,
    package_dir={PACKAGE_NAME: PACKAGE_DIR},
    scripts=SCRIPTS,
    platforms=PLATFORMS,
    classifiers=CLASSIFIERS,
    install_requires=DEPENDENCIES,
    tests_require=TEST_REQUIREMENTS,
    setup_requires=['pytest-runner'],
    include_package_data=True,
    zip_safe=False,
    package_data={"conkit.misc.ext.c_bandwidth": ["c_bandwidth.h", "c_bandwidth.pyx", "c_bandwidth.c"]}

)
