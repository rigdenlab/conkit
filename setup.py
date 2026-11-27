"""Cython extension definitions — all other metadata lives in pyproject.toml"""
from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy as np

extensions = [
    Extension(
        ext.replace('/', '.').rsplit('.', 1)[0],
        [ext],
        include_dirs=[np.get_include()],
    )
    for ext in [
        "conkit/core/ext/c_contactmap.pyx",
        "conkit/core/ext/c_sequencefile.pyx",
        "conkit/misc/ext/c_bandwidth.pyx",
    ]
]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=extensions,
)
