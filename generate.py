#!/usr/bin/env python

from distutils.extension import Extension
from Cython.Build import cythonize

ext = Extension(
    "graphcontraction",
    sources=["py/graphcontraction.pyx", "src/eigenmap.cpp"],
    include_dirs=["py/","src/","include/gc/"],
    language="c++"
)

ext2 = Extension(
    "libeigenmap",
    sources=["src/eigenmap.cpp", "py/libeigenmap.pyx"],
    include_dirs=["py/","src/","include/gc/"],
    language="c++"
)

cythonize([ext])

