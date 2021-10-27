#!/usr/bin/python
# -*- coding: utf-8 -*-
import setuptools
from setuptools import setup, Extension

# Need to install cython and numpy for extensions to work.
# See here: https://luminousmen.com/post/resolve-cython-and-numpy-dependencies
# and here: https://github.com/pypa/pip/issues/5761
from setuptools import dist
dist.Distribution().fetch_build_eggs(['Cython>=0.15.1', 'numpy>=1.10'])

import numpy
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()

cythonize("gotoh/*.pyx", include_path=[numpy.get_include()])

setup(
    install_requires=[
        'numpy',
    ],
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        'wheel',
        'setuptools>=18.0',
        'cython',
    ],
    ext_modules=[Extension(
        'gotoh.gotoh', 
        ['gotoh/gotoh.c'],
        include_dirs=[numpy.get_include()]),
    ],
    name="gotoh",
    version="0.1.1",
    author="Robert Turnbull",
    author_email="rob@robturnbull.com",
    description="Sequence Alignment with different penalties for opening gaps and extending them.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rbturnbull/gotoh",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)