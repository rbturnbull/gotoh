#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy
import setuptools
from setuptools import setup, Extension
from Cython.Build import cythonize

with open("README.md", "r") as fh:
    long_description = fh.read()

cythonize("gotoh_counts/*.pyx", include_path=[numpy.get_include()])

setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        'setuptools>=18.0',
        'cython',
    ],
    ext_modules=[Extension(
        'gotoh_counts.gotoh_counts', 
        ['gotoh_counts/gotoh_counts.c'],
        include_dirs=[numpy.get_include()]),
        ],
    name="gotoh_counts-rbturnbull",
    version="0.0.2",
    author="Robert Turnbull",
    author_email="rob@robturnbull.com",
    description="Aligns two sequences and returns the number of characters which match, mismatch, open gaps or extend gaps.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rbturnbull/gotoh_counts",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)