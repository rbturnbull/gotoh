# adapted from https://blagovdaryu.hashnode.dev/tremendously-speed-up-python-code-with-cython-and-package-it-with-poetry
import os
import shutil
from distutils.core import Distribution, Extension

from Cython.Build import build_ext, cythonize
import numpy as np


cython_dir = "gotoh"
extension = Extension(
    "gotoh.gotoh",
    [
        os.path.join("gotoh", "gotoh.pyx"),
    ],
    include_dirs= [np.get_include()],
)

ext_modules = cythonize([extension], include_path=[cython_dir])
dist = Distribution({"ext_modules": ext_modules})
cmd = build_ext(dist)
cmd.ensure_finalized()
cmd.run()

for output in cmd.get_outputs():
    relative_extension = os.path.relpath(output, cmd.build_lib)
    shutil.copyfile(output, relative_extension)
