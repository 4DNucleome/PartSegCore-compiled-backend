# PartSegCore_compiled_backend

![Contributions](https://img.shields.io/badge/Contributions-Welcome-brightgreen.svg)
[![Wheels](https://github.com/4DNucleome/PartSegCore-compiled-backend/actions/workflows/wheels.yml/badge.svg?branch=master)](https://github.com/4DNucleome/PartSegCore-compiled-backend/actions/workflows/wheels.yml)
[![PyPI version](https://badge.fury.io/py/PartSegCore_compiled_backend.svg)](https://pypi.org/project/PartSegCore-compiled-backend/)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/partsegcore-compiled-backend/badges/version.svg)](https://anaconda.org/conda-forge/partsegcore-compiled-backend)
[![Python Version](https://img.shields.io/pypi/pyversions/partsegcore-compiled-backend.svg)](https://pypi.org/project/partsegcore-compiled-backend)
[![Documentation Status](https://readthedocs.org/projects/partsegcore-compiled-backend/badge/?version=stable)](https://partsegcore-compiled-backend.readthedocs.io/en/stable/)
[![Publication DOI](https://img.shields.io/badge/Publication%20DOI-10.1186%2Fs12859--021--03984--1-blue)](https://doi.org/10.1186/s12859-021-03984-1)

This is a package with all cython/c++ backend for ParSeg

Current release extracted from PartSegCore module to avoid building multiple wheels
for the main package and speedup tests.

This package requires `libomp` for build.
On linux it can be installed with `apt install libomp-dev` or `yum install libomp-devel`.
On macOS it can be installed with `brew install libomp`.

As `clang` version created on macOS do not have native openmp `brew` do
not set `libomp` as default visible.
If you install `libomp` with `brew` and want to build package
from source you need to set `OMP` variable to `libomp` location.

For example:

```bash
export OMP="$(brew --prefix libomp)"
pip install .
```

As currently there is no wheel for macOS ARM64
it is required to build package from source.

## Cite as

Bokota, G., Sroka, J., Basu, S. et al. PartSeg: a tool for quantitative feature extraction
from 3D microscopy images for dummies. BMC Bioinformatics 22, 72 (2021).
[https://doi.org/10.1186/s12859-021-03984-1](https://doi.org/10.1186/s12859-021-03984-1)
