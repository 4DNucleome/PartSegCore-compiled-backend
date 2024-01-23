# PartSegCore_compiled_backend

This is package with all cython/c++ backends for ParSeg

Current release extracted from PartSegCore module to avoid build multiple wheels
for main package and speedup tests.

This package requires `libomp` for build.
On linux it can be installed with `apt install libomp-dev` or `yum install libomp-devel`.
On mac it can be installed with `brew install libomp`.

As `clang` version created on macos do not have native openmp `brew` do
not set `libomp` as default visible.
If you install `libomp` with `brew` and want to build package
from source you need to set `OMP` variable to `libomp` location.

For example:

```bash
export OMP="$(brew --prefix libomp)"
pip install .
```

As currently there are no wheel for macos arm
it is required to build package from source.
