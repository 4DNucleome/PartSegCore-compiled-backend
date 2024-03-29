#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
build_dir=${DIR}/libs_build

git clone --depth 1 --branch llvmorg-17.0.6 https://github.com/llvm/llvm-project
pushd llvm-project/openmp
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
make
make install

popd
rm -rf llvm-project
