#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
build_dir=${DIR}/libs_build

echo MACOSX_DEPLOYMENT_TARGET $MACOSX_DEPLOYMENT_TARGET

git clone --depth 1 --branch llvmorg-17.0.6 https://github.com/llvm/llvm-project
pushd llvm-project/openmp
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
make
sudo make install

popd
rm -rf llvm-project
