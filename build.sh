#!/bin/bash
mkdir Binaries
mkdir -p Build
cd Build
rm -rf *
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build . --config debug
