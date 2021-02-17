mkdir Binaries
mkdir Build
pushd Build
rm -rf *
cmake -G "Visual Studio 15 2017" ..
cmake --build . --config debug
popd