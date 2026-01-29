if [[ "$(uname)" == "Darwin" ]]; then
	export CC="$(brew --prefix llvm)/bin/clang"
	export CXX="$(brew --prefix llvm)/bin/clang++"
fi

rm -rf build
mkdir build
cd build
cmake .. -DMONAD_USE_OPENMP=ON -DMONAD_BUILD_APPS=ON -DMONAD_BUILD_TESTS=ON
make -j8
