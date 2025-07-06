#!/usr/bin/env bash

# Apply post-3.3.1 patch TRIQS/triqs@bb75cd9af40203295c3d11478f7905e69f009450
patch -d $PREFIX -p1 < "${RECIPE_DIR}/fix-triqs-mesh-std-apply-cpp23.patch"

mkdir build
cd build

export CXXFLAGS="$CXXFLAGS -D_LIBCPP_DISABLE_AVAILABILITY"
source $PREFIX/share/triqs/triqsvars.sh

# Build SOM
cmake ${CMAKE_ARGS} \
    -DCMAKE_CXX_COMPILER=${BUILD_PREFIX}/bin/$(basename ${CXX}) \
    -DCMAKE_C_COMPILER=${BUILD_PREFIX}/bin/$(basename ${CC}) \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    ..

make -j2 VERBOSE=1

# Run unit tests
PATH="$(pwd)/test/c++:$PATH" ctest --output-on-failure

# Install SOM
make install

# Set correct paths in som-targets.cmake
if [[ "$target_platform" == "osx-arm64" ]]; then
  sed "s|$BUILD_PREFIX|$PREFIX|g" ${PREFIX}/lib/cmake/som/som-targets.cmake > \
    tmp_file
  cp tmp_file ${PREFIX}/lib/cmake/som/som-targets.cmake
fi
