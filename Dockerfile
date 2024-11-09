FROM flatironinstitute/triqs:3.3.0 as base
LABEL maintainer="Igor Krivenko"
LABEL description="Stochastic Optimization Method for Analytic Continuation"
ARG APPNAME=som

USER root
RUN useradd -m -s /bin/bash -u 999 build && echo "build:build" | chpasswd
RUN apt-get update && \
    apt-get install -y --no-install-recommends make g++-12 apt-utils file
RUN sh -c 'echo -e "\nrmaps_base_oversubscribe = 1" >> \
          /etc/openmpi/openmpi-mca-params.conf'

COPY requirements.txt /src/$APPNAME/requirements.txt
RUN pip3 install -r /src/$APPNAME/requirements.txt

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build

ENV CC gcc-12
ENV CXX g++-12
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} \
          -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_SHARED_LIBS=ON \
          -DBuild_Documentation=OFF \
          -DBUILD_DEBIAN_PACKAGE=ON \
    && make -j4 || make -j1 VERBOSE=1 \
    && ctest --output-on-failure \
    && cpack

USER root
RUN make install \
    && mkdir -p /build/repo \
    && mv som-*.deb /build/repo
