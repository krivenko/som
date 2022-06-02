FROM flatironinstitute/triqs:3.1.0 as base
LABEL maintainer="Igor Krivenko <igor.s.krivenko@gmail.com>"
LABEL description="Stochastic Optimization Method for Analytic Continuation"
ARG APPNAME=som

USER root
RUN useradd -m -s /bin/bash -u 999 build && echo "build:build" | chpasswd
RUN apt-get update && \
    apt-get install -y --no-install-recommends make g++-10 apt-utils

COPY requirements.txt /src/$APPNAME/requirements.txt
RUN pip3 install -r /src/$APPNAME/requirements.txt

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build

ENV CC gcc-10
ENV CXX g++-10
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} \
          -DCMAKE_BUILD_TYPE=Release \
          -DBUILD_SHARED_LIBS=ON \
          -DBuild_Documentation=OFF \
          -DBUILD_DEBIAN_PACKAGE=ON \
    && make -j4 || make -j1 VERBOSE=1 \
    && cpack

USER root
RUN make install \
    && mkdir -p /build/repo \
    && mv som-*.deb /build/repo
