ARG PARENT_IMAGE_NAME

FROM ${PARENT_IMAGE_NAME}:latest

ARG LIBINT_VERSION

RUN export INSTALL_PATH=`pwd`/install \
    && wget https://github.com/evaleev/libint/releases/download/v${LIBINT_VERSION}/libint-${LIBINT_VERSION}.tgz \
    && tar -zxf libint-${LIBINT_VERSION}.tgz \
    && cd libint-${LIBINT_VERSION} \
    && cmake -GNinja -H. -Bbuild -DCMAKE_INSTALL_PREFIX=${INSTALL_PATH} -DBUILD_SHARED_LIBS=ON \
    && cmake --build build --target install
