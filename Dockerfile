FROM ubuntu:latest
ARG github_token
ARG username
ARG password

# This Dockerfile aim to reproduce the setup in the Github Actions cloud
# instances as closely as possible. This is useful when trying to debug
# failures during the Actions CI runs.
#
# Note:
# In this Dockerfile we use the developer's personal GITHUB_TOKEN to access
# private repositories. For this to work the GITHUB_TOKEN must be passed to
# to Docker when Docker is invoked. To pass the GITHUB_TOKEN to the Docker
# instance invoke Docker as:
#
#    docker build . \
#           --build-arg github_token=$CPP_GITHUB_TOKEN \
#           --build-arg username=<username> \
#           --build-arg password=<password>
#
# Where $GITHUB_TOKEN returns the value of the GITHUB_TOKEN environment
# variable (which should be set to your personal GITHUB_TOKEN). The 
# GITHUB_TOKEN before the "=" is the name of the variable that can be used
# in the Dockerfile.

RUN mkdir -p /Integrals
WORKDIR /Integrals
# We will pull everything in over the network and not use local directories
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y git curl wget
#
# BEWARE: this Python stuff installs all sorts of GCC-7 crap that breaks
#         the alternatives set up. So this Python stuff needs to be installed
#         before installing GCC-8 and running the alternatives set up
#         (whether Python still works after that? Who knows?).
#
RUN apt-get install -y python3-pip python3-dev python3-virtualenv
RUN pip3 install gcovr
#
RUN apt-get install -y gcc-8 g++-8
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 95 --slave /usr/bin/g++ g++ /usr/bin/g++-8 --slave /usr/bin/gcov gcov /usr/bin/gcov-8 --slave /usr/bin/gcc-ar gcc-ar /usr/bin/gcc-ar-8 --slave /usr/bin/gcc-nm gcc-nm /usr/bin/gcc-nm-8 --slave /usr/bin/gcc-ranlib gcc-ranlib /usr/bin/gcc-ranlib-8 --slave /usr/bin/gcov-dump gcov-dump /usr/bin/gcov-dump-8 --slave /usr/bin/gcov-tool gcov-tool /usr/bin/gcov-tool-8
RUN update-alternatives --install /usr/bin/cpp cpp /usr/bin/cpp-8 95
RUN wget https://github.com/Kitware/CMake/releases/download/v3.16.3/cmake-3.16.3-Linux-x86_64.sh
RUN yes | /bin/sh cmake-3.16.3-Linux-x86_64.sh
RUN apt-get install -y vim
RUN apt-get install -y libgmp-dev
RUN apt-get install -y libblas-dev liblapack-dev
RUN apt-get install -y libboost-dev libboost-all-dev
RUN apt-get install -y libeigen3-dev
#RUN apt-get install -y openmpi-bin libopenmpi-dev
RUN apt-get install -y mpich
RUN apt-get install -y autotools-dev autoconf
RUN echo "GCC: " `g++ --version`
RUN echo "GCC: " `gcc --version`
RUN git clone https://github.com/evaleev/libint
RUN cd libint; ./autogen.sh; mkdir build; cd build; ../configure; make export
RUN tar -zxf libint/build/libint-*.tgz
RUN cd libint-*; ../cmake-3.16.3-Linux-x86_64/bin/cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=/Integrals/install -DBUILD_SHARED_LIBS=ON; cd build; make; make install
RUN git clone https://github.com/evaleev/tiledarray
RUN export CXX=`which mpicxx`; export CC=`which mpicc`; cd tiledarray; ../cmake-3.16.3-Linux-x86_64/bin/cmake -H. -Bbuild -DCMAKE_INSTALL_PREFIX=/Integrals/install -DBUILD_SHARED_LIBS=ON -DCMAKE_CXX_FLAGS="-std=c++17" -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_C_COMPILER=${CC}; cd build; make; make install
RUN git config --global credential.helper store
# The git clone command will store the Github username and password (in plain text, ***shivers***)
#RUN git clone https://$username:$password@github.com/NWChemEx-Project/Integrals
RUN git clone https://$username:$password@github.com/hjjvandam/Integrals
WORKDIR /Integrals/Integrals
RUN git branch github-actions
RUN git checkout github-actions
#RUN git pull https://github.com/NWChemEx-Project/Integrals github-actions
RUN git pull https://github.com/hjjvandam/Integrals github-actions
RUN export CXX=`which mpicxx`; export CC=`which mpicc`; /Integrals/cmake-3.16.3-Linux-x86_64/bin/cmake -H. -Bbuild -DCMAKE_PREFIX_PATH=/Integrals/install -DCMAKE_INSTALL_PREFIX=/Integrals/install -DCMAKE_BUILD_TYPE=DEBUG -DBUILD_TESTING=TRUE -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_FLAGS="-std=c++17 -O0 -fPIC --coverage" -DCMAKE_C_FLAGS="-O0 -fPIC --coverage" -DCMAKE_EXE_LINKER_FLAGS="-O0 -fPIC -fprofile-arcs" -DCPP_GITHUB_TOKEN=$github_token -DBUILD_SHARED_LIBS=ON
#RUN /Integrals/cmake-3.16.3-Linux-x86_64/bin/ctest --build-and-test build build --build-nocmake --build-generator "Unix Makefiles" --test-command make test
RUN cd build; make; make test
WORKDIR /Integrals
RUN gcovr --root Integrals --filter Integrals --exclude Integrals/tests --xml Integrals/coverage.xml
