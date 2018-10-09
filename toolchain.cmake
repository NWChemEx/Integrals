#Pass to CMake with variable CMAKE_TOOLCHAIN_FILE
set(PYTHON_EXECUTABLE /opt/rh/rh-python36/root/usr/bin/python3.6)
set(CMAKE_CXX_COMPILER /opt/rh/devtoolset-7/root/usr/bin/g++)
set(CMAKE_C_COMPILER /opt/rh/devtoolset-7/root/usr/bin/gcc)
list(APPEND CMAKE_PREFIX_PATH /home/jboschen/from_scratch/install)
list(APPEND CMAKE_PREFIX_PATH /usr/local/libint/2.4.2)
#list(APPEND CMAKE_PREFIX_PATH /usr/local/libint/2.5.0-beta.2)
