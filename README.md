<!--
  ~ Copyright 2022 NWChemEx-Project
  ~
  ~ Licensed under the Apache License, Version 2.0 (the "License");
  ~ you may not use this file except in compliance with the License.
  ~ You may obtain a copy of the License at
  ~
  ~ http://www.apache.org/licenses/LICENSE-2.0
  ~
  ~ Unless required by applicable law or agreed to in writing, software
  ~ distributed under the License is distributed on an "AS IS" BASIS,
  ~ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ~ See the License for the specific language governing permissions and
  ~ limitations under the License.
-->

Integrals
===========

[![Nightly Workflow](https://github.com/NWChemEx/Integrals/actions/workflows/nightly.yaml/badge.svg)](https://github.com/NWChemEx/Integrals/actions/workflows/nightly.yaml)

Repository for the generic integral API and implementations for specific
integral libraries. Integrals is a PluginPlay plugin; its modules are written
in C++ and are also usable from Python.

Installing a Pre-Built Wheel
----------------------------

If you just want to use Integrals from Python, install the published wheel.
The wheels are built for Python 3.12:

```.sh
python3.12 -m venv venv
source venv/bin/activate
pip install nwchemex-integrals
```

The rest of this file describes building Integrals from source.

Prerequisites
-------------

- A C++20 compiler. CI uses GCC 14 or Clang 18 on Linux and AppleClang on
  macOS. On macOS do **not** use GCC: the published wheels of Integrals'
  dependencies are built with AppleClang (libc++), and GCC (libstdc++) cannot
  link against them. The build stops with an error if you try.
- CMake 3.14 or newer. Ninja is optional but recommended.
- Python 3 with its development headers, and pip. Use Python 3.12 if you want
  the build to reuse the published wheels of the dependencies (see
  [How dependencies are found](#how-dependencies-are-found)).
- MPI, Boost (headers only) and Eigen3. The build can not obtain these for you.

On Ubuntu:

```.sh
sudo apt-get install g++-14 cmake ninja-build libopenmpi-dev openmpi-bin \
                     libboost-dev libeigen3-dev
```

On macOS, with Homebrew:

```.sh
brew install cmake ninja open-mpi boost eigen
```

How dependencies are found
--------------------------

Integrals depends on [SimDE](https://github.com/NWChemEx/SimDE) (and through it
the rest of the NWChemEx stack: PluginPlay, Chemist, TensorWrapper,
ParallelZone, Utilities, ...) and on [Libint](https://github.com/evaleev/libint).
The build system comes from [NWXCMake](https://github.com/NWChemEx/NWXCMake).
You do not need to install any of these yourself:

- **NWXCMake** is taken from the `nwchemex-nwxcmake` pip package if it is
  installed in the Python CMake finds; otherwise it is downloaded from GitHub.
- **NWChemEx dependencies** (SimDE and everything below it) are taken from
  their `nwchemex-*` wheels if they are installed in the active Python
  environment. Otherwise they are downloaded from the `master` branch on GitHub
  and built as part of this build, which takes considerably longer.
- **Libint** is taken from an existing installation if CMake can find one
  (`find_package(Libint2)`, e.g., via `CMAKE_PREFIX_PATH`). Otherwise the
  v2.11.0 release is downloaded and built, which takes a while.

The quickest builds come from a Python 3.12 virtual environment that has the
dependency wheels installed, which is what the instructions below set up.

Building with pip
-----------------

This builds Integrals and installs it, including the Python module, into the
active virtual environment:

```.sh
git clone https://github.com/NWChemEx/Integrals.git
cd Integrals
python3.12 -m venv venv
source venv/bin/activate
pip install scikit-build-core setuptools-scm nwchemex-nwxcmake nwchemex-simde
pip install . --no-build-isolation
```

The third `pip install` puts Integrals' build requirements (listed under
`[build-system]` in `pyproject.toml`) into the environment. Together with
`--no-build-isolation`, this makes Integrals link against the dependency wheels
that stay installed, rather than copies in a temporary build environment that
pip deletes afterwards.

For development, do an editable install instead:

```.sh
pip install -e ".[dev]" --no-build-isolation
```

This turns on `DEVELOPER_SETUP` (see the table below), so it also builds the
tests, uses a Debug build and installs the pre-commit hooks. The build tree is
`build/`.

Building with CMake
-------------------

Integrals can also be built with CMake directly. Activating the virtual
environment from the previous section first is optional, but if the dependency
wheels are installed there, CMake reuses them:

```.sh
git clone https://github.com/NWChemEx/Integrals.git
cd Integrals
cmake -B build -G Ninja \
      -DPython_EXECUTABLE="$(command -v python3)" \
      -DBUILD_TESTING=ON \
      -DCMAKE_INSTALL_PREFIX=<where/you/want/to/install/Integrals>
cmake --build build --parallel
cmake --install build
```

Setting `Python_EXECUTABLE` makes sure CMake uses the Python of your active
environment (and therefore its wheels) rather than some other Python it finds
first. Commonly used options are:

| Option                               | Default | Description                                                                                         |
|--------------------------------------|---------|-----------------------------------------------------------------------------------------------------|
| `BUILD_TESTING`                      | `OFF`   | Build the C++ unit tests and register the Python tests with CTest.                                  |
| `BUILD_PYBIND11_BINDINGS`            | `ON`    | Build the `integrals` Python module.                                                                |
| `ENABLE_SIGMA`                       | `OFF`   | Use [Sigma](https://github.com/QCUncertainty/sigma) for uncertainty quantification (UQ modules).   |
| `INTEGRALS_ENABLE_OPENMP`            | `OFF`   | Thread the integral evaluation with OpenMP. On macOS, install `libomp` and add its prefix to `CMAKE_PREFIX_PATH`. |
| `DEVELOPER_SETUP`                    | `OFF`   | Debug build, turns on `BUILD_TESTING`, and installs the pre-commit hooks.                           |
| `ENABLE_XHOST`                       | `ON`    | Tune dependencies for this machine's CPU. Turn off if the build must run on other machines.          |
| `NWX_ECOSYSTEM_FROM_SOURCE`          | `OFF`   | Ignore installed `nwchemex-*` wheels and build every NWChemEx dependency from GitHub.                |
| `FETCHCONTENT_SOURCE_DIR_<NAME>`     | unset   | Use a local checkout of dependency `<NAME>` (e.g. `FETCHCONTENT_SOURCE_DIR_SIMDE=/path/to/SimDE`).  |

If you change which Python or which dependency wheels the build should use,
delete `build/` and configure again; earlier choices are cached.

Running the tests
-----------------

After a CMake build with `-DBUILD_TESTING=ON` (or an editable pip install),
run the C++ and Python tests with:

```.sh
ctest --test-dir build --output-on-failure
```

After a pip install, the Python tests can also be run directly:

```.sh
pip install pytest
pytest tests/python
```

Examples
--------

The [examples](examples) directory contains small programs showing how to use
Integrals. `examples/primitive_error_models` is built along with the tests;
the Python examples have their own `README.md`.
