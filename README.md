[![Actions](https://github.com/hjjvandam/Integrals/workflows/C_C++_CI/badge.svg)](https://github.com/NWChemEx-Project/Integrals)

[![Codecov](https://codecov.io/github/hjjvandam/Integrals/branch/master/graphs/sunburst.svg?token=PO7KR1oCuH)](https://codecov.io/github/NWChemEx-Project/Integrals/branch/master)

Integrals
===========

Repository for the generic integral API and implementations for specific integral libraries.

Building Integrals
------------------

Integrals is built using 
[CPP](https://github.com/CMakePackagingProject/CMakePackagingProject.git).
Assuming you have already installed CPP, that you are on a sane Unix-like 
computer, and you are willing to let Integrals build all dependencies, then 
the following will suffice to build Integrals:

```
git clone https://github.com/NWChemEx-Project/Integrals.git
cd SDE
cmake -H. -Bbuild -DCPP_GITHUB_TOKEN=<your super-secret token> \
                  -DMPI_ROOT=<path/to/your/mpi/installation> \
                  -DCMAKE_PREFIX_PATH=<where/you/installed/CPP> \                  
                  -DCMAKE_INSTALL_PREFIX=<where/you/want/to/install/Integrals>
cd build
cmake --build .
#May need to run as an admin depending on where you are installing
cmake --build . --target install  
```
The build process is not capable of building MPI so you will have to provide a
path to a known installation via the `MPI_ROOT` variable. The GitHub token is
necessary because, at the moment, LibChemist, TAMM, SDE, and Utilities are 
private repositories (instructions for generating a token are 
[here](https://help.github.com/articles/creating-a-personal-access-token-for-the-command-line/)).

For finer-grained control over the build we direct the reader to the more 
thorough CPP build instructions located 
[here](https://cmakepackagingproject.readthedocs.io/en/latest/end_user/quick_start.html)

and note that Integrals depends on several other projects:

- [utilities](https://github.com/NWChemEx-Project/Utilities)
- [SDE](https://github.com/NWChemEx-Project/SDE)
  - [bphash](https://github.com/bennybp/BPHash)
  - [pybind11](https://github.com/pybind/pybind11)
    - Requires a development version of Python
  - [cereal](https://github.com/USCiLab/cereal)
- [LibChemist](https://github.com/NWChemEx-Project/LibChemist)  
- [TAMM](https://github.com/NWChemEx-Project/TAMM)
  - [HPTT](https://github.com/ajaypanyala/hptt)
  - [MSGSL](https://github.com/Microsoft/GSL)
  - [LAPACKE](http://www.netlib.org/lapack/)
  - [Eigen3](https://github.com/eigenteam/eigen-git-mirror)
  - [GlobalArrays](https://github.com/GlobalArrays/ga)
  - MPI
- [LibInt](https://github.com/evaleev/libint)     
- (For testing only)[Catch2](https://github.com/catchorg/Catch2)

We recommend building Libint separately for the moment because a proper 
autotools build (the only autotools variant supported by CPP) builds the 
generator and many more integrals than are necessary.
