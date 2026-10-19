<!--
  ~ Copyright 2026 NWChemEx-Project
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

ERI Error Explorer
==================

`eri_error_explorer.py` explores how the basis set parameters, the separation
between centers, and the screening tolerance affect the error in the
four-center electron repulsion integrals (ERIs). See the docstring at the top
of the script for details.

Setting up
----------

1. Build Integrals, with its Python bindings, as described in the
   [top-level README](../../README.md). Building with pip in a Python 3.12
   virtual environment is the simplest option.
2. In the same environment, install ChemCache (the script takes its STO-3G
   basis set from ChemCache) and NumPy:

   ```.sh
   pip install nwchemex-chemcache numpy
   ```

   ChemCache and Integrals must share one copy of SimDE and the libraries
   below it. That is the case if Integrals was built against the dependency
   wheels installed in this environment, which the top-level README's
   instructions do. If you built Integrals with `NWX_ECOSYSTEM_FROM_SOURCE=ON`,
   or in an environment without the wheels, build ChemCache from source the
   same way instead of installing its wheel.
3. If you built Integrals with CMake and did not install it into the
   environment, add the build directory to `PYTHONPATH` so Python can find the
   `integrals` module:

   ```.sh
   export PYTHONPATH=/path/to/Integrals/build:$PYTHONPATH
   ```

   This step is not needed after a pip install.

Running
-------

From the root of the Integrals repository:

```.sh
python examples/eri_error_explorer/eri_error_explorer.py
```

With the default settings this takes a few seconds. The script prints the
sampled shells, one line per shell quartet and tolerance, and a summary per
separation and tolerance.

Changing what is explored
-------------------------

The parameters are hard-coded near the top of the script. Edit them and run
the script again:

| Parameter               | Default            | Meaning                                                                           |
|-------------------------|--------------------|-----------------------------------------------------------------------------------|
| `SHELL_SAMPLES`         | `[(1, 0), (8, 2)]` | `(Z, i)` pairs: shell `i` of the STO-3G basis set for atomic number `Z` (H 1s, O 2p). |
| `R_MIN`, `R_MAX`        | `1.0`, `5.0`       | Range of center separations, in bohr.                                             |
| `N_SEPARATIONS`         | `3`                | Number of evenly spaced separations in that range.                                |
| `TOLERANCES`            | `[1e-6, 1e-8]`     | Screening tolerances to study.                                                    |
| `REFERENCE_TOLERANCE`   | `1e-25`            | Tolerance used to compute the reference ("exact") ERIs.                           |

Each sampled shell is placed on all four corners of a regular tetrahedron
whose edge length is the separation, so the number of shells, and therefore of
shell quartets, grows quickly with the number of samples.
