<!--
  ~ Copyright 2025 NWChemEx-Project
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

Building the Integrals Documentation
=====================================

This directory contains the source files for generating the Sphinx
documentation for `Integrals`. General instructions for building
documentation found throughout the NWChemEx project are available at:
https://github.com/NWChemEx/NWChemEx/blob/master/docs/README.md

Obtaining the Documentation's Dependencies
-------------------------------------------

The documentation's dependencies can be installed via Python's `pip`
command. Commands are assumed to be run from this directory (the same
directory as this README file).

~~~.sh
# These first two steps are strongly recommended, but not required
python3 -m venv venv
. venv/bin/activate
pip3 install -r requirements.txt
~~~

Building the Documentation
----------------------------

With the dependencies installed, build the documentation with:

~~~.sh
make html BUILDDIR=${BUILDDIR}
~~~

where `${BUILDDIR}` is the directory where you want the resulting HTML to be
placed (*e.g.* `build`).

Viewing the Documentation Locally
------------------------------------

After building, the main index of the resulting HTML will be located at
`${BUILDDIR}/html/index.html` and can be viewed by pointing your web browser
of choice at that file, either by opening it directly or by using
`file:///path/to/${BUILDDIR}/html/index.html` as the URL. Alternatively, for
a locally served copy that resembles how the docs are hosted online, run a
simple HTTP server from the build directory:

~~~.sh
python3 -m http.server --directory ${BUILDDIR}/html
~~~

and then navigate to `http://localhost:8000` in your browser.
