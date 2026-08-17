.. Copyright 2025 NWChemEx-Project
..
.. Licensed under the Apache License, Version 2.0 (the "License");
.. you may not use this file except in compliance with the License.
.. You may obtain a copy of the License at
..
.. http://www.apache.org/licenses/LICENSE-2.0
..
.. Unless required by applicable law or agreed to in writing, software
.. distributed under the License is distributed on an "AS IS" BASIS,
.. WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
.. See the License for the specific language governing permissions and
.. limitations under the License.

#####################################
Electron-Nuclear Attraction Integrals
#####################################

Over a set of :math:`N_b` atomic orbitals (AOs),
:math:`\lbrace\phi_\mu\left(\vec{r}\right)\rbrace`, and a set of nuclei with
charges :math:`Z_A` located at :math:`\vec{R}_A`, the electron-nuclear
attraction integrals are:

.. math::

   V^{en}_{\mu\nu} = -\sum_A Z_A \int
     \phi_\mu\left(\vec{r}\right)
     \frac{1}{\left|\vec{r}-\vec{R}_A\right|}
     \phi_\nu\left(\vec{r}\right) d\vec{r}.
