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

################
Four-Center ERIs
################

.. |eri4| replace:: :math:`\left(\mu\nu|\lambda\sigma\right)`

Electron repulsion integrals (ERIs) over a set of :math:`N_b` atomic orbitals
(AOs), :math:`\lbrace\phi_\mu\left(\vec{r}\right)\rbrace`, are the
four-index quantities:

.. math::
   :label: eri4

   \left(\mu\nu|\lambda\sigma\right) = \int\int
     \phi_\mu\left(\vec{r}_1\right)\phi_\nu\left(\vec{r}_1\right)
     \frac{1}{r_{12}}
     \phi_\lambda\left(\vec{r}_2\right)\phi_\sigma\left(\vec{r}_2\right)
     d\vec{r}_1 d\vec{r}_2.

There are formally :math:`O\left(N_b^4\right)` unique |eri4| and, even with
permutational and screening-based reductions, computing, storing, and
contracting them is the dominant cost of most Hartree-Fock and
post-Hartree-Fock methods. :doc:`density_fitting` discusses an approximation
that avoids forming |eri4| explicitly.
