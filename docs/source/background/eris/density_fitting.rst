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

###############
Density Fitting
###############

.. |eri4| replace:: :math:`\left(\mu\nu|\lambda\sigma\right)`
.. |eri3| replace:: :math:`\left(P|\kappa\lambda\right)`
.. |eri2| replace:: :math:`\left(P|Q\right)`
.. |Bcal| replace:: :math:`B_{P,\mu\nu}`

Density fitting (DF), also known as the resolution-of-the-identity (RI)
approximation, reduces the cost of computing the four-center ERI |eri4|
(see :doc:`eris` for its definition and notation) by approximately
factorizing it into a product of three-index quantities. This page derives
the DF approximation using the Coulomb metric, which is the variant
implemented in the ``Integrals`` package (see
:ref:`df-implementation-notes` below).

**************************
The Density-Fitting Ansatz
**************************

DF/RI introduces an auxiliary basis set,
:math:`\lbrace\chi_P\left(\vec{r}\right)\rbrace_{P=1}^{N_{aux}}`, and
approximates each AO charge distribution (product density)
:math:`\rho_{\mu\nu}\left(\vec{r}\right) \equiv
\phi_\mu\left(\vec{r}\right)\phi_\nu\left(\vec{r}\right)` as a linear
expansion in that auxiliary basis :cite:`whitten1973,dunlap1979`:

.. math::

   \rho_{\mu\nu}\left(\vec{r}\right) \approx
     \tilde{\rho}_{\mu\nu}\left(\vec{r}\right) =
     \sum_P^{N_{aux}} c^P_{\mu\nu}\chi_P\left(\vec{r}\right).

Typically :math:`N_{aux}` is chosen to be a modest multiple of :math:`N_b`
(rather than scaling as :math:`N_b^2`, the number of AO pairs), which is
what makes the approximation useful.

*************************************
Fitting Coefficients: Coulomb Metric
*************************************

The expansion coefficients :math:`c^P_{\mu\nu}` are determined
variationally, by minimizing the residual self-repulsion energy of the
error density :math:`\Delta\rho_{\mu\nu} =
\rho_{\mu\nu} - \tilde\rho_{\mu\nu}` under the Coulomb metric
:cite:`whitten1973,dunlap1979`:

.. math::

   E^{err}_{\mu\nu} = \int\int
     \Delta\rho_{\mu\nu}\left(\vec{r}_1\right)
     \frac{1}{r_{12}}
     \Delta\rho_{\mu\nu}\left(\vec{r}_2\right)
     d\vec{r}_1 d\vec{r}_2 \geq 0.

Defining the two-center Coulomb-metric matrix over the auxiliary basis,

.. math::

   \left(P|Q\right) = \int\int
     \chi_P\left(\vec{r}_1\right)\frac{1}{r_{12}}
     \chi_Q\left(\vec{r}_2\right) d\vec{r}_1 d\vec{r}_2,

and the three-center integral between the auxiliary basis and an AO pair,

.. math::

   \left(P|\mu\nu\right) = \int\int
     \chi_P\left(\vec{r}_1\right)\frac{1}{r_{12}}
     \rho_{\mu\nu}\left(\vec{r}_2\right) d\vec{r}_1 d\vec{r}_2,

setting :math:`\partial E^{err}_{\mu\nu} / \partial c^P_{\mu\nu} = 0` gives
the well-known result:

.. math::
   :label: fit-coeffs

   c^P_{\mu\nu} = \sum_Q^{N_{aux}} \left(P|Q\right)^{-1}\left(Q|\mu\nu\right).

Substituting the fitted densities back into the definition of the
four-center ERI, Eq. :eq:`eri4`, yields the DF approximation
:cite:`vahtras1993,feyereisen1993`:

.. math::
   :label: df-eri4

   \left(\mu\nu|\lambda\sigma\right) \approx
     \sum_P^{N_{aux}}\sum_Q^{N_{aux}}
     \left(\mu\nu|P\right)\left(P|Q\right)^{-1}\left(Q|\lambda\sigma\right).

Because the metric matrix :math:`\left(P|Q\right)` is symmetric
positive-definite, Eq. :eq:`df-eri4` can always be written as a contraction
of two three-index tensors,

.. math::
   :label: df-three-index

   \left(\mu\nu|\lambda\sigma\right) \approx
     \sum_P^{N_{aux}} B_{P,\mu\nu} B_{P,\lambda\sigma},

for *any* factorization :math:`\left(P|Q\right)^{-1} = \sum_R M_{PR}M_{QR}`
of the inverse metric, with :math:`B_{P,\mu\nu} = \sum_Q M_{PQ}\left(Q|\mu\nu\right)`.
The two common choices for :math:`M` are the symmetric inverse square root,
:math:`M = \left(P|Q\right)^{-1/2}`, and a Cholesky factor of the *inverse*
metric; the latter is used in this package and is discussed next.

.. _df-implementation-notes:

*******************************************
Implementation: Cholesky-Factored Metric
*******************************************

Rather than diagonalizing :math:`\left(P|Q\right)` to form the symmetric
inverse square root, the ``Integrals`` package Cholesky-factorizes the
Coulomb metric,

.. math::

   \left(P|Q\right) = \sum_R L_{PR}L_{QR} = \left(LL^T\right)_{PQ},

with :math:`L` lower triangular, and forms :math:`B_{P,\mu\nu}` using
:math:`M = L^{-1}`:

.. math::
   :label: cholesky-fit

   B_{P,\mu\nu} = \sum_Q \left(L^{-1}\right)_{PQ}\left(Q|\mu\nu\right).

This is equivalent to Eq. :eq:`df-three-index` since
:math:`\sum_R\left(L^{-1}\right)_{PR}\left(L^{-1}\right)_{QR} =
\left[\left(LL^T\right)^{-1}\right]_{PQ} = \left(P|Q\right)^{-1}`, but avoids
an explicit eigendecomposition of the metric; the Cholesky-decomposition
route to the fitted three-index tensor is discussed and benchmarked against
the symmetric-inverse-square-root (RI) route by
:cite:t:`weigend2009`. Concretely:

- The ``CoulombMetric`` module computes the two-center integral |eri2| (via
  its "Two-center ERI" submodule, which defaults to the raw ``ERI2``
  integral driver), Cholesky-decomposes it, and returns :math:`M = L^{-1}`.
- The ``DFIntegral`` module computes the raw three-center integral |eri3|
  (via its "Three-center ERI" submodule, defaulting to the raw ``ERI3``
  driver) and its "Coulomb Metric" submodule (defaulting to
  ``CoulombMetric``), then contracts them to form |Bcal| following
  Eq. :eq:`cholesky-fit`.

**********************
Density-Fitted J and K
**********************

In self-consistent-field methods the Coulomb (J) and exchange (K) matrices
are built by contracting |eri4| with an AO density matrix :math:`D_{\kappa\lambda}`:

.. math::

   J_{\mu\nu} = \sum_{\kappa\lambda}\left(\mu\nu|\kappa\lambda\right)D_{\kappa\lambda},
   \qquad
   K_{\mu\nu} = \sum_{\kappa\lambda}\left(\mu\lambda|\kappa\nu\right)D_{\kappa\lambda}.

Substituting the DF approximation, Eq. :eq:`df-three-index`, for |eri4|
avoids ever forming the four-index quantities explicitly.

**Density-Fitted J** (``JDensityFitted``). Contract the density into the
fitted tensor to form an intermediate over only the auxiliary index,

.. math::

   j_P = \sum_{\kappa\lambda} D_{\kappa\lambda}B_{P,\kappa\lambda},

then contract back onto the AO pair:

.. math::

   J_{\mu\nu} = \sum_P^{N_{aux}} j_P B_{P,\mu\nu}.

This costs :math:`O\left(N_{aux}N_b^2\right)` to form :math:`j_P` and
:math:`O\left(N_{aux}N_b^2\right)` to form :math:`J_{\mu\nu}`, i.e. no step
scales worse than cubically in the system size.

**Density-Fitted K** (``KDensityFitted``). Exchange cannot be reduced to a
single index-:math:`P` intermediate because the density matrix and the two
fitted tensors do not share a common pair of AO indices. Instead, an
intermediate retaining the auxiliary index and one AO index is formed,

.. math::

   K^{aux}_{P,\mu,\lambda} = \sum_\kappa D_{\kappa\lambda} B_{P,\mu\kappa},

and then contracted with the second fitted tensor:

.. math::

   K_{\mu\nu} = \sum_P^{N_{aux}}\sum_\lambda K^{aux}_{P,\mu,\lambda}B_{P,\nu\lambda}.

Because the intermediate :math:`K^{aux}_{P,\mu,\lambda}` retains three
indices (rather than the single auxiliary index :math:`j_P` used for J),
forming density-fitted K is more expensive than density-fitted J, scaling
as :math:`O\left(N_{aux}N_b^3\right)`.

**********
References
**********

.. bibliography::
   :style: unsrt
   :filter: docname in docnames
   :labelprefix: DF
