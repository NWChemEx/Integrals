.. Copyright 2026 NWChemEx-Project
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

##########################
Primitive Screening Models
##########################

.. |eri4| replace:: :math:`\left(\mu\nu|\lambda\sigma\right)`
.. |Kij| replace:: :math:`K_{ij}`
.. |Qij| replace:: :math:`Q_{ij}`

Most of the |eri4| defined in :doc:`eris` are negligible. Screening is the
practice of predicting, before an integral is computed, that its contribution
falls below a threshold and may be skipped. This page collects the screening
models used by the ``Integrals`` package, states each one precisely, and --
in :ref:`scr-rigor` -- records which of them are actually rigorous bounds.
That last question turns out to matter a great deal, because a model that is
adequate for *deciding what to skip* is not automatically adequate for
*quantifying what skipping cost*. :doc:`error_bounds` takes up the latter
problem.

**************************************
Notation and the Primitive Expansion
**************************************

Each contracted AO :math:`\phi_\mu` is a fixed linear combination of primitive
Gaussians. Writing the expansion out inside Eq. :eq:`eri4` turns a contracted
ERI into a sum over quartets of primitives:

.. math::
   :label: scr-prim-expansion

   \left(\mu\nu|\lambda\sigma\right) =
     \sum_{i\in\mu}\sum_{j\in\nu}\sum_{k\in\lambda}\sum_{l\in\sigma}
     c_i c_j c_k c_l \left[ij|kl\right],

where square brackets denote an integral over primitives and :math:`c_i` is
the contraction coefficient of primitive :math:`i`. Throughout, the
:math:`c_i` are the *renormalized* coefficients returned by the
``Normalize`` property type -- the same values ``PrimitiveContractor``
consumes through its ``Primitive Normalization`` submodule. These absorb both
the per-primitive normalization factor and the unit-norm scaling of the
contracted shell, so Eq. :eq:`scr-prim-expansion` holds as written, with no
further normalization factors.

Screening operates on Eq. :eq:`scr-prim-expansion` one quartet at a time, so
the quantities of interest are properties of *pairs*. For two primitives with
exponents :math:`\alpha_i`, :math:`\alpha_j` centered at :math:`\vec{A}`,
:math:`\vec{B}`, the Gaussian product theorem gives

.. math::
   :label: scr-gpt

   \gamma_{ij} = \alpha_i + \alpha_j, \qquad
   \rho_{ij} = \frac{\alpha_i\alpha_j}{\gamma_{ij}}, \qquad
   \vec{P}_{ij} = \frac{\alpha_i\vec{A} + \alpha_j\vec{B}}{\gamma_{ij}},

and the product of the two radial parts is a single Gaussian of exponent
:math:`\gamma_{ij}` centered at :math:`\vec{P}_{ij}`, scaled by the overlap
factor

.. math::
   :label: scr-kij

   K_{ij} = \exp\left(-\rho_{ij}\left|\vec{A}-\vec{B}\right|^2\right).

|Kij| decays as the square of the distance between the two centers, and it is
the origin of essentially all screening: a pair of primitives on distant
centers has an exponentially small product density, so any integral involving
that pair is exponentially small.

*******************************
The Primitive [ss|ss] Integral
*******************************

For four s-type primitives the integral in Eq. :eq:`scr-prim-expansion` has a
closed form. With :math:`\vec{P}` and :math:`\vec{Q}` the product centers of
the bra and ket pairs,

.. math::
   :label: scr-ssss

   \left[ij|kl\right] =
     \frac{2\pi^{5/2}}
          {\gamma_{ij}\gamma_{kl}\sqrt{\gamma_{ij}+\gamma_{kl}}}
     K_{ij}K_{kl}\,F_0\left(T\right), \qquad
   T = \frac{\gamma_{ij}\gamma_{kl}}{\gamma_{ij}+\gamma_{kl}}
       \left|\vec{P}-\vec{Q}\right|^2,

where :math:`F_0` is the zeroth-order Boys function :cite:`boys1950`,

.. math::
   :label: scr-boys

   F_0\left(T\right) = \int_0^1 \exp\left(-Tu^2\right)du
     = \frac{1}{2}\sqrt{\frac{\pi}{T}}\,\mathrm{erf}\left(\sqrt{T}\right).

Two properties of :math:`F_0` are used repeatedly below. It is bounded,
:math:`0 < F_0\left(T\right) \leq F_0\left(0\right) = 1`, and for large
:math:`T` it decays as
:math:`F_0\left(T\right) \rightarrow \frac{1}{2}\sqrt{\pi/T}`, which is the
classical :math:`1/\left|\vec{P}-\vec{Q}\right|` Coulomb tail. Equation
:eq:`scr-ssss` also makes plain that the s-type primitive ERI is *strictly
positive*, a fact :doc:`error_bounds` relies on.

**********************
Coarse Pair Estimates
**********************

The cheapest useful estimate discards everything in Eq. :eq:`scr-ssss` except
the exponential factors. ``BlackBoxPrimitiveEstimator`` computes

.. math::
   :label: scr-coarse

   \bar{K}_{ij} = \left|c_i\right|\left|c_j\right|K_{ij},

which is the quantity libint uses for coarse screening and which the
``Integrals`` package builds with the ``coarse_k_ij`` helper. It requires only
exponents, centers, and coefficients -- no integrals at all -- so it can be
tabulated once per basis-set pair and reused for every quartet.

******************
The Fine Metric
******************

The next refinement restores the exponent-dependent prefactor. libint's
``ShellPair`` data carries, and the ``fine_k_ij`` helper reproduces,

.. math::
   :label: scr-fine-pair

   Q_{ij} = \bar{K}_{ij}\,\frac{\sqrt{2}\,\pi^{5/4}}{\gamma_{ij}},

from which the quartet-level fine screening quantity is formed as
:math:`\left|Q_{ij}Q_{kl}\right|/\sqrt{\gamma_{ij}+\gamma_{kl}}`. Multiplying
Eq. :eq:`scr-fine-pair` out for the bra and ket pairs gives

.. math::
   :label: scr-fine-identity

   \frac{Q_{ij}Q_{kl}}{\sqrt{\gamma_{ij}+\gamma_{kl}}} =
     \left|c_ic_jc_kc_l\right|
     \frac{2\pi^{5/2}}
          {\gamma_{ij}\gamma_{kl}\sqrt{\gamma_{ij}+\gamma_{kl}}}
     K_{ij}K_{kl},

since :math:`\sqrt{2}\cdot\sqrt{2} = 2` and
:math:`\pi^{5/4}\cdot\pi^{5/4} = \pi^{5/2}`. Comparing with Eq.
:eq:`scr-ssss`, the fine metric is *exactly* the contracted contribution of an
s-type primitive quartet with :math:`F_0\left(T\right)` replaced by one. It is
the :math:`\left|\vec{P}-\vec{Q}\right|\rightarrow 0` limit of the true
integral.

This identity explains both the strength and the weakness of the fine metric.
Because :math:`F_0 \leq 1`, it is an upper bound on the s-type contribution,
and a cheap one. But it is a *distance-blind* upper bound: it discards
:math:`F_0\left(T\right)` entirely, and :math:`F_0` is small precisely for the
well-separated quartets that screening is most interested in. The fine metric
is therefore loosest exactly where it is used most.

********************
Cauchy-Schwarz
********************

A qualitatively different estimate follows from the observation that the
electron-repulsion operator defines a positive-definite bilinear form on
charge distributions. The Cauchy-Schwarz inequality then applies directly
:cite:`haser1989`:

.. math::
   :label: scr-cauchy-schwarz

   \left|\left[ij|kl\right]\right| \leq
     \sqrt{\left[ij|ij\right]}\sqrt{\left[kl|kl\right]}.

Unlike Eqs. :eq:`scr-coarse` and :eq:`scr-fine-pair`, this holds for *any*
angular momentum, Cartesian or solid-harmonic, with no restriction whatsoever
-- it is a property of the operator, not of the basis functions. The cost is
that it requires the hyper-diagonal integrals :math:`\left[ij|ij\right]`,
which must actually be computed, though only :math:`O(N^2)` of them.

``CauchySchwarzPrimitiveEstimator`` computes this family of quantities. Its
implementation is described in :ref:`scr-cs-implementation` below, since the
precise quantity it returns matters when the estimate is reused as an error
bound.

Cauchy-Schwarz shares the fine metric's distance blindness: the right-hand
side of Eq. :eq:`scr-cauchy-schwarz` does not depend on
:math:`\left|\vec{P}-\vec{Q}\right|` at all, so it cannot capture the
:math:`1/\left|\vec{P}-\vec{Q}\right|` decay of a well-separated quartet.
Recovering that decay while preserving rigor is the subject of the
distance-including estimates, which the package does not currently implement;
:doc:`error_bounds` surveys them and explains why they are the most valuable
direction for future work.

.. _scr-cs-implementation:

Implementation Notes
====================

``CauchySchwarzPrimitiveEstimator`` decontracts both basis sets, computes the
primitive quartet tensor for the hyper-diagonal, and reduces each
:math:`\left(ij|ij\right)` shell block to a single number with
``utils::rank2_shell_norm``. Two details are worth stating explicitly, because
both affect how the result may be used:

#. The reduction is a **Frobenius norm over the shell block**, not a maximum
   over its elements. Writing :math:`M_{ij}` for the matrix of
   :math:`\left(ij|ij\right)` values within the block, the module forms
   :math:`\left\|M_{ij}\right\|_F`. Since :math:`M_{ij}` is a Gram matrix, its
   Frobenius norm is at least its largest diagonal element, so the result is a
   valid bound on every element of the block, but a loose one -- loose by up to
   the square root of the block size.

#. The returned value is scaled by the contraction coefficients,
   :math:`\left|c_ic_j\right|\left\|M_{ij}\right\|_F`, and is **not**
   square-rooted. This is the natural form for the ``PairScreener`` use case,
   where a per-pair magnitude is compared against a tolerance directly. It is
   *not* the form that appears in Eq. :eq:`scr-cauchy-schwarz`: assembling a
   Cauchy-Schwarz product from it requires taking the square root of the norm
   and being careful not to apply the coefficient factors twice.
   :doc:`error_bounds` writes the bound in terms of the unscaled norm to keep
   this unambiguous.

There is no double counting of normalization between the two factors. The
decontracted basis carries unit contraction coefficients but is converted to
libint with normalization embedded, so :math:`M_{ij}` is built from normalized
primitives, while :math:`c_i` supplies the contraction coefficient exactly
once.

.. _scr-rigor:

**************************
Rigor of These Estimates
**************************

The three models above are not interchangeable, and the difference is not
merely one of tightness.

Equations :eq:`scr-coarse` and :eq:`scr-fine-pair` are derived from the
s-type closed form, Eq. :eq:`scr-ssss`. Neither retains any part of the
angular prefactor that a primitive with :math:`l > 0` contributes. For an
s-type quartet the fine metric is a rigorous upper bound, by the identity
:eq:`scr-fine-identity` together with :math:`F_0 \leq 1`. For :math:`l > 0`
no such guarantee exists, and the true integral may exceed the estimate. This
is not a defect of the implementation: libint's ``ScreeningMethod::Original``,
which ``PrimitiveContractor`` reproduces, is a screening heuristic chosen for
its cost, and it is documented as such. Used at a tight threshold it discards
only genuinely negligible terms. Used as the basis of an *error model*,
however, it produces estimates that the true error can and does exceed.

Equation :eq:`scr-cauchy-schwarz` has no such restriction. It is a rigorous
bound at every angular momentum. This is the reason :doc:`error_bounds` builds
on Cauchy-Schwarz rather than on the metric that the screening gates
themselves use.

******************************
The Gates as Implemented
******************************

``PrimitiveContractor`` evaluates Eq. :eq:`scr-prim-expansion` directly, and
applies four tests in sequence. Writing :math:`t` for the threshold, a quartet
is skipped if any of the following holds:

#. :math:`\bar{K}_{ij} < t`, discarding the bra pair outright;
#. :math:`\bar{K}_{kl} < t`, discarding the ket pair outright;
#. :math:`\bar{K}_{ij}\bar{K}_{kl} \leq t`;
#. :math:`\left|Q_{ij}Q_{kl}\right|/\sqrt{\gamma_{ij}+\gamma_{kl}} < t`.

The first three are the coarse tests of Eq. :eq:`scr-coarse` and the last is
the fine test of Eq. :eq:`scr-fine-pair`. The threshold is exposed as the
``Screening Threshold`` input and defaults to :math:`10^{-16}`; the sequence
matches libint's ``ScreeningMethod::Original``.

A separate module, ``ScreenPrimitivePairs``, applies a single tolerance test
to a matrix of pair estimates supplied by any ``PrimitivePairEstimator``
submodule, and returns the surviving pair indices. This is the composable
entry point: it is agnostic to which of the estimates above is used.

****************************
From Screening to Error
****************************

Screening and error estimation ask different questions of the same
quantities. Screening asks whether a contribution is small enough to skip, and
a good heuristic answer suffices -- if the estimate is occasionally too small,
the consequence is a slightly larger error than the threshold nominally
promises. Error estimation asks what the accumulated cost of every skipped
contribution actually was, and here a heuristic will not do: an estimate that
the true value can exceed provides no guarantee at all.

Nothing requires the two to use the same model. The screening gates may
remain the cheap, distance-blind heuristics documented above while the error
model is built on the rigorous inequality of Eq. :eq:`scr-cauchy-schwarz`.
:doc:`error_bounds` develops that error model.

**********
References
**********

.. bibliography::
   :style: unsrt
   :filter: docname in docnames
   :labelprefix: SCR
