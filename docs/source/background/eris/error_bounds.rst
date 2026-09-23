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

#############################
Bounding the Screening Error
#############################

.. |eri4| replace:: :math:`\left(\mu\nu|\lambda\sigma\right)`
.. |delta| replace:: :math:`\delta_{\mu\nu\lambda\sigma}`
.. |Bq| replace:: :math:`B_q`

:doc:`screening` describes how contributions to |eri4| are discarded. This
page asks what the discarding cost. The answer developed here is a rigorous
two-sided interval,
:math:`\left[\delta_{\min},\delta_{\max}\right]`, that is guaranteed to
contain the error in every screened integral.

The construction has two ingredients. A magnitude bound, taken from the
Cauchy-Schwarz inequality, brackets each discarded contribution symmetrically
about zero. A *sign oracle* then collapses that symmetric interval to a
half-interval for those contributions whose sign can be established a priori,
tightening the result and -- more usefully -- making it asymmetric. The second
ingredient is optional: the bound is valid with no sign information at all,
and improves monotonically as more signs are determined.

*****************************
The Exact Screening Error
*****************************

Let :math:`t` be the screening threshold and :math:`S\left(t\right)` the set
of primitive quartets that the gates of :doc:`screening` discard. Write a
quartet as :math:`q = \left(ijkl\right)` and its exact contribution to Eq.
:eq:`scr-prim-expansion` as

.. math::
   :label: err-tq

   t_q = c_ic_jc_kc_l\left[ij|kl\right].

The screened integral is the sum of Eq. :eq:`err-tq` over the quartets that
survive. Subtracting the exact result, which sums over all quartets, leaves

.. math::
   :label: err-delta

   \delta_{\mu\nu\lambda\sigma} \equiv
     \left(\mu\nu|\lambda\sigma\right)^{t} -
     \left(\mu\nu|\lambda\sigma\right)^{\rm exact} =
     -\sum_{q\in S\left(t\right)} t_q.

The sign convention in Eq. :eq:`err-delta` is *computed minus exact*, and it
is used without exception below. A positive |delta| means the screened value
is too large. The leading minus sign is easy to lose: discarding a positive
contribution makes the computed integral too *small*.

*****************************************
A Bracket from Per-Quartet Intervals
*****************************************

Suppose that for each discarded quartet we can produce an interval
:math:`t_q \in \left[\ell_q, u_q\right]`. Because Eq. :eq:`err-delta` is a
plain sum, the intervals add, and negating flips the endpoints:

.. math::
   :label: err-generic-bracket

   \delta_{\mu\nu\lambda\sigma} \in
     \left[-\sum_{q\in S\left(t\right)}u_q,\;
           -\sum_{q\in S\left(t\right)}\ell_q\right].

This is the entire structure of the model. Everything that follows is a
matter of supplying :math:`\ell_q` and :math:`u_q`. Nothing in Eq.
:eq:`err-generic-bracket` requires the intervals to be tight, symmetric, or
sign-resolved; it requires only that each one genuinely contain its
:math:`t_q`.

**********************
The Magnitude Bound
**********************

The default interval comes from Cauchy-Schwarz, Eq. :eq:`scr-cauchy-schwarz`.
Applying it to Eq. :eq:`err-tq`,

.. math::
   :label: err-bq

   \left|t_q\right| \leq B_q \equiv
     \left|c_ic_jc_kc_l\right|
     \sqrt{\left[ij|ij\right]}\sqrt{\left[kl|kl\right]},

so that :math:`t_q \in \left[-B_q, +B_q\right]` unconditionally.

The choice of Cauchy-Schwarz here is deliberate and is the single most
important design decision on this page. As :ref:`scr-rigor` records, the
coarse and fine metrics that drive the screening gates are derived from the
s-type closed form and retain no angular prefactor; for :math:`l > 0` the true
integral may exceed them. An error model built on the fine metric therefore
produces intervals that the true error can escape, which defeats the purpose.
Equation :eq:`scr-cauchy-schwarz` is a rigorous bound at every angular
momentum, so Eq. :eq:`err-bq` cannot be violated.

In practice :math:`\sqrt{\left[ij|ij\right]}` is evaluated per shell block
rather than per AO. Writing :math:`M_{ij}` for the block of
:math:`\left(ij|ij\right)` values, any quantity at least as large as the
largest diagonal element of :math:`M_{ij}` may be substituted; the Frobenius
norm :math:`\left\|M_{ij}\right\|_F` computed by
``CauchySchwarzPrimitiveEstimator`` qualifies, giving the usable form

.. math::
   :label: err-bq-blocked

   B_q = \left|c_ic_jc_kc_l\right|
     \sqrt{\left\|M_{ij}\right\|_F}\sqrt{\left\|M_{kl}\right\|_F}.

Note the square roots, and note that the coefficient factors appear once.
:ref:`scr-cs-implementation` describes what that module returns in its raw
form; the two differ, and the difference matters.

****************************************
Tightening with Sign Information
****************************************

The sign of a discarded contribution factorizes,

.. math::
   :label: err-sign-factor

   \mathrm{sign}\left(t_q\right) =
     \mathrm{sign}\left(c_ic_jc_kc_l\right)\cdot
     \mathrm{sign}\left(\left[ij|kl\right]\right),

and the first factor is free. ``PrimitiveNormalization`` emits one
renormalized coefficient per primitive, shared by all AO components of that
primitive, so the coefficient product is available wherever the contraction
itself is. It is also genuinely informative: standard contracted basis sets
carry coefficients of both signs, so this factor is not a constant.

The second factor requires an argument. If the primitive charge distributions
:math:`\rho_{ij} = g_ig_j` and :math:`\rho_{kl} = g_kg_l` are each of one sign
throughout space, then :math:`\left[ij|kl\right]` is the Coulomb interaction
of two same-signed distributions and is therefore positive. This gives a
*sufficient* -- not necessary -- criterion, so the oracle must be permitted to
answer "unknown."

Define :math:`\sigma\left(q\right) \in \left\{+1,-1,0\right\}`, with
:math:`0` denoting an indeterminate sign:

Tier 0: all-s quartets
======================

If all four primitives have :math:`l = 0`, both distributions are strictly
positive -- Eq. :eq:`scr-ssss` shows this directly -- and

.. math::
   :label: err-tier0

   \sigma\left(q\right) = \mathrm{sign}\left(c_ic_jc_kc_l\right).

This tier is available in any basis and requires no information beyond the
angular momenta and the coefficients.

Tier 1: one-signed Cartesian products
=====================================

For Cartesian primitives, :math:`\rho_{ij}` is a product of three independent
factors, one per axis. Along :math:`x` the factor is
:math:`\left(x-A_x\right)^{a_x}\left(x-B_x\right)^{b_x}`, which changes sign
at :math:`A_x` when :math:`a_x` is odd and at :math:`B_x` when :math:`b_x` is
odd. It is therefore one-signed if and only if

* :math:`a_x` and :math:`b_x` are both even, **or**
* :math:`A_x = B_x` and :math:`a_x + b_x` is even,

and likewise for :math:`y` and :math:`z`. When this holds for both the bra and
ket pairs, Eq. :eq:`err-tier0` again applies.

The second clause carries most of the coverage in practice. A same-center
:math:`p_x \cdot p_x` pair gives :math:`\left(x-A_x\right)^2`, which is
one-signed, and the pairs that matter for screening are heavily concentrated
on or near a single center. The first clause alone would admit only functions
with all-even Cartesian exponents.

This tier is unavailable for ``pure`` shells, which is a per-shell property in
this package. Real solid harmonics with :math:`l > 0` have angular nodes, so
their products are not one-signed and the criterion cannot fire.

.. _err-cartesian-caveat:

Cartesian does not guarantee a sign
===================================

It is worth being explicit that working in a Cartesian basis does *not* make
every sign determinable. Consider :math:`\left(p_xp_x|ss\right)` with the two
:math:`p` functions on distinct centers. The bra factor is
:math:`\left(x-A_x\right)\left(x-B_x\right)`, which is negative between the
centers and positive outside, so Tier 1 declines it and
:math:`\sigma\left(q\right) = 0`.

Because one-signedness is sufficient rather than necessary, such a quartet may
still have a definite sign; the oracle simply cannot certify it cheaply. This
is a limit on coverage, not a failure of the bound. What a Cartesian basis
buys is far *more* coverage than a solid-harmonic one, where only
:math:`l = 0` qualifies at all.

**********************
Assembling the Bracket
**********************

Partition :math:`S\left(t\right)` by the oracle into :math:`S_+`,
:math:`S_-`, and :math:`S_0`. A determined sign collapses
:math:`\left[-B_q, +B_q\right]` to a half-interval anchored at zero:
:math:`t_q \in \left[0, B_q\right]` on :math:`S_+` and
:math:`t_q \in \left[-B_q, 0\right]` on :math:`S_-`. Substituting into Eq.
:eq:`err-generic-bracket`,

.. math::
   :label: err-bracket

   \delta_{\max} = \sum_{q\in S_-}B_q + \sum_{q\in S_0}B_q,
   \qquad
   \delta_{\min} = -\sum_{q\in S_+}B_q - \sum_{q\in S_0}B_q.

Four properties follow immediately and are worth stating together.

**It cannot be violated.** The magnitude bound is rigorous at every angular
momentum and the oracle declines to guess. Containment is guaranteed by
construction, not by empirical calibration.

**It is valid with no signs at all.** If :math:`S_0 = S\left(t\right)`, Eq.
:eq:`err-bracket` reduces to the symmetric envelope
:math:`\pm\sum_q B_q`. The sign oracle can therefore be introduced
incrementally, one tier at a time, with no risk of regression: adding sign
coverage can only shrink the interval.

**It can become one-sided.** If every discarded quartet has a determined
positive sign, then :math:`\delta_{\max} = 0` and the bracket is
:math:`\left[-\sum_q B_q, 0\right]` -- the screened integral is *guaranteed*
to be an underestimate, with a rigorous bound on by how much. A symmetric
envelope can never produce a statement of this kind at any threshold.

**It composes with interval arithmetic.** Equation :eq:`err-bracket` is an
asymmetric interval, so its midpoint is a natural best estimate and its
half-width a natural uncertainty, in the form the package's uncertainty
machinery already consumes.

*****************************
How Much the Sign Buys
*****************************

It is easy to overstate the value of the sign oracle, so it is worth
quantifying. The width of the bracket is

.. math::
   :label: err-width

   W = \delta_{\max} - \delta_{\min}
     = \sum_{q\in S_0}2B_q + \sum_{q\in S_+\cup S_-}B_q
     = \left(1+f_0\right)\sum_{q\in S\left(t\right)}B_q,

where

.. math::
   :label: err-f0

   f_0 = \frac{\sum_{q\in S_0}B_q}{\sum_{q\in S\left(t\right)}B_q}

is the :math:`B`-weighted fraction of indeterminate quartets. Relative to the
symmetric envelope, which is the :math:`f_0 = 1` case, sign determination
narrows the bracket by a factor of :math:`2/\left(1+f_0\right)`. **The best
possible improvement in width is therefore a factor of two**, attained only
when every sign is known.

A factor of two is real but bounded, and it is not the main reason to do this.
The qualitative payoff is asymmetry. A symmetric interval always contains
zero, so it can only ever support the statement "the error is no larger than
this." An asymmetric one may exclude zero, supporting the strictly stronger
statement "the computed value is too small, by no more than this." Direction
is information that no amount of tightening a symmetric envelope will produce.

Two further consequences of Eq. :eq:`err-width` are worth noting. Tightening
:math:`B_q` itself -- for instance via the distance-including bounds discussed
below -- is unbounded in its potential benefit, and so is ultimately the more
valuable direction. And because :math:`f_0` rises with angular momentum, the
sign oracle contributes most in exactly the s- and p-dominated regimes where
the integrals are cheapest to begin with.

***********************************************
Cartesian Screening and the Pure-Basis Bracket
***********************************************

Tier 1 requires Cartesian primitives, but results are usually wanted in a
solid-harmonic basis. These are reconciled without needing the sign of the
transformed integral at all, because the transformation is linear and
*intervals propagate through linear maps exactly*.

Let :math:`w_c` denote the product of the four transformation coefficients
relating a solid-harmonic quartet to the Cartesian quartet :math:`c`, and let
:math:`\left[\ell_c, u_c\right]` be the Cartesian bracket. Then

.. math::
   :label: err-transform

   \ell_{\rm pure} = \sum_c \begin{cases}
       w_c\ell_c & w_c > 0\\ w_cu_c & w_c < 0
     \end{cases},
   \qquad
   u_{\rm pure} = \sum_c \begin{cases}
       w_cu_c & w_c > 0\\ w_c\ell_c & w_c < 0
     \end{cases}.

Each endpoint is obtained by pairing every coefficient with whichever end of
its interval extremizes the sum. The result is rigorous, and one-sidedness
survives wherever cancellation among the :math:`w_c` does not destroy it.

The widening this introduces is worth naming. Equation :eq:`err-transform`
treats the Cartesian contributions as independent when they are in fact
correlated, so the resulting interval is wider than the true attainable range.
The effect is mild for :math:`d` shells, where only a few Cartesian components
contribute to each solid-harmonic function, and grows with :math:`l`.

The transformation coefficients themselves need not be derived: libint exposes
them in sparse row-compressed form, together with a closed-form generator. The
only new machinery required is an application of those coefficients that
carries a pair of endpoints rather than a single value, since the routines
libint provides transform values.

***************
Diagnostics
***************

Three cheap quantities characterize how well the model is doing on a given
system, and are more informative than the interval width alone:

* :math:`f_0` from Eq. :eq:`err-f0`, the :math:`B`-weighted fraction of
  indeterminate quartets. This measures sign *coverage* and, by Eq.
  :eq:`err-width`, fixes the width improvement exactly.
* The asymmetry ratio
  :math:`\left|\delta_{\max}+\delta_{\min}\right|/
  \left(\delta_{\max}-\delta_{\min}\right)`, which is zero for a symmetric
  envelope and one for a fully one-sided bracket.
* Whether the bracket excludes zero, per element. This is the binary form of
  the previous item and the one that matters for reporting: it is the
  condition under which a *direction* can be asserted.

.. _err-tightness:

*****************************
Practical Tightness Caveat
*****************************

The dominant source of looseness in Eq. :eq:`err-bracket` is not the sign
oracle but the magnitude bound. Two effects compound.

First, the shell-block Frobenius norm of :ref:`scr-cs-implementation`
overestimates the largest element of the block by up to the square root of the
block size. Replacing it with a per-element maximum over the block would
recover this without sacrificing rigor.

Second, and more importantly, Cauchy-Schwarz is distance-blind. The right-hand
side of Eq. :eq:`err-bq` does not depend on
:math:`\left|\vec{P}-\vec{Q}\right|`, whereas the true integral decays as
:math:`1/\left|\vec{P}-\vec{Q}\right|` by Eq. :eq:`scr-ssss`. Since screening
discards well-separated quartets preferentially, the bound is loosest over
precisely the set it is summed across.

***********************
Refinements Deferred
***********************

Two directions were considered and are recorded here so the choice is
traceable.

**Distance-including bounds.** The :math:`1/\left|\vec{P}-\vec{Q}\right|`
decay can be restored while preserving rigor, which addresses the second
effect in :ref:`err-tightness` directly and, per Eq. :eq:`err-width`, offers
unbounded improvement rather than the sign oracle's factor of two. Rigorous
and non-rigorous variants are both well documented in the literature
:cite:`gill1994,maurer2012,thompson2017`.

**A Boys-damped signed estimate.** Restoring :math:`F_0\left(T\right)` to the
fine metric and keeping the coefficient signs yields, by Eq.
:eq:`scr-fine-identity`, an estimate that is *exact* for s-type quartets and
cheaper than Cauchy-Schwarz, requiring no hyper-diagonal integrals. It was
rejected as the starting point for two reasons: it is a prediction rather than
a bound, so it supports no guarantee; and it inherits the fine metric's
failure to bound :math:`l > 0` contributions, which is the very behavior that
motivated moving off the fine metric. It remains attractive as a *bias
estimate* to accompany the bracket, and could be added without disturbing the
bound.

*******************************************
Relation to the Existing Modules
*******************************************

``PrimitiveErrorModel`` implements the unsigned predecessors of this model. It
walks the same decontracted quartet loop as ``PrimitiveContractor`` and
applies the same gates, accumulating one of three per-quartet estimates for
each skipped quartet. In the language of this page, each is a choice of
:math:`B_q` with :math:`S_0 = S\left(t\right)` -- that is, a symmetric
envelope with no sign information:

* ``"Tolerance"`` takes :math:`B_q = t`, which is not a bound at all but a
  nominal per-quartet charge;
* ``"Coarse"`` takes :math:`B_q = \bar{K}_{ij}\bar{K}_{kl}`;
* ``"Fine"`` takes
  :math:`B_q = \left|Q_{ij}Q_{kl}\right|/\sqrt{\gamma_{ij}+\gamma_{kl}}`,
  subject to the :math:`l > 0` caveat of :ref:`scr-rigor`.

Validation
==========

``AnalyticError`` computes the true signed error directly, as the difference
between a screened result and a tight-threshold benchmark. It therefore
supplies exactly the quantity Eq. :eq:`err-delta` defines, element by element,
and the model can be checked against it without any statistical machinery.

The acceptance criterion is containment, and it is absolute: every element of
the true error must lie within :math:`\left[\delta_{\min},\delta_{\max}\right]`,
with **zero** violations across the test set. Only once containment holds is
tightness meaningful, at which point the diagnostics above quantify it. A
model that is tight but occasionally violated is strictly less useful than one
that is loose but never violated, since only the latter supports a guarantee.
The ``primitive_error_models`` example provides a driver for these
comparisons.

***************
Limitations
***************

Sign coverage degrades as angular momentum rises. Tier 0 requires an all-s
quartet; Tier 1 requires one-signed Cartesian products and is unavailable for
``pure`` shells. For basis sets dominated by high-:math:`l` solid-harmonic
functions, :math:`f_0 \rightarrow 1` and Eq. :eq:`err-bracket` degenerates to
the symmetric envelope. The model does not fail in that limit -- it simply
stops adding value over what a magnitude bound alone provides.

The bound also says nothing about how the error propagates into derived
quantities. Screening errors in |eri4| enter the Fock matrix and the total
energy through contractions with the density, and bounding those requires
tracking correlations between elements that Eq. :eq:`err-bracket` treats
independently. That analysis is outside the scope of this page.

**********
References
**********

.. bibliography::
   :style: unsrt
   :filter: docname in docnames
   :labelprefix: ERR
