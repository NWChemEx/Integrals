/*
 * Copyright 2026 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/** @file primitive_quartet_screener.hpp
 *  @brief The primitive-quartet screening decisions libint makes for ERI4s.
 */

#pragma once
#include "primitive_pair_estimators.hpp"
#include <cmath>
#include <cstddef>
#include <simde/simde.hpp>
#include <vector>

namespace integrals::libint::detail_ {

/** @brief Why a primitive quartet was screened out.
 *
 *  The gates are checked in the order listed. A quartet that fails several of
 *  them is labeled by the first one it fails.
 */
enum class ScreeningReason : int {
    /// The quartet passed every gate and is part of the integral.
    kept = 0,
    /// The coarse estimate of the bra pair is below the threshold.
    bra_pair = 1,
    /// The coarse estimate of the ket pair is below the threshold.
    ket_pair = 2,
    /// The product of the coarse bra and ket estimates is not above the
    /// threshold.
    coarse = 3,
    /// The fine estimate of the quartet is below the threshold.
    fine = 4
};

/** @brief Decides which primitive quartets libint screens out of an ERI4.
 *
 *  This class reproduces libint's `ScreeningMethod::Original` for four-center
 *  ERIs. It precomputes the coarse estimates (`coarse_k_ij`), the fine
 *  estimates (`fine_k_ij`), and the exponent sums (`gamma_ij`) for the bra
 *  and the ket. Then it applies the same gates as libint, one primitive
 *  quartet at a time.
 *
 *  Primitive indices are primitive-shell indices, i.e., they run over the
 *  primitives of each basis set (see `build_prim_ao_to_prim_shell_map`).
 */
class PrimitiveQuartetScreener {
public:
    /// Type of the basis sets
    using basis_type = simde::type::ao_basis_set;

    /// Type of a primitive index
    using size_type = std::size_t;

    /** @brief Precomputes the pair quantities for the bra and the ket.
     *
     *  @param[in] bs0 The basis set for the first bra index.
     *  @param[in] bs1 The basis set for the second bra index.
     *  @param[in] bs2 The basis set for the first ket index.
     *  @param[in] bs3 The basis set for the second ket index.
     *  @param[in] thresh The screening threshold.
     *
     *  @throw std::bad_alloc if there is a problem allocating the matrices.
     *                        Strong throw guarantee.
     */
    PrimitiveQuartetScreener(const basis_type& bs0, const basis_type& bs1,
                             const basis_type& bs2, const basis_type& bs3,
                             double thresh) :
      m_thresh_(thresh),
      m_K_bra_(coarse_k_ij(bs0, bs1)),
      m_K_ket_(coarse_k_ij(bs2, bs3)),
      m_Q_bra_(fine_k_ij(bs0, bs1)),
      m_Q_ket_(fine_k_ij(bs2, bs3)),
      m_gamma_bra_(gamma_ij(bs0, bs1)),
      m_gamma_ket_(gamma_ij(bs2, bs3)) {}

    /** @brief Is the bra pair (pi, pj) screened out on its own?
     *
     *  If it is, every quartet with this bra pair is screened out, which lets
     *  callers skip the ket loop.
     *
     *  @param[in] pi The primitive index for the first bra basis set.
     *  @param[in] pj The primitive index for the second bra basis set.
     *
     *  @return True if the bra pair is screened out and false otherwise.
     *
     *  @throw None No throw guarantee.
     */
    bool bra_pair_screened(size_type pi, size_type pj) const noexcept {
        return m_K_bra_[pi][pj] < m_thresh_;
    }

    /** @brief Determines whether the quartet (pi pj|pk pl) is screened out.
     *
     *  @param[in] pi The primitive index for the first bra basis set.
     *  @param[in] pj The primitive index for the second bra basis set.
     *  @param[in] pk The primitive index for the first ket basis set.
     *  @param[in] pl The primitive index for the second ket basis set.
     *
     *  @return The first gate the quartet fails, or ScreeningReason::kept if
     *          it passes them all.
     *
     *  @throw None No throw guarantee.
     */
    ScreeningReason screen(size_type pi, size_type pj, size_type pk,
                           size_type pl) const noexcept {
        if(bra_pair_screened(pi, pj)) return ScreeningReason::bra_pair;
        if(m_K_ket_[pk][pl] < m_thresh_) return ScreeningReason::ket_pair;
        if(coarse_estimate(pi, pj, pk, pl) <= m_thresh_)
            return ScreeningReason::coarse;
        if(fine_estimate(pi, pj, pk, pl) < m_thresh_)
            return ScreeningReason::fine;
        return ScreeningReason::kept;
    }

    /** @brief The coarse estimate @f$K_{ij}K_{kl}@f$ of a quartet.
     *
     *  @param[in] pi The primitive index for the first bra basis set.
     *  @param[in] pj The primitive index for the second bra basis set.
     *  @param[in] pk The primitive index for the first ket basis set.
     *  @param[in] pl The primitive index for the second ket basis set.
     *
     *  @return The coarse estimate of the quartet.
     *
     *  @throw None No throw guarantee.
     */
    double coarse_estimate(size_type pi, size_type pj, size_type pk,
                           size_type pl) const noexcept {
        return m_K_bra_[pi][pj] * m_K_ket_[pk][pl];
    }

    /** @brief The fine estimate of a quartet.
     *
     *  The fine estimate is
     *  @f$|Q_{ij}Q_{kl}| / \sqrt{\gamma_{ij} + \gamma_{kl}}@f$.
     *
     *  @param[in] pi The primitive index for the first bra basis set.
     *  @param[in] pj The primitive index for the second bra basis set.
     *  @param[in] pk The primitive index for the first ket basis set.
     *  @param[in] pl The primitive index for the second ket basis set.
     *
     *  @return The fine estimate of the quartet.
     *
     *  @throw None No throw guarantee.
     */
    double fine_estimate(size_type pi, size_type pj, size_type pk,
                         size_type pl) const noexcept {
        const auto Q_ijkl     = m_Q_bra_[pi][pj] * m_Q_ket_[pk][pl];
        const auto gamma_ijkl = m_gamma_bra_[pi][pj] + m_gamma_ket_[pk][pl];
        return std::abs(Q_ijkl / std::sqrt(gamma_ijkl));
    }

private:
    /// Type of the precomputed pair matrices
    using matrix_type = std::vector<std::vector<double>>;

    /// The screening threshold
    double m_thresh_;

    /// Coarse estimates for the bra and ket pairs
    matrix_type m_K_bra_;
    matrix_type m_K_ket_;

    /// Fine estimates for the bra and ket pairs
    matrix_type m_Q_bra_;
    matrix_type m_Q_ket_;

    /// Exponent sums for the bra and ket pairs
    matrix_type m_gamma_bra_;
    matrix_type m_gamma_ket_;
};

} // namespace integrals::libint::detail_
