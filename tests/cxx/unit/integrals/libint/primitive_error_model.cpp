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

#include <integrals/libint/detail_/primitive_pair_estimators.hpp>
#undef DEPRECATED
#include "../testing/testing.hpp"
#include <array>
#include <cmath>
#include <integrals/integrals.hpp>
#include <integrals/utils/primitive_index_helpers.hpp>

using eri4_pt = simde::ERI4;
using pt      = integrals::property_types::Uncertainty<eri4_pt>;

using namespace integrals::testing;

namespace {

using integrals::libint::detail_::coarse_k_ij;
using integrals::libint::detail_::fine_k_ij;
using integrals::libint::detail_::gamma_ij;

bool quartet_skipped(double K_ij, double K_kl, double Q_ij, double Q_kl,
                     double gamma_ij, double gamma_kl, double thresh) {
    if(K_ij < thresh) return true;
    if(K_kl < thresh) return true;
    if(K_ij * K_kl <= thresh) return true;
    const double g    = gamma_ij + gamma_kl;
    const double pfac = std::abs(Q_ij * Q_kl / std::sqrt(g));
    return pfac < thresh;
}

struct SkipTotals {
    std::size_t n_skip;
    double sum_coarse;
    double sum_fine;
};

SkipTotals reference_skip_totals(const simde::type::ao_basis_set& bs0,
                                 const simde::type::ao_basis_set& bs1,
                                 const simde::type::ao_basis_set& bs2,
                                 const simde::type::ao_basis_set& bs3,
                                 double thresh) {
    const auto K_bra     = coarse_k_ij(bs0, bs1);
    const auto K_ket     = coarse_k_ij(bs2, bs3);
    const auto gamma_bra = gamma_ij(bs0, bs1);
    const auto gamma_ket = gamma_ij(bs2, bs3);
    const auto Q_bra     = fine_k_ij(bs0, bs1);
    const auto Q_ket     = fine_k_ij(bs2, bs3);

    auto map0  = integrals::utils::build_prim_ao_to_cgto_map(bs0);
    auto map1  = integrals::utils::build_prim_ao_to_cgto_map(bs1);
    auto map2  = integrals::utils::build_prim_ao_to_cgto_map(bs2);
    auto map3  = integrals::utils::build_prim_ao_to_cgto_map(bs3);
    auto pmap0 = integrals::utils::build_prim_ao_to_prim_shell_map(bs0);
    auto pmap1 = integrals::utils::build_prim_ao_to_prim_shell_map(bs1);
    auto pmap2 = integrals::utils::build_prim_ao_to_prim_shell_map(bs2);
    auto pmap3 = integrals::utils::build_prim_ao_to_prim_shell_map(bs3);

    std::array<std::size_t, 4> nprims{map0.size(), map1.size(), map2.size(),
                                      map3.size()};
    SkipTotals out{0, 0.0, 0.0};

    for(std::size_t i = 0; i < nprims[0]; ++i) {
        const auto pi = pmap0[i];
        for(std::size_t j = 0; j < nprims[1]; ++j) {
            const auto pj     = pmap1[j];
            const double K_ij = K_bra[pi][pj];
            const double g_ij = gamma_bra[pi][pj];
            const double Q_ij = Q_bra[pi][pj];
            for(std::size_t k = 0; k < nprims[2]; ++k) {
                const auto pk = pmap2[k];
                for(std::size_t l = 0; l < nprims[3]; ++l) {
                    const auto pl     = pmap3[l];
                    const double K_kl = K_ket[pk][pl];
                    const double g_kl = gamma_ket[pk][pl];
                    const double Q_kl = Q_ket[pk][pl];
                    if(!quartet_skipped(K_ij, K_kl, Q_ij, Q_kl, g_ij, g_kl,
                                        thresh)) {
                        continue;
                    }
                    ++out.n_skip;
                    out.sum_coarse += K_ij * K_kl;
                    out.sum_fine +=
                      std::abs(Q_ij * Q_kl / std::sqrt(g_ij + g_kl));
                }
            }
        }
    }
    return out;
}

double buffer_sum(const simde::type::tensor& t) {
    using tensorwrapper::buffer::get_raw_data;
    const auto buf = get_raw_data<double>(t.buffer());
    double s       = 0.0;
    for(double x : buf) s += x;
    return s;
}

} // namespace

TEST_CASE("PrimitiveErrorModel") {
    auto mm   = initialize_integrals();
    auto& mod = mm.at("Primitive Error Model");
    simde::type::v_ee_type v_ee{};

    auto run_mode = [&](auto&& bk, double tol, const std::string& mode) {
        auto copy = mod.unlocked_copy();
        copy.change_input("Error estimate", mode);
        return copy.run_as<pt>(bk, tol);
    };

    SECTION("aggregate sums match reference skip totals (H2 STO-3G)") {
        auto aobs = h2_sto3g_basis_set();
        simde::type::aos_squared bra(aobs, aobs);
        simde::type::aos_squared ket(aobs, aobs);
        chemist::braket::BraKet mnls(bra, v_ee, ket);

        const double tol = 1e-12;
        auto ref         = reference_skip_totals(aobs, aobs, aobs, aobs, tol);

        auto t_tol    = run_mode(mnls, tol, "Tolerance");
        auto t_coarse = run_mode(mnls, tol, "Coarse");
        auto t_fine   = run_mode(mnls, tol, "Fine");

        REQUIRE(
          buffer_sum(t_tol) ==
          Catch::Approx(static_cast<double>(ref.n_skip) * tol).margin(1e-20));
        REQUIRE(buffer_sum(t_coarse) ==
                Catch::Approx(ref.sum_coarse).margin(1e-10));
        REQUIRE(buffer_sum(t_fine) ==
                Catch::Approx(ref.sum_fine).margin(1e-10));
    }

    SECTION("invalid Error estimate throws") {
        auto aobs = h2_sto3g_basis_set();
        simde::type::aos_squared bra(aobs, aobs);
        simde::type::aos_squared ket(aobs, aobs);
        chemist::braket::BraKet mnls(bra, v_ee, ket);

        auto copy = mod.unlocked_copy();
        copy.change_input("Error estimate", std::string("not_a_mode"));
        REQUIRE_THROWS_AS(copy.run_as<pt>(mnls, 1e-10), std::invalid_argument);
    }
}
