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

#include "../testing/testing.hpp"
#include <integrals/utils/primitive_index_helpers.hpp>
#include <string>
#include <utility>
#include <vector>

using namespace integrals::testing;
using tensorwrapper::operations::approximately_equal;

namespace {

using tensor_type = simde::type::tensor;

/* Adds up the primitive AO contributions into contracted AOs. If kept_only is
 * true, only quartets the screening keeps are added. Returns the contracted
 * ERIs and the number of primitive quartets that were screened out.
 */
auto contract(const simde::type::ao_basis_set& bs,
              const tensor_type& contributions, const tensor_type& reasons,
              bool kept_only) {
    using tensorwrapper::buffer::get_raw_data;
    const auto pcontrib = get_raw_data<double>(contributions.buffer());
    const auto preasons = get_raw_data<double>(reasons.buffer());
    const auto n_aos    = bs.n_aos();
    const auto n_prims  = bs.n_primitives();
    const auto to_cgto  = integrals::utils::build_prim_ao_to_cgto_map(bs);
    const auto to_prim  = integrals::utils::build_prim_ao_to_prim_shell_map(bs);
    const auto n_pao    = to_cgto.size();
    std::size_t n_dropped = 0;
    for(std::size_t i = 0; i < preasons.size(); ++i)
        if(preasons[i] != 0.0) ++n_dropped;

    std::vector<double> ao_data(n_aos * n_aos * n_aos * n_aos, 0.0);
    std::size_t ijkl = 0;
    for(std::size_t i = 0; i < n_pao; ++i)
        for(std::size_t j = 0; j < n_pao; ++j)
            for(std::size_t k = 0; k < n_pao; ++k)
                for(std::size_t l = 0; l < n_pao; ++l, ++ijkl) {
                    const auto pqrs =
                      ((to_prim[i] * n_prims + to_prim[j]) * n_prims +
                       to_prim[k]) *
                        n_prims +
                      to_prim[l];
                    if(kept_only && preasons[pqrs] != 0.0) continue;
                    const auto mnls =
                      ((to_cgto[i] * n_aos + to_cgto[j]) * n_aos + to_cgto[k]) *
                        n_aos +
                      to_cgto[l];
                    ao_data[mnls] += pcontrib[ijkl];
                }

    tensorwrapper::shape::Smooth shape({n_aos, n_aos, n_aos, n_aos});
    tensorwrapper::buffer::Contiguous buffer(std::move(ao_data), shape);
    return std::make_pair(tensor_type(shape, std::move(buffer)), n_dropped);
}

} // namespace

TEST_CASE("Primitive Quartet Screening ERI4") {
    using test_pt = integrals::property_types::PrimitiveQuartetScreening;
    using eri4_pt = simde::ERI4;

    auto mm = initialize_integrals();
    REQUIRE(mm.count("Primitive Quartet Screening ERI4"));
    auto& test_mod = mm.at("Primitive Quartet Screening ERI4");
    auto& corr_mod = mm.at("ERI4");

    auto aobs = water_sto3g_basis_set();
    simde::type::aos aos(aobs);
    simde::type::aos_squared aos_squared(aos, aos);
    simde::type::v_ee_type op{};
    chemist::braket::BraKet braket(aos_squared, op, aos_squared);

    auto eri4_at = [&](double thresh) {
        auto copy = corr_mod.unlocked_copy();
        copy.change_input("Threshold", thresh);
        return copy.run_as<eri4_pt>(braket);
    };

    const auto n_prims    = aobs.n_primitives();
    const auto n_quartets = n_prims * n_prims * n_prims * n_prims;

    // N.b. std::to_string prints 1e-16 and 1e-8 the same way, which would
    // give the sections the same name, so the names are spelled out.
    using pair_type = std::pair<const char*, double>;
    for(const auto& [name, thresh] :
        {pair_type{"1e-16", 1e-16}, pair_type{"1e-8", 1e-8},
         pair_type{"1e-6", 1e-6}}) {
        SECTION(std::string("Tolerance: ") + name) {
            const auto& [contributions, reasons, coarse, fine] =
              test_mod.run_as<test_pt>(braket, thresh);

            using tensorwrapper::buffer::get_raw_data;
            REQUIRE(get_raw_data<double>(reasons.buffer()).size() ==
                    n_quartets);
            REQUIRE(get_raw_data<double>(coarse.buffer()).size() == n_quartets);
            REQUIRE(get_raw_data<double>(fine.buffer()).size() == n_quartets);

            // The kept quartets reproduce libint's screened ERIs
            auto [kept, n_dropped] =
              contract(aobs, contributions, reasons, true);
            REQUIRE(approximately_equal(kept, eri4_at(thresh), 1.0e-14));

            // All quartets reproduce the unscreened ERIs
            auto all = contract(aobs, contributions, reasons, false).first;
            REQUIRE(approximately_equal(all, eri4_at(1.0e-16), 1.0e-14));

            // At a loose tolerance something must be screened out, otherwise
            // the "kept" check above can't tell kept from dropped quartets.
            if(thresh == 1e-6) REQUIRE(n_dropped > 0);
        }
    }
}
