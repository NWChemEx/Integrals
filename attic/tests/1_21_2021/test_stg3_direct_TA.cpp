/*
 * Copyright 2022 NWChemEx-Project
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

#include "test_common_TA.hpp"
#include <integrals/integralsmm.hpp>
#include <integrals/nwx_direct/stg_direct_type.hpp>
#include "H2O_STO3G_STG[1]_3C.hpp"

TEST_CASE("STG3CDirect") {
    using integral_type = property_types::STG3CDirect<double>;

    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [molecule, bs] = make_molecule();
    auto stg_exponent = 1.0;
    mm.at("STG3Direct").change_input("Tile size", std::vector<std::size_t>{7});
    mm.at("STG3Direct").change_input("Screening Threshold", 0.000001);
    auto [X] = mm.at("STG3Direct").run_as<integral_type>(bs, bs, bs, std::size_t{0}, stg_exponent);

    TensorType real_X(X.world(), X.trange());
    real_X("k, l, m") = X("k, l, m");

    REQUIRE(chemist::ta_helpers::allclose(real_X, TensorType(real_X.world(), real_X.trange(), corr)));
}