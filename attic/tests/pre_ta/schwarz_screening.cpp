/*
 * Copyright 2023 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

#include "test_common.hpp"
#include <integrals/schwarz_screening.hpp>

using namespace integrals::libint::detail_;

std::vector<std::vector<double>> corr {
{2.1874792352627042, 0.3699640325144978, 0.1564525878919850, 0.0606886000006391, 0.0606886000006391},
{0.3699640325144978, 0.9039946468459082, 0.4248745604347619, 0.3222057218175323, 0.3222057218175323},
{0.1564525878919850, 0.4248745604347619, 0.9381679451862797, 0.2740491521872362, 0.2740491521872362},
{0.0606886000006391, 0.3222057218175323, 0.2740491521872362, 0.8801170058122373, 0.1336470421634229},
{0.0606886000006391, 0.3222057218175323, 0.2740491521872362, 0.1336470421634229, 0.8801170058122373},
};

const double eps  = 10000 * std::numeric_limits<double>::epsilon();
const double marg = 100 * std::numeric_limits<double>::epsilon();

TEST_CASE("Testing Schwarz Screening functions") {

    auto[molecule, bs] = make_molecule();
    auto scr = schwarz_screening<libint2::Operator::coulomb>(bs, bs);

    for (std::size_t i = 0; i < 5; ++i) {
        for (std::size_t j = 0; j < 5; ++j) {
            REQUIRE(scr(i,j == Approx(corr[i][j]).epsilon(eps).margin(marg)));
        }
    }

    REQUIRE(schwarz_estimate<4>(scr, {0,0,0,0}, 0.0) == false);
    REQUIRE(schwarz_estimate<4>(scr, {0,0,0,0}, 12.0) == true);
}