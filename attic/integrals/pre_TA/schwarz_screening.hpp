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

#pragma once
#include "integrals/libint_functor.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"

namespace integrals::libint::detail_ {

template<libint2::Operator op>
auto schwarz_screening(const chemist::AOBasisSet& bs1,
                       const chemist::AOBasisSet& bs2) {
    using matrix_type =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    const auto nsh1          = bs1.nshells();
    const auto nsh2          = bs2.nshells();
    const bool bs1_equiv_bs2 = (bs1 == bs2);
    matrix_type rv           = matrix_type::Zero(nsh1, nsh2);
    LibintFunctor<4> fxn;

    fxn.bs[0] = nwx_libint::make_basis(bs1);
    fxn.bs[1] = nwx_libint::make_basis(bs2);
    fxn.bs[2] = fxn.bs[0];
    fxn.bs[3] = fxn.bs[1];
    auto max_prims =
            std::max(fxn.bs[0].max_nprim(fxn.bs[0]), fxn.bs[1].max_nprim(fxn.bs[1]));
    auto max_l =
            std::max(fxn.bs[0].max_l(fxn.bs[0]), fxn.bs[1].max_l(fxn.bs[1]));

    fxn.engine =
            make_engine<op, 4>(chemist::Molecule{}, max_prims, max_l, 0.0, 0);

    const auto& buf = fxn.engine.results();

    for(std::size_t s1 = 0, s12 = 0; s1 != nsh1; ++s1) {
        std::size_t n1 =
                bs1[s1].size(); // number of basis functions in this shell

        std::size_t s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
        for(std::size_t s2 = 0; s2 <= s2_max; ++s2, ++s12) {
            std::size_t n2  = bs2[s2].size();
            std::size_t n12 = n1 * n2;
            std::array<std::size_t, 4> shells{s1, s2, s1, s2};
            fxn(shells);

            Eigen::Map<const matrix_type> buf_mat(buf[0], n12, n12);
            auto norm2 = buf_mat.lpNorm<Eigen::Infinity>();
            rv(s1, s2) = std::sqrt(norm2);
            if(bs1_equiv_bs2) rv(s2, s1) = rv(s1, s2);
        }
    }
    return rv;
}

template<std::size_t NBases, typename element_type = double>
bool schwarz_estimate(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat,
                      const std::array<std::size_t, NBases>& shells,
                      const element_type threshold) {
    if constexpr(NBases == 3) {
        return threshold > mat(shells[1], shells[2]);
    } else if constexpr(NBases == 4) {
        return threshold >
               mat(shells[0], shells[1]) * mat(shells[2], shells[3]);
    } else {
        return false;
    }
}

} // end namespace