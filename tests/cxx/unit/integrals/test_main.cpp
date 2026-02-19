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

#include <iomanip>
#include <libint2.hpp>
#include <vector>

namespace {

/// Makes 4 s shells each with the same three primitives
auto make_shells() {
    using float_type         = double;
    using vector_type        = libint2::svector<float_type>;
    using contraction_type   = libint2::Shell::Contraction;
    using contraction_vector = libint2::svector<contraction_type>;
    using shell_type         = libint2::Shell;
    using atom_bases_type    = std::vector<shell_type>;

    const int l        = 0; // s shell
    const bool is_pure = false;
    vector_type coeffs{0.1543289673, 0.5353281423, 0.4446345422};
    contraction_type contraction{l, is_pure, coeffs};
    contraction_vector conts(1, contraction);

    vector_type exponents{3.425250914, 0.6239137298, 0.1688554040};
    atom_bases_type rv;
    std::array<float_type, 3> r0{0.0, -0.143222342980786, 0.0};
    std::array<float_type, 3> r1{4.638033502034240, 4.136556880358410, 3.0};
    std::array<float_type, 3> r2{1.638033502034240, 1.136556880358410, 0.0};
    std::array<float_type, 3> r3{3.0, 2.856777657019214, 3.0};

    rv.push_back(shell_type(exponents, conts, r0));
    rv.push_back(shell_type(exponents, conts, r1));
    rv.push_back(shell_type(exponents, conts, r2));
    rv.push_back(shell_type(exponents, conts, r3));

    return rv;
}

/// Same basis set, but decontracted
auto make_decontracted_shells() {
    using float_type         = double;
    using vector_type        = libint2::svector<float_type>;
    using contraction_type   = libint2::Shell::Contraction;
    using contraction_vector = libint2::svector<contraction_type>;
    using shell_type         = libint2::Shell;
    using atom_bases_type    = std::vector<shell_type>;

    const int l        = 0; // s shell
    const bool is_pure = false;

    vector_type coeffs{1.0};
    contraction_type contraction{l, is_pure, coeffs};
    contraction_vector conts(1, contraction);

    vector_type e0{3.425250914};
    vector_type e1{0.6239137298};
    vector_type e2{0.1688554040};

    std::vector<atom_bases_type> rv(4);
    std::array<float_type, 3> r0{0.0, -0.143222342980786, 0.0};
    std::array<float_type, 3> r1{4.638033502034240, 4.136556880358410, 3.0};
    std::array<float_type, 3> r2{1.638033502034240, 1.136556880358410, 0.0};
    std::array<float_type, 3> r3{3.0, 2.856777657019214, 3.0};

    rv[0].push_back(shell_type(e0, conts, r0));
    rv[0].push_back(shell_type(e1, conts, r0));
    rv[0].push_back(shell_type(e2, conts, r0));

    rv[1].push_back(shell_type(e0, conts, r1));
    rv[1].push_back(shell_type(e1, conts, r1));
    rv[1].push_back(shell_type(e2, conts, r1));

    rv[2].push_back(shell_type(e0, conts, r2));
    rv[2].push_back(shell_type(e1, conts, r2));
    rv[2].push_back(shell_type(e2, conts, r2));

    rv[3].push_back(shell_type(e0, conts, r3));
    rv[3].push_back(shell_type(e1, conts, r3));
    rv[3].push_back(shell_type(e2, conts, r3));

    return rv;
}

auto make_k(std::size_t prim_i, std::size_t prim_j, std::size_t center_i,
            std::size_t center_j) {
    using float_type   = double;
    using vector_type  = std::vector<float_type>;
    using point_type   = std::array<float_type, 3>;
    using point_vector = std::vector<point_type>;

    vector_type exponents{3.425250914, 0.6239137298, 0.1688554040};
    point_vector centers{point_type{0.0, -0.143222342980786, 0.0},
                         point_type{4.638033502034240, 4.136556880358410, 3.0},
                         point_type{1.638033502034240, 1.136556880358410, 0.0},
                         point_type{3.0, 2.856777657019214, 3.0}};

    const auto ei = exponents[prim_i];
    const auto ej = exponents[prim_j];
    const auto ri = centers[center_i];
    const auto rj = centers[center_j];

    const auto eij = (ei * ej) / (ei + ej);
    const auto dx  = ri[0] - rj[0];
    const auto dy  = ri[1] - rj[1];
    const auto dz  = ri[2] - rj[2];
    const auto dr  = std::sqrt(dx * dx + dy * dy + dz * dz);
    return std::exp(-eij * dr * dr);
}

} // namespace

int main(int argc, char* argv[]) {
    using float_type = double;

    libint2::initialize();
    libint2::Shell::do_enforce_unit_normalization(false);

    auto shells = make_shells();

    float_type tol = 1.0E-10;

    libint2::Engine engine(libint2::Operator::coulomb, 3, 0, 0, 1.0E-16);
    const auto& buf_vec = engine.results();
    libint2::Engine screener(libint2::Operator::coulomb, 3, 0, 0, tol);
    const auto& screen_buf_vec = screener.results();

    float_type corr_value          = 0.0;
    float_type corr_screened_value = 0.0;
    engine.compute(shells[0], shells[1], shells[2], shells[3]);
    corr_value = buf_vec[0][0];
    screener.compute(shells[0], shells[1], shells[2], shells[3]);
    corr_screened_value = screen_buf_vec[0][0];

    std::vector<float_type> coeffs{0.1543289673, 0.5353281423, 0.4446345422};

    auto decon_shells    = make_decontracted_shells();
    double prim_value    = 0.0;
    double prim_screened = 0.0;
    double error         = 0.0;
    for(std::size_t i = 0; i < 3; ++i) {
        const auto c0p    = coeffs[i];
        const auto c0     = shells[0].contr[0].coeff[i];
        const auto& prim0 = decon_shells[0][i];

        for(std::size_t j = 0; j < 3; ++j) {
            const auto c1p          = coeffs[j];
            const auto c1           = shells[1].contr[0].coeff[j];
            const auto& prim1       = decon_shells[1][j];
            const auto K01          = c0 * c1 * make_k(i, j, 0, 1);
            const bool K01_is_small = std::abs(K01) < tol;

            for(std::size_t k = 0; k < 3; ++k) {
                const auto c2p    = coeffs[k];
                const auto c2     = shells[2].contr[0].coeff[k];
                const auto& prim2 = decon_shells[2][k];

                for(std::size_t l = 0; l < 3; ++l) {
                    const auto c3p           = coeffs[l];
                    const auto c3            = shells[3].contr[0].coeff[l];
                    const auto& prim3        = decon_shells[3][l];
                    const auto K23           = c2 * c3 * make_k(k, l, 2, 3);
                    const bool K23_is_small  = std::abs(K23) < tol;
                    const bool prod_is_small = std::abs(K01 * K23) < tol;

                    engine.compute(prim0, prim1, prim2, prim3);
                    if(buf_vec[0] == nullptr) continue;

                    auto p_ijkl = c0p * c1p * c2p * c3p * buf_vec[0][0];

                    if(K01_is_small || K23_is_small || prod_is_small) {
                        error += p_ijkl;
                    } else {
                        prim_screened += p_ijkl;
                    }
                    prim_value += p_ijkl;
                }
            }
        }
    }

    std::cout << std::fixed << std::setprecision(16);

    auto diff_value = std::abs(corr_value - prim_value);
    std::cout << "Correct value:       " << corr_value << std::endl;
    std::cout << "My contracted value: " << prim_value << std::endl;
    std::cout << "Difference:          " << diff_value << std::endl;

    auto diff_screen = std::abs(corr_screened_value - prim_screened);
    std::cout << "Correct screened value: " << corr_screened_value << std::endl;
    std::cout << "My screened value:      " << prim_screened << std::endl;
    std::cout << "Difference:             " << diff_screen << std::endl;

    auto corr_error = std::abs(corr_screened_value - corr_value);
    auto diff_error = std::abs(corr_error - error);
    std::cout << "Correct error: " << corr_error << std::endl;
    std::cout << "My error:      " << error << std::endl;
    std::cout << "Difference:    " << diff_error << std::endl;

    libint2::finalize();

    return 0;
}
