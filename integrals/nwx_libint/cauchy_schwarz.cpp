#include "integrals/nwx_libint/cauchy_schwarz.hpp"
#include "integrals/nwx_libint/nwx_libint.hpp"

namespace nwx_libint {

    template<std::size_t NBases, libint2::Operator op>
    void CauchySchwarz<NBases, op>::initialize(const basis_vec& basis_sets, factory_type& factory) {

        // Handle different potential cases
        if constexpr (NBases == 2) {
            cs_mat1 = make_mat(basis_sets[0], factory);

            if (basis_sets[0] == basis_sets[1]) {
                cs_mat2 = cs_mat1; // don't waste time if these are the same
            } else {
                cs_mat2 = make_mat(basis_sets[1], factory);
            }
        } else if constexpr (NBases == 3) {
            cs_mat1 = make_mat(basis_sets[0], factory);
            cs_mat2 = make_mat(basis_sets[1], basis_sets[2], factory);
        } else if constexpr (NBases == 4) {
            cs_mat1 = make_mat(basis_sets[0], basis_sets[1], factory);

            if ((basis_sets[0] == basis_sets[1]) && (basis_sets[2] == basis_sets[3])) {
                cs_mat2 = cs_mat1; // Same as before
            } else {
                cs_mat2 = make_mat(basis_sets[2], basis_sets[3], factory);
            }
        }
    }

    template<std::size_t NBases, libint2::Operator op>
    bool CauchySchwarz<NBases, op>::tile(const basis_vec& basis_sets,
                                         const TiledArray::Range& range,
                                         double cs_thresh) {
        // Shell lists for current tile
        std::vector<size_vec> shell_list;

        // Fill in the above vectors
        for (int i = 0; i < NBases; ++i) {
            shell_list.push_back(nwx_libint::aos2shells(basis_sets[i], range.lobound()[i], range.upbound()[i]));
        }

        // Find the largest coefficient in each submatrix and use them for the approximation
        if constexpr (NBases == 2) {
            double maxVal1 = 0.0, maxVal2 = 0.0;

            for (auto i : shell_list[0]) { maxVal1 = std::max(maxVal1, cs_mat1[0][i]); }
            for (auto i : shell_list[1]) { maxVal1 = std::max(maxVal1, cs_mat2[0][i]); }

            return (maxVal1 * maxVal2) < cs_thresh;
        } else if constexpr (NBases == 3) {
            double maxVal1 = 0.0, maxVal2 = 0.0;

            for (auto i : shell_list[0]) { maxVal1 = std::max(maxVal1, cs_mat1[0][i]); }
            for (auto i : shell_list[1]) {
                for (auto j : shell_list[2]) {
                    maxVal2 = std::max(maxVal2, cs_mat2[i][j]);
                }
            }

            return (maxVal1 * maxVal2) < cs_thresh;
        } else if constexpr (NBases == 4) {
            double maxVal1 = 0.0, maxVal2 = 0.0;

            for (auto i : shell_list[0]) {
                for (auto j : shell_list[1]) {
                    maxVal1 = std::max(maxVal1, cs_mat1[i][j]);
                }
            }
            for (auto i : shell_list[2]) {
                for (auto j : shell_list[3]) {
                    maxVal2 = std::max(maxVal2, cs_mat2[i][j]);
                }
            }

            return (maxVal1 * maxVal2) < cs_thresh;
        }
    }

    template<std::size_t NBases, libint2::Operator op>
    bool CauchySchwarz<NBases, op>::shellset(size_vec shells, double cs_thresh) {
        // Check approximation product vs provided threshold value
        if constexpr (NBases == 2) {
            return (cs_mat1[0][shells[0]] * cs_mat2[0][shells[1]]) < cs_thresh;
        } else if constexpr (NBases == 3) {
            return (cs_mat1[0][shells[0]] * cs_mat2[shells[1]][shells[2]]) < cs_thresh;
        } else if constexpr (NBases == 4) {
            return (cs_mat1[shells[0]][shells[1]] * cs_mat2[shells[2]][shells[3]]) < cs_thresh;
        }
        return false;
    }

    template<std::size_t NBases, libint2::Operator op>
    double CauchySchwarz<NBases, op>::cs_approx(const shell_vec& shells, libint2::Engine engine) {
        // Buffer of results
        const auto& buf = engine.results();

        // Appropriate call to engine to compute values
        if (shells.size() == 1) {
            engine.set(libint2::BraKet::xs_xs);
            engine.compute(shells[0], shells[0]);
        } else if (shells.size() == 2) {
            engine.set(libint2::BraKet::xx_xx);
            engine.compute(shells[0], shells[1], shells[0], shells[1]);
        }
        auto vals = buf[0];

        // Determine the number of compute values
        std::size_t nvals = 1;
        for (const auto& shell : shells) { nvals *= (shell.size() * shell.size()); }

        // Find the norm and take the square root
        double infinity_norm = 0.0;
        for (int i = 0; i < nvals; ++i) { infinity_norm = std::max(infinity_norm, std::abs(vals[i])); }
        infinity_norm = std::sqrt(infinity_norm);

        return infinity_norm;
    }

    template<std::size_t NBases, libint2::Operator op>
    auto CauchySchwarz<NBases, op>::make_mat(const basis_type& bs, factory_type& factory) {
        auto& world = TA::get_default_world();
        if (not libint2::initialized()) { libint2::initialize(); } // in case it was finalized

        // Place for values to go
        approx_vec mat(1, double_vec(bs.size(), 0.0));

        // Lambda to fill in the values
        auto into_mat = [&] (int i) { mat[0][i] = cs_approx({bs[i]}, factory(NBases, op)); };

        // Calculate values
        for (int i = 0; i < bs.size(); ++i) { world.taskq.add(into_mat, i); }
        world.gop.fence();

        return mat;
    }

    template<std::size_t NBases, libint2::Operator op>
    auto CauchySchwarz<NBases, op>::make_mat(const basis_type& bs1, const basis_type& bs2, factory_type& factory) {
        auto& world = TA::get_default_world();
        if (not libint2::initialized()) { libint2::initialize(); } // in case it was finalized

        // Place for values to go
        approx_vec mat(bs1.size(), double_vec(bs2.size(), 0.0));

        // Check if the basis sets are the same
        bool same_bs = (bs1 == bs2);

        // Lambda to fill in the values
        auto into_mat = [&] (int i, int j) {
            mat[i][j] = cs_approx({bs1[i], bs2[j]}, factory(NBases, op));
            if (same_bs && (i!=j)) { mat[j][i] = mat[i][j]; } // cut down on work
        };

        // Calculate values
        for (int i = 0; i < bs1.size(); ++i) {
            auto len = (same_bs) ? i : bs2.size() - 1; // only do lower triangle, since it's mirrored
            for (int j = 0; j <= len; ++j) {
                world.taskq.add(into_mat, i, j);
            }
        }
        world.gop.fence();

        return mat;
    }

    template<std::size_t NBases, libint2::Operator op>
    void CauchySchwarz<NBases, op>::set_sub_screen(const basis_vec& basis_sets, const TiledArray::Range& range,
                                                   approx_vec& mat1, approx_vec& mat2) const {
        if (cs_mat1.empty() || cs_mat2.empty()) {
            // Default out if cs_mats aren't set to anything
            mat1 = {}, mat2 = {};
        } else {
            // Shell lists for current tile
            std::vector<size_vec> shell_list;

            // Fill in the above vectors
            for (int i = 0; i < NBases; ++i) {
                shell_list.push_back(nwx_libint::aos2shells(basis_sets[i], range.lobound()[i], range.upbound()[i]));
            }

            if constexpr (NBases == 2) {
                mat1 = approx_vec(1, double_vec(shell_list[0].size(), 0.0));
                mat2 = approx_vec(1, double_vec(shell_list[1].size(), 0.0));

                for (int i = 0; i < shell_list[0].size(); ++i) { mat1[0][i] = cs_mat1[0][shell_list[0][i]]; }
                for (int i = 0; i < shell_list[1].size(); ++i) { mat2[0][i] = cs_mat2[0][shell_list[1][i]]; }

            } else if constexpr (NBases == 3) {
                mat1 = approx_vec(1, double_vec(shell_list[0].size(), 0.0));
                mat2 = approx_vec(shell_list[1].size(), double_vec(shell_list[2].size(), 0.0));
                for (int i = 0; i < shell_list[0].size(); ++i) { mat1[0][i] = cs_mat1[0][shell_list[0][i]]; }
                for (int i = 0; i < shell_list[1].size(); ++i) {
                    for (int j = 0; j < shell_list[2].size(); ++j) {
                        mat2[i][j] = cs_mat2[shell_list[1][i]][shell_list[2][j]];
                    }
                }

            } else if constexpr (NBases == 4) {
                mat1 = approx_vec(shell_list[0].size(), double_vec(shell_list[1].size(), 0.0));
                mat2 = approx_vec(shell_list[2].size(), double_vec(shell_list[3].size(), 0.0));
                for (int i = 0; i < shell_list[0].size(); ++i) {
                    for (int j = 0; j < shell_list[1].size(); ++j) {
                        mat1[i][j] = cs_mat1[shell_list[0][i]][shell_list[1][j]];
                    }
                }
                for (int i = 0; i < shell_list[2].size(); ++i) {
                    for (int j = 0; j < shell_list[3].size(); ++j) {
                        mat2[i][j] = cs_mat2[shell_list[2][i]][shell_list[3][j]];
                    }
                }
            }
        }
    }

    template class CauchySchwarz<2, libint2::Operator::overlap>;
    template class CauchySchwarz<2, libint2::Operator::kinetic>;
    template class CauchySchwarz<2, libint2::Operator::nuclear>;
    template class CauchySchwarz<2, libint2::Operator::coulomb>;
    template class CauchySchwarz<3, libint2::Operator::coulomb>;
    template class CauchySchwarz<4, libint2::Operator::coulomb>;
    template class CauchySchwarz<2, libint2::Operator::stg>;
    template class CauchySchwarz<3, libint2::Operator::stg>;
    template class CauchySchwarz<4, libint2::Operator::stg>;
    template class CauchySchwarz<2, libint2::Operator::yukawa>;
    template class CauchySchwarz<3, libint2::Operator::yukawa>;
    template class CauchySchwarz<4, libint2::Operator::yukawa>;
    template class CauchySchwarz<2, libint2::Operator::emultipole1>;
    template class CauchySchwarz<2, libint2::Operator::emultipole2>;
    template class CauchySchwarz<2, libint2::Operator::emultipole3>;
    template class CauchySchwarz<4, libint2::Operator::delta>;

} // namespace nwx_libint