#include "integrals/nwx_libint/cauchy_schwarz.hpp"

namespace nwx_libint {

    template<size_type NBases, libint2::Operator op>
    CauchySchwarz<NBases, op>::CauchySchwarz(nwx_libint::LibintFactory<NBases, op>& factory, basis_vec& bs_vec) :
            basis_sets(bs_vec), factory(factory) {

        // Handle different potential cases
        if constexpr (NBases == 2) {
            cs_mat1 = make_mat(basis_sets[0]);

            if (basis_sets[0] == basis_sets[1]) {
                cs_mat2 = cs_mat1; // don't waste time if these are the same
            } else {
                cs_mat2 = make_mat(basis_sets[1]);
            }
        } else if constexpr (NBases == 3) {
            cs_mat1 = make_mat(basis_sets[0]);
            cs_mat2 = make_mat(basis_sets[1], basis_sets[2]);
        } else if constexpr (NBases == 4) {
            cs_mat1 = make_mat(basis_sets[0], basis_sets[1]);

            if ((basis_sets[0] == basis_sets[1]) && (basis_sets[2] == basis_sets[3])) {
                cs_mat2 = cs_mat1; // Same as before
            } else {
                cs_mat2 = make_mat(basis_sets[2], basis_sets[3]);
            }
        }
    }

    template<size_type NBases, libint2::Operator op>
    bool CauchySchwarz<NBases, op>::tile(const TiledArray::Range& range, double cs_thresh){
        // The corners indices for the submatrix
        size_vec corners;
        // The length of the submatrix in the different dimensions
        size_vec lengths;

        // Fill in the above vectors
        for (int i = 0; i < NBases; ++i) {
            auto shells = nwx_TA::aos2shells(basis_sets[i], range.lobound()[i], range.upbound()[i]);
            corners.push_back(shells[0]);
            lengths.push_back(shells.size());
        }

        // Find the largest coefficient in each submatrix and use them for the approximation
        if constexpr (NBases == 2) {
            auto maxCoeff1 = cs_mat1.block(corners[0], 0, lengths[0], 1).maxCoeff();
            auto maxCoeff2 = cs_mat2.block(corners[1], 0, lengths[1], 1).maxCoeff();

            return (maxCoeff1 * maxCoeff2) < cs_thresh;
        } else if constexpr (NBases == 3) {
            auto maxCoeff1 = cs_mat1.block(corners[0], 0, lengths[0], 1).maxCoeff();
            auto maxCoeff2 = cs_mat2.block(corners[1], corners[2], lengths[1], lengths[2]).maxCoeff();

            return (maxCoeff1 * maxCoeff2) < cs_thresh;
        } else if constexpr (NBases == 4) {
            auto maxCoeff1 = cs_mat1.block(corners[0], corners[1], lengths[0], lengths[1]).maxCoeff();
            auto maxCoeff2 = cs_mat2.block(corners[2], corners[3], lengths[2], lengths[3]).maxCoeff();

            return (maxCoeff1 * maxCoeff2) < cs_thresh;
        }
    }

    template<size_type NBases, libint2::Operator op>
    bool CauchySchwarz<NBases, op>::shellset(size_vec shells, double cs_thresh) {
        // Check approximation product vs provided threshold value
        if constexpr (NBases == 2) {
            return (cs_mat1(shells[0], 0) * cs_mat2(shells[1], 0)) < cs_thresh;
        } else if constexpr (NBases == 3) {
            return (cs_mat1(shells[0], 0) * cs_mat2(shells[1], shells[2])) < cs_thresh;
        } else if constexpr (NBases == 4) {
            return (cs_mat1(shells[0], shells[1]) * cs_mat2(shells[2], shells[3])) < cs_thresh;
        }
        return false;
    }

    template<size_type NBases, libint2::Operator op>
    double CauchySchwarz<NBases, op>::cs_approx(const shell_vec& shells) {
        // Get an engine from the factory
        auto engine = factory();
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

    template<size_type NBases, libint2::Operator op>
    Eigen::MatrixXd CauchySchwarz<NBases, op>::make_mat(const basis_type& bs) {
        auto& world = TA::get_default_world();
        if (not libint2::initialized()) { libint2::initialize(); } // in case it was finalized

        // Place for values to go
        Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(bs.size(), 1);

        // Lambda to fill in the values
        auto into_mat = [&] (int i) { mat(i, 0) = cs_approx({bs[i]}); };

        // Calculate values
        for (int i = 0; i < bs.size(); ++i) { world.taskq.add(into_mat, i); }
        world.gop.fence();

        return mat;
    }

    template<size_type NBases, libint2::Operator op>
    Eigen::MatrixXd CauchySchwarz<NBases, op>::make_mat(const basis_type& bs1, const basis_type& bs2) {
        auto& world = TA::get_default_world();
        if (not libint2::initialized()) { libint2::initialize(); } // in case it was finalized

        // Place for values to go
        Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(bs1.size(), bs2.size());

        // Check if the basis sets are the same
        bool same_bs = (bs1 == bs2);

        // Lambda to fill in the values
        auto into_mat = [&] (int i, int j) {
            mat(i, j) = cs_approx({bs1[i], bs2[j]});
            if (same_bs && (i!=j)) { mat(j, i) = mat(i, j); } // cut down on work
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