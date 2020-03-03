#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"

namespace nwx_libint {

    // Class to handle the Cauchy-Schwarz approximations for an integral
    template<size_type NBases, libint2::Operator op>
    struct CauchySchwarz {
        using shell_vec = std::vector<libint2::Shell>;
        using basis_type = libint2::BasisSet;
        using basis_vec = std::vector<basis_type>;
        using size_vec = std::vector<std::size_t>;

        // Matrices to hold the approximations for the two sides of the integral
        Eigen::MatrixXd cs_mat1;
        Eigen::MatrixXd cs_mat2;

        // References to the basis sets and the integral factory
        basis_vec& basis_sets;
        nwx_libint::LibintFactory<NBases, op>& factory;

        CauchySchwarz(nwx_libint::LibintFactory<NBases, op>& factory, basis_vec& bs_vec);

        /** @brief Check if a whole tile passes screening
         *
         *  @param range The range of the tile
         *  @param cs_thresh The screening threshold
         *  @returns true if the tile is screened out
         */
        bool tile(const TiledArray::Range& range, double cs_thresh);

        /** @brief Check if a shell set passes screening
         *
         *  @param shells The indices of the shells
         *  @param cs_thresh The screening threshold
         *  @returns true if the tile is screened out
         */
        bool shellset(size_vec shells, double cs_thresh);

        /** @brief Calculate the Cauchy-Schwarz approximation for a shell or pair of shells
         *
         *  @param shells The LibInt2 shell(s)
         *  @returns The square root of the norm of the approximate integral
         */
        double cs_approx(const shell_vec& shells);

        /** @brief Make matrix with approximation values for two-index integral
         *
         *  @param bs The LibInt2 basis set
         *  @returns The matrix with the approximation values
         */
        Eigen::MatrixXd make_mat(const basis_type& bs);

        /** @brief Make matrix with approximation values for four-index integral
         *
         *  @param bs1 One of the LibInt2 basis sets
         *  @param bs2 The other LibInt2 basis set
         *  @returns The matrix with the approximation values
         */
        Eigen::MatrixXd make_mat(const basis_type& bs1, const basis_type& bs2);

    }; // Class CauchySchwarz

    extern template class CauchySchwarz<2, libint2::Operator::overlap>;
    extern template class CauchySchwarz<2, libint2::Operator::kinetic>;
    extern template class CauchySchwarz<2, libint2::Operator::nuclear>;
    extern template class CauchySchwarz<2, libint2::Operator::coulomb>;
    extern template class CauchySchwarz<3, libint2::Operator::coulomb>;
    extern template class CauchySchwarz<4, libint2::Operator::coulomb>;
    extern template class CauchySchwarz<2, libint2::Operator::stg>;
    extern template class CauchySchwarz<3, libint2::Operator::stg>;
    extern template class CauchySchwarz<4, libint2::Operator::stg>;
    extern template class CauchySchwarz<2, libint2::Operator::yukawa>;
    extern template class CauchySchwarz<3, libint2::Operator::yukawa>;
    extern template class CauchySchwarz<4, libint2::Operator::yukawa>;
    extern template class CauchySchwarz<2, libint2::Operator::emultipole1>;
    extern template class CauchySchwarz<2, libint2::Operator::emultipole2>;
    extern template class CauchySchwarz<2, libint2::Operator::emultipole3>;
    extern template class CauchySchwarz<4, libint2::Operator::delta>;

} // namespace nwx_libint