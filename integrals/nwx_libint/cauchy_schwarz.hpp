#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/nwx_libint/nwx_libint_factory.hpp"

namespace nwx_libint {

    // Class to handle the Cauchy-Schwarz approximations for an integral
    template<std::size_t NBases, libint2::Operator op>
    struct CauchySchwarz {
        using size_type = integrals::type::size;
        using shell_vec = std::vector<libint2::Shell>;
        using basis_type = libint2::BasisSet;
        using basis_vec = std::vector<basis_type>;
        using size_vec = std::vector<size_type >;
        using factory_type = nwx_libint::LibintFactory;
        using double_vec = std::vector<double>;
        using approx_vec = std::vector<double_vec>;

        // Matrices to hold the approximations for the two sides of the integral
        approx_vec cs_mat1 = {};
        approx_vec cs_mat2 = {};

        /** @brief Fill in screening matrices
         *
         *  @param basis_sets The basis sets of the integral
         *  @param factory A factory that produces appropriate LibInt engines
         */
        void initialize(const basis_vec& basis_sets, factory_type& factory);

        /** @brief Check if a whole tile passes screening
         *
         *  @param basis_sets The basis sets of the integral
         *  @param range The range of the tile
         *  @param cs_thresh The screening threshold
         *  @returns true if the tile is screened out
         */
        bool tile(const basis_vec& basis_sets, const TiledArray::Range& range, double cs_thresh);

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
         *  @param engine A LibInt2 engine
         *  @returns The square root of the norm of the approximate integral
         */
        double cs_approx(const shell_vec& shells, libint2::Engine engine);

        /** @brief Make matrix with approximation values for two-index integral
         *
         *  @param bs The LibInt2 basis set
         *  @param factory A factory that produces appropriate LibInt engines
         *  @returns The matrix with the approximation values
         */
        auto make_mat(const basis_type& bs, factory_type& factory);

        /** @brief Make matrix with approximation values for four-index integral
         *
         *  @param bs1 One of the LibInt2 basis sets
         *  @param bs2 The other LibInt2 basis set
         *  @param factory A factory that produces appropriate LibInt engines
         *  @returns The matrix with the approximation values
         */
        auto make_mat(const basis_type& bs1, const basis_type& bs2, factory_type& factory);

        /** @brief Set input vectors as tile specific subset of this objects screening vectors
         *
         *  @param basis_sets The basis sets of the integral
         *  @param range The range of the tile
         *  @param mat1 The vector that is a subset of cs_mat1
         *  @param mat2 The vector that is a subset of cs_mat2
         */
        void set_sub_screen(const basis_vec& basis_sets, const TiledArray::Range& range,
                            approx_vec& mat1, approx_vec& mat2) const;

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