#pragma once
#include <libint2.hpp>
#include <tiledarray.h>
#include "integrals/nwx_libint/nwx_libint_factory.hpp"
#include "integrals/nwx_TA/nwx_TA_utils.hpp"

namespace nwx_TA {

    template<typename val_type, libint2::Operator op, std::size_t NBases>
    struct FillNDFunctor {

        using basis_vec = std::vector<libint2::BasisSet>;
        using size_vec = std::vector<std::size_t>;

        // The collected LibInt2 basis sets needed for the integral
        basis_vec LIBasis_sets;

        // The factory that produces the appropriate LibInt2 engines
        nwx_libint::LibintFactory<NBases, op> factory;

        // Number of arrays returned by operator
        std::size_t nopers = libint2::operator_traits<op>::nopers;

        // Initialize and finalize LibInt2
        FillNDFunctor() { libint2::initialize(); }
        ~FillNDFunctor() { libint2::finalize(); }

        // Complies with the TA API for these functions
        float operator()(val_type& tile, const TiledArray::Range& range);

    private:

        /** @brief The top level function that starts the recursive calls of the other functions.
         *         Gets the tile @p tile and range @p range from TA, then initializes the tile and
         *         declares the vectors for the offsets and shells.
         *
         *  @param[in] tile The tile to be filled
         *  @param[in] range The range of the tile
         *  @returns The norm of the filled tile
         */
        float _fill(val_type& tile, const TiledArray::Range& range);

        /** @brief Recursive function that transverses all of the dimensions of the current tile
         *         and finds the shells that need to be computed to fill the tile.
         *
         *  @param tile The tile to be filled
         *  @param range The range of the tile
         *  @param tile_engine The LibInt2 engine that computes integrals
         *  @param offsets Vector containing the coordinate offset based on the current AOs
         *  @param shells Vector tracking the indices of the current shells
         *  @param depth The current dimension of the desired tensor
         */
        void _index_shells(val_type& tile,
                           const TiledArray::Range& range,
                           libint2::Engine& tile_engine,
                           size_vec& offsets,
                           size_vec& shells,
                           int depth);

        /** @brief Recursive function that transverses all of the dimensions of the current tile
         *         and fills the LibInt2 results into the correct coordinate position in the tile
         *
         *  @param tile The tile to be filled
         *  @param offsets Vector containing the coordinate offset based on the current AOs
         *  @param shells Vector tracking the indices of the current shells
         *  @param int_vals Pointer to the beginning of the integral values
         *  @param indexer Index for placing values into array
         *  @param int_i Index for accessing the current integral value to be filled
         *  @param depth The current dimension of the desired tensor
         */
        void _fill_from_libint(val_type& tile,
                               const size_vec& offsets,
                               const size_vec& shells,
                               const double* int_vals,
                               size_vec& indexer,
                               int& int_i,
                               int depth);

        /** @brief Wrap the call of LibInt2 engine so it can take a variable number of shell inputs.
         *
         * @tparam Is A variadic parameter pack of integers from [0,NBases) to expand.
         * @param tile_engine The LibInt2 engine that computes integrals
         * @param shells The index of the requested shell block
         * @return An std::vector filled with the requested block per operator component
         */
        template<std::size_t... Is>
        void _call_libint(libint2::Engine& tile_engine,
                          size_vec shells,
                          std::index_sequence<Is ...>);
    };

    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::overlap, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::kinetic, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::nuclear, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::coulomb, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::coulomb, 3>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::coulomb, 4>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::stg, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::stg, 3>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::stg, 4>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::yukawa, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::yukawa, 3>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::yukawa, 4>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole1, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole2, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole3, 2>;
    extern template class FillNDFunctor<TA::TensorD, libint2::Operator::delta, 4>;

} // namespace nwx_TA