#pragma once
#include "basis_set_serialize.hpp"
#include "cauchy_schwarz_screener.hpp"
#include "integrals/types.hpp"
#include "nwx_libint_factory.hpp"
#include <libint2.hpp>
#include <tiledarray.h>

namespace nwx_TA {

template<typename val_type, libint2::Operator op, std::size_t NBases>
struct FillNDFunctor {
    using size_type    = integrals::type::size;
    using basis_vec    = std::vector<libint2::BasisSet>;
    using size_vec     = std::vector<size_type>;
    using element_type = typename val_type::numeric_type;

    // The collected LibInt2 basis sets needed for the integral
    basis_vec LIBasis_sets;

    // The factory that produces the appropriate LibInt2 engines
    nwx_libint::LibintFactory factory = nwx_libint::LibintFactory();

    // Cauchy-Schwarz Screening Threshold
    double cs_thresh = 0.0;
    nwx_libint::CauchySchwarzScreener<NBases> screen =
      nwx_libint::CauchySchwarzScreener<NBases>();

    // Number of arrays returned by operator
    size_type nopers = libint2::operator_traits<op>::nopers;

    void initialize(const basis_vec& sets, size_type deriv, element_type thresh,
                    element_type cs_thresh);

    /** @brief The top level function that starts the recursive calls of the
     * other functions. Core version.
     *
     *  Gets the tile @p tile and range @p range from TA, then initializes the
     * tile and declares the vectors for the offsets and shells. Complies with
     * the TA API for these functions
     *
     *  @param[in] tile The tile to be filled
     *  @param[in] range The range of the tile
     *  @returns The norm of the filled tile
     */
    float operator()(val_type& tile, const TiledArray::Range& range);

    /** @brief The top level function that starts the recursive calls of the
     * other functions. Direct version.
     *
     *  Assumes that LIBasis_sets only contains the shells needed for the
     * current tile.
     *
     *  @param[in] range The range of the tile
     *  @returns The filled tile
     */
    val_type operator()(const TiledArray::Range& range);

    /** @brief Serialize this object.
     *
     *  @tparam Archive the archive type
     *  @param ar The archive
     */
    template<typename Archive>
    void serialize(Archive& ar) {
        ar& nopers;
        ar& cs_thresh;

        ar& screen.cs_mat1;
        ar& screen.cs_mat2;

        ar& factory.max_nprims;
        ar& factory.max_l;
        ar& factory.thresh;
        ar& factory.deriv;
        ar& factory.stg_exponent;
        ar& factory.origin;
        ar& factory.qs;

        serialize_basis_sets(ar, LIBasis_sets);
    }

private:
    /** @brief Recursive function that transverses all of the dimensions of the
     * current tile and finds the shells that need to be computed to fill the
     * tile.
     *
     *  @param tile The tile to be filled
     *  @param range The range of the tile
     *  @param tile_engine The LibInt2 engine that computes integrals
     *  @param offsets Vector containing the coordinate offset based on the
     * current AOs
     *  @param shells Vector tracking the indices of the current shells
     *  @param depth The current dimension of the desired tensor
     */
    void _index_shells(val_type& tile, const TiledArray::Range& range,
                       libint2::Engine& tile_engine, size_vec& offsets,
                       size_vec& shells, std::vector<size_vec>& tile_shells,
                       int depth);

    /** @brief Recursive function that transverses all of the dimensions of the
     * current tile and fills the LibInt2 results.
     *
     *  @param tile The tile to be filled
     *  @param offsets Vector containing the coordinate offset based on the
     * current AOs
     *  @param shells Vector tracking the indices of the current shells
     *  @param int_vals Pointer to the beginning of the integral values
     *  @param indexer Index for placing values into array
     *  @param int_i Index for accessing the current integral value to be filled
     *  @param depth The current dimension of the desired tensor
     */
    void _fill_from_libint(val_type& tile, const size_vec& offsets,
                           const size_vec& shells, const double* int_vals,
                           size_vec& indexer, int& int_i, int depth);

    /** @brief Wrap the call of LibInt2 engine so it can take a variable number
     * of shell inputs.
     *
     * @tparam Is A variadic parameter pack of integers from [0,NBases) to
     * expand.
     * @param tile_engine The LibInt2 engine that computes integrals
     * @param shells The index of the requested shell block
     * @return An std::vector filled with the requested block per operator
     * component
     */
    template<std::size_t... Is>
    void _call_libint(libint2::Engine& tile_engine, size_vec shells,
                      std::index_sequence<Is...>);
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
extern template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole1,
                                    2>;
extern template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole2,
                                    2>;
extern template class FillNDFunctor<TA::TensorD, libint2::Operator::emultipole3,
                                    2>;
extern template class FillNDFunctor<TA::TensorD, libint2::Operator::delta, 4>;

} // namespace nwx_TA