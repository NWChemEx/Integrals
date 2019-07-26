#pragma once
#include "integrals/integralpimpl.hpp"
#include "integrals/libint_integral.hpp"
#include <libchemist/basis_set_map.hpp>

namespace integrals::libint::detail_ {

/// Functor that handles the work of putting the results from the integral
/// engine calculation into the appropriate places in a TAMM tensor.
template<libint2::Operator op, std::size_t NBases, typename element_type>
struct TAMMIntFunctor {
    using base_type        = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO         = typename base_type::tiled_AO;
    using basis_array_type = typename base_type::basis_array_type;
    using size_type        = std::size_t;
    using fxn_type         = typename base_type::fxn_type;
    using matrix_type =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    constexpr static size_type nopers = libint2::operator_traits<op>::nopers;

    // tensor level information
    fxn_type fxn; // libint functor
    const tiled_AO tAO; // tiled index spaces
    const std::array<std::vector<size_type>, NBases> atom_blocks; // atom indices for the start of each tile
    const basis_array_type bases; // basis sets
    std::array<libchemist::BasisSetMap, NBases> maps;
    size_type iopers; // index to libint buffer for multicomponent (i.e. x,y,z) integrals
    const element_type schwarz_thresh; // integral screening theshold
    matrix_type Scr; // Scr(i,j) = <ij|op|ij>

    // block level information
    std::array<size_type, NBases> idx; // tensor blockid with multicomponent index removed
    typename fxn_type::shell_index shells;
    std::array<size_type, NBases> ao_off; //
    std::array<typename libchemist::BasisSetMap::range, NBases> ao_ranges; // The corresponding AO ranges for shells
    size_type x_libint; // libint buffer index
    tamm::span<element_type> tamm_buf; // the tamm::Tensor buffer for a block

    TAMMIntFunctor(
      const tiled_AO& tAO,
      const std::array<std::vector<size_type>, NBases>& atom_blocks,
      const basis_array_type& bases, fxn_type&& libint_fxn,
      const element_type schwarz_thresh);


    void operator()(const tamm::IndexVector& blockid,
                    tamm::span<element_type> buff);

    void fxn_call(size_type depth);

    template <bool zero>
    void fill(size_type x_tamm, size_type depth);
};

extern template class TAMMIntFunctor<libint2::Operator::overlap, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::kinetic, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::nuclear, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::coulomb, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::coulomb, 3, double>;
extern template class TAMMIntFunctor<libint2::Operator::coulomb, 4, double>;
extern template class TAMMIntFunctor<libint2::Operator::stg, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::stg, 3, double>;
extern template class TAMMIntFunctor<libint2::Operator::stg, 4, double>;
extern template class TAMMIntFunctor<libint2::Operator::yukawa, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::yukawa, 3, double>;
extern template class TAMMIntFunctor<libint2::Operator::yukawa, 4, double>;
extern template class TAMMIntFunctor<libint2::Operator::emultipole1, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::emultipole2, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::emultipole3, 2, double>;

} // namespace integrals::libint::detail_
