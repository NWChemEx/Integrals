#pragma once
#include "Integrals/LibintIntegral.hpp"
#include "Integrals/IntegralPIMPL.hpp"
#include <LibChemist/BasisSetMap.hpp>

namespace Integrals::Libint::detail_ {

/// Functor that handles the work of putting the results from the integral engine calculation into the appropriate
/// places in a TAMM tensor.
template<libint2::Operator op, std::size_t NBases, typename element_type>
struct TAMMIntFunctor {
    using base_type = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO = typename base_type::tiled_AO;
    using basis_array_type = typename base_type::basis_array_type;
    using size_type = std::size_t;
    using fxn_type = typename base_type::fxn_type;

    using matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using range_array = std::array<typename LibChemist::BasisSetMap::range,NBases>;

    constexpr static size_type nopers = libint2::operator_traits<op>::nopers;

    const tiled_AO tAO;
    const std::array<std::vector<size_type>, NBases> atom_blocks;
    const basis_array_type bases;
    fxn_type fxn;
    bool screen = false;
    matrix Scr;

    TAMMIntFunctor(const tiled_AO &tAO,
                   const std::array<std::vector<size_type>, NBases> &atom_blocks,
                   const basis_array_type &bases,
                   fxn_type &&fxn);

    void operator()(const tamm::IndexVector& blockid, tamm::span<element_type> buff);

    void fxn_call(typename fxn_type::shell_index& shells,
                  range_array& ao_ranges,
                  size_type depth,
                  std::array<LibChemist::BasisSetMap, NBases>& maps,
                  std::array<size_type, NBases>& idx,
                  size_type off_nopers,
                  size_type iopers,
                  std::array<size_type, NBases>& ao_off,
                  tamm::span<element_type> tamm_buf);

    void fill(size_type x_tamm,
              size_type x_libint,
              size_type depth,
              size_type off,
              std::array<size_type, NBases>& idx,
              range_array& ao_ranges,
              std::array<size_type, NBases>& ao_off,
              tamm::span<element_type> tamm_buf,
              std::vector<double>& libint_buf);
};

extern template class TAMMIntFunctor<libint2::Operator::overlap, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::kinetic, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::nuclear, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::coulomb, 2, double>;
//extern template class TAMMIntFunctor<libint2::Operator::coulomb, 3, double>;
//extern template class TAMMIntFunctor<libint2::Operator::coulomb, 4, double>;
extern template class TAMMIntFunctor<libint2::Operator::emultipole1, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::emultipole2, 2, double>;
extern template class TAMMIntFunctor<libint2::Operator::emultipole3, 2, double>;

}
