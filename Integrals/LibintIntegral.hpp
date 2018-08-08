#pragma once
#include "Integrals/nwx_libint/nwx_libint.hpp"
#include <LibChemist/Defaults/PropertyTypes.hpp>
#include <memory>

/** @file LibintIntegral.hpp
 *
 *  This file contains the API for a module that builds a TAMM tensor filled
 *  with integrals generated from Libint.
 *
 *  The actual guts of the module are controlled and implemented *via* the PIMPL
 *  idiom.  Note tha even though our public API is templated this is possible so
 *  long as we forward declare all necessary instantiations.
 *
 */
namespace Integrals {
namespace detail_ {

template<libint2::Operator op, std::size_t NBases,typename element_type = double>
class LibintIntegralPIMPL;

/**
 * @brief The class that actually serves as the module for computing integrals
 *        for NWChemEx.
 *
 *
 * @tparam op The libint2 operator enumeration involved in the integral.
 * @tparam NBases The total number of AO bases in the bra and ket
 * @tparam Deriv The derivative order to compute
 * @tparam element_type The literal type of the elements in the tensor
 */
template<libint2::Operator op, std::size_t NBases,typename element_type = double>
struct LibintIntegral : LibChemist::AOIntegral<NBases, element_type> {
    /// Typedef of base class
    using base_type = LibChemist::AOIntegral<NBases, element_type>;

    /// Pull typdefs from baseclass into scope.
    ///@{
    using tensor_type = typename base_type::tensor_type;
    using molecule_type = typename base_type::molecule_type;
    using basis_array_type = typename base_type::basis_array_type;
    using size_type        = typename base_type::size_type;
    ///@}

    LibintIntegral();
    ~LibintIntegral() noexcept;


    /**
     * @brief Computes the tensor representation of an operator in the
     *        requested AO basis sets.
     *
     * @param mol The molecular system associated with the basis sets.
     * @param bases An array of the basis sets in the integral from left to
     *        right across the braket.
     * @return The tensor representation of the operator.
     * @throw ??? If engine construction throws.  Strong throw guarantee.
     *
     * @par Data Races:
     * The content of both @p mol and @p bases will be accessed, if concurrent
     * modifications occur data races will ensue.
     */
    tensor_type run(const molecule_type& mol,
                    const basis_array_type& bases,
                    size_type deriv=0) override;

private:
    std::unique_ptr<LibintIntegralPIMPL<op, NBases, element_type>> pimpl_;
};

extern template class LibintIntegral<libint2::Operator::overlap, 2, double>;


} // namespace detail_

///Typedef for AO overlap matrix
using LibIntOverlap =
  detail_::LibintIntegral<libint2::Operator::overlap, 2, double>;

/////Typedef for kinetic energy of electrons
//template<typename element_type=double>
//using LibIntKinetic =
//  detail_::LibIntIntegral<libint2::Operator::kinetic, 2, element_type>;
//
/////Typedef for the density-fitting metric
//template<typename element_type=double>
//using LibIntMetric =
//    detail_::LibIntIntegral<libint2::Operator::coulomb, 2, element_type>;
//
///// Nucleus-electron attraction
//template<typename element_type=double>
//using LibIntNucleusElectronAttraction =
//    detail_::LibIntIntegral<libint2::Operator::nuclear, 2, element_type>;
//
/////Typedef for the density-fitting metric
//template<typename element_type=double>
//using LibIntDF3C2E =
//    detail_::LibIntIntegral<libint2::Operator::coulomb, 3, element_type>;
//
///// Standard 4 electron integral
//template<typename element_type=double>
//using LibIntERI =
//    detail_::LibIntIntegral<libint2::Operator::coulomb, 4, element_type>;

} // namespace Integrals
