#pragma once
#include "Integrals/nwx_libint/nwx_libint.hpp"
#include <LibChemist/Defaults/PropertyTypes.hpp>
#include <SDE/ModuleBase.hpp>
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
///Namespace for classes needed to implement integrals using Libint
namespace Libint {
///Namespace for classes that are considered implementation details
namespace detail_ {

///Forward declaration of the implementation of the LibintIntegral class
template<libint2::Operator op, std::size_t NBases,
        typename element_type = double>
class IntegralPIMPL;

///Enum containing the possible implementations
enum class implementation_type {direct, core};

/**
 * @brief The class that actually serves as the module for computing integrals
 *        for NWChemEx.
 *
 * All of the integrals needed by NWX are obtained from this class.  In order to
 * store the definitions in the source file (and avoid template overhead) all
 * explicit instantiations must be included in this header file (see right below
 * this class for the actual instantiations).  This class is largely code
 * factorization as the process of creating a set of integrals is basically
 * identical regardless of what kind of integral is being formed.  Users of the
 * integrals library should deal exclusively with the typedefs outside the
 * detail_ namespace (*e.g.*, Overlap, and Kinetic).
 *
 * @tparam op The libint2 operator enumeration involved in the integral.
 * @tparam NBases The total number of AO bases in the bra and ket
 * @tparam element_type The literal type of the elements in the tensor
 */
template<libint2::Operator op, std::size_t NBases, typename element_type = double>
//struct Integral : public ModuleBase<Integral<op, NBases, element_type>> {
struct Integral : public SDE::ModuleBase {
    /// Typedef of base class
    using base_type = LibChemist::AOIntegral<NBases, element_type>;

    /// Pull typdefs from baseclass into scope.
    ///@{
    using tensor_type      = typename base_type::tensor_type;
    using molecule_type    = typename base_type::molecule_type;
    using basis_array_type = typename base_type::basis_array_type;
    using size_type        = typename base_type::size_type;
    ///@}

    ///Initializes the buffers libint will need
    Integral(implementation_type impl = implementation_type::direct) ;
    ///Frees the buffers libint uses
    ~Integral() noexcept;

    /**
     * @brief Computes the tensor representation of an operator in the
     *        requested AO basis sets.
     *
     * @param mol The molecular system associated with the basis sets.
     * @param bases An array of the basis sets in the integral from left to
     *        right across the braket.
     * @param deriv What derivative of the integral should we compute?
     *        Default is 0 (*i.e.* just the integral).
     *
     * @return The tensor representation of the operator.
     * @throw ??? If engine construction throws.  Strong throw guarantee.
     *
     * @par Data Races:
     * The content of both @p mol and @p bases will be accessed, if concurrent
     * modifications occur data races will ensue.
     */
    /*tensor_type run(const molecule_type& mol,
                    const basis_array_type& bases,
                    size_type deriv=0) override;*/

private:
    ///The object that actually implements this class
    std::unique_ptr<IntegralPIMPL<op, NBases, element_type>> pimpl_;
};


// Explicit instantiations of the Integral class
extern template class Integral<libint2::Operator::overlap, 2, double>;
extern template class Integral<libint2::Operator::kinetic, 2, double>;
extern template class Integral<libint2::Operator::nuclear, 2, double>;
extern template class Integral<libint2::Operator::coulomb, 2, double>;
extern template class Integral<libint2::Operator::coulomb, 3, double>;
extern template class Integral<libint2::Operator::coulomb, 4, double>;
extern template class Integral<libint2::Operator::emultipole1, 2, double>;
extern template class Integral<libint2::Operator::emultipole2, 2, double>;
extern template class Integral<libint2::Operator::emultipole3, 2, double>;

} // namespace detail_

///@defgroup Libint Integral classes
///@{
///The matrix containing the overlap of the AO basis set
using Overlap = detail_::Integral<libint2::Operator::overlap, 2, double>;

///Kinetic energy of electrons
using Kinetic = detail_::Integral<libint2::Operator::kinetic, 2, double>;

/// Nucleus-electron attraction
using Nuclear = detail_::Integral<libint2::Operator::nuclear, 2, double>;

///Typedef for the density-fitting metric
using Metric = detail_::Integral<libint2::Operator::coulomb, 2, double>;

///Typedef for the density-fit ERI
using DFERI = detail_::Integral<libint2::Operator::coulomb, 3, double>;

///Typedef of the canonical ERI
using ERI = detail_::Integral<libint2::Operator::coulomb, 4, double>;

///Electric dipole
using EDipole = detail_::Integral<libint2::Operator::emultipole1, 2, double>;

///Electric quadrupole
using EQuadrupole = detail_::Integral<libint2::Operator::emultipole2, 2, double>;

///Electric octopole
using EOctopole = detail_::Integral<libint2::Operator::emultipole3, 2, double>;

///@}

} //namespace Libint
} // namespace Integrals
