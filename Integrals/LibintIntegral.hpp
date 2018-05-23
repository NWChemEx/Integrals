#pragma once
#include "Integrals/nwx_libint/nwx_libint.hpp"
#include <SDE/Defaults/DefaultModuleTypes.hpp>

/** @file LibIntIntegral.hpp
 *
 *  For the most part the set-up of a tensor containing a lambda function to
 *  libint is the same regardless of which integral we are actually computing.
 *  This header file contains the machinery necessary to encapsulate that logic.
 *
 *  Typical usage would be something like:
 *  ```.cpp
 *  LibIntKinetic<> TBuilder;
 *  auto T = TBuilder.run(molecule, array<BasisSet,2>{bs1, bs2});
 *  ```
 *  That is generally the input is just the molecular system and the basis sets
 *  for the bra and ket (the template parameters, are the derivative order and
 *  the element type within the tensor, we're using the defaults in this
 *  example).
 *
 *  Ultimately all of the integrals are just typedefs of the LibIntIntegral
 *  class.  The inner workings of the LibIntIntegral class rely on two classes:
 *  LibIntFunctor, which is the functor that gets passed to the Tamm tensor and
 *  MakeEngine, which factors out the unique set-up of libint for select
 *  integral types.
 *
 */
namespace Integrals {
/// namespace for implementation details
namespace detail_ {

/**
 * @brief Primary template for building a LibInt engine.
 *
 * Exactly how to set-up the LibInt engine depends on the operator, the
 * number of basis sets, and the requested derivative order.  This class handles
 * the usual set-up, but can be specialized for any particular combination of
 * the above parameters if need be.
 *
 * @tparam op The libint enumeration for the operator
 * @tparam NBases The number of AO indices in the bra plus the number in the ket
 * @tparam Deriv What order derivative should we compute
 */
template<libint2::Operator op, std::size_t, std::size_t Deriv>
struct MakeEngine {
    /**
     * @brief Function that actually makes the engine.
     *
     * @param max_prims The maximum number of primitives possessed by any AO in
     *        the integral.
     * @param max_l The maximum angular momentum possessed by any AO in the
     *        integral.
     * @param thresh The cut-off for declaring an integral to be zero.
     *
     * @return The engine consistent with all of the input.
     *
     * @throws ??? If engine creation throws.  Strong throw guarantee.
     */
    static auto engine(const LibChemist::Molecule&, std::size_t max_prims,
                       std::size_t max_l, double thresh) {
        return libint2::Engine(op, max_prims, max_l, Deriv, thresh);
    }
};

/// MakeEngine specialization to the electron-nucleus attraction
template<std::size_t Deriv>
struct MakeEngine<libint2::Operator::nuclear, 2, Deriv> {
    /// Format LibInt wants the nuclear charges in.
    using charge_array = std::vector<std::pair<double,std::array<double,3>>>;

    ///@copydoc MakeEngine::engine
    static auto engine(const LibChemist::Molecule& molecule,
                       std::size_t max_prims, std::size_t max_l,
                       double thresh) {
        charge_array qs;
        for(const auto& ai: molecule.atoms)
            qs.push_back({ai.properties.at(LibChemist::Atom::Property::charge),
                          {ai.coords[0],ai.coords[1],ai.coords[2]}});
        libint2::Engine engine(libint2::Operator::nuclear, max_prims,
                               max_l, Deriv, thresh);
        engine.set_params(qs);
        return engine;
    }
};

/// MakeEngine specialization to the DF-metric integrals
template<std::size_t Deriv>
struct MakeEngine<libint2::Operator::coulomb, 2, Deriv> {
    ///@copydoc MakeEngine::engine
    static auto engine(const LibChemist::Molecule& molecule,
                       std::size_t max_prims, std::size_t max_l,
                       double thresh) {
        libint2::Engine engine(libint2::Operator::nuclear, max_prims,
                               max_l, Deriv, thresh);
        engine.set_braket(libint2::BraKet::xs_xs);
        return engine;
    }
};

/// MakeEngine specialization to the DF-Coulomb integrals
template<std::size_t Deriv>
struct MakeEngine<libint2::Operator::coulomb, 3, Deriv> {
    ///@copydoc MakeEngine::engine
    static auto engine(const LibChemist::Molecule& molecule,
                       std::size_t max_prims, std::size_t max_l,
                       double thresh) {
        libint2::Engine engine(libint2::Operator::nuclear, max_prims,
                               max_l, Deriv, thresh);
        engine.set_braket(libint2::BraKet::xs_xx);
        return engine;
    }
};

/**
 * @brief This class is the functor we use to fill a tensor
 *
 * @tparam NBases The total number of AO basis sets involved in the bra and the
 *         ket.
 */
template<std::size_t NBases>
struct LibIntFunctor {
    /// The type of a shell block returned by this functor
    using shell_block_type = std::vector<const double*>;

    /// The type of the index to a shell block
    using shell_index = std::array<std::size_t, NBases>;

    /// The object that actually computes the integrals, made by LibIntIntegral
    libint2::Engine engine;

    /// The basis sets, from left to right, in the integral
    std::array<libint2::BasisSet, NBases> bs;

    /**
     * @brief Initializes the functor by initializing libint.
     *
     * @throws ??? if libint2::initialize() throws.  Strong throw guarantee.
     */
    LibIntFunctor(){libint2::initialize();}

    /**
     *  @brief Tells libint we're done with it.
     *
     *  @throws ??? If libint2::finalize throws.  Its a destructor so no
     *  guarantees.
     */
    ~LibIntFunctor(){libint2::finalize();}

    /**
     * @brief The function responsible for actually computing a block of
     * integrals.
     *
     *
     * @param shells The index of the shell block we want.  Index 0 is assumed
     *        to map to bs[0], index 1 to bs[1], etc.
     * @return The requested block of integrals.  For things like dipole
     *         moments multiple blocks are returned corresponding to the
     *         various components.
     *
     * @throw ??? if libint2::Engine throws.  Strong throw guarantee.
     * @throw std::bad_alloc if there is insufficient memory to copy libint's
     *        return over (the pointers, not the actual integrals).
     *
     * @par Data Races:
     * Calls to this function modify the internal state and data races may occur
     * if multiple threads call this function concurrently.
     */
    shell_block_type operator()(shell_index shells)  {
        return call_libint_(shells, std::make_index_sequence<NBases>());
    }

private:
    /**
     * @brief The function that operator() dispatches to.
     *
     * We need to unpack the indices given to operator() into distinct arguments
     * and not an array.  That's what this function does via the usual
     * std::index_sequence trick.
     *
     * @tparam I A variadic parameter pack of integers from [0,NBases) to
     *         expand.
     * @param shells The index of the requested shell block
     * @return An std::vector filled with the requested block per operator
     *         component
     * @throws std::bad_alloc if there is insufficient memory to copy the
     * pointers over.  Strong throw guarantee.
     *
     * @par Data Races:
     * Calls to this function modify the internal state and data races may occur
     * if multiple threads call this function concurrently.
     */
    template<std::size_t...I>
    shell_block_type call_libint_(shell_index shells,
                                  std::index_sequence<I...>) {
        const auto& buf_vec=engine.results();
        engine.compute((bs[I][shells[I]])...);
        std::vector<const double*> rv(buf_vec.begin(),buf_vec.end());
        return rv;
    }

}; // Class LibIntFunctor

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
template<libint2::Operator op, std::size_t NBases,
         std::size_t Deriv=0, typename element_type = double>
struct LibIntIntegral : SDE::Integral<NBases, Deriv, element_type> {
    /// Typedef of base class
    using base_type = SDE::Integral<NBases, Deriv, element_type>;

    /// Pull typdefs from baseclass into scope.
    ///@{
    using tensor_type = typename base_type::tensor_type;
    using molecule_type = typename base_type::molecule_type;
    using basis_array_type = typename base_type::basis_array_type;
    ///@}

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
    tensor_type run(const molecule_type& mol, const basis_array_type& bases)
    override {
        const double thresh = 1.0E-16; // should come from parameters
        std::array<tamm::IndexSpace, NBases> AOs; //AO spaces per mode
        std::array<tamm::TiledIndexSpace, NBases> tAOs; //tiled version of AOs
        detail_::LibIntFunctor<NBases> fxn; // Call-back to give to TAMM
        size_t max_prims = 0; // max primitives in any basis set
        int max_l = 0; // max angular momentum in any basis set

        for(std::size_t basis_i = 0; basis_i < NBases; ++basis_i) {
            const auto& basis = bases[basis_i];
            const auto nshells = basis.types.size();

            // Make index spaces
            AOs[basis_i] = tamm::IndexSpace{tamm::range(0, basis.size())};
            std::vector<unsigned int> tiling(nshells);
            for(std::size_t shell_i = 0; shell_i < nshells; ++shell_i)
                tiling[basis_i] = basis.shellsize(shell_i);
            tAOs[basis_i] = tamm::TiledIndexSpace{AOs[basis_i], tiling};

            // update functor's state
            fxn.bs[basis_i] = nwx_libint::make_basis(basis);
            auto& LIbasis = fxn.bs[basis_i];
            max_prims = std::max(max_prims, LIbasis.max_nprim(LIbasis));
            max_l = std::max(max_l, LIbasis.max_l(LIbasis));
        }

        // Make engine and return (typedef so next line doesn't exceed 80 chars)
        using engine_maker = detail_::MakeEngine<op, NBases, Deriv>;
        fxn.engine = engine_maker::engine(mol, max_prims, max_l, thresh);

        return run_(std::move(AOs), std::move(fxn),
                    std::make_index_sequence<NBases>());
    }

private:
    template<typename AO_type, typename FunctorType, std::size_t...I>
    auto run_(AO_type tAOs, FunctorType fxn, std::index_sequence<I...>) {
        return tensor_type(std::move(tAOs[I])..., std::move(fxn));
    }
};
} // namespace detail_

///Typedef for AO overlap matrix
template<std::size_t Deriv=0, typename element_type=double>
using LibIntOverlap =
  detail_::LibIntIntegral<libint2::Operator::overlap, 2, Deriv, element_type>;

///Typedef for kinetic energy of electrons
template<std::size_t Deriv=0, typename element_type=double>
using LibIntKinetic =
  detail_::LibIntIntegral<libint2::Operator::kinetic, 2, Deriv, element_type>;

///Typedef for the density-fitting metric
template<std::size_t Deriv=0, typename element_type=double>
using LibIntMetric =
    detail_::LibIntIntegral<libint2::Operator::coulomb, 2, Deriv, element_type>;

/// Nucleus-electron attraction
template<std::size_t Deriv=0, typename element_type=double>
using LibIntNucleusElectronAttraction =
    detail_::LibIntIntegral<libint2::Operator::nuclear, 2, Deriv, element_type>;

///Typedef for the density-fitting metric
template<std::size_t Deriv=0, typename element_type=double>
using LibIntDF3C2E =
    detail_::LibIntIntegral<libint2::Operator::coulomb, 3, Deriv, element_type>;

/// Standard 4 electron integral
template<std::size_t Deriv=0, typename element_type=double>
using LibIntERI =
    detail_::LibIntIntegral<libint2::Operator::coulomb, 4, Deriv, element_type>;

} // namespace Integrals
