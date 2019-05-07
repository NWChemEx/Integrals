#pragma once
#include <LibChemist/Molecule.hpp>
#include <libint2.hpp>

namespace Integrals::Libint::detail_ {

//Functor that wraps the call to libint in a uniform API, used by PIMPLs
template<std::size_t NBases>
struct LibintFunctor {
    using size_type = std::size_t;
    // The type of a shell block returned by this functor
    using shell_block_type = std::vector<const double*>;

    // The type of the index to a shell block
    using shell_index = std::array<size_type, NBases>;

    // The object that actually computes the integrals, made by LibIntIntegral
    libint2::Engine engine;

    // The basis sets, from left to right, in the integral
    std::array<libint2::BasisSet, NBases> bs;

    //initializes libint
    LibintFunctor(){libint2::initialize();}
    //finalizes libint
    ~LibintFunctor(){libint2::finalize();}

    //Takes a set of shell indices returns the shell block
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
     * @tparam Is A variadic parameter pack of integers from [0,NBases) to
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
    template<size_type...Is>
    shell_block_type call_libint_(shell_index shells,
                                  std::index_sequence<Is...>) {
        const auto& buf_vec=engine.results();
        engine.compute((bs[Is][shells[Is]])...);
        std::vector<const double*> rv(buf_vec.begin(),buf_vec.end());
        return rv;
    }

}; // Class LibIntFunctor

//Factors out the building of a Libint2 engine.
template<libint2::Operator op, std::size_t NBases>
static auto make_engine(const LibChemist::Molecule& molecule,
                        const typename LibintFunctor<NBases>::size_type max_prims,
                        const typename LibintFunctor<NBases>::size_type max_l,
                        const double thresh,
                        const typename LibintFunctor<NBases>::size_type deriv) {
    libint2::Engine engine(op, max_prims, max_l, deriv, thresh);
    //Take care of any special set-up
    if constexpr (NBases==2 && op == libint2::Operator::nuclear){
        std::vector<std::pair<double,std::array<double,3>>> qs;
        for(const auto& ai: molecule)
            qs.push_back({static_cast<const double&>(ai.Z()), ai.coords()});
        engine.set_params(qs);
    }
    else if constexpr (NBases == 2 && op ==libint2::Operator::coulomb) {
        engine.set(libint2::BraKet::xs_xs);
    }
    else if constexpr (NBases == 3 && op == libint2::Operator::coulomb) {
        engine.set(libint2::BraKet::xs_xx);
    }
    return engine;
}

}
