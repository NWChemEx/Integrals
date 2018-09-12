#include "Integrals/LibintIntegral.hpp"

namespace Integrals::Libint::detail_ {

using size_type = std::size_t;
using molecule_type = LibChemist::Molecule;

static auto get_execution_context() {
    tamm::ProcGroup pg{GA_MPI_Comm()};
    auto *pMM = tamm::MemoryManagerLocal::create_coll(pg);
    tamm::Distribution_NW dist;
    return tamm::ExecutionContext(pg, &dist, pMM);
}

//Functor that wraps the call to libint in a uniform API, used by PIMPLs
template<size_type NBases>
struct LibIntFunctor {
    // The type of a shell block returned by this functor
    using shell_block_type = std::vector<const double*>;

    // The type of the index to a shell block
    using shell_index = std::array<size_type, NBases>;

    // The object that actually computes the integrals, made by LibIntIntegral
    libint2::Engine engine;

    // The basis sets, from left to right, in the integral
    std::array<libint2::BasisSet, NBases> bs;

    //initializes libint
    LibIntFunctor(){libint2::initialize();}
    //finalizes libint
    ~LibIntFunctor(){libint2::finalize();}

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

//Defines the API for the LibintIntegral PIMPL (move to header if needed)
template<libint2::Operator op, size_type NBases,typename element_type>
struct IntegralPIMPL {
    // Typedef of the main class
    using main_type = Integral<op, NBases, element_type>;
    using tensor_type = typename main_type::tensor_type;
    using basis_array_type = typename main_type::basis_array_type;
    using tiled_AO = std::array<tamm::TiledIndexSpace, NBases>;
    using fxn_type = LibIntFunctor<NBases>;

    //Public API to PIMPL
    virtual tensor_type run_impl(const tiled_AO& tAOs,
                                 const basis_array_type& bases,
                                 fxn_type&& fxn) {
        return run_impl_(tAOs, bases, std::move(fxn));
    }

private:
    //Implemented by derived class
    virtual tensor_type run_impl_(const tiled_AO& tAOs,
                                  const basis_array_type& bases,
                                  fxn_type&& fxn) = 0;
};

/*******************************************************************************
 *   PIMPL Implementations.
 *
 *   These are the various ways that one can build a tensor filled with
 *   integrals.
 ******************************************************************************/

///Builder that assumes you can get the integrals in core memory
template<libint2::Operator op, size_type NBases, typename element_type>
struct CoreIntegrals : IntegralPIMPL<op, NBases, element_type> {
public:
    using base_type = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO = typename base_type::tiled_AO;
    using tensor_type = typename base_type::tensor_type;
    using basis_array_type = typename base_type::basis_array_type;
    using fxn_type = typename base_type::fxn_type;
private:
    ///Wraps forwarding the tiled index space into the ctor
    template<size_type...Is>
    tensor_type make_tensor(const tiled_AO& tAO, std::index_sequence<Is...>){
        return tensor_type{tAO[Is]...};
    }

    template<size_type depth>
    void fill(const basis_array_type& bases,
              tensor_type& A,
              std::array<size_type, NBases> idx,
              fxn_type&& fxn)
    {
        if constexpr (depth == NBases) {
            auto buffer = fxn(idx);
            std::size_t nbfs = 1ul;
            for(std::size_t i=0; i<NBases; ++i)
                nbfs *= bases[i][idx[i]].size();

            std::vector<double> buffer2(buffer[0], buffer[0]+nbfs);
            tamm::IndexVector temp(idx.begin(), idx.end());
            A.put(temp, tamm::span<double>(buffer2.data(), buffer2.size()));
        }
        else {
            for(size_type i = 0; i< bases[depth].nshells(); ++i){
                idx[depth] = i;
                //Reset indices after this one
                for(size_type j=depth+1; j < NBases; ++j)idx[j]=0;
                fill<depth+1>(bases, A, idx, std::move(fxn));
            }
        }

    }

    tensor_type run_impl_(const tiled_AO& tAO,
                          const basis_array_type& bases,
                          fxn_type&& fxn) {
        tensor_type A = std::move(make_tensor(tAO,
                std::make_index_sequence<NBases>()));
        tamm::ProcGroup pg{GA_MPI_Comm()};
        auto *pMM = tamm::MemoryManagerLocal::create_coll(pg);
        tamm::Distribution_NW dist;
        tamm::ExecutionContext ec(pg, &dist, pMM);
        tamm::Tensor<double>::allocate(&ec, A);
        fill<0>(bases, A, {}, std::move(fxn));
        return A;
    }
};


/*******************************************************************************
 *  Implementation and instantiations of the Integral class.
 ******************************************************************************/

//Factors out the building of a Libint2 engine.
template<libint2::Operator op, size_type NBases>
static auto make_engine(const molecule_type& molecule, size_type max_prims,
                        size_type max_l, double thresh, size_type deriv) {
    libint2::Engine engine(op, max_prims, max_l, deriv, thresh);
    //Take care of any special set-up
    if constexpr (NBases==2 && op == libint2::Operator::nuclear){
        std::vector<std::pair<double,std::array<double,3>>> qs;
        for(const auto& ai: molecule)
            qs.push_back({static_cast<const double&>(ai.Z()), ai.coords()});
        engine.set_params(qs);
    }
    else if constexpr (NBases == 2 && op ==libint2::Operator::coulomb) {
        engine.set_braket(libint2::BraKet::xs_xs);
    }
    else if constexpr (NBases == 3 && op == libint2::Operator::coulomb) {
        engine.set_braket(libint2::BraKet::xs_xx);
    }
    return engine;
}

template<libint2::Operator op, size_type NBases, typename element_type>
typename Integral<op, NBases, element_type>::tensor_type
Integral<op, NBases, element_type>::run(
        const Integral<op, NBases, element_type>::molecule_type& mol,
        const Integral<op, NBases, element_type>::basis_array_type & bases,
        Integral<op, NBases, element_type>::size_type deriv) {
    const double thresh = 1.0E-16; // should come from parameters
    std::array<tamm::IndexSpace, NBases> AOs; //AO spaces per mode
    std::array<tamm::TiledIndexSpace, NBases> tAOs; //tiled version of AOs
    size_t max_prims = 0; // max primitives in any basis set
    int max_l = 0; // max angular momentum in any basis set
    LibIntFunctor<NBases> fxn;
    for(size_type basis_i = 0; basis_i < NBases; ++basis_i) {
        const auto& basis = bases[basis_i];
        const auto nshells = basis.nshells();

        // Make index spaces
        AOs[basis_i] = tamm::IndexSpace{tamm::range(0, basis.nbf())};
        std::vector<unsigned int> tiling;
        for(const auto& shelli : basis) tiling.push_back(shelli.size());
        tAOs[basis_i] = tamm::TiledIndexSpace{AOs[basis_i], tiling};

        // update functor's state
        fxn.bs[basis_i] = nwx_libint::make_basis(basis);
        auto& LIbasis = fxn.bs[basis_i];
        max_prims = std::max(max_prims, LIbasis.max_nprim(LIbasis));
        max_l = std::max(max_l, LIbasis.max_l(LIbasis));
    }

    fxn.engine = make_engine<op, NBases>(mol, max_prims, max_l, thresh, deriv);
    return pimpl_->run_impl(tAOs, bases, std::move(fxn));
}

template<libint2::Operator op, size_type NBases, typename element_type>
Integral<op, NBases, element_type>::Integral() :
pimpl_(std::make_unique<CoreIntegrals<op, NBases, element_type>>()) {}

template<libint2::Operator op, size_type NBases, typename element_type>
Integral<op, NBases, element_type>::~Integral() noexcept = default;

template class Integral<libint2::Operator::overlap, 2, double>;
template class Integral<libint2::Operator::kinetic, 2, double>;
template class Integral<libint2::Operator::nuclear, 2, double>;
template class Integral<libint2::Operator::coulomb, 2, double>;
template class Integral<libint2::Operator::coulomb, 3, double>;
template class Integral<libint2::Operator::coulomb, 4, double>;

} //namespace Integrals::detail_