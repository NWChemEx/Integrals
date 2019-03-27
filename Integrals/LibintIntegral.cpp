#include "Integrals/LibintIntegral.hpp"
#include <SDE/ModuleBase.hpp>
#include <LibChemist/BasisSetMap.hpp>
#include <Utilities/IterTools.hpp>
#include <tamm/tamm_utils.hpp>

namespace Integrals::Libint::detail_ {

using size_type = std::size_t;
using molecule_type = LibChemist::Molecule;

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
    // For integrals with multiple components
    constexpr static size_type extra = (libint2::operator_traits<op>::nopers > 1) ? 1 : 0;
    // Typedef of the main class
    using main_type = Integral<op, NBases, element_type>;
    using tensor_type = typename main_type::tensor_type;
    using basis_array_type = typename main_type::basis_array_type;
    using tiled_AO = std::vector<tamm::TiledIndexSpace>;
    using fxn_type = LibIntFunctor<NBases>;

    //Public API to PIMPL
    virtual tensor_type run_impl(const tiled_AO& tAOs,
                                 const std::array<std::vector<size_type>, NBases>& atom_blocks,
                                 const basis_array_type& bases,
                                 fxn_type&& fxn) {
        return run_impl_(tAOs, atom_blocks, bases, std::move(fxn));
    }

private:
    //Implemented by derived class
    virtual tensor_type run_impl_(const tiled_AO& tAOs,
                                  const std::array<std::vector<size_type>, NBases>& atom_blocks,
                                  const basis_array_type& bases,
                                  fxn_type&& fxn) = 0;
};

/*******************************************************************************
 *   PIMPL Implementations.
 *
 *   These are the various ways that one can build a tensor filled with
 *   integrals.
 ******************************************************************************/

/// Functor that handles the work of putting the results from the integral engine calculation into the appropriate
/// places in a TAMM tensor.
template<libint2::Operator op, size_type NBases, typename element_type>
struct IntegralLambda {
    using base_type = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO = typename base_type::tiled_AO;
    using basis_array_type = typename base_type::basis_array_type;
    using fxn_type = typename base_type::fxn_type;

    constexpr static size_type nopers = libint2::operator_traits<op>::nopers;

    const tiled_AO tAO;
    const std::array<std::vector<size_type>, NBases> atom_blocks;
    const basis_array_type bases;
    fxn_type fxn;

    IntegralLambda(const tiled_AO& tAO,
                   const std::array<std::vector<size_type>, NBases>& atom_blocks,
                   const basis_array_type& bases,
                   fxn_type&& fxn) : tAO{std::move(tAO)}, atom_blocks{std::move(atom_blocks)},
                                     bases{std::move(bases)}, fxn{std::move(fxn)} {}

    ~IntegralLambda(){}

    void operator()(const tamm::IndexVector& blockid, tamm::span<element_type> buff) {
        std::array<size_type, NBases> idx;
        std::array<LibChemist::BasisSetMap, NBases> maps;
        std::array<size_type, NBases> ao_off;
        const size_type off = (nopers > 1) ? 1 : 0;
        for (size_type i = off; i < NBases + off; ++i) {
            idx[i-off] = blockid[i];
            maps[i-off] = LibChemist::BasisSetMap{bases[i-off]};
            ao_off[i-off] = maps[i-off].atom_to_ao(atom_blocks[i-off][(idx[i-off])]).first;
        }
        const size_type index = (nopers > 1) ? blockid[0] : 0;

        typename fxn_type::shell_index shells;
        using range_array = std::array<typename LibChemist::BasisSetMap::range,NBases>;
        range_array ao_ranges;

        std::function<void(typename fxn_type::shell_index&, range_array&, size_type)> lambda_impl;
        lambda_impl = [&](typename fxn_type::shell_index& shells,
                          range_array& ao_ranges,
                          size_type depth) mutable -> void {

            auto first_shell = maps[depth].atom_to_shell((atom_blocks[depth])[(idx[depth])]).first;
            auto second_shell = maps[depth].atom_to_shell(((atom_blocks[depth])[(idx[depth]) + 1]) - 1).second;

            for (size_type si = first_shell; si < second_shell; ++si) {
                shells[depth] = si;
                ao_ranges[depth] = maps[depth].shell_to_ao(si);
                if (depth == NBases - 1) {
                    std::size_t nbfs = 1ul;
                    for(std::size_t i=0; i<NBases; ++i)
                        nbfs *= bases[i][shells[i]].size();

                    auto buffer = fxn(shells);

                    std::vector<double> libint_buffer;
                    if (buffer[index] == nullptr) {
                        libint_buffer = std::vector<double>(nbfs, 0.0);
                    } else {
                        libint_buffer = std::vector<double>(buffer[index], buffer[index] + nbfs);
                    }

                    size_type x_libint = 0;
                    std::function<void(size_type,size_type)> fill_recursive;
                    fill_recursive = [&](size_type x_tamm, size_type depth) -> void {

                        const auto first_ao = ao_ranges[depth].first;
                        const auto second_ao = ao_ranges[depth].second;
                        const auto tsize = tAO[depth+off].tile_size(idx[depth]);

                        x_tamm = x_tamm*tsize + (first_ao - ao_off[depth]);

                        for (size_type mui = first_ao; mui < second_ao; ++mui) {
                            if (depth == NBases - 1) {
                                buff[x_tamm] = libint_buffer[x_libint++];
                            } else {
                                fill_recursive(x_tamm,depth+1);
                            }
                            x_tamm++;
                        }
                    };
                    fill_recursive(0,0);
                } else {
                    lambda_impl(shells, ao_ranges, depth + 1);
                }
            }
        };
        lambda_impl(shells, ao_ranges, 0);
    }
};

///Builder that assumes you can get the integrals in core memory, with custom tile sizes
template<libint2::Operator op, size_type NBases, typename element_type>
struct CoreIntegrals : IntegralPIMPL<op, NBases, element_type> {
public:
    using base_type = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO = typename base_type::tiled_AO;
    using tensor_type = typename base_type::tensor_type;
    using basis_array_type = typename base_type::basis_array_type;
    using fxn_type = typename base_type::fxn_type;
private:
    constexpr static size_type nopers = libint2::operator_traits<op>::nopers;

    tensor_type run_impl_(const tiled_AO& tAO,
                          const std::array<std::vector<size_type>, NBases>& atom_blocks,
                          const basis_array_type& bases,
                          fxn_type&& fxn) {

        tamm::ProcGroup pg{GA_MPI_Comm()};
        auto *pMM = tamm::MemoryManagerGA::create_coll(pg);
        tamm::Distribution_NW dist;
        tamm::ExecutionContext ec(pg, &dist, pMM);

        IntegralLambda<op,NBases,element_type> lambda{tAO,atom_blocks,bases,std::move(fxn)};

        tensor_type A{tAO};
        tamm::Tensor<element_type>::allocate(&ec, A);
        tamm::fill_tensor<element_type>(ec,A,lambda);
        return A;
    }
};

///Builder that constructs a tensor using a lambda function for direct integrals, with custom tile sizes
template<libint2::Operator op, size_type NBases, typename element_type>
struct DirectIntegrals : IntegralPIMPL<op, NBases, element_type> {
public:
    using base_type = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO = typename base_type::tiled_AO;
    using tensor_type = typename base_type::tensor_type;
    using basis_array_type = typename base_type::basis_array_type;
    using fxn_type = typename base_type::fxn_type;
private:
    constexpr static size_type nopers = libint2::operator_traits<op>::nopers;

    tensor_type run_impl_(const tiled_AO& tAO,
                          const std::array<std::vector<size_type>, NBases>& atom_blocks,
                          const basis_array_type& bases,
                          fxn_type&& fxn) {

        IntegralLambda<op,NBases,element_type> lambda{tAO,atom_blocks,bases,std::move(fxn)};
        tensor_type A{tAO,lambda};
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
        engine.set(libint2::BraKet::xs_xs);
    }
    else if constexpr (NBases == 3 && op == libint2::Operator::coulomb) {
        engine.set(libint2::BraKet::xs_xx);
    }
    return engine;
}

template<libint2::Operator op, size_type NBases, typename element_type>
SDE::type::result_map Integral<op, NBases, element_type>::run_(SDE::type::input_map inputs,
                                                           SDE::type::submodule_map submods) const {

    const auto [mol, bases, deriv] = LibChemist::AOIntegral<NBases, element_type>::unwrap_inputs(inputs);
    const auto thresh = inputs.at("Threshold").value<double>();
    const auto tile_size = inputs.at("Tile Size").value<size_type>();

    std::array<tamm::IndexSpace, NBases> AOs; //AO spaces per mode
    std::array<std::vector<size_type>, NBases> atom_blocks; //AO spaces per mode
    constexpr static size_type nopers = libint2::operator_traits<op>::nopers; // for integrals with multiple components
    constexpr size_type extra = (nopers > 1) ? 1 : 0; // increase size of tAOs if multiple components
    std::vector<tamm::TiledIndexSpace> tAOs(NBases+extra); //tiled version of AOs
    size_t max_prims = 0; // max primitives in any basis set
    int max_l = 0; // max angular momentum in any basis set
    LibIntFunctor<NBases> fxn;


    if (nopers > 1) {
        std::vector<unsigned int> tiling(nopers,1);
        tAOs[0] = tamm::TiledIndexSpace{tamm::IndexSpace{tamm::range(0, nopers)}, tiling};
    }

    for(size_type basis_i = 0; basis_i < NBases; ++basis_i) {
        const auto& basis = bases[basis_i];
        LibChemist::BasisSetMap map(basis);

        AOs[basis_i] = tamm::IndexSpace{tamm::range(0, basis.nbf())};

        std::vector<unsigned int> tiling;
        atom_blocks[basis_i].push_back(0);
        unsigned int count = 0;
        for (size_type ai = 0; ai < mol.size(); ++ai) {
            auto shell_range = map.atom_to_shell(ai);
            for (size_type si = shell_range.first; si < shell_range.second; ++si) {
                count += basis[si].size();
            }
            if (count >= tile_size) {
                tiling.push_back(count);
                atom_blocks[basis_i].push_back(ai+1);
                count = 0;
            }
        }
        if (count != 0) {
            if (count < (tile_size / 5) && !tiling.empty()) {
                tiling.back() += count;
                atom_blocks[basis_i].back() = mol.size();
            } else {
                tiling.push_back(count);
                atom_blocks[basis_i].push_back(mol.size());
            }
        }

        tAOs[basis_i+extra] = tamm::TiledIndexSpace{AOs[basis_i], tiling};

        // update functor's state
        fxn.bs[basis_i] = nwx_libint::make_basis(basis);
        auto& LIbasis = fxn.bs[basis_i];
        max_prims = std::max(max_prims, LIbasis.max_nprim(LIbasis));
        max_l = std::max(max_l, LIbasis.max_l(LIbasis));

    }

    fxn.engine = make_engine<op, NBases>(mol, max_prims, max_l, thresh, deriv);
    auto result = results();
    return LibChemist::AOIntegral<NBases, element_type>::wrap_results(result,
                                                                  pimpl_->run_impl(tAOs, atom_blocks, bases, std::move(fxn)));
}

template<libint2::Operator op, size_type NBases, typename element_type>
Integral<op, NBases, element_type>::Integral(implementation_type impl) : SDE::ModuleBase(this) {

    satisfies_property_type<LibChemist::AOIntegral<NBases, element_type>>();
    description("Computes integrals of many-body operators with Libint");
    citation("Libint: A library for the evaluation of molecular integrals of "
      "many-body operators over Gaussian functions, Version 2.4.2 Edward F. Valeev, "
      "http://libint.valeyev.net/");

    add_input<double>("Threshold")
        .set_description("Convergence threshold of integrals")
        .set_default(1.0E-16);

    add_input<size_type>("Tile Size")
        .set_description("Size threshold for tiling tensors by atom blocks")
        .set_default(size_type{180});

    switch (impl) {
        case implementation_type::direct:
            pimpl_ = std::make_unique<DirectIntegrals<op, NBases, element_type>>();
            break;
        case implementation_type::core:
            pimpl_ = std::make_unique<CoreIntegrals<op, NBases, element_type>>();
            break;
    }
}

template<libint2::Operator op, size_type NBases, typename element_type>
Integral<op, NBases, element_type>::~Integral() noexcept = default;

template class Integral<libint2::Operator::overlap, 2, double>;
template class Integral<libint2::Operator::kinetic, 2, double>;
template class Integral<libint2::Operator::nuclear, 2, double>;
template class Integral<libint2::Operator::coulomb, 2, double>;
template class Integral<libint2::Operator::coulomb, 3, double>;
template class Integral<libint2::Operator::coulomb, 4, double>;
template class Integral<libint2::Operator::emultipole1, 2, double>;
template class Integral<libint2::Operator::emultipole2, 2, double>;
template class Integral<libint2::Operator::emultipole3, 2, double>;

} //namespace Integrals::detail_
