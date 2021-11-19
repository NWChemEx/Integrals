#include "integrals/integralpimpl.hpp"
#include "integrals/libint_functor.hpp"
#include "integrals/tamm_int_functor.hpp"
#include <chemist/basis_set_map.hpp>
#include <sde/module_base.hpp>
#include <tamm/tamm_utils.hpp>
#include <utilities/iter_tools.hpp>

namespace integrals::libint::detail_ {

using size_type     = std::size_t;
using molecule_type = chemist::Molecule;

/*******************************************************************************
 *   PIMPL Implementations.
 *
 *   These are the various ways that one can build a tensor filled with
 *   integrals.
 ******************************************************************************/

/// Builder that assumes you can get the integrals in core memory, with custom
/// tile sizes
template<libint2::Operator op, size_type NBases, typename element_type>
struct CoreIntegrals : IntegralPIMPL<op, NBases, element_type> {
public:
    using base_type        = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO         = typename base_type::tiled_AO;
    using tensor_type      = typename base_type::tensor_type;
    using basis_array_type = typename base_type::basis_array_type;
    using fxn_type         = typename base_type::fxn_type;

private:
    constexpr static size_type nopers = libint2::operator_traits<op>::nopers;

    tensor_type run_impl_(
      const tiled_AO& tAO,
      const std::array<std::vector<size_type>, NBases>& atom_blocks,
      const basis_array_type& bases, fxn_type&& fxn,
      const element_type schwarz_thresh) {
        tamm::ProcGroup pg{GA_MPI_Comm()};
        auto* pMM = tamm::MemoryManagerGA::create_coll(pg);
        tamm::Distribution_NW dist;
        tamm::ExecutionContext ec(pg, &dist, pMM);

        TAMMIntFunctor<op, NBases, element_type> lambda{
          tAO, atom_blocks, bases, std::move(fxn), schwarz_thresh};

        tensor_type A{tAO};
        tamm::Tensor<element_type>::allocate(&ec, A);
        tamm::fill_tensor<element_type>(A, lambda);
        return A;
    }
};

/// Builder that constructs a tensor using a lambda function for direct
/// integrals, with custom tile sizes
template<libint2::Operator op, size_type NBases, typename element_type>
struct DirectIntegrals : IntegralPIMPL<op, NBases, element_type> {
public:
    using base_type        = IntegralPIMPL<op, NBases, element_type>;
    using tiled_AO         = typename base_type::tiled_AO;
    using tensor_type      = typename base_type::tensor_type;
    using basis_array_type = typename base_type::basis_array_type;
    using fxn_type         = typename base_type::fxn_type;

private:
    constexpr static size_type nopers = libint2::operator_traits<op>::nopers;

    tensor_type run_impl_(
      const tiled_AO& tAO,
      const std::array<std::vector<size_type>, NBases>& atom_blocks,
      const basis_array_type& bases, fxn_type&& fxn,
      const element_type schwarz_thresh) {
        TAMMIntFunctor<op, NBases, element_type> lambda{
          tAO, atom_blocks, bases, std::move(fxn), schwarz_thresh};
        tensor_type A{tAO, lambda};
        return A;
    }
};

/*******************************************************************************
 *  Implementation and instantiations of the Integral class.
 ******************************************************************************/

template<libint2::Operator op, size_type NBases, typename element_type>
sde::type::result_map Integral<op, NBases, element_type>::run_(
  sde::type::input_map inputs, sde::type::submodule_map submods) const {
    const auto[mol, bases, deriv] =
      property_types::AOIntegral<NBases, element_type>::unwrap_inputs(inputs);
    const auto thresh = inputs.at("Threshold").value<element_type>();
    const auto schwarz_thresh =
      inputs.at("Screening Threshold").value<element_type>();
    const auto tile_size = inputs.at("Tile Size").value<size_type>();
    const auto op_params = inputs.at("Operator Parameters").value<double>();

    std::array<tamm::IndexSpace, NBases> AOs; // AO spaces per mode
    std::array<std::vector<size_type>, NBases> atom_blocks; // AO spaces per
                                                            // mode
    constexpr static size_type nopers =
      libint2::operator_traits<op>::nopers; // for integrals with multiple
                                            // components
    constexpr size_type extra =
      (nopers > 1) ? 1 : 0; // increase size of tAOs if multiple components
    std::vector<tamm::TiledIndexSpace> tAOs(NBases +
                                            extra); // tiled version of AOs
    size_t max_prims = 0; // max primitives in any basis set
    int max_l        = 0; // max angular momentum in any basis set
    LibintFunctor<NBases> fxn;

    if(nopers > 1) {
        std::vector<unsigned int> tiling(nopers, 1);
        tAOs[0] = tamm::TiledIndexSpace{
          tamm::IndexSpace{tamm::range(0, nopers)}, tiling};
    }

    for(size_type basis_i = 0; basis_i < NBases; ++basis_i) {
        const auto& basis = bases[basis_i];
        chemist::BasisSetMap map(basis);

        AOs[basis_i] = tamm::IndexSpace{tamm::range(0, basis.nbf())};

        std::vector<unsigned int> tiling;
        atom_blocks[basis_i].push_back(0);
        unsigned int count = 0;
        for(size_type ai = 0; ai < mol.size(); ++ai) {
            auto shell_range = map.atom_to_shell(ai);
            for(size_type si = shell_range.first; si < shell_range.second;
                ++si) {
                count += basis[si].size();
            }
            if(count >= tile_size) {
                tiling.push_back(count);
                atom_blocks[basis_i].push_back(ai + 1);
                count = 0;
            }
        }
        if(count != 0) {
            if(count < (tile_size / 5) && !tiling.empty()) {
                tiling.back() += count;
                atom_blocks[basis_i].back() = mol.size();
            } else {
                tiling.push_back(count);
                atom_blocks[basis_i].push_back(mol.size());
            }
        }

        tAOs[basis_i + extra] = tamm::TiledIndexSpace{AOs[basis_i], tiling};

        // update functor's state
        fxn.bs[basis_i] = nwx_libint::make_basis(basis);
        auto& LIbasis   = fxn.bs[basis_i];
        max_prims       = std::max(max_prims, LIbasis.max_nprim(LIbasis));
        max_l           = std::max(max_l, LIbasis.max_l(LIbasis));
    }

    fxn.engine  = make_engine<op, NBases>(mol, max_prims, max_l, thresh, deriv, sde::type::any{op_params});
    auto result = results();
    return property_types::AOIntegral<NBases, element_type>::wrap_results(
      result, pimpl_->run_impl(tAOs, atom_blocks, bases, std::move(fxn),
                               schwarz_thresh));
}

template<libint2::Operator op, size_type NBases, typename element_type>
Integral<op, NBases, element_type>::Integral(implementation_type impl) :
  sde::ModuleBase(this) {
    satisfies_property_type<property_types::AOIntegral<NBases, element_type>>();
    description("Computes integrals of many-body operators with Libint");
    citation("Libint: A library for the evaluation of molecular integrals of "
             "many-body operators over Gaussian functions, Version 2.4.2 "
             "Edward F. Valeev, "
             "http://libint.valeyev.net/");

    add_input<element_type>("Threshold")
      .set_description("Convergence threshold of integrals")
      .set_default(1.0E-16);

    add_input<element_type>("Screening Threshold")
      .set_description("Threshold for Cauchy-Schwarz screening")
      .set_default(0.0);

    add_input<size_type>("Tile Size")
      .set_description("Size threshold for tiling tensors by atom blocks")
      .set_default(size_type{180});

    // really should be sde::type::any but not
    add_input<double>("Operator Parameters")
      .set_description("Operator parameters (only required for some operators)")
      .set_default(std::numeric_limits<double>::quiet_NaN());

  switch(impl) {
        case implementation_type::direct:
            pimpl_ =
              std::make_unique<DirectIntegrals<op, NBases, element_type>>();
            break;
        case implementation_type::core:
            pimpl_ =
              std::make_unique<CoreIntegrals<op, NBases, element_type>>();
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
template class Integral<libint2::Operator::stg, 2, double>;
template class Integral<libint2::Operator::stg, 3, double>;
template class Integral<libint2::Operator::stg, 4, double>;
template class Integral<libint2::Operator::yukawa, 2, double>;
template class Integral<libint2::Operator::yukawa, 3, double>;
template class Integral<libint2::Operator::yukawa, 4, double>;
template class Integral<libint2::Operator::emultipole1, 2, double>;
template class Integral<libint2::Operator::emultipole2, 2, double>;
template class Integral<libint2::Operator::emultipole3, 2, double>;

} // namespace integrals::libint::detail_
