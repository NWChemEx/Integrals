#include "cs_screened_integrals.hpp"
#include "detail_/fill_ND_functor.hpp"
#include "detail_/nwx_TA_utils.hpp"
#include "detail_/nwx_libint.hpp"
#include "detail_/special_setup.hpp"
#include "detail_/type_traits.hpp"
#include "detail_/unpack_bases.hpp"
#include <simde/cauchy_schwarz_approximation.hpp>

namespace integrals {

using size_type   = std::size_t;
using size_vector = std::vector<size_type>;
using pair_vector = std::vector<std::pair<size_type, size_type>>;

template<typename PropType>
TEMPLATED_MODULE_CTOR(CauchySchwarzScreened, PropType) {
    using element_type   = double; // TODO: Get from PropType
    using cs_approx_type = simde::ShellNorms;

    description("Computes an in-core integral with libint");
    satisfies_property_type<PropType>();

    add_input<double>("Threshold").set_default(1.0E-16);
    add_input<size_vector>("Tile Size").set_default(size_vector{180});
    add_input<double>("Screening Threshold").set_default(0.0);
    add_input<pair_vector>("Atom Tile Groups").set_default(pair_vector{});

    add_submodule<cs_approx_type>("Shell Norms")
      .set_description(
        "Computes the Cauchy-Schwarz Matrix for a pair of basis sets");
}

template<typename PropType>
TEMPLATED_MODULE_RUN(CauchySchwarzScreened, PropType) {
    using element_type   = double; // TODO: Get from PropType
    using tensor_type    = TA::TSpArrayD;
    using value_type     = typename tensor_type::value_type;
    using cs_approx_type = simde::ShellNorms;
    using aospace        = simde::type::ao_space;

    auto& world      = TA::get_default_world(); // TODO: Get from runtime
    auto thresh      = inputs.at("Threshold").value<element_type>();
    auto tile_size   = inputs.at("Tile Size").value<size_vector>();
    auto cs_thresh   = inputs.at("Screening Threshold").value<element_type>();
    auto atom_ranges = inputs.at("Atom Tile Groups").value<pair_vector>();

    constexpr auto n_centers = simde::n_centers_v<PropType>;
    using op_type            = simde::operator_t<PropType>;
    auto op_str              = op_type{}.as_string();
    auto op = inputs.at(op_str).template value<const op_type&>();
    constexpr auto libint_op = op_v<op_type>;

    auto bs     = detail_::unpack_bases<n_centers>(inputs);
    auto trange = nwx_TA::select_tiling(bs, tile_size, atom_ranges);
    auto fill   = nwx_TA::FillNDFunctor<value_type, libint_op, n_centers>();
    fill.initialize(nwx_libint::make_basis_sets(bs), 0, thresh, cs_thresh);

    SpecialSetup<op_type>::setup(fill, op);

    if(cs_thresh > 0.0) {
        if constexpr(n_centers == 4) {
            auto bra1      = inputs.at("bra 1").value<aospace>();
            auto bra2      = inputs.at("bra 1").value<aospace>();
            auto [cs_mat1] = submods.at("Shell Norms")
                               .run_as<cs_approx_type>(bra1, bra2);
            fill.screen.cs_mat1 = cs_mat1;
        }

        auto ket1 = inputs.at("ket 1").value<aospace>();
        auto ket2 = inputs.at("ket 1").value<aospace>();
        auto [cs_mat2] =
          submods.at("Shell Norms").run_as<cs_approx_type>(ket1, ket2);
        fill.screen.cs_mat2 = cs_mat2;
    }

    auto I  = TiledArray::make_array<tensor_type>(world, trange, fill);
    auto rv = results();
    return PropType::wrap_results(rv, simde::type::tensor(I));
}

template class CauchySchwarzScreened<simde::ERI3>;
template class CauchySchwarzScreened<simde::ERI4>;
template class CauchySchwarzScreened<simde::STG3>;
template class CauchySchwarzScreened<simde::STG4>;
template class CauchySchwarzScreened<simde::Yukawa3>;
template class CauchySchwarzScreened<simde::Yukawa4>;

} // namespace integrals
