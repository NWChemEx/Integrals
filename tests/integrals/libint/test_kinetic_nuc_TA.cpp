#include "../../../src/integrals/libint/detail_/bases_helper.hpp"
#include "../../../src/integrals/libint/detail_/make_shape.hpp"
#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <simde/simde.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

// using namespace mokup;
using namespace integrals::detail_;

TEST_CASE("Kinetic Nuclear") {
    using size_vector = std::vector<std::size_t>;

    auto& world = TA::get_default_world();

    auto op_i = integrals::property::kinetic;
    auto op_m = mokup::property::kinetic_nuc_int;
    auto name = mokup::molecule::h2o;
    auto bs   = mokup::basis_set::sto3g;
    auto aos  = mokup::get_bases(name, bs);
    auto aob  = aos.basis_set();
    auto aol  = make_libint_basis_set(aob); // libint version of the basis set
    std::vector bases{bs, bs};
    std::vector basao{aol, aol};
    auto corr = mokup::get_ao_data(name, bases, op_m);
    auto mol  = mokup::get_molecule(name);
    auto shl  = aob.shells();
    auto nsh  = aob.n_shells();
    auto nbf  = aob.n_aos();
    auto nat  = mol.size();
    auto ncrd = 3ul;

    // Make a dense tensor (this is just for testing)

    std::vector<std::size_t> leading_extents{nat, ncrd};

    // Lambda to calculate gradients
    auto l = [&](const auto& lo, const auto& up, auto* tile) {
        size_t n_element = 1;
        for(size_t idx = 0; idx < lo.size(); idx++) {
            n_element *= up[idx] - lo[idx];
        }
        for(size_t idx = 0; idx < n_element; idx++) tile[idx] = 0.0;

        // Make table from shells to atoms
        std::vector<size_t> shell2atom{};
        auto atom2shell = aob.shell_offsets();
        size_t idx      = 0;
        for(size_t iatm = 0; iatm < nat; iatm++) {
            for(; idx < atom2shell[iatm + 1]; idx++) shell2atom.push_back(iatm);
        }

        auto integral_factory = integrals::Factory(op_i, aob, aob);
        size_t ish            = 0;
        size_t ibf_begin      = 0;
        for(const auto& ishell : shl) {
            auto ibf_size    = ishell.size();
            size_t iatm      = shell2atom[ish];
            size_t jsh       = 0;
            size_t jbf_begin = 0;
            for(const auto& jshell : shl) {
                auto res      = integral_factory.compute(ish, jsh);
                size_t jatm   = shell2atom[jsh];
                auto jbf_size = jshell.size();
                for(size_t icrd = 0; icrd < ncrd; icrd++) {
                    auto resc       = res[icrd];
                    size_t idx      = 0;
                    size_t ibf_tile = ibf_begin;
                    for(size_t ibf = 0; ibf < ibf_size; ibf++) {
                        size_t jbf_tile = jbf_begin;
                        for(size_t jbf = 0; jbf < jbf_size; jbf++) {
                            auto idx_tile =
                              jbf_tile +
                              nbf * (ibf_tile + nbf * (icrd + ncrd * iatm));
                            tile[idx_tile] += resc[idx];
                            idx++;
                            jbf_tile++;
                        } // for-jbf
                        ibf_tile++;
                    } // for-ibf
                }     // for-icrd
                for(size_t icrd = 0; icrd < ncrd; icrd++) {
                    auto resc       = res[ncrd + icrd];
                    size_t idx      = 0;
                    size_t ibf_tile = ibf_begin;
                    for(size_t ibf = 0; ibf < ibf_size; ibf++) {
                        size_t jbf_tile = jbf_begin;
                        for(size_t jbf = 0; jbf < jbf_size; jbf++) {
                            auto idx_tile =
                              jbf_tile +
                              nbf * (ibf_tile + nbf * (icrd + ncrd * jatm));
                            tile[idx_tile] += resc[idx];
                            idx++;
                            jbf_tile++;
                        } // for-jbf
                        ibf_tile++;
                    } // for-ibf
                }     // for-icrd
                jbf_begin += jbf_size;
                jsh++;
            };
            ibf_begin += ibf_size;
            ish++;
        }
    };
    using field_t = typename simde::type::tensor::field_type;
    simde::type::tensor simde_ekin(
      l, integrals::detail_::make_shape(basao, leading_extents),
      tensorwrapper::tensor::default_allocator<field_t>());
    simde::type::tensor simde_corr = simde::type::tensor(corr);
    REQUIRE(tensorwrapper::tensor::allclose(simde_ekin, simde_corr));
}
