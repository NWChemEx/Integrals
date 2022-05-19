#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <simde/simde.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

// using namespace mokup;

TEST_CASE("Kinetic Nuclear") {
    // using op_type       = simde::type::el_kinetic_nuc;
    // using integral_type = simde::AO_NuclearTensorRepresentation<2, op_type>;
    using size_vector = std::vector<std::size_t>;

    auto& world = TA::get_default_world();
    // pluginplay::ModuleManager mm;
    // integrals::load_modules(mm);

    auto op_i = integrals::property::kinetic;
    auto op_m = mokup::property::kinetic_nuc_int;
    auto name = mokup::molecule::h2o;
    auto bs   = mokup::basis_set::sto3g;
    auto aos  = mokup::get_bases(name, bs);
    auto aob  = aos.basis_set();
    std::vector bases{bs, bs};
    auto corr = mokup::get_ao_data(name, bases, op_m, world);
    auto mol  = mokup::get_molecule(name);
    auto shl  = aob.shells();
    auto nsh  = aob.n_shells();
    auto nbf  = aob.n_aos();
    auto nat  = mol.size();

    // Make table from shells to atoms

    size_t* shell2atom = new size_t[nsh];
    auto atom2shell    = aob.shell_offsets();
    size_t idx         = 0;
    for(size_t iatm = 0; iatm < nat; iatm++) {
        for(; idx < atom2shell[iatm + 1]; idx++) shell2atom[idx] = iatm;
    }

    // Make a dense tensor (this is just for testing)

    TA::TiledRange1 TR_at{0, nat};
    TA::TiledRange1 TR_crd{0, 3};
    TA::TiledRange1 TR_bf{0, nbf};
    TA::TiledRange trange{TR_at, TR_crd, TR_bf, TR_bf};
    TA::Tensor<float> tile_norm(trange.tiles_range());
    tile_norm[0] = 1.0;
    TA::SparseShape<float> shape(world, tile_norm, trange);
    TA::TSpArrayD ekin(world, trange, shape);
    TA::TSpArrayD::value_type tile(
      ekin.trange().make_tile_range(ekin.begin().index()));
    for(size_t idx = 0; idx < nat * 3 * nbf * nbf; idx++) tile[idx] = 0.0;

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
            size_t idx    = 0;
            for(size_t icrd = 0; icrd < 3; icrd++) {
                size_t ibf_tile = ibf_begin;
                for(size_t ibf = 0; ibf < ibf_size; ibf++) {
                    size_t jbf_tile = jbf_begin;
                    for(size_t jbf = 0; jbf < jbf_size; jbf++) {
                        auto idx_tile =
                          jbf_tile + nbf * (ibf_tile + nbf * (icrd + 3 * iatm));
                        tile[idx_tile] += res[idx];
                        idx++;
                        jbf_tile++;
                    } // for-jbf
                    ibf_tile++;
                } // for-ibf
            } // for-icrd
            for(size_t icrd = 0; icrd < 3; icrd++) {
                size_t ibf_tile = ibf_begin;
                for(size_t ibf = 0; ibf < ibf_size; ibf++) {
                    size_t jbf_tile = jbf_begin;
                    for(size_t jbf = 0; jbf < jbf_size; jbf++) {
                        auto idx_tile =
                          jbf_tile + nbf * (ibf_tile + nbf * (icrd + 3 * jatm));
                        tile[idx_tile] += res[idx];
                        idx++;
                        jbf_tile++;
                    } // for-jbf
                    ibf_tile++;
                } // for-ibf
            } // for-icrd
            jbf_begin += jbf_size;
            jsh++;
        };
        ibf_begin += ibf_size;
        ish++;
    }
    *(ekin.begin()) = tile;
    delete[] shell2atom;
    simde::type::tensor simde_ekin = simde::type::tensor(ekin);
    simde::type::tensor simde_corr = simde::type::tensor(corr);
    REQUIRE(tensorwrapper::tensor::allclose(simde_ekin, simde_corr));
}
