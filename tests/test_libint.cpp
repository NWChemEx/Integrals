#include "TestCommon.hpp"
#include <integrals/integralsmm.hpp>
#include <property_types/aointegral.hpp>
#include <integrals/nwx_libint/nwx_libint.hpp>
#include <chrono>
#include <libchemist/molecule_stream_parser.hpp>

TEST_CASE("Testing libint itself") {

    std::string file = "/home/jboschen/molecules/h2o_hexamer.xyz";
    std::string basis = "cc-pvdz";
    std::string fit_basis = "cc-pvdz-jkfit";

    std::ifstream fs;
    fs.open(file);

    libchemist::PeriodicTable pt;
    auto mol = libchemist::parse_molecule_stream(fs, libchemist::XYZParser(), pt);
    auto bs = libchemist::apply_basis(basis,mol);
    auto fit_bs = libchemist::apply_basis(fit_basis,mol);

    libint2::initialize();
    auto lint_bs = nwx_libint::make_basis(bs);
    auto lint_fit = nwx_libint::make_basis(fit_bs);

    auto max_nprim = std::max(lint_bs.max_nprim(lint_bs),lint_fit.max_nprim(lint_fit));
    auto max_l = std::max(lint_bs.max_l(lint_bs),lint_fit.max_l(lint_fit));

    libint2::Engine engine(libint2::Operator::coulomb, max_nprim, max_l);
    engine.set(libint2::BraKet::xs_xx);
    const auto& buf_vec = engine.results();

    auto int_st = std::chrono::high_resolution_clock::now();

    for (auto s1=0; s1!=lint_fit.size(); ++s1) {
        for (auto s2=0; s2!=lint_bs.size(); ++s2) {
            for (auto s3=0; s3!=lint_bs.size(); ++s3) {
                engine.compute(lint_fit[s1],lint_bs[s2],lint_bs[s3]);
            }
        }
    }
    auto int_en = std::chrono::high_resolution_clock::now();
    std::cout << "Time ERI3 = " << std::chrono::duration<double>( int_en - int_st ).count() << std::endl;

    libint2::finalize();

    sde::ModuleManager mm;
    integrals::libint::load_modules(mm,700,integrals::libint::detail_::implementation_type::core);
    std::array<libchemist::AOBasisSet, 3> bases = {fit_bs, bs, bs};
    int_st = std::chrono::high_resolution_clock::now();
    auto[Ints] =
    mm.at("ERI3").run_as<property_types::AOIntegral<3,double>>(mol, bases, std::size_t{0});
    int_en = std::chrono::high_resolution_clock::now();
    std::cout << "Time ERI3 Module = " << std::chrono::duration<double>( int_en - int_st ).count() << std::endl;

}