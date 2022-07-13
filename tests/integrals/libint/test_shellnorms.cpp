#include "integrals/integrals.hpp"
#include <catch2/catch.hpp>
#include <mokup/mokup.hpp>
#include <simde/cauchy_schwarz_approximation.hpp>

using namespace mokup;

std::vector<std::vector<double>> eri_check{
  {2.1874792352627042, 0.3699640325144978, 0.1564525878919850,
   0.0606886000006391, 0.0606886000006391},
  {0.3699640325144978, 0.9039946468459082, 0.4248745604347619,
   0.3222057218175323, 0.3222057218175323},
  {0.1564525878919850, 0.4248745604347619, 0.9381679451862797,
   0.2740491521872362, 0.2740491521872362},
  {0.0606886000006391, 0.3222057218175323, 0.2740491521872362,
   0.8801170058122373, 0.1336470421634229},
  {0.0606886000006391, 0.3222057218175323, 0.2740491521872362,
   0.1336470421634229, 0.8801170058122373},
};

std::vector<std::vector<double>> stg_check{
  {0.8712382481019230, 0.1833574527284384, 0.0499799084050583,
   0.0298548303674855, 0.0298548303674855},
  {0.1833574527284384, 0.5049372580585069, 0.2429078076543889,
   0.1760787877695699, 0.1760787877695699},
  {0.0499799084050583, 0.2429078076543889, 0.5178471413265817,
   0.1563617001398491, 0.1563617001398491},
  {0.0298548303674855, 0.1760787877695699, 0.1563617001398491,
   0.4834953496227799, 0.0689823476694450},
  {0.0298548303674855, 0.1760787877695699, 0.1563617001398491,
   0.0689823476694450, 0.4834953496227799},
};

std::vector<std::vector<double>> yuk_check{
  {1.9779965679341631, 0.3052904918032274, 0.1515386376453260,
   0.0502970255598737, 0.0502970255598737},
  {0.3052904918032274, 0.5375014882330088, 0.3623472474047580,
   0.1781992287690466, 0.1781992287690466},
  {0.1515386376453260, 0.3623472474047580, 0.5862346971968084,
   0.1769470889613180, 0.1769470889613180},
  {0.0502970255598737, 0.1781992287690466, 0.1769470889613180,
   0.5201331598096510, 0.0647323127080033},
  {0.0502970255598737, 0.1781992287690466, 0.1769470889613180,
   0.0647323127080033, 0.5201331598096510},
};

const double eps  = 100000000 * std::numeric_limits<double>::epsilon();
const double marg = 1000000 * std::numeric_limits<double>::epsilon();

TEST_CASE("Cauchy-Schwarz") {
    pluginplay::ModuleManager mm;
    integrals::load_modules(mm);

    auto name = molecule::h2o;
    auto bs   = basis_set::sto3g;
    auto mol  = get_molecule(name);
    auto aos  = get_bases(name, bs);

    SECTION("Coulomb") {
        using op_t  = simde::type::el_el_coulomb;
        using cs_pt = simde::ShellNorms<op_t>;
        op_t op;
        auto [rv] = mm.at("Shell Norms Coulomb").run_as<cs_pt>(aos, op, aos);
        for(int i = 0; i < 5; ++i) {
            for(int j = 0; j < 5; ++j) {
                REQUIRE(rv[i][j] ==
                        Approx(eri_check[i][j]).epsilon(eps).margin(marg));
            }
        }
    }

    SECTION("STG") {
        using op_t  = simde::type::el_el_stg;
        using cs_pt = simde::ShellNorms<op_t>;
        op_t op(chemist::operators::STG(1.0, 1.0));
        auto [rv] = mm.at("Shell Norms STG").run_as<cs_pt>(aos, op, aos);
        for(int i = 0; i < 5; ++i) {
            for(int j = 0; j < 5; ++j) {
                REQUIRE(rv[i][j] ==
                        Approx(stg_check[i][j]).epsilon(eps).margin(marg));
            }
        }
    }

    SECTION("Yukawa") {
        using op_t  = simde::type::el_el_yukawa;
        using cs_pt = simde::ShellNorms<op_t>;
        chemist::Electron e;
        op_t op(chemist::operators::STG(1.0, 1.0), e, e);
        auto [rv] = mm.at("Shell Norms Yukawa").run_as<cs_pt>(aos, op, aos);
        for(int i = 0; i < 5; ++i) {
            for(int j = 0; j < 5; ++j) {
                REQUIRE(rv[i][j] ==
                        Approx(yuk_check[i][j]).epsilon(eps).margin(marg));
            }
        }
    }
}
