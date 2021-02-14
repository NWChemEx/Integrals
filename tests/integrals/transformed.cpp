#include "../test_common_TA.hpp"
#include "integrals/integrals.hpp"
#include <libchemist/ta_helpers/ta_helpers.hpp>

template<typename T>
using tensor_type = integrals::type::tensor<T>;

namespace {
template<typename ElementType>
auto make_transform(const integrals::type::ao_space_t<ElementType>& bs,
                    TA::World& world, ElementType value) {
    std::initializer_list<ElementType> row{value, value, value, value,
                                           value, value, value};
    tensor_type<ElementType> C(world, {row, row, row, row, row, row, row});
    return libchemist::orbital_space::DerivedSpace<ElementType>(C, bs);
}

} // namespace

TEST_CASE("Two Center Integrals") {
    using element_type = double;
    using tensor_type  = tensor_type<element_type>;
    auto& world        = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [mol, bs] = make_molecule();

    tensor_type corr;

    SECTION("Transform Dipole") {
        using base_pt   = integrals::pt::edipole<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        std::array<double, 3> origin{0, 0, 0};
        const std::string key("EDipole");
        auto [X_ao] = mm.at(key).run_as<base_pt>(origin, bs, bs);

        const auto tkey = std::string("Transformed ") + key;
        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2          = make_transform(bs, world, 2.0);
            corr("q, i, mu") = C2.C()("nu, i") * X_ao("q, nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, bs, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2          = make_transform(bs, world, 2.0);
            corr("q, nu, i") = C2.C()("mu, i") * X_ao("q, nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("q, i, j") =
              C2.C()("nu, i") * X_ao("q, nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("Transform Quadrupole") {
        using base_pt   = integrals::pt::equadrupole<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        std::array<double, 3> origin{0, 0, 0};
        const std::string key("EQuadrupole");
        auto [X_ao] = mm.at(key).run_as<base_pt>(origin, bs, bs);

        const auto tkey = std::string("Transformed ") + key;
        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2          = make_transform(bs, world, 2.0);
            corr("q, i, mu") = C2.C()("nu, i") * X_ao("q, nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, bs, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2          = make_transform(bs, world, 2.0);
            corr("q, nu, i") = C2.C()("mu, i") * X_ao("q, nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("q, i, j") =
              C2.C()("nu, i") * X_ao("q, nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("Transform Octopole") {
        using base_pt   = integrals::pt::eoctopole<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        std::array<double, 3> origin{0, 0, 0};
        const std::string key("EOctopole");
        auto [X_ao] = mm.at(key).run_as<base_pt>(origin, bs, bs);

        const auto tkey = std::string("Transformed ") + key;
        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2          = make_transform(bs, world, 2.0);
            corr("q, i, mu") = C2.C()("nu, i") * X_ao("q, nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, bs, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2          = make_transform(bs, world, 2.0);
            corr("q, nu, i") = C2.C()("mu, i") * X_ao("q, nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("q, i, j") =
              C2.C()("nu, i") * X_ao("q, nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, origin, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("ERI2C") {
        using base_pt   = integrals::pt::eri2c<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("ERI2");
        auto [X_ao] = mm.at(key).run_as<base_pt>(bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("i, mu") = C2.C()("nu, i") * X_ao("nu, mu");
            auto [X]      = mm.at(tkey).run_as<prop_type>(C2, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("nu, i") = C2.C()("mu, i") * X_ao("nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2      = make_transform(bs, world, 2.0);
            auto C3      = make_transform(bs, world, 3.0);
            corr("i, j") = C2.C()("nu, i") * X_ao("nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("Kinetic") {
        using base_pt   = integrals::pt::kinetic<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("Kinetic");
        auto [X_ao] = mm.at(key).run_as<base_pt>(bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("i, mu") = C2.C()("nu, i") * X_ao("nu, mu");
            auto [X]      = mm.at(tkey).run_as<prop_type>(C2, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("nu, i") = C2.C()("mu, i") * X_ao("nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2      = make_transform(bs, world, 2.0);
            auto C3      = make_transform(bs, world, 3.0);
            corr("i, j") = C2.C()("nu, i") * X_ao("nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("Nuclear") {
        using base_pt   = integrals::pt::nuclear<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("Nuclear");
        auto [X_ao] = mm.at(key).run_as<base_pt>(mol, bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, mol, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("i, mu") = C2.C()("nu, i") * X_ao("nu, mu");
            auto [X]      = mm.at(tkey).run_as<prop_type>(C2, bs, mol, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("nu, i") = C2.C()("mu, i") * X_ao("nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, mol, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2      = make_transform(bs, world, 2.0);
            auto C3      = make_transform(bs, world, 3.0);
            corr("i, j") = C2.C()("nu, i") * X_ao("nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, mol, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("Overlap") {
        using base_pt   = integrals::pt::overlap<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("Overlap");
        auto [X_ao] = mm.at(key).run_as<base_pt>(bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("i, mu") = C2.C()("nu, i") * X_ao("nu, mu");
            auto [X]      = mm.at(tkey).run_as<prop_type>(C2, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("nu, i") = C2.C()("mu, i") * X_ao("nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2      = make_transform(bs, world, 2.0);
            auto C3      = make_transform(bs, world, 3.0);
            corr("i, j") = C2.C()("nu, i") * X_ao("nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("STG2") {
        using base_pt   = integrals::pt::stg2c<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("STG2");
        const element_type gamma = 1.2;
        auto [X_ao]              = mm.at(key).run_as<base_pt>(gamma, bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("i, mu") = C2.C()("nu, i") * X_ao("nu, mu");
            auto [X] = mm.at(tkey).run_as<prop_type>(C2, bs, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("nu, i") = C2.C()("mu, i") * X_ao("nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2      = make_transform(bs, world, 2.0);
            auto C3      = make_transform(bs, world, 3.0);
            corr("i, j") = C2.C()("nu, i") * X_ao("nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("Yukawa2") {
        using base_pt   = integrals::pt::yukawa2c<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("Yukawa2");
        const element_type gamma = 1.2;
        auto [X_ao]              = mm.at(key).run_as<base_pt>(gamma, bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("i, mu") = C2.C()("nu, i") * X_ao("nu, mu");
            auto [X] = mm.at(tkey).run_as<prop_type>(C2, bs, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2       = make_transform(bs, world, 2.0);
            corr("nu, i") = C2.C()("mu, i") * X_ao("nu, mu");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2      = make_transform(bs, world, 2.0);
            auto C3      = make_transform(bs, world, 3.0);
            corr("i, j") = C2.C()("nu, i") * X_ao("nu, mu") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, gamma, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }
}

TEST_CASE("Three Center Integrals") {
    using element_type = double;
    using tensor_type  = tensor_type<element_type>;
    auto& world        = TA::get_default_world();
    sde::ModuleManager mm;
    integrals::load_modules(mm);
    auto [mol, bs] = make_molecule();

    tensor_type corr;

    SECTION("ERI3") {
        using base_pt   = integrals::pt::eri3c<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("ERI3");
        auto [X_ao] = mm.at(key).run_as<base_pt>(bs, bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, bs, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2               = make_transform(bs, world, 2.0);
            corr("i, mu, lambda") = C2.C()("nu, i") * X_ao("nu, mu, lambda");
            auto [X] = mm.at(tkey).run_as<prop_type>(C2, bs, bs, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2               = make_transform(bs, world, 2.0);
            corr("nu, i, lambda") = C2.C()("mu, i") * X_ao("nu, mu, lambda");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, bs, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 2") {
            auto C2           = make_transform(bs, world, 2.0);
            corr("nu, mu, i") = C2.C()("lambda, i") * X_ao("nu, mu, lambda");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, bs, C2, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 0 and 1") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("i, j, lambda") =
              C2.C()("nu, i") * X_ao("nu, mu, lambda") * C3.C()("mu,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, bs, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 0 and 2") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("i, mu, j") =
              C2.C()("nu, i") * X_ao("nu, mu, lambda") * C3.C()("lambda,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, bs, C3, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 1 and 2") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("nu, i, j") =
              C2.C()("mu, i") * X_ao("nu, mu, lambda") * C3.C()("lambda,j");

            auto [X] = mm.at(tkey).run_as<prop_type>(bs, C2, C3, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2         = make_transform(bs, world, 2.0);
            auto C3         = make_transform(bs, world, 3.0);
            auto C4         = make_transform(bs, world, 4.0);
            corr("i, j, k") = C2.C()("nu, i") * X_ao("nu, mu, lambda") *
                              C3.C()("mu,j") * C4.C()("lambda,k");

            auto [X] = mm.at(tkey).run_as<prop_type>(C2, C3, C4, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("STG3") {
        using base_pt   = integrals::pt::stg3c<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("STG3");
        const element_type gamma = 1.2;
        auto [X_ao] = mm.at(key).run_as<base_pt>(gamma, bs, bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, bs, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2               = make_transform(bs, world, 2.0);
            corr("i, mu, lambda") = C2.C()("nu, i") * X_ao("nu, mu, lambda");
            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, bs, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2               = make_transform(bs, world, 2.0);
            corr("nu, i, lambda") = C2.C()("mu, i") * X_ao("nu, mu, lambda");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, C2, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 2") {
            auto C2           = make_transform(bs, world, 2.0);
            corr("nu, mu, i") = C2.C()("lambda, i") * X_ao("nu, mu, lambda");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, bs, C2, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 0 and 1") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("i, j, lambda") =
              C2.C()("nu, i") * X_ao("nu, mu, lambda") * C3.C()("mu,j");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, C3, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 0 and 2") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("i, mu, j") =
              C2.C()("nu, i") * X_ao("nu, mu, lambda") * C3.C()("lambda,j");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, bs, C3, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 1 and 2") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("nu, i, j") =
              C2.C()("mu, i") * X_ao("nu, mu, lambda") * C3.C()("lambda,j");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, C2, C3, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2         = make_transform(bs, world, 2.0);
            auto C3         = make_transform(bs, world, 3.0);
            auto C4         = make_transform(bs, world, 4.0);
            corr("i, j, k") = C2.C()("nu, i") * X_ao("nu, mu, lambda") *
                              C3.C()("mu,j") * C4.C()("lambda,k");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, C3, C4, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }

    SECTION("Yukawa3") {
        using base_pt   = integrals::pt::yukawa3c<element_type>;
        using prop_type = integrals::pt::transformed<base_pt>;

        const std::string key("Yukawa3");
        const element_type gamma = 1.2;
        auto [X_ao] = mm.at(key).run_as<base_pt>(gamma, bs, bs, bs);

        const auto tkey = std::string("Transformed ") + key;

        SECTION("No transform") {
            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, bs, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, X_ao));
        }
        SECTION("Mode 0") {
            auto C2               = make_transform(bs, world, 2.0);
            corr("i, mu, lambda") = C2.C()("nu, i") * X_ao("nu, mu, lambda");
            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, bs, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 1") {
            auto C2               = make_transform(bs, world, 2.0);
            corr("nu, i, lambda") = C2.C()("mu, i") * X_ao("nu, mu, lambda");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, C2, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Mode 2") {
            auto C2           = make_transform(bs, world, 2.0);
            corr("nu, mu, i") = C2.C()("lambda, i") * X_ao("nu, mu, lambda");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, bs, C2, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 0 and 1") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("i, j, lambda") =
              C2.C()("nu, i") * X_ao("nu, mu, lambda") * C3.C()("mu,j");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, C3, bs, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 0 and 2") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("i, mu, j") =
              C2.C()("nu, i") * X_ao("nu, mu, lambda") * C3.C()("lambda,j");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, bs, C3, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("Modes 1 and 2") {
            auto C2 = make_transform(bs, world, 2.0);
            auto C3 = make_transform(bs, world, 3.0);
            corr("nu, i, j") =
              C2.C()("mu, i") * X_ao("nu, mu, lambda") * C3.C()("lambda,j");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(bs, C2, C3, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
        SECTION("All modes") {
            auto C2         = make_transform(bs, world, 2.0);
            auto C3         = make_transform(bs, world, 3.0);
            auto C4         = make_transform(bs, world, 4.0);
            corr("i, j, k") = C2.C()("nu, i") * X_ao("nu, mu, lambda") *
                              C3.C()("mu,j") * C4.C()("lambda,k");

            auto [X] =
              mm.at(tkey).run_as<prop_type>(C2, C3, C4, gamma, bs, bs, bs);
            REQUIRE(libchemist::ta_helpers::allclose(X, corr));
        }
    }
}