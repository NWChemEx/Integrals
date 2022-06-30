// #include "integrals/transforms/detail_/run_fundamental.hpp"
// #include <catch2/catch.hpp>
// #include <mokup/mokup.hpp>

// using namespace integrals::detail_;

// using map_type = simde::space_map_t<simde::type::ao_space>;

// TEST_CASE("run_fundamental") {
//     simde::type::el_el_coulomb r12;
//     simde::type::ao_space aos;

//     map_type mode2ao;
//     mode2ao.emplace(0, aos);
//     mode2ao.emplace(1, aos);

//     simde::type::tensor t;

//     auto mod2 =
//       pluginplay::make_lambda<simde::ERI2>([&](auto&& b, auto&& op, auto&& k)
//       {
//           REQUIRE(b == aos);
//           REQUIRE(op == r12);
//           REQUIRE(k == aos);
//           return t;
//       });

//     SECTION("Two bases") {
//         auto rv = run_fundamental<2>(mod2, mode2ao, r12);
//         REQUIRE(t == rv);
//     }

//     SECTION("Three bases") {
//         auto mod3 = pluginplay::make_lambda<simde::ERI3>(
//           [&](auto&& b, auto&& op, auto&& k0, auto&& k1) {
//               REQUIRE(b == aos);
//               REQUIRE(op == r12);
//               REQUIRE(k0 == aos);
//               REQUIRE(k1 == aos);
//               return t;
//           });
//         mode2ao.emplace(2, aos);
//         auto rv = run_fundamental<3>(mod3, mode2ao, r12);
//     }

//     SECTION("Four bases") {
//         auto mod4 = pluginplay::make_lambda<simde::ERI4>(
//           [&](auto&& b0, auto&& b1, auto&& op, auto&& k0, auto&& k1) {
//               REQUIRE(b0 == aos);
//               REQUIRE(b1 == aos);
//               REQUIRE(op == r12);
//               REQUIRE(k0 == aos);
//               REQUIRE(k1 == aos);
//               return t;
//           });
//         mode2ao.emplace(2, aos);
//         mode2ao.emplace(3, aos);
//         auto rv = run_fundamental<4>(mod4, mode2ao, r12);
//     }

//     SECTION("Raises an error if N != mode2ao.size()") {
//         REQUIRE_THROWS_AS(run_fundamental<3>(mod2, mode2ao, r12),
//                           std::runtime_error);
//     }

//     // These shouldn't compile (uncomment to check)
//     // run_fundamental<1>(mod2, mode2ao, r12);
//     // run_fundamental<5>(mod2, mode2ao, r12);
// }
