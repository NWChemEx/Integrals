/*
 * Copyright 2026 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "../testing/testing.hpp"
#include <integrals/property_types.hpp>

using namespace integrals::testing;
using tensorwrapper::operations::approximately_equal;

using pt = integrals::property_types::PrimitivePairEstimator;

namespace {

/* Correct values for H2 come from a MWE I made that recreates the integral
 * screening using libint directly. Since the MWE exactly matches libint, I
 * strongly believe these values are also correct.
 */

template<typename FloatType>
auto corr_k01() {
    std::vector<FloatType> buffer{
      3.693001008106247e-38, 4.768233946869884e-13, 8.942504834854372e-06,
      4.768233946869884e-13, 1.739771026598759e-08, 3.399368340300339e-05,
      8.942504834854372e-06, 3.399368340300339e-05, 1.129153135303198e-04};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

template<typename FloatType>
auto corr_k23() {
    std::vector<FloatType> buffer{
      4.074187000872577e-12, 5.057097771784614e-05, 2.503166969635399e-03,
      5.057097771784614e-05, 9.642627760653152e-04, 3.565848769689949e-03,
      2.503166969635399e-03, 3.565848769689949e-03, 2.170621430835504e-03};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

/* Correct values for the l=2 test computed using libint's renorm() convention,
 * which applies both per-primitive normalization and contracted-shell
 * unit-norm scaling.
 */
template<typename FloatType>
auto corr_k01_l2() {
    std::vector<FloatType> buffer{
      // Primitive 0 (zeta=3.425250914): 6 Cartesian AOs (xx,xy,xz,yy,yz,zz)
      3.257464233125241e-37, 4.205888788821956e-12, 7.887863986537421e-05,
      3.257464233125241e-37, 4.205888788821956e-12, 7.887863986537421e-05,
      3.257464233125241e-37, 4.205888788821956e-12, 7.887863986537421e-05,
      3.257464233125241e-37, 4.205888788821956e-12, 7.887863986537421e-05,
      3.257464233125241e-37, 4.205888788821956e-12, 7.887863986537421e-05,
      3.257464233125241e-37, 4.205888788821956e-12, 7.887863986537421e-05,
      // Primitive 1 (zeta=0.6239137298): 6 Cartesian AOs
      7.661078931859855e-13, 2.795274583136903e-08, 5.461734777212891e-05,
      7.661078931859855e-13, 2.795274583136903e-08, 5.461734777212891e-05,
      7.661078931859855e-13, 2.795274583136903e-08, 5.461734777212891e-05,
      7.661078931859855e-13, 2.795274583136903e-08, 5.461734777212891e-05,
      7.661078931859855e-13, 2.795274583136903e-08, 5.461734777212891e-05,
      7.661078931859855e-13, 2.795274583136903e-08, 5.461734777212891e-05,
      // Primitive 2 (zeta=0.1688554040): 6 Cartesian AOs
      3.888498955506963e-06, 1.478158579140043e-05, 4.909933926030022e-05,
      3.888498955506963e-06, 1.478158579140043e-05, 4.909933926030022e-05,
      3.888498955506963e-06, 1.478158579140043e-05, 4.909933926030022e-05,
      3.888498955506963e-06, 1.478158579140043e-05, 4.909933926030022e-05,
      3.888498955506963e-06, 1.478158579140043e-05, 4.909933926030022e-05,
      3.888498955506963e-06, 1.478158579140043e-05, 4.909933926030022e-05};
    tensorwrapper::shape::Smooth shape({18, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

/* Correct values for water computed using libint's renorm() convention. */
template<typename FloatType>
auto corr_k01_water() {
    std::vector<FloatType> buffer{
      1.807902168137022e+01, 1.748523965357513e+01, 5.449386318678551e+00,
      1.748523965357513e+01, 1.691095962664857e+01, 5.270408290134114e+00,
      5.449386318678551e+00, 5.270408290134114e+00, 1.642556316020210e+00};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

template<typename FloatType>
auto corr_k23_water() {
    std::vector<FloatType> buffer{
      7.980105744079762e-10, 2.571668526586013e-04, 4.110063898356164e-03,
      2.571668526586013e-04, 2.521623116668593e-03, 5.370401777554983e-03,
      4.110063898356164e-03, 5.370401777554983e-03, 2.815605205301997e-03};
    tensorwrapper::shape::Smooth shape({3, 3});
    tensorwrapper::buffer::Contiguous cont(std::move(buffer), shape);
    return simde::type::tensor(std::move(cont), shape);
}

} // namespace

TEST_CASE("BlackBoxPrimitiveEstimator") {
    auto mm          = initialize_integrals();
    auto& mod        = mm.at("Black Box Primitive Pair Estimator");
    using float_type = double;

    SECTION("H2 Dimer/STO-3G (03|12)") {
        auto&& [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();

        SECTION("K(bra0, bra1)") {
            auto K01  = mod.run_as<pt>(bra0, bra1);
            auto corr = corr_k01<float_type>();
            REQUIRE(approximately_equal(K01, corr, 1E-8));
        }

        SECTION("K(ket0, ket1)") {
            auto K23  = mod.run_as<pt>(ket0, ket1);
            auto corr = corr_k23<float_type>();
            REQUIRE(approximately_equal(K23, corr, 1E-8));
        }
    }

    SECTION("H2 Dimer/STO-3G l=2 (03|12)") {
        auto [bra0, bra1, ket0, ket1] = get_h2_dimer_0312_bases();
        bra0.shell(0).l()             = 2;
        bra0.shell(0).pure()          = chemist::ShellType::cartesian;

        auto K    = mod.run_as<pt>(bra0, bra1);
        auto corr = corr_k01_l2<float_type>();
        REQUIRE(approximately_equal(K, corr, 1E-8));
    }

    SECTION("H2O/STO-3G (00|34)") {
        auto&& [bra0, bra1, ket0, ket1] = get_h2o_0034_bases();

        SECTION("K(bra0, bra1)") {
            auto K01  = mod.run_as<pt>(bra0, bra1);
            auto corr = corr_k01_water<float_type>();
            REQUIRE(approximately_equal(K01, corr, 1E-8));
        }

        SECTION("K(ket0, ket1)") {
            auto K23  = mod.run_as<pt>(ket0, ket1);
            auto corr = corr_k23_water<float_type>();
            REQUIRE(approximately_equal(K23, corr, 1E-8));
        }
    }
}
