#include "../testing/testing.hpp"
#include <catch2/catch_test_macros.hpp>
#include <integrals/utils/ao_quartet_iterator.hpp>
#include <tuple>

using namespace integrals;
using ao_type             = simde::type::ao_basis_set;
using const_ao_type       = const simde::type::ao_basis_set;
using iterator_type       = utils::AOQuartetIterator<ao_type>;
using const_iterator_type = utils::AOQuartetIterator<const_ao_type>;
using offset_type         = iterator_type::offset_value;

namespace {

template<typename IteratorType>
void test_out_of_range(IteratorType&& aqi) {
    REQUIRE_THROWS_AS(aqi.shell_offsets(), std::out_of_range);
    REQUIRE_THROWS_AS(aqi.shell_quartet(), std::out_of_range);
    REQUIRE_THROWS_AS(aqi.absolute_ao_offsets(), std::out_of_range);
    REQUIRE_THROWS_AS(aqi.relative_ao_offsets(), std::out_of_range);
}

template<typename IteratorType, typename ShellTypes>
void test_quartet(IteratorType& aqi, offset_type corr_shell_offsets,
                  ShellTypes&& corr_shells, offset_type corr_rel_ao_offsets,
                  offset_type corr_abs_ao_offsets) {
    REQUIRE(aqi);

    auto shell_offsets  = aqi.shell_offsets();
    auto shell_quartet  = aqi.shell_quartet();
    auto abs_ao_offsets = aqi.absolute_ao_offsets();
    auto rel_ao_offsets = aqi.relative_ao_offsets();
    for(std::size_t i = 0; i < 4; ++i) {
        REQUIRE(shell_offsets[i] == corr_shell_offsets[i]);
        REQUIRE(rel_ao_offsets[i] == corr_rel_ao_offsets[i]);
        REQUIRE(abs_ao_offsets[i] == corr_abs_ao_offsets[i]);
        REQUIRE(shell_quartet[i] == corr_shells[i]);
    }
}

template<typename IteratorType, typename AOTypes>
void test_h2_sto3g_shell_quartet(IteratorType& aqi, AOTypes&& aos) {
    offset_type ao_off{0, 0, 0, 0};
    using shell_itr = std::decay_t<IteratorType>::shell_quartet_iterator_type;
    shell_itr sqi(aos, aos, aos, aos);

    // To make our life easy, we will assume that all shells have the same
    // number of AOs. Thus the absolute AO offset is just the shell offset
    // multiplied by the number of AOs per shell
    auto n_aos = aos.shell(0).size();
    assert(n_aos == aos.shell(1).size());

    while(sqi) {
        auto corr_shell_offset = sqi.shell_offsets();

        // These are the offsets for the first AOs in the quartet
        offset_type corr_abs_ao_offset(corr_shell_offset);
        for(std::size_t i = 0; i < 4; ++i) { corr_abs_ao_offset[i] *= n_aos; }

        for(ao_off[0] = 0; ao_off[0] < n_aos; ++ao_off[0]) {
            for(ao_off[1] = 0; ao_off[1] < n_aos; ++ao_off[1]) {
                for(ao_off[2] = 0; ao_off[2] < n_aos; ++ao_off[2]) {
                    for(ao_off[3] = 0; ao_off[3] < n_aos; ++ao_off[3]) {
                        // Correct indices are abs offsets + rel offsets
                        offset_type corr_abs_index(corr_abs_ao_offset);
                        for(std::size_t i = 0; i < 4; ++i) {
                            corr_abs_index[i] += ao_off[i];
                        }

                        test_quartet(aqi, corr_shell_offset,
                                     sqi.current_quartet(), ao_off,
                                     corr_abs_index);
                        ++aqi;
                    }
                }
            }
        }

        ++sqi;
    }
    REQUIRE_FALSE(aqi);
    test_out_of_range(aqi);
}
} // namespace

using itr_types = std::tuple<iterator_type, const_iterator_type>;

TEMPLATE_LIST_TEST_CASE("AOQuartetIterator", "", itr_types) {
    SECTION("Constructing with empty basis sets") {
        ao_type bra0, bra1, ket0, ket1;
        TestType aqi(bra0, bra1, ket0, ket1);
        test_out_of_range(aqi);
        REQUIRE_FALSE(aqi);
    }

    SECTION("Iterating over a single shell quartet") {
        auto&& [bra0, bra1, ket0, ket1] = testing::get_h2_dimer_0312_bases();
        offset_type corr_shell_offset{0, 0, 0, 0};
        std::array corr_shells = {bra0.shell(0), bra1.shell(0), ket0.shell(0),
                                  ket1.shell(0)};

        SECTION("(ss|ss)") {
            offset_type corr_rel_ao_offset{0, 0, 0, 0};
            offset_type corr_abs_ao_offset{0, 0, 0, 0};

            TestType aqi(bra0, bra1, ket0, ket1);
            test_quartet(aqi, corr_shell_offset, corr_shells,
                         corr_rel_ao_offset, corr_abs_ao_offset);
            ++aqi;
            REQUIRE_FALSE(aqi);
            test_out_of_range(aqi);
        }
        SECTION("(ps|ss)") {
            offset_type corr_rel_ao_offset{0, 0, 0, 0};
            offset_type corr_abs_ao_offset{0, 0, 0, 0};
            bra0.shell(0).l() = 1;
            corr_shells[0]    = bra0.shell(0);

            TestType aqi(bra0, bra1, ket0, ket1);
            for(std::size_t i = 0; i < 3; ++i) {
                corr_rel_ao_offset[0] = i;
                corr_abs_ao_offset[0] = i;
                test_quartet(aqi, corr_shell_offset, corr_shells,
                             corr_rel_ao_offset, corr_abs_ao_offset);
                ++aqi;
            }
            REQUIRE_FALSE(aqi);
            test_out_of_range(aqi);
        }
        SECTION("(sd|ss)") {
            offset_type corr_rel_ao_offset{0, 0, 0, 0};
            offset_type corr_abs_ao_offset{0, 0, 0, 0};
            bra1.shell(0).l() = 2;
            corr_shells[1]    = bra1.shell(0);

            TestType aqi(bra0, bra1, ket0, ket1);
            for(std::size_t i = 0; i < 6; ++i) {
                corr_rel_ao_offset[1] = i;
                corr_abs_ao_offset[1] = i;
                test_quartet(aqi, corr_shell_offset, corr_shells,
                             corr_rel_ao_offset, corr_abs_ao_offset);
                ++aqi;
            }
            REQUIRE_FALSE(aqi);
            test_out_of_range(aqi);
        }
        SECTION("(sd|ps)") {
            offset_type corr_rel_ao_offset{0, 0, 0, 0};
            offset_type corr_abs_ao_offset{0, 0, 0, 0};
            ket0.shell(0).l() = 1;
            bra1.shell(0).l() = 2;
            corr_shells[1]    = bra1.shell(0);
            corr_shells[2]    = ket0.shell(0);

            TestType aqi(bra0, bra1, ket0, ket1);
            for(std::size_t i = 0; i < 6; ++i) {
                corr_rel_ao_offset[1] = i;
                corr_abs_ao_offset[1] = i;
                for(std::size_t j = 0; j < 3; ++j) {
                    corr_rel_ao_offset[2] = j;
                    corr_abs_ao_offset[2] = j;
                    test_quartet(aqi, corr_shell_offset, corr_shells,
                                 corr_rel_ao_offset, corr_abs_ao_offset);
                    ++aqi;
                }
            }
            REQUIRE_FALSE(aqi);
            test_out_of_range(aqi);
        }
    }

    SECTION("Iterating over multiple quartets") {
        auto aos = testing::h2_sto3g_basis_set();

        // n.b. our test function requires each shell to have the same number of
        // AOs
        SECTION("All s shells") {
            TestType aqi(aos, aos, aos, aos);
            test_h2_sto3g_shell_quartet(aqi, aos);
        }

        SECTION("All p shells") {
            aos.shell(0).l() = 1;
            aos.shell(1).l() = 1;
            TestType aqi(aos, aos, aos, aos);
            test_h2_sto3g_shell_quartet(aqi, aos);
        }
    }
}