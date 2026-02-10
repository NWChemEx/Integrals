#include "../testing/testing.hpp"
#include <integrals/utils/shell_quartet_iterator.hpp>

using namespace integrals;
using ao_type             = simde::type::ao_basis_set;
using const_ao_type       = const simde::type::ao_basis_set;
using iterator_type       = utils::ShellQuartetIterator<ao_type>;
using const_iterator_type = utils::ShellQuartetIterator<const_ao_type>;
using offset_type         = iterator_type::offset_value;
using quartet_reference   = iterator_type::quartet_reference;

namespace {
template<typename IteratorType, typename ShellTypes>
void test_quartet(IteratorType& sqi, offset_type corr_offsets,
                  ShellTypes&& corr_shells) {
    REQUIRE(sqi);

    auto offsets = sqi.shell_offsets();
    auto quartet = sqi.current_quartet();
    for(std::size_t i = 0; i < 4; ++i) {
        REQUIRE(offsets[i] == corr_offsets[i]);
        REQUIRE(quartet[i] == corr_shells[i]);
    }
}

template<typename IteratorType, typename AOTypes>
void test_h2_sto3g_shell_quartet(IteratorType& sqi, AOTypes&& aos) {
    for(std::size_t i = 0; i < 2; ++i) {
        for(std::size_t j = 0; j < 2; ++j) {
            for(std::size_t k = 0; k < 2; ++k) {
                for(std::size_t l = 0; l < 2; ++l) {
                    offset_type corr_offset{i, j, k, l};
                    quartet_reference corr_quartet{aos.shell(i), aos.shell(j),
                                                   aos.shell(k), aos.shell(l)};
                    test_quartet(sqi, corr_offset, corr_quartet);
                    ++sqi;
                }
            }
        }
    }
    REQUIRE_FALSE(sqi);

    REQUIRE_THROWS_AS(sqi.shell_offsets(), std::out_of_range);
    REQUIRE_THROWS_AS(sqi.current_quartet(), std::out_of_range);
}
} // namespace

TEST_CASE("ShellQuartetIterator") {
    SECTION("Constructing with empty basis sets") {
        ao_type bra0, bra1, ket0, ket1;
        iterator_type sqi(bra0, bra1, ket0, ket1);
        REQUIRE_THROWS_AS(sqi.shell_offsets(), std::out_of_range);
        REQUIRE_THROWS_AS(sqi.current_quartet(), std::out_of_range);
        REQUIRE_FALSE(sqi);

        const_iterator_type csqi(bra0, bra1, ket0, ket1);
        REQUIRE_THROWS_AS(csqi.shell_offsets(), std::out_of_range);
        REQUIRE_THROWS_AS(csqi.current_quartet(), std::out_of_range);
        REQUIRE_FALSE(csqi);
    }

    SECTION("Iterating over a single quartet") {
        auto&& [bra0, bra1, ket0, ket1] = testing::get_h2_dimer_0312_bases();
        offset_type corr_offset{0, 0, 0, 0};

        SECTION("Mutable shells") {
            quartet_reference corr_quartet{bra0.shell(0), bra1.shell(0),
                                           ket0.shell(0), ket1.shell(0)};

            iterator_type sqi(bra0, bra1, ket0, ket1);
            test_quartet(sqi, corr_offset, corr_quartet);
            ++sqi;
            REQUIRE_FALSE(sqi);
            REQUIRE_THROWS_AS(sqi.shell_offsets(), std::out_of_range);
            REQUIRE_THROWS_AS(sqi.current_quartet(), std::out_of_range);
        }
        SECTION("Immutable shells") {
            std::array corr_quartet{
              std::as_const(bra0).shell(0), std::as_const(bra1).shell(0),
              std::as_const(ket0).shell(0), std::as_const(ket1).shell(0)};

            const_iterator_type sqi(bra0, bra1, ket0, ket1);
            test_quartet(sqi, corr_offset, corr_quartet);
            ++sqi;
            REQUIRE_FALSE(sqi);
            REQUIRE_THROWS_AS(sqi.shell_offsets(), std::out_of_range);
            REQUIRE_THROWS_AS(sqi.current_quartet(), std::out_of_range);
        }
    }

    SECTION("Iterating over multiple quartets") {
        auto aos = testing::h2_sto3g_basis_set();
        iterator_type sqi(aos, aos, aos, aos);
        test_h2_sto3g_shell_quartet(sqi, aos);

        const_iterator_type csqi(aos, aos, aos, aos);
        test_h2_sto3g_shell_quartet(csqi, aos);
    }
}