#pragma once
#include <property_types/types.hpp>
#include <sde/types.hpp>
#include <random>

namespace integrals::type {

using namespace property_types::type;

} // namespace integrals::type

// For sdeAny to wrap a more general TA::DistArray (i.e. DirectTile case)
namespace TiledArray {
    template<typename T1, typename T2, typename U1, typename U2>
    bool operator==(const TA::DistArray<T1, T2> &lhs, const TA::DistArray<U1, U2> &rhs) {
        return false;
    }

    template<typename T1, typename T2>
    inline void hash_object(const TA::DistArray<T1, T2> &t, sde::type::hasher &h) {
        std::mt19937 rng;
        rng.seed(std::random_device()());
        std::uniform_real_distribution<double> dist;
        h(dist(rng), dist(rng), dist(rng), dist(rng));
    }
} // namespace TiledArray
