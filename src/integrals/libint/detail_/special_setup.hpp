#pragma once
#include <simde/types.hpp>

namespace integrals {

template<typename Op>
struct SpecialSetup {
    template<typename FillerType>
    static void setup(FillerType&, const Op&) {}
};

template<>
struct SpecialSetup<simde::type::el_nuc_coulomb> {
    template<typename FillerType>
    static auto setup(FillerType& fill, const simde::type::el_nuc_coulomb& op) {
        const auto& nuclei = op.at<1>();
        std::vector<std::pair<double, std::array<double, 3>>> qs;
        for(const auto& ai : nuclei)
            qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());
        fill.factory.qs = qs;
    }
};

template<>
struct SpecialSetup<simde::type::el_el_stg> {
    template<typename FillerType>
    static auto setup(FillerType& fill, const simde::type::el_el_stg& op) {
        const auto& stg           = op.at<0>();
        fill.factory.stg_exponent = stg.exponent;
    }
};

template<>
struct SpecialSetup<simde::type::el_el_yukawa> {
    template<typename FillerType>
    static auto setup(FillerType& fill, const simde::type::el_el_yukawa& op) {
        const auto& stg           = op.at<0>();
        fill.factory.stg_exponent = stg.exponent;
    }
};

} // namespace integrals
