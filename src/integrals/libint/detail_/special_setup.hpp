#pragma once

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
        const auto& nuclei = op.get<1>();
        std::vector<std::pair<double, std::array<double, 3>>> qs;
        for(const auto& ai : mol)
            qs.emplace_back(static_cast<const double&>(ai.Z()), ai.coords());
        fill.factory.qs = qs;
    }
};

} // namespace integrals
