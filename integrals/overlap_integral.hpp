#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/overlap.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    class OverlapInt : public sde::ModuleBase {

        using overlap_type = property_types::OverlapIntegral<element_type>;
        using size_type = std::size_t;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        OverlapInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };


    extern template class OverlapInt<double>;

    using Overlap = OverlapInt<double>;
} // namespace integrals