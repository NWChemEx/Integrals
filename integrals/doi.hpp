#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/doi.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    class DOInt : public sde::ModuleBase {

        using doi_type = property_types::DOI<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        DOInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    extern template class DOInt<double>;

    using DOI = DOInt<double>;
} // namespace integrals