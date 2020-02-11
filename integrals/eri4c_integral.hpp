#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    class ERI4CInt : public sde::ModuleBase {

        using eri4c_type = property_types::ERI4CIntegral<element_type>;
        using size_type = std::size_t;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        ERI4CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };


    extern template class ERI4CInt<double>;

    using ERI4 = ERI4CInt<double>;
} // namespace integrals