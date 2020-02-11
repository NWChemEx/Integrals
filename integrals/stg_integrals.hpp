#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/stg.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    class STG2CInt : public sde::ModuleBase {

        using stg2c_type = property_types::STG2CIntegral<element_type>;
        using size_type = std::size_t;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        STG2CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class STG3CInt : public sde::ModuleBase {

        using stg3c_type = property_types::STG3CIntegral<element_type>;
        using size_type = std::size_t;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        STG3CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class STG4CInt : public sde::ModuleBase {

        using stg4c_type = property_types::STG4CIntegral<element_type>;
        using size_type = std::size_t;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        STG4CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    extern template class STG2CInt<double>;
    extern template class STG3CInt<double>;
    extern template class STG4CInt<double>;

    using STG2 = STG2CInt<double>;
    using STG3 = STG3CInt<double>;
    using STG4 = STG4CInt<double>;
} // namespace integrals