#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/electron_repulsion.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    class ERI2CInt : public sde::ModuleBase {

        using eri2c_type = property_types::ERI2CIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        ERI2CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class ERI3CInt : public sde::ModuleBase {

        using eri3c_type = property_types::ERI3CIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        ERI3CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class ERI4CInt : public sde::ModuleBase {

        using eri4c_type = property_types::ERI4CIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        ERI4CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    extern template class ERI2CInt<double>;
    extern template class ERI3CInt<double>;
    extern template class ERI4CInt<double>;

    using ERI2 = ERI2CInt<double>;
    using ERI3 = ERI3CInt<double>;
    using ERI4 = ERI4CInt<double>;
} // namespace integrals