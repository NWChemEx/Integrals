#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/yukawa.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    class Yukawa2CInt : public sde::ModuleBase {

        using yukawa2c_type = property_types::Yukawa2CIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        Yukawa2CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class Yukawa3CInt : public sde::ModuleBase {

        using yukawa3c_type = property_types::Yukawa3CIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        Yukawa3CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class Yukawa4CInt : public sde::ModuleBase {

        using yukawa4c_type = property_types::Yukawa4CIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        Yukawa4CInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    extern template class Yukawa2CInt<double>;
    extern template class Yukawa3CInt<double>;
    extern template class Yukawa4CInt<double>;

    using Yukawa2 = Yukawa2CInt<double>;
    using Yukawa3 = Yukawa3CInt<double>;
    using Yukawa4 = Yukawa4CInt<double>;
} // namespace integrals