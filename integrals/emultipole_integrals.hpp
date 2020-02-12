#pragma once
#include <sde/module_base.hpp>
#include <property_types/ao_integrals/emultipole.hpp>
#include "integrals/types.hpp"

namespace integrals {

    template<typename element_type = double>
    class EDipoleInt : public sde::ModuleBase {

        using eDipole_type = property_types::EDipoleIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        EDipoleInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class EQuadrupoleInt : public sde::ModuleBase {

        using eQuadrupole_type = property_types::EQuadrupoleIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        EQuadrupoleInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    template<typename element_type = double>
    class EOctopoleInt : public sde::ModuleBase {

        using eOctopole_type = property_types::EOctopoleIntegral<element_type>;
        using size_type = std::size_t;
        using size_vec = std::vector<size_type>;
        using tensor = typename integrals::type::tensor<element_type>;
        using basis_set = integrals::type::basis_set<element_type>;

    public:
        EOctopoleInt();

    private:
        sde::type::result_map run_(sde::type::input_map inputs,
                                   sde::type::submodule_map submods) const override;
    };

    extern template class EDipoleInt<double>;
    extern template class EQuadrupoleInt<double>;
    extern template class EOctopoleInt<double>;

    using EDipole = EDipoleInt<double>;
    using EQuadrupole = EQuadrupoleInt<double>;
    using EOctopole = EOctopoleInt<double>;
} // namespace integrals