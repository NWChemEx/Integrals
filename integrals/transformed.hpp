#pragma once

namespace integrals {

/** @brief Property type returning a
 *
 *  Goal is to get a syntax like
 *
 *  mod.run_as<Transformed<FockMatrix>>(bs1, bs2, occ, ao1, ao2)
 */
template<typename WrappedType>
class Transformed : sde::PropertyType<Transformed<WrappedType>, WrappedType> {};

template<typename TransPropType>
struct GenericTransform : {
    using base_type = TransPropType::WrappedType;
    set_property_type<TransPropType>();
    add_submodule<base_type>();
}

} // namespace integrals