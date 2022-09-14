#pragma once
#include <simde/types.hpp>
#include <tensorwrapper/tensor/allclose.hpp>

/// Roundabout comparison because of lazy eval and allclose
/// TODO: remove when allclose can handle direct tensors
inline bool direct_allclose(simde::type::tensor direct_tensor,
                            simde::type::tensor ref) {
    simde::type::tensor data_tensor, diff;
    auto idx = ref.make_annotation();

    diff(idx)        = direct_tensor(idx) - ref(idx);
    data_tensor(idx) = ref(idx) + diff(idx);

    return tensorwrapper::tensor::allclose(data_tensor, ref);
}