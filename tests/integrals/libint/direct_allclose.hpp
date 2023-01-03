/*
 * Copyright 2023 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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