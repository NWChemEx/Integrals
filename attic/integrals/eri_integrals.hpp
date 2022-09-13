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
#include <sde/module_base.hpp>

namespace integrals {

template<typename element_type = double>
DECLARE_MODULE(ERI2CInt);

template<typename element_type = double>
DECLARE_MODULE(ERI3CInt);

template<typename element_type = double>
DECLARE_MODULE(ERI4CInt);

extern template class ERI2CInt<double>;
extern template class ERI3CInt<double>;
extern template class ERI4CInt<double>;

using ERI2 = ERI2CInt<double>;
using ERI3 = ERI3CInt<double>;
using ERI4 = ERI4CInt<double>;
} // namespace integrals