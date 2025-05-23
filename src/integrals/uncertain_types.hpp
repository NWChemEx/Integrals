/*
 * Copyright 2025 NWChemEx-Project
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
#include <stdexcept>
#ifdef ENABLE_SIGMA
#include <sigma/sigma.hpp>
#endif

namespace integrals::type {
#ifdef ENABLE_SIGMA

static constexpr bool has_sigma() { return true; }

using uncertain_float  = sigma::UFloat;
using uncertain_double = sigma::UDouble;
#else

static constexpr bool has_sigma() { return false; }

using uncertain_float  = float;
using uncertain_double = double;
#endif

} // namespace integrals::type