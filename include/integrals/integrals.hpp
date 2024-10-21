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

/** @file integrals.hpp
 *
 *  This is a convenience header which brings in the entire public API of the
 *  integrals library. It should not be included anywhere in the integrals
 *  source files or public header files (it's okay to use it in the tests if
 *  your unit test also needs most of the headers included by it).
 */
#pragma once
#include "integrals/integrals_mm.hpp"

/** @namespace integrals
 * 
 *  @brief The primary namespace for the Integrals library
 */
namespace integrals {} // end namespace integrals
