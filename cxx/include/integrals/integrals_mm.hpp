/*
 * Copyright 2024 NWChemEx-Project
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
#include <pluginplay/plugin/plugin.hpp>

namespace integrals {

/** @brief The Integrals plugin.
 *
 *  Load it into a ModuleManager with `mm.load_plugin<integrals::Integrals>()`
 *  (`mm.load_plugin(integrals.Integrals())` from Python). The plugin's name is
 *  "integrals". All of the plugin's default submodules are provided by the
 *  plugin itself, so it does not need other plugins to be loaded.
 */
DECLARE_PLUGIN(Integrals);

} // end namespace integrals
