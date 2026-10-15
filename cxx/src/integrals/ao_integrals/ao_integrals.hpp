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
#include <simde/simde.hpp>

namespace integrals::ao_integrals {

DECLARE_MODULE(AOIntegralsDriver);
DECLARE_MODULE(JFourCenter);
DECLARE_MODULE(KFourCenter);
DECLARE_MODULE(JDensityFitted);
DECLARE_MODULE(KDensityFitted);
DECLARE_MODULE(DFIntegral);
DECLARE_MODULE(CoulombMetric);
DECLARE_MODULE(UQDriver);
DECLARE_MODULE(UQAtomSymmBlockedDriver);
DECLARE_MODULE(UQInitializer);

inline void set_defaults(pluginplay::ModuleMap& modules) {
    modules.change_submod("AO integral driver", "Coulomb matrix",
                          "Four center J builder");
    modules.change_submod("AO integral driver", "Exchange matrix",
                          "Four center K builder");
    modules.change_submod("Density Fitted J builder", "DF ERI",
                          "Density Fitting Integral");
    modules.change_submod("Density Fitted K builder", "DF ERI",
                          "Density Fitting Integral");
    modules.change_submod("Density Fitting Integral", "Coulomb Metric",
                          "Coulomb Metric");
    modules.change_submod("UQ Driver", "ERIs", "ERI4");
    modules.change_submod("UQ Driver", "ERI Error", "Primitive Error Model");
    modules.change_submod("UQ Atom Symm Blocked Driver", "ERIs", "ERI4");
    modules.change_submod("UQ Atom Symm Blocked Driver", "ERI Error",
                          "Primitive Error Model");
    modules.change_submod("UQ Atom Symm Blocked Driver", "UQ Initializer",
                          "UQ Initializer");
}

inline void load_modules(pluginplay::ModuleMap& modules) {
    modules.add_module<AOIntegralsDriver>("AO integral driver");
    modules.add_module<JFourCenter>("Four center J builder");
    modules.add_module<KFourCenter>("Four center K builder");
    modules.add_module<JDensityFitted>("Density Fitted J builder");
    modules.add_module<KDensityFitted>("Density Fitted K builder");
    modules.add_module<DFIntegral>("Density Fitting Integral");
    modules.add_module<CoulombMetric>("Coulomb Metric");
    modules.add_module<UQDriver>("UQ Driver");
    modules.add_module<UQAtomSymmBlockedDriver>("UQ Atom Symm Blocked Driver");
    modules.add_module<UQInitializer>("UQ Initializer");
}

} // namespace integrals::ao_integrals
