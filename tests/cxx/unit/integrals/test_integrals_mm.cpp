/*
 * Copyright 2026 NWChemEx-Project
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

#include "testing/testing.hpp"

TEST_CASE("Integrals plugin") {
    integrals::Integrals plugin;
    const auto& modules = plugin.get_modules();

    SECTION("Identity") {
        REQUIRE(plugin.name() == "integrals");
        // Every default comes from the plugin itself
        REQUIRE(plugin.preferred_providers().empty());
    }

    SECTION("Modules from each subdirectory") {
        REQUIRE(modules.count("AO integral driver"));
        REQUIRE(modules.count("ERI4"));
        REQUIRE(modules.count("Decontract Basis Set"));
    }

    SECTION("Defaults set across subdirectories") {
        REQUIRE(modules.has_submod("AO integral driver", "Kinetic"));
        REQUIRE(modules.has_submod("CauchySchwarz Estimator",
                                   "Decontract Basis Set"));
        REQUIRE(modules.has_submod("Screen Primitive Pairs",
                                   "Primitive Pair Estimator"));
    }

    SECTION("Benchmark ERI4 is a tighter copy of ERI4") {
        const auto& inputs = modules.at("Benchmark ERI4")->inputs();
        REQUIRE(inputs.at("Threshold").value<double>() == 1.0E-16);
        const auto& submod =
          modules.at("CauchySchwarz Estimator")->submods().at("ERI4");
        REQUIRE(submod.value() == *modules.at("Benchmark ERI4"));
    }

    SECTION("Analytic Benchmark ERI4 is a tighter copy of ERI4") {
        const auto& inputs = modules.at("Analytic Benchmark ERI4")->inputs();
        const auto& eri4   = modules.at("ERI4")->inputs();
        const auto thresh  = inputs.at("Threshold").value<double>();
        REQUIRE(thresh < eri4.at("Threshold").value<double>());
    }

    SECTION("Integrals share the UQ Initializers") {
        // Returns the module bound to the "UQ Initializer" submodule of key
        auto uq_initializer = [&](const std::string& key) {
            return &modules.at(key)->submods().at("UQ Initializer").value();
        };
        const auto* uq    = modules.at("UQ Initializer").get();
        const auto* no_uq = modules.at("No UQ Initializer").get();
        for(const auto* key :
            {"Evaluate 2-Index BraKet", "Evaluate 3-Index BraKet",
             "Evaluate 4-Index BraKet", "Overlap", "Kinetic", "Nuclear"})
            REQUIRE(uq_initializer(key) == uq);
        for(const auto* key : {"ERI2", "ERI3", "ERI4", "Benchmark ERI4",
                               "Analytic Benchmark ERI4"})
            REQUIRE(uq_initializer(key) == no_uq);
        REQUIRE(uq_initializer("UQ Atom Symm Blocked Driver") == uq);
    }

    SECTION("Can be loaded into a ModuleManager") {
        pluginplay::ModuleManager mm;
        mm.load_plugin<integrals::Integrals>();
        REQUIRE(mm.has_plugin("integrals"));
        REQUIRE(mm.size() == modules.size());
        REQUIRE(mm.at("AO integral driver").submods().at("Kinetic").value() ==
                mm.at("Kinetic"));
    }
}
