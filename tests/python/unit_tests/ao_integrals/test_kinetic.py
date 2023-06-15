#
# Copyright 2023 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pluginplay
import simde
import integrals
import unittest


class TestKinetic(unittest.TestCase):

    def test_integral(self):

        mm = pluginplay.ModuleManager()
        integrals.load_modules(mm)
        mod = mm.at('Kinetic')

        # Uncomment when pt and tensor class is exposed
        # pt = simde.AOTensorRepresentation2T()
        # ints = mod.run_as(pt, aos, t, aos)
        # self.AssertEqual(ints, corr)
