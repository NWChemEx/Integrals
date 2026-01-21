#
# Copyright 2026 NWChemEx-Project
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

import unittest

import parallelzone as pz
import pluginplay as pp
from chemist import PointD, ShellType
from chemist.basis_set import AOBasisSetD, AtomicBasisSetD, ContractedGaussianD

import integrals


class TestDecontractBasisSet(unittest.TestCase):
    def test_temp(self):
        origin = PointD(0.0, 0.0, 0.0)
        coefs = [0.1543289673, 0.5353281423, 0.4446345422]
        exps = [3.425250914, 0.6239137298, 0.1688554040]
        cartesian = ShellType.cartesian

        # Inputs
        aobs = AOBasisSetD()
        h = AtomicBasisSetD("STO-3G", 1, origin)
        h_cg = ContractedGaussianD(coefs, exps, origin)
        h.add_shell(cartesian, 0, h_cg)
        aobs.add_center(h)

        # Correct Results
        correct_result = AOBasisSetD()
        h2 = AtomicBasisSetD("STO-3G", 1, origin)
        for exp in exps:
            dcg = ContractedGaussianD([1.0], [exp], origin)
            h2.add_shell(cartesian, 0, dcg)
        correct_result.add_center(h2)

        # Test module
        mod_name = "Decontract Basis Set"
        pt = integrals.DecontractBasisSet()
        result = self.mm.run_as(pt, mod_name, aobs)
        self.assertEqual(result, correct_result)

    def setUp(self):
        self.mm = pp.ModuleManager(pz.runtime.RuntimeView())
        integrals.load_modules(self.mm)
