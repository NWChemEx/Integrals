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
