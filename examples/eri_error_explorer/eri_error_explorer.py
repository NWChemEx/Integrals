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

"""Explores how screening introduces error into four-center ERIs.

This example looks at how three knobs impact the error in the four-center
electron repulsion integrals (ERIs):

1. The basis set parameters, i.e., which STO-3G shells are used.
2. The separation between the centers the shells sit on.
3. The screening tolerance (the "Threshold" input of the "ERI4" module).

The sampled shells are placed on every vertex of a regular tetrahedron whose
edge length is the center separation, i.e., every pair of centers is the same
distance apart. The ERIs are then computed one shell quartet at a time, only
for symmetry-unique quartets. Every quartet is computed once at
``REFERENCE_TOLERANCE``, which defines the "exact" ERIs, and once per
screening tolerance. The error is the difference between the screened and the
reference ERIs.

Exploring a different region of parameter space is done by changing the
hard-coded parameters below and re-running the module.
"""

import math

import chemcache
import integrals
import numpy as np
import parallelzone as pz
import pluginplay as pp
import simde
from chemist import Electron, PointD
from chemist.basis_set import AOBasisSetD, AtomicBasisSetD, ShellD
from chemist.braket import BraKet, _BraKet_ERI4
from chemist.qm_operator import Coulomb
from chemist.wavefunction import AOs

# -----------------------------------------------------------------------------
# Parameters controlling the run
# -----------------------------------------------------------------------------

#: Shells to sample. Each entry is (Z, i) where i is the index of the shell in
#: ChemCache's STO-3G atomic basis set for the atom with atomic number Z. The
#: default is the H 1s and the O 2p shell.
SHELL_SAMPLES = [(1, 0), (8, 2)]

#: The center separations (in bohr) are ``N_SEPARATIONS`` evenly spaced values
#: from ``R_MIN`` to ``R_MAX`` (inclusive).
R_MIN = 1.0
R_MAX = 5.0
N_SEPARATIONS = 3

#: Screening tolerances to study.
TOLERANCES = [1.0e-6, 1.0e-8]

#: Screening tolerance used to compute the reference ("exact") ERIs.
REFERENCE_TOLERANCE = 1.0e-25

# -----------------------------------------------------------------------------
# Implementation
# -----------------------------------------------------------------------------

_L_LABELS = "spdfghik"

#: Number of centers, i.e., vertices of a tetrahedron.
_N_CENTERS = 4


class SampledShell:
    """The parameters of a shell taken from the STO-3G basis set.

    The center of a shell is not part of the shell's parameters, so that the
    same parameters can be placed on any center.

    :param label: A human-readable name for the shell, e.g., "O[2] p".
    :type label: str
    :param Z: The atomic number of the atom the shell came from.
    :type Z: int
    :param am: The angular momentum of the shell.
    :type am: int
    :param pure: Whether the shell is pure or Cartesian.
    :type pure: chemist.ShellType
    :param coefs: The contraction coefficients of the shell.
    :type coefs: list(float)
    :param exps: The exponents of the primitives in the shell.
    :type exps: list(float)
    """

    def __init__(self, label, Z, am, pure, coefs, exps):
        self.label = label
        self.Z = Z
        self.am = am
        self.pure = pure
        self.coefs = coefs
        self.exps = exps

    def at(self, center, center_label):
        """Places the shell on ``center``.

        :param center: The (x, y, z) coordinates of the center.
        :type center: tuple(float)
        :param center_label: A name for the center, used for printing.
        :type center_label: str

        :return: A single-shell AO space and a label for it.
        :rtype: tuple(chemist.wavefunction.AOs, str)
        """
        x, y, z = center
        shell = ShellD(self.pure, self.am, self.coefs, self.exps, x, y, z)
        atom = AtomicBasisSetD("sto-3g", self.Z, PointD(x, y, z), [shell])
        bs = AOBasisSetD()
        bs.add_center(atom)
        return AOs(bs), f"{self.label}@{center_label}"


def get_sto3g_shell(mm, Z, i):
    """Pulls the parameters of shell ``i`` of atom ``Z`` out of ChemCache.

    :param mm: A ModuleManager with ChemCache loaded.
    :type mm: pluginplay.ModuleManager
    :param Z: The atomic number of the atom.
    :type Z: int
    :param i: The index of the shell in the atom's STO-3G basis set.
    :type i: int

    :return: The parameters of the requested shell.
    :rtype: SampledShell
    """
    atom = mm.run_as(simde.AtomicBasisSetFromZ(), "sto-3g atomic basis", Z)
    shell = atom.at(i)
    prims = [shell.primitive(p) for p in range(shell.n_primitives())]
    coefs = [prim.coefficient for prim in prims]
    exps = [prim.exponent for prim in prims]
    symbol = mm.run_as(simde.SymbolFromZ(), "Symbol from Z", Z)
    label = f"{symbol}[{i}] {_L_LABELS[shell.l]}"
    return SampledShell(label, Z, shell.l, shell.pure, coefs, exps)


def tetrahedron_centers(r):
    """The vertices of a regular tetrahedron with edge length ``r``.

    :param r: The edge length of the tetrahedron.
    :type r: float

    :return: The (x, y, z) coordinates of the four vertices.
    :rtype: list(tuple(float))
    """
    scale = r / (2.0 * math.sqrt(2.0))
    signs = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    return [tuple(scale * s for s in sign) for sign in signs]


def unique_quartets(n):
    """Yields the symmetry-unique shell quartets for ``n`` shells.

    A quartet (i, j, k, l) is symmetry unique if i >= j, k >= l, and the
    pair index of (i, j) is not less than that of (k, l).

    :param n: The number of shells.
    :type n: int

    :return: A generator over the unique quartets.
    :rtype: generator(tuple(int))
    """
    pairs = [(i, j) for i in range(n) for j in range(i + 1)]
    for ij, (i, j) in enumerate(pairs):
        for k, m in pairs[: ij + 1]:
            yield i, j, k, m


def make_eri4_modules(mm):
    """Makes one copy of the "ERI4" module per screening tolerance.

    Modules are locked once they run, so rather than changing the threshold
    of "ERI4" between runs, each tolerance gets its own copy.

    :param mm: A ModuleManager with Integrals loaded.
    :type mm: pluginplay.ModuleManager

    :return: The key of the reference module and a dictionary from each
        tolerance to the key of its module.
    :rtype: tuple(str, dict(float, str))
    """
    ref_key = "ERI4 reference"
    mm.copy_module("ERI4", ref_key)
    mm.change_input(ref_key, "Threshold", REFERENCE_TOLERANCE)

    tol_keys = {}
    for tol in TOLERANCES:
        tol_keys[tol] = f"ERI4 tol={tol:.1e}"
        mm.copy_module("ERI4", tol_keys[tol])
        mm.change_input(tol_keys[tol], "Threshold", tol)
    return ref_key, tol_keys


def compute_eri4(mm, key, aos):
    """Computes the ERIs for one shell quartet.

    :param mm: A ModuleManager with Integrals loaded.
    :type mm: pluginplay.ModuleManager
    :param key: The key of the module to use.
    :type key: str
    :param aos: The single-shell AO spaces of the quartet.
    :type aos: tuple(chemist.wavefunction.AOs)

    :return: The (n_i, n_j, n_k, n_l) block of ERIs.
    :rtype: numpy.ndarray
    """
    v_ee = Coulomb(Electron(), Electron())
    braket = BraKet((aos[0], aos[1]), v_ee, (aos[2], aos[3]))
    # N.B. C++ modules currently fail to convert the erased BraKet to the
    # type ERI4 wants, so we hand them the concrete type.
    t = mm.run_as(simde.ERI4(), key, _BraKet_ERI4(braket))
    return np.array(t)


def print_shells(shells):
    """Prints the parameters of the sampled shells.

    :param shells: The sampled shells.
    :type shells: list(SampledShell)
    """
    print("Sampled STO-3G shells:")
    for shell in shells:
        print(f"  {shell.label} (Z = {shell.Z}, l = {shell.am}, {shell.pure})")
        print("    exponents    :", " ".join(f"{e:.10e}" for e in shell.exps))
        print("    coefficients :", " ".join(f"{c:.10e}" for c in shell.coefs))
    print()


def print_summary(r, tol, records):
    """Summarizes the errors for one separation and tolerance.

    :param r: The center separation.
    :type r: float
    :param tol: The screening tolerance.
    :type tol: float
    :param records: For each quartet, a tuple of the quartet's label, the
        number of distinct centers, and the absolute error block.
    :type records: list(tuple(str, int, numpy.ndarray))
    """
    print(f"Summary for R = {r:.4f}, tol = {tol:.1e}")
    header = (
        f"  {'centers':>7} {'quartets':>8} {'max|err|':>12} {'sum|err|':>12}"
        f"  worst quartet"
    )
    print(header)
    for n_centers in range(1, _N_CENTERS + 1):
        subset = [rec for rec in records if rec[1] == n_centers]
        if not subset:
            continue
        worst = max(subset, key=lambda rec: rec[2].max())
        total = sum(rec[2].sum() for rec in subset)
        print(
            f"  {n_centers:>7} {len(subset):>8} {worst[2].max():12.4e}"
            f" {total:12.4e}  {worst[0]}"
        )
    print()


def main():
    rv = pz.runtime.RuntimeView()

    # N.b. passing None for the cache turns memoization off
    mm = pp.ModuleManager(rv, None)
    mm.load_plugin(chemcache.ChemCache())
    mm.load_plugin(integrals.Integrals())

    ref_key, tol_keys = make_eri4_modules(mm)
    shells = [get_sto3g_shell(mm, Z, i) for Z, i in SHELL_SAMPLES]
    print_shells(shells)

    row = "{:>8} {:>8} {:>3} {:>12} {:>12} {:>12}  {}"
    header = (
        "R",
        "tol",
        "nc",
        "max|ERI|",
        "max|err|",
        "max|err|/tol",
        "quartet",
    )
    for r in np.linspace(R_MIN, R_MAX, N_SEPARATIONS):
        centers = tetrahedron_centers(r)

        # Each element is (AOs, label, center index)
        basis = []
        for c, center in enumerate(centers):
            for shell in shells:
                aos, label = shell.at(center, f"C{c}")
                basis.append((aos, label, c))

        print(f"R = {r:.4f} bohr ({len(basis)} shells)")
        print(row.format(*header))

        records = {tol: [] for tol in TOLERANCES}
        for quartet in unique_quartets(len(basis)):
            aos = tuple(basis[s][0] for s in quartet)
            labels = [basis[s][1] for s in quartet]
            label = f"({labels[0]} {labels[1]}|{labels[2]} {labels[3]})"
            n_centers = len({basis[s][2] for s in quartet})

            ref = compute_eri4(mm, ref_key, aos)
            max_eri = np.abs(ref).max()
            for tol, key in tol_keys.items():
                error = np.abs(compute_eri4(mm, key, aos) - ref)
                records[tol].append((label, n_centers, error))
                print(
                    row.format(
                        f"{r:.4f}",
                        f"{tol:.1e}",
                        n_centers,
                        f"{max_eri:.4e}",
                        f"{error.max():.4e}",
                        f"{error.max() / tol:.4e}",
                        label,
                    )
                )
        print()

        for tol in TOLERANCES:
            print_summary(r, tol, records[tol])


if __name__ == "__main__":
    main()
