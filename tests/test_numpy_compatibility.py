"""Run from the repository root with: python -m unittest discover -s tests -v.

The energy derivative check exposes the remaining NumPy 2.5 complex-output
incompatibility; it is intentionally not marked as an expected failure.
"""

import importlib
from pathlib import Path
import unittest

import numpy as np

from colradpy.energy_balance import energy_balance
from colradpy.read_adf04_py3_class import read_adf04
from colradpy.solve_matrix_exponential import (
    eval_matrix_exponential_solution,
    solve_matrix_exponential,
    solve_matrix_exponential_source,
)


cr_module = importlib.import_module("colradpy.colradpy_class")
INPUT = Path(__file__).resolve().parents[1] / "examples" / "input"
EXAMPLES = ("cpb03_ls#be0.dat", "cpb03_ls#be1.dat", "be2_adf04", "be3_adf04")


class NumpyCompatibilityTests(unittest.TestCase):
    def test_membership_preserves_duplicates_and_inversion(self):
        values = np.array([3, 1, 3, 2])
        np.testing.assert_array_equal(
            cr_module._isin(values, [1, 2]), [False, True, False, True]
        )
        np.testing.assert_array_equal(
            cr_module._isin(values, [1, 2], invert=True),
            [True, False, True, False],
        )

    def test_quadrature_on_nonuniform_grids(self):
        x = np.array([[0.0, 1.0, 3.0], [0.0, 2.0, 3.0]])
        # The integral of 2*x + 1 from 0 to 3 is 12 on either grid.
        np.testing.assert_allclose(
            cr_module._trapezoid(2 * x + 1, x=x, axis=1), [12.0, 12.0]
        )

    def test_bundled_models_on_grid_and_paired_inputs(self):
        for name in EXAMPLES:
            path = INPUT / name
            raw = read_adf04(str(path))
            temperatures = raw["input_file"]["temp_grid"][1:3] / 11604.5
            for paired in (False, True):
                with self.subTest(file=name, paired=paired):
                    model = cr_module.colradpy(
                        str(path),
                        metas=np.array([0]),
                        temp_grid=temperatures,
                        electron_den=np.array([1e12, 1e13]),
                        temp_dens_pair=paired,
                    )
                    # Degenerate levels in bundled data produce a pre-existing
                    # zero-energy wavelength warning, outside this check's scope.
                    with np.errstate(divide="ignore"):
                        model.solve_cr()
                    processed = model.data["processed"]
                    for key in ("pops", "pecs", "scd", "acd", "qcd", "xcd"):
                        self.assertTrue(np.isfinite(processed[key]).all(), key)
                    expected_shape = (2,) if paired else (2, 2)
                    self.assertEqual(processed["pops"].shape[2:], expected_shape)

    def test_nist_configuration_string_operations(self):
        model = cr_module.colradpy.__new__(cr_module.colradpy)
        model.data = {
            "atomic": {"config": np.array(["2s2", "2s1.2p1"])},
            "nist": {"levels": [{"conf": "1s2.2s2"}]},
        }
        model.format_config_to_nist()
        np.testing.assert_array_equal(
            model.data["atomic"]["nist_conf_form"],
            ["1s2.2s2", "1s2.2s1.2p1"],
        )

    def test_two_state_population_evolution(self):
        matrix = np.array([[-1.0, 0.0], [1.0, 0.0]])[:, :, None, None]
        initial = np.array([1.0, 0.0])
        times = np.array([0.0, 1.0, 2.0])
        expected = np.array([np.exp(-times), 1 - np.exp(-times)])
        for with_source in (False, True):
            with self.subTest(with_source=with_source):
                if with_source:
                    pops, _, _ = solve_matrix_exponential_source(
                        matrix, initial, np.zeros(2), times
                    )
                else:
                    pops, _, _ = solve_matrix_exponential(matrix, initial, times)
                np.testing.assert_allclose(pops[:, :, 0, 0], expected, atol=1e-14)
                np.testing.assert_allclose(pops.sum(axis=0), 1.0, atol=1e-14)

    def test_energy_derivative_accepts_solver_populations(self):
        matrix = np.array([[-1.0, 0.0], [1.0, 0.0]])[:, :, None, None]
        initial = np.array([0.5, 0.5])
        _, eigenvalues, eigenvectors = solve_matrix_exponential(
            matrix, initial, np.array([0.0, 1.0])
        )

        def populations(time):
            return eval_matrix_exponential_solution(
                time, initial, eigenvalues, eigenvectors
            )

        # Isolate the energy RHS from file parsing and constructor unit issues.
        state = energy_balance.__new__(energy_balance)
        state.data = {"user": {
            "mass": 12.0, "ion_mass": 2.0, "ion_dens": 1e13,
            "polarizability": 12.0, "ion_charge_number": 1, "ion_temp": 10.0,
        }}
        state.energy_matrix = matrix[:, :, 0, 0]
        state.debye_length = np.array([1e-5])
        derivative = state._evolve_energy(
            0.0, np.array([0.75, 0.75]), populations
        )
        self.assertEqual(derivative.shape, (2,))
        self.assertTrue(np.isfinite(derivative).all())


if __name__ == "__main__":
    unittest.main()
