"""Run from the repository root with: python -m unittest discover -s tests -v.

Checks reading the bundled NIST level files (conf,term,J,energy per line),
including lines with commas inside the configuration or the J column.
"""

from pathlib import Path
import os
import re
import unittest

import numpy as np

from colradpy import colradpy
from colradpy.nist_read_txt import _split_level, get_nist_txt


ROOT = Path(__file__).resolve().parents[1]
NIST = str(ROOT / "atomic" / "nist_energies") + "/"


class NistReadTxtTests(unittest.TestCase):
    def test_well_formed_levels(self):
        neutral = get_nist_txt(NIST, "Be", 1)
        self.assertEqual(neutral[0], {"conf": "1s2.2s2", "term": "1S", "j_val": "0",
                                      "energy": 0.0})
        self.assertEqual(neutral[2], {"conf": "1s2.2s.2p", "term": "3P", "j_val": "1",
                                      "energy": 21978.925})
        # Fractional J values, compared with str(Fraction(j)) when shifting energies.
        ion = get_nist_txt(NIST, "Be", 2)
        self.assertEqual((ion[2]["term"], ion[2]["j_val"]), ("2P", "3/2"))
        for level in neutral + ion:
            for key in ("conf", "term", "j_val"):
                self.assertIs(type(level[key]), str)
            self.assertIs(type(level["energy"]), float)

    def test_commas_inside_fields(self):
        # Comma in the configuration.
        self.assertEqual(_split_level("2,3.sp-,3P,0,2428700."),
                         ("2,3.sp-", "3P", "0", 2428700.))
        # Several J values of an unresolved level.
        self.assertEqual(_split_level("2s.2p,3P,0,1,2,470310"),
                         ("2s.2p", "3P", "0,1,2", 470310.))
        self.assertEqual(_split_level("2s.2p3,3D,1/2 or 3/2,262180"),
                         ("2s.2p3", "3D", "1/2 or 3/2", 262180.))
        # Comma in the configuration with empty term and J.
        self.assertEqual(_split_level("5s<1/2,F=4>.5g,,,134512.23"),
                         ("5s<1/2,F=4>.5g", "", "", 134512.23))
        # A numeric term of a four-field line stays the term.
        self.assertEqual(_split_level("4f13.6s2,21,1,500.5"), ("4f13.6s2", "21", "1", 500.5))

    def test_files_with_commas_inside_fields(self):
        levels = get_nist_txt(NIST, "Be", 3)
        self.assertIn({"conf": "2,3.sp+", "term": "1P", "j_val": "1", "energy": 2441800.},
                      levels)
        levels = get_nist_txt(NIST, "He", 0)
        self.assertIn({"conf": "2p.3d", "term": "3D", "j_val": "1,2,3", "energy": 514473.5},
                      levels)
        levels = get_nist_txt(NIST, "In", 2)
        self.assertIn({"conf": "5s<1/2,F=5>.14h", "term": "", "j_val": "", "energy": 149959.83},
                      levels)

    def test_all_bundled_files(self):
        lines = 0
        for name in os.listdir(NIST):
            match = re.fullmatch(r"([a-z]+)_(\d+)", name)
            if match:
                with self.subTest(file=name):
                    levels = get_nist_txt(NIST, *match.groups())
                    self.assertTrue(levels)
                    lines += len(levels)
        self.assertGreater(lines, 100000)

    def test_models_read_bundled_levels(self):
        for name, conf in (("cpb03_ls#be0.dat", "1s2.2s2"), ("be2_adf04", "1s2")):
            with self.subTest(file=name):
                model = colradpy(str(ROOT / "examples" / "input" / name), metas=np.array([0]),
                                 temp_grid=np.array([10.]), electron_den=np.array([1e13]))
                model.get_nist_levels_txt()
                self.assertEqual(model.data["nist"]["levels"][0]["conf"], conf)


if __name__ == "__main__":
    unittest.main()
