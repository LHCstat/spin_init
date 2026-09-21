"""Format conversion checks that do not require VASP or symlinks."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


EXTXYZ_TWO_FRAMES = """2
Lattice="2 0 0 0 3 0 0 0 4" Properties=species:S:1:pos:R:3:spin_length:R:1:initial_magmoms:R:3:spin_forces_vert:R:3:forces:R:3 energy=-1 stress="1 0 0 0 2 0 0 0 3" pbc="T T T"
Ni 0 0 0 2 0 0 1 0.1 0.2 0.3 1 2 3
O 1 2 3 0 0 0 0 0 0 0 -1 -2 -3
2
Lattice="2 0 0 0 3 0 0 0 4" Properties=species:S:1:pos:R:3:spin_length:R:1:initial_magmoms:R:3:spin_forces_vert:R:3:forces:R:3 energy=-2 stress="0 0 0 0 0 0 0 0 0" pbc="T T T"
Ni 0.5 0 0 3 1 0 0 0.4 0.5 0.6 4 5 6
O 1 2 2 0 0 0 0 0 0 0 -4 -5 -6
"""


class TestOut2Npy(unittest.TestCase):
    def test_writes_all_pdf_raw_fields_for_each_frame(self):
        """Catches a missing frame, field, column, or virial sign."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "data.extxyz"
            output = root / "output"
            source.write_text(EXTXYZ_TWO_FRAMES)

            result = subprocess.run(
                [sys.executable, "-m", "dpgen.data.out2npy", str(source), str(output)],
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                (output / "type_map.raw").read_text().splitlines(), ["Ni", "O"]
            )
            np.testing.assert_array_equal(np.loadtxt(output / "type.raw"), [0, 1])
            expected = {
                "box.raw": [[2, 0, 0, 0, 3, 0, 0, 0, 4]] * 2,
                "coord.raw": [[0, 0, 0, 1, 2, 3], [0.5, 0, 0, 1, 2, 2]],
                "energy.raw": [[-1], [-2]],
                "force.raw": [[1, 2, 3, -1, -2, -3], [4, 5, 6, -4, -5, -6]],
                "force_mag.raw": [
                    [0.1, 0.2, 0.3, 0, 0, 0],
                    [0.4, 0.5, 0.6, 0, 0, 0],
                ],
                "spin.raw": [[0, 0, 2, 0, 0, 0], [3, 0, 0, 0, 0, 0]],
                "virial.raw": [
                    [-24, 0, 0, 0, -48, 0, 0, 0, -72],
                    [0] * 9,
                ],
            }
            for name, rows in expected.items():
                with self.subTest(name=name):
                    np.testing.assert_allclose(np.loadtxt(output / name, ndmin=2), rows)
            expected_shapes = {
                "box": (2, 9),
                "coord": (2, 6),
                "energy": (2,),
                "force": (2, 6),
                "force_mag": (2, 6),
                "spin": (2, 6),
                "virial": (2, 9),
            }
            for name, shape in expected_shapes.items():
                with self.subTest(name=name):
                    actual = np.load(output / "set" / (name + ".npy"))
                    self.assertEqual(actual.shape, shape)
                    self.assertEqual(actual.dtype, np.dtype("float64"))
                    np.testing.assert_allclose(
                        actual.reshape(2, -1),
                        np.loadtxt(output / (name + ".raw"), ndmin=2),
                    )
            self.assertTrue((output / "type.raw").is_file())
            self.assertTrue((output / "type_map.raw").is_file())

    def test_existing_output_with_extxyz_is_supported_and_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            output.mkdir()
            source = output / "data.extxyz"
            source.write_text(EXTXYZ_TWO_FRAMES)

            result = subprocess.run(
                [sys.executable, "-m", "dpgen.data.out2npy", str(source), str(output)],
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(source.read_text(), EXTXYZ_TWO_FRAMES)
            self.assertEqual(np.load(output / "set/energy.npy").shape, (2,))

    def test_existing_raw_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            output.mkdir()
            source = output / "data.extxyz"
            source.write_text(EXTXYZ_TWO_FRAMES)
            (output / "energy.raw").write_text("keep me")

            result = subprocess.run(
                [sys.executable, "-m", "dpgen.data.out2npy", str(source), str(output)],
                capture_output=True,
                text=True,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertEqual((output / "energy.raw").read_text(), "keep me")
            self.assertFalse((output / "set").exists())

    def test_single_frame_keeps_energy_one_dimensional(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "data.extxyz"
            output = root / "output"
            source.write_text(EXTXYZ_TWO_FRAMES.split("\n2\n", 1)[0] + "\n")

            result = subprocess.run(
                [sys.executable, "-m", "dpgen.data.out2npy", str(source), str(output)],
                capture_output=True,
                text=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(np.load(output / "set/energy.npy").shape, (1,))
            self.assertEqual(np.load(output / "set/coord.npy").shape, (1, 6))

    def test_incomplete_frame_does_not_publish_partial_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            output.mkdir()
            source = output / "data.extxyz"
            source.write_text(EXTXYZ_TWO_FRAMES.rsplit("\n", 2)[0] + "\n")

            result = subprocess.run(
                [sys.executable, "-m", "dpgen.data.out2npy", str(source), str(output)],
                capture_output=True,
                text=True,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertTrue(source.is_file())
            self.assertFalse((output / "energy.raw").exists())
            self.assertFalse((output / "set").exists())


if __name__ == "__main__":
    unittest.main()
