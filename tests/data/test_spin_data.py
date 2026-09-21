"""Windows-safe Stage 5 checks; Linux/VASP validation remains separate."""

import argparse
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import dpdata
import numpy as np

from dpgen.data import spin_data, spin_init
from tests.data.test_spin_init import POSCAR_TEXT


def _magnetization_outcar(moment):
    def block(axis, value):
        return (
            f" magnetization ({axis})\n"
            " # of ion       s       p       d       tot\n"
            " ------------------------------------------\n"
            f" 1  0.0  0.0  0.0  {value:.8f}\n"
            " 2  0.0  0.0  0.0  0.00000000\n"
            " tot 0.0 0.0 0.0 0.0\n"
        )

    return (
        " vasp.6.4.2  test\n"
        + block("x", 0.0)
        + block("y", 0.0)
        + block("z", 0.5)
        + block("x", moment[0])
        + block("y", moment[1])
        + block("z", moment[2])
    )


OSZICAR_TEXT = """ 1 MW_int
 1 0.0 0.0 1.0
 2 0.0 0.0 0.0
 1 lambda*MW_perp
 1 0.0 0.0 0.1
 2 0.0 0.0 0.0
 2 MW_int
 1 0.0 0.0 2.0
 2 0.0 0.0 0.0
 2 lambda*MW_perp
 1 0.1 0.2 0.3
"""

INCAR_TEXT = """LNONCOLLINEAR = .TRUE.
NSW = 0
MAGMOM = 0 0 1 0 0 0
M_CONSTR = 0 0 1 0 0 0
"""


def _labeled_frame():
    return dpdata.LabeledSystem(
        data={
            "atom_names": ["Fe"],
            "atom_numbs": [2],
            "atom_types": np.array([0, 0]),
            "orig": np.zeros(3),
            "cells": np.array([np.eye(3) * 2.0]),
            "coords": np.array([[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]]),
            "energies": np.array([-2.0]),
            "forces": np.array([[[0.1, 0.2, 0.3], [0.0, 0.0, 0.0]]]),
        }
    )


class TestSpinMagneticParser(unittest.TestCase):
    def test_reads_last_outcar_xyz_tables(self):
        with tempfile.TemporaryDirectory() as temporary:
            outcar = Path(temporary) / "OUTCAR"
            outcar.write_text(_magnetization_outcar((0.0, 0.0, 1.003)))
            major, moments = spin_data.read_final_magnetization(outcar, 2)
            outcar.write_text(
                _magnetization_outcar((0.0, 0.0, 1.003)).replace(
                    "vasp.6.4.2", "vasp.5.4.4"
                )
            )
            legacy_major, legacy_moments = spin_data.read_final_magnetization(
                outcar, 2
            )

        self.assertEqual(major, 6)
        self.assertEqual(legacy_major, 5)
        np.testing.assert_allclose(moments, [[0.0, 0.0, 1.003], [0.0, 0.0, 0.0]])
        np.testing.assert_allclose(legacy_moments, moments)

    def test_vasp5_doubles_perpendicular_force_and_vasp6_does_not(self):
        with tempfile.TemporaryDirectory() as temporary:
            oszicar = Path(temporary) / "OSZICAR"
            oszicar.write_text(OSZICAR_TEXT)
            v5 = spin_data.read_final_spin_force(oszicar, 2, 5)
            v6 = spin_data.read_final_spin_force(oszicar, 2, 6)

        np.testing.assert_allclose(v5, [[0.4, 0.8, 1.2], [0, 0, 0]])
        np.testing.assert_allclose(v6, [[0.2, 0.4, 0.6], [0, 0, 0]])

    def test_rmse_ignores_initial_zero_but_counts_final_zero(self):
        initial = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])
        final = np.array([[0.0, 0.0, 1.003], [100.0, 0.0, 0.0]])
        self.assertAlmostEqual(spin_data.magnetic_rmse(initial, final), 0.003)
        final[0] = 0.0
        self.assertAlmostEqual(spin_data.magnetic_rmse(initial, final), 1.0)

    def test_missing_final_axis_and_empty_oszicar_vectors_are_errors(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outcar = root / "OUTCAR"
            outcar.write_text(
                _magnetization_outcar((0.0, 0.0, 1.0)).replace(
                    "magnetization (z)", "magnetization (q)"
                )
            )
            with self.assertRaisesRegex(ValueError, "x/y/z"):
                spin_data.read_final_magnetization(outcar, 2)
            oszicar = root / "OSZICAR"
            oszicar.write_text("1 MW_int\n2 lambda*MW_perp\n")
            with self.assertRaisesRegex(ValueError, "no atom rows"):
                spin_data.read_final_spin_force(oszicar, 2, 6)


class TestSpinDataCollector(unittest.TestCase):
    def test_real_outcar_fixture_has_aligned_dpdata_final_frame(self):
        outcar = (
            Path(__file__).parent
            / "out_data_02_md/02.md/sys-004/scale-1.000/000000/OUTCAR"
        )
        frame, rotation, index = spin_data._last_dpdata_frame(outcar)
        self.assertGreaterEqual(index, 0)
        self.assertGreater(frame.get_natoms(), 0)
        np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-6)

    def test_filters_high_rmse_and_writes_aligned_raw_npy_and_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tasks = [
                "scale-1.000/000000/00/origin/000000",
                "scale-1.000/000000/00/Rotation-000/R1",
                "scale-1.000/000000/00/Canting-001/C1",
            ]
            stage = root / "03.spin"
            stage.mkdir()
            (stage / "tasks.json").write_text(json.dumps({"tasks": tasks}))
            for task, magnitude in zip(tasks, (1.003, 1.1, 0.997)):
                folder = stage / task
                folder.mkdir(parents=True)
                (folder / "POSCAR").write_text(POSCAR_TEXT)
                (folder / "INCAR").write_text(INCAR_TEXT)
                (folder / "OUTCAR").write_text(
                    _magnetization_outcar((magnitude, 0.0, 0.0))
                )
                (folder / "OSZICAR").write_text(OSZICAR_TEXT)

            rotation = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
            with mock.patch.object(spin_data, "check_spin_results", return_value=3):
                with mock.patch.object(
                    spin_data,
                    "_last_dpdata_frame",
                    side_effect=lambda _: (_labeled_frame(), rotation, 4),
                ):
                    count = spin_data.collect_spin_data({"out_dir": str(root)})

            self.assertEqual(count, 2)
            deepmd = root / "04.data/deepmd"
            selected = json.loads((deepmd / "selection.json").read_text())
            self.assertEqual(selected["selected"][0]["task"], tasks[0])
            self.assertEqual(selected["selected"][1]["task"], tasks[2])
            self.assertEqual(selected["rejected"][0]["task"], tasks[1])
            system = deepmd / selected["selected"][0]["system"]
            np.testing.assert_allclose(
                np.load(system / "set.000/spin.npy"),
                [[0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]],
            )
            np.testing.assert_allclose(
                np.load(system / "set.000/spin_force.npy"),
                [[-0.4, 0.2, 0.6, 0, 0, 0], [-0.4, 0.2, 0.6, 0, 0, 0]],
            )
            np.testing.assert_allclose(
                np.load(system / "set.000/spin_length.npy"),
                [[1.003, 0], [0.997, 0]],
            )
            self.assertTrue((system / "spin.raw").is_file())
            self.assertTrue((system / "spin_force.raw").is_file())
            self.assertTrue((system / "spin_length.raw").is_file())
            self.assertTrue((system / "set.000/coord.npy").is_file())
            frames = json.loads((system / "frames.json").read_text())
            self.assertEqual(frames[0]["outcar_frame"], 4)
            self.assertEqual([row["task"] for row in frames], [tasks[0], tasks[2]])

    def test_existing_output_stage_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage = root / "03.spin"
            stage.mkdir()
            (stage / "tasks.json").write_text(
                json.dumps({"tasks": ["scale-1.000/000000/00/origin/000000"]})
            )
            output = root / "04.data"
            output.mkdir()
            sentinel = output / "keep.txt"
            sentinel.write_text("existing")
            with self.assertRaisesRegex(RuntimeError, "already exists"):
                spin_data.collect_spin_data({"out_dir": str(root)})
            self.assertEqual(sentinel.read_text(), "existing")

    def test_cli_stage_five_requires_convert_machine_but_not_md_incar_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            param = root / "param.json"
            param.write_text(
                json.dumps(
                    {
                        "stages": [5],
                        "out_dir": str(root),
                        "from_poscar_path": "not-read",
                        "super_cell": [1, 1, 1],
                        "scale": [1.0],
                        "pert_numb": 0,
                        "pert_box": 0.0,
                        "pert_atom": 0.0,
                        "md_incar": "not-read",
                        "md_nstep": 1,
                        "potcars": ["not-read"],
                    }
                )
            )
            with self.assertRaisesRegex(ValueError, "convert-data"):
                spin_init.gen_spin_init(
                    argparse.Namespace(PARAM=str(param), MACHINE=None)
                )
