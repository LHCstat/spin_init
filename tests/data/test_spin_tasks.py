"""Stage 4 contracts, without requiring a VASP executable."""

import argparse
import json
import logging
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from dpgen.data import spin_init
from tests.data.test_spin_init import POSCAR_TEXT, _skip_without_symlink_privilege

TEMPLATE = """SYSTEM = magnetic_test
ENCUT = 400
LNONCOLLINEAR = .TRUE.
NSW = 0
IBRION = -1
I_CONSTRAINED_M = 2
LAMBDA = 10
RWIGS = 1.2
MAGMOM = \\
0 0 2 \\
3*0.0
M_CONSTR = 0 0 2 3*0.0
"""


class TestSpinTasks(unittest.TestCase):
    def setUp(self):
        # tests.generator disables all logging at import during discovery.
        # Restore its global state afterward, while keeping real log assertions.
        self.addCleanup(logging.disable, logging.root.manager.disable)
        logging.disable(logging.NOTSET)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.incar = self.root / "INCAR_scf_tmp"
        self.incar.write_text(TEMPLATE)
        self.potcar = self.root / "POTCAR_Fe"
        self.potcar.write_bytes(b"test potential\n")
        self.kpoints = self.root / "KPOINTS"
        self.kpoints.write_text("Gamma\n0\nGamma\n1 1 1\n0 0 0\n")
        for frame in ("00", "01"):
            path = self.root / "02.disp/scale-1.000/000000" / frame
            path.mkdir(parents=True)
            (path / "POSCAR").write_text(POSCAR_TEXT)
        self.jdata = {
            "stages": [4],
            "out_dir": str(self.root),
            "from_poscar_path": "unused",
            "super_cell": [1, 1, 1],
            "scale": [1.0],
            "pert_numb": 0,
            "pert_box": 0.0,
            "pert_atom": 0.0,
            "md_incar": "missing-md-incar",
            "md_nstep": 3,
            "spin_incar": str(self.incar),
            "spin_pert_numb": 0,
            "potcars": [str(self.potcar)],
        }

    def module(self):
        from dpgen.data import spin_tasks

        return spin_tasks

    def test_parse_continuations_repetitions_and_scientific_notation(self):
        module = self.module()
        for value, expected in [("2", 2.0), ("1e-3", 0.001), ("-2D+0", -2.0)]:
            self.incar.write_text(TEMPLATE.replace("0 0 2", "0 0 " + value))
            incar, moments = module.read_spin_incar(self.incar, 2)
            np.testing.assert_allclose(moments, [[0, 0, expected], [0, 0, 0]])
            self.assertEqual(incar["ENCUT"], 400)

    def test_bad_templates_fail_with_file_context(self):
        module = self.module()
        bad = [
            "",
            TEMPLATE.replace("MAGMOM", "MAMGOM"),
            TEMPLATE.replace("3*0.0", "2*0.0"),
            TEMPLATE.replace("M_CONSTR = 0 0 2", "M_CONSTR = 0 0 3"),
            TEMPLATE.replace("0 0 2", "0 0 nan"),
            TEMPLATE.replace(".TRUE.", ".FALSE."),
            TEMPLATE.replace("NSW = 0", "NSW = -1"),
            TEMPLATE + "MAGMOM = 6*0\n",
        ]
        for text in bad:
            with self.subTest(text=text):
                self.incar.write_text(text)
                with self.assertRaisesRegex(ValueError, "INCAR_scf_tmp"):
                    module.read_spin_incar(self.incar, 2)

    def test_write_syncs_vectors_and_preserves_other_parameters(self):
        module = self.module()
        incar, _ = module.read_spin_incar(self.incar, 2)
        target = self.root / "generated.INCAR"
        module.write_spin_incar(incar, [[1, 0, 2], [0, 0, 0]], target)
        actual, vectors = module.read_spin_incar(target, 2)
        np.testing.assert_array_equal(vectors, [[1, 0, 2], [0, 0, 0]])
        self.assertEqual(actual["ENCUT"], 400)
        self.assertEqual(actual["I_CONSTRAINED_M"], 2)
        self.assertEqual(self.incar.read_text(), TEMPLATE)

    def test_relaxation_templates_preserve_ionic_settings(self):
        module = self.module()
        for ibrion in (1, 2, 3):
            with self.subTest(ibrion=ibrion):
                text = (
                    TEMPLATE.replace("NSW = 0", "NSW = 50").replace(
                        "IBRION = -1", f"IBRION = {ibrion}"
                    )
                    + "ISIF = 3\nEDIFFG = -0.02\n"
                )
                self.incar.write_text(text)
                incar, _ = module.read_spin_incar(self.incar, 2)
                target = self.root / "relax.INCAR"
                module.write_spin_incar(incar, [[2, 0, 0], [0, 0, 0]], target)
                actual, moments = module.read_spin_incar(target, 2)
                self.assertEqual(actual["NSW"], 50)
                self.assertEqual(actual["IBRION"], ibrion)
                self.assertEqual(actual["ISIF"], 3)
                self.assertEqual(actual["EDIFFG"], -0.02)
                np.testing.assert_array_equal(moments, [[2, 0, 0], [0, 0, 0]])
                self.assertEqual(self.incar.read_text(), text)

    def test_zero_steps_does_not_require_or_inject_ibrion(self):
        module = self.module()
        for suffix in ("IBRION = 2\n", ""):
            self.incar.write_text(TEMPLATE.replace("IBRION = -1\n", suffix))
            incar, moments = module.read_spin_incar(self.incar, 2)
            target = self.root / "static.INCAR"
            module.write_spin_incar(incar, moments, target)
            actual, _ = module.read_spin_incar(target, 2)
            self.assertEqual("IBRION" in actual, bool(suffix))

    def test_all_snapshots_receive_named_configs_without_writing(self):
        module = self.module()
        self.jdata["spin_pert_numb"] = 1

        def provider(moments, count):
            return {"R1": [[2, 0, 0], [0, 0, 0]]}

        tasks = module.plan_spin_tasks(self.jdata, {}, perturb=provider)
        self.assertEqual(
            [t["task"] for t in tasks],
            [
                "scale-1.000/000000/00/origin/000000",
                "scale-1.000/000000/00/R1",
                "scale-1.000/000000/01/origin/000000",
                "scale-1.000/000000/01/R1",
            ],
        )
        self.assertFalse((self.root / "03.spin").exists())
        np.testing.assert_array_equal(tasks[1]["moments"], [[2, 0, 0], [0, 0, 0]])

    def test_multiple_spin_incars_add_indexed_layer_and_keep_initial_moments(self):
        module = self.module()
        second_incar = self.root / "INCAR_second"
        second_incar.write_text(TEMPLATE.replace("0 0 2", "0 2 0"))
        self.jdata["spin_incar"] = [str(self.incar), str(second_incar)]

        tasks = module.plan_spin_tasks(self.jdata, {})

        self.assertEqual(
            [task["task"] for task in tasks],
            [
                "scale-1.000/000000/00/incar-000/origin/000000",
                "scale-1.000/000000/00/incar-001/origin/000000",
                "scale-1.000/000000/01/incar-000/origin/000000",
                "scale-1.000/000000/01/incar-001/origin/000000",
            ],
        )
        np.testing.assert_array_equal(tasks[0]["moments"], [[0, 0, 2], [0, 0, 0]])
        np.testing.assert_array_equal(tasks[1]["moments"], [[0, 2, 0], [0, 0, 0]])
        self.assertFalse((self.root / "03.spin").exists())

    def test_multiple_spin_incars_share_one_reproducible_rng_sequence(self):
        from dpgen.data.spin_perturb import build_spin_perturbation

        module = self.module()
        second_incar = self.root / "INCAR_same_moments"
        second_incar.write_text(TEMPLATE)
        self.jdata["spin_incar"] = [str(self.incar), str(second_incar)]
        parameters = [{"Random": {"num": 1, "seed": 2468}}]
        count_a, provider_a = build_spin_perturbation(parameters)
        count_b, provider_b = build_spin_perturbation(parameters)
        self.jdata["spin_pert_numb"] = count_a

        tasks_a = module.plan_spin_tasks(self.jdata, {}, perturb=provider_a)
        tasks_b = module.plan_spin_tasks(self.jdata, {}, perturb=provider_b)
        random_a = [
            task["moments"] for task in tasks_a if task["task"].endswith("Rand1")
        ]
        random_b = [
            task["moments"] for task in tasks_b if task["task"].endswith("Rand1")
        ]

        self.assertEqual(len(random_a), 4)
        for moments_a, moments_b in zip(random_a, random_b):
            np.testing.assert_array_equal(moments_a, moments_b)
        self.assertFalse(np.array_equal(random_a[0], random_a[1]))

    def test_multiple_spin_incar_input_rejects_empty_duplicate_and_bad_entries(self):
        module = self.module()
        invalid_inputs = [
            ([], "nonempty"),
            ([str(self.incar), str(self.incar)], "duplicate"),
            ([str(self.incar), 3], r"spin_incar\[1\]"),
        ]
        for spin_incar, message in invalid_inputs:
            self.jdata["spin_incar"] = spin_incar
            with (
                self.subTest(spin_incar=spin_incar),
                self.assertRaisesRegex(ValueError, message),
            ):
                module.plan_spin_tasks(self.jdata, {})
        self.assertFalse((self.root / "03.spin").exists())

    def test_bad_second_spin_incar_reports_index_path_and_snapshot(self):
        module = self.module()
        bad_incar = self.root / "INCAR_bad"
        bad_incar.write_text(TEMPLATE.replace("MAGMOM", "MAMGOM"))
        self.jdata["spin_incar"] = [str(self.incar), str(bad_incar)]

        with self.assertRaisesRegex(
            ValueError, r"scale-1\.000/000000/00.*incar-001.*INCAR_bad"
        ):
            module.plan_spin_tasks(self.jdata, {})
        self.assertFalse((self.root / "03.spin").exists())

    def test_missing_second_spin_incar_reports_index_path_and_snapshot(self):
        module = self.module()
        missing_incar = self.root / "INCAR_missing"
        self.jdata["spin_incar"] = [str(self.incar), str(missing_incar)]

        with self.assertRaisesRegex(
            FileNotFoundError,
            r"scale-1\.000/000000/00.*incar-001.*INCAR_missing",
        ):
            module.plan_spin_tasks(self.jdata, {})
        self.assertFalse((self.root / "03.spin").exists())

    def test_independent_groups_plan_one_baseline_and_each_operation(self):
        from dpgen.data.spin_perturb import build_spin_perturbation

        module = self.module()
        count, provider = build_spin_perturbation(
            [
                {"Rotation": {"angle": 90, "axis": [0, 1, 0]}},
                {"Scale": {"pert": 0.1, "pert_step": 0.1}},
            ]
        )
        self.jdata["spin_pert_numb"] = count

        tasks = module.plan_spin_tasks(self.jdata, {}, perturb=provider)

        self.assertEqual(
            [task["task"] for task in tasks],
            [
                "scale-1.000/000000/00/origin/000000",
                "scale-1.000/000000/00/Rotation-000/R1",
                "scale-1.000/000000/00/Scale-001/S1",
                "scale-1.000/000000/00/Scale-001/S2",
                "scale-1.000/000000/01/origin/000000",
                "scale-1.000/000000/01/Rotation-000/R1",
                "scale-1.000/000000/01/Scale-001/S1",
                "scale-1.000/000000/01/Scale-001/S2",
            ],
        )
        np.testing.assert_array_equal(tasks[0]["moments"], [[0, 0, 2], [0, 0, 0]])
        np.testing.assert_allclose(
            tasks[1]["moments"], [[2, 0, 0], [0, 0, 0]], atol=1e-14
        )
        np.testing.assert_allclose(tasks[2]["moments"], [[0, 0, 1.8], [0, 0, 0]])
        np.testing.assert_allclose(tasks[3]["moments"], [[0, 0, 2.2], [0, 0, 0]])
        self.assertFalse((self.root / "03.spin").exists())

    def test_multiple_incars_each_have_baseline_and_independent_groups(self):
        from dpgen.data.spin_perturb import build_spin_perturbation

        second_incar = self.root / "INCAR_second"
        second_incar.write_text(TEMPLATE.replace("0 0 2", "0 2 0"))
        self.jdata["spin_incar"] = [str(self.incar), str(second_incar)]
        count, provider = build_spin_perturbation(
            [
                {"Rotation": {"angle": 90, "axis": [1, 0, 0]}},
                {"Scale": {"pert": 0.1, "pert_step": 0.1}},
            ]
        )
        self.jdata["spin_pert_numb"] = count

        tasks = self.module().plan_spin_tasks(self.jdata, {}, perturb=provider)

        self.assertEqual(len(tasks), 16)
        self.assertEqual(
            [task["task"] for task in tasks[:8]],
            [
                "scale-1.000/000000/00/incar-000/origin/000000",
                "scale-1.000/000000/00/incar-000/Rotation-000/R1",
                "scale-1.000/000000/00/incar-000/Scale-001/S1",
                "scale-1.000/000000/00/incar-000/Scale-001/S2",
                "scale-1.000/000000/00/incar-001/origin/000000",
                "scale-1.000/000000/00/incar-001/Rotation-000/R1",
                "scale-1.000/000000/00/incar-001/Scale-001/S1",
                "scale-1.000/000000/00/incar-001/Scale-001/S2",
            ],
        )
        np.testing.assert_allclose(tasks[2]["moments"], [[0, 0, 1.8], [0, 0, 0]])
        np.testing.assert_allclose(tasks[6]["moments"], [[0, 1.8, 0], [0, 0, 0]])

    def test_make_writes_independent_incars_and_manifest(self):
        from dpgen.data.spin_perturb import build_spin_perturbation

        module = self.module()
        second_incar = self.root / "INCAR_second"
        second_incar.write_text(TEMPLATE.replace("0 0 2", "0 2 0"))
        self.jdata["spin_incar"] = [str(self.incar), str(second_incar)]
        count, provider = build_spin_perturbation(
            [
                {"Rotation": {"angle": 90, "axis": [0, 1, 0]}},
                {"Scale": {"pert": 0.1, "pert_step": 0.1}},
            ]
        )
        self.jdata["spin_pert_numb"] = count
        # Only the OS symlink operation is disabled here. The real link
        # contract has separate capability-gated Linux tests below.
        with mock.patch.object(module, "_make_relative_symlink"):
            paths = module.make_spin_tasks(self.jdata, {}, perturb=provider)

        stage = self.root / "03.spin"
        self.assertEqual(len(paths), 16)
        self.assertEqual(json.loads((stage / "tasks.json").read_text())["tasks"], paths)
        self.assertEqual(module._saved_tasks(self.jdata)[1], paths)
        parent = stage / "scale-1.000/000000/00/incar-000"
        for relative, expected in (
            ("origin/000000", [[0, 0, 2], [0, 0, 0]]),
            ("Rotation-000/R1", [[2, 0, 0], [0, 0, 0]]),
            ("Scale-001/S1", [[0, 0, 1.8], [0, 0, 0]]),
            ("Scale-001/S2", [[0, 0, 2.2], [0, 0, 0]]),
        ):
            incar, moments = module.read_spin_incar(parent / relative / "INCAR", 2)
            np.testing.assert_allclose(moments, expected, atol=1e-14)
            self.assertEqual(incar["ENCUT"], 400)
        self.assertEqual(self.incar.read_text(), TEMPLATE)
        for frame in ("00", "01"):
            for label in ("incar-000", "incar-001"):
                parent = stage / "scale-1.000/000000" / frame / label
                self.assertTrue((parent / "origin/000000/INCAR").is_file())
                self.assertFalse((parent / "000000").exists())

    def test_make_single_incar_writes_baseline_inside_origin(self):
        module = self.module()
        with mock.patch.object(module, "_make_relative_symlink"):
            paths = module.make_spin_tasks(self.jdata, {})
        expected = [
            "scale-1.000/000000/00/origin/000000",
            "scale-1.000/000000/01/origin/000000",
        ]
        self.assertEqual(paths, expected)
        stage = self.root / "03.spin"
        self.assertEqual(
            json.loads((stage / "tasks.json").read_text())["tasks"], expected
        )
        self.assertEqual(module._saved_tasks(self.jdata)[1], expected)
        for task in expected:
            incar, moments = module.read_spin_incar(stage / task / "INCAR", 2)
            np.testing.assert_array_equal(moments, [[0, 0, 2], [0, 0, 0]])
            self.assertEqual(incar["ENCUT"], 400)
            self.assertFalse((stage / task).parent.parent.joinpath("000000").exists())

    def test_discovers_extra_snapshot_parents_and_sorts_frames_numerically(self):
        module = self.module()
        for parent, frames in {
            "scale-1.000/000000": ("11", "100"),
            "scale-1.000/000001": ("00",),
        }.items():
            for frame in frames:
                path = self.root / "02.disp" / parent / frame
                path.mkdir(parents=True)
                (path / "POSCAR").write_text(POSCAR_TEXT)

        tasks = module.plan_spin_tasks(self.jdata, {})

        self.assertEqual(
            [task["task"] for task in tasks],
            [
                "scale-1.000/000000/00/origin/000000",
                "scale-1.000/000000/01/origin/000000",
                "scale-1.000/000000/11/origin/000000",
                "scale-1.000/000000/100/origin/000000",
                "scale-1.000/000001/00/origin/000000",
            ],
        )

    def test_no_placeholder_perturbations_or_path_escape(self):
        module = self.module()
        self.jdata["spin_pert_numb"] = 1
        with self.assertRaises(NotImplementedError):
            module.plan_spin_tasks(self.jdata, {})
        for name in (
            "../outside",
            "000000",
            "origin",
            "origin/000000",
            "R/1",
            "Rotation-000/../outside",
            "Unknown-000/R1",
            "Rotation-000/R1/extra",
        ):
            with self.subTest(name=name), self.assertRaises(ValueError):
                module.plan_spin_tasks(self.jdata, {}, perturb=lambda m, n: {name: m})
        self.assertFalse((self.root / "03.spin").exists())

    def test_deprecated_canting_parameters_are_rejected_before_task_planning(self):
        from dpgen.data.spin_perturb import build_spin_perturbation

        with self.assertRaisesRegex(ValueError, "Rcut"):
            build_spin_perturbation([{"Canting": {"angle": 30, "Rcut": 0.5}}])
        self.assertFalse((self.root / "03.spin").exists())

    def test_missing_snapshot_and_reserved_forward_file_fail_before_writing(self):
        module = self.module()
        with self.assertRaisesRegex(ValueError, "INCAR"):
            module.plan_spin_tasks(
                self.jdata,
                {"fp_user_forward_files": [str(self.incar.parent / "INCAR")]},
            )
        (self.root / "02.disp/scale-1.000/000000/01/POSCAR").unlink()
        with self.assertRaisesRegex(FileNotFoundError, "000000.*01"):
            module.plan_spin_tasks(self.jdata, {})

    def test_real_symlinks_and_per_task_incar(self):
        _skip_without_symlink_privilege()
        module = self.module()
        paths = module.make_spin_tasks(self.jdata, {})
        module.materialize_spin_forward_files(
            self.root / "03.spin",
            paths,
            {"KPOINTS": self.kpoints.resolve()},
        )
        for task in paths:
            path = self.root / "03.spin" / task
            self.assertTrue((path / "POSCAR").is_symlink())
            self.assertTrue((path / "POTCAR").is_symlink())
            self.assertTrue((path / "KPOINTS").is_symlink())
            self.assertFalse((path / "INCAR").is_symlink())
            self.assertEqual((path / "POSCAR").read_text(), POSCAR_TEXT)
            self.assertEqual((path / "POTCAR").read_bytes(), self.potcar.read_bytes())
            module.read_spin_incar(path / "INCAR", 2)
        with self.assertRaisesRegex(RuntimeError, "already exists"):
            module.make_spin_tasks(self.jdata, {})

    def test_multiple_spin_incars_create_real_symlinks_in_indexed_layers(self):
        _skip_without_symlink_privilege()
        module = self.module()
        second_incar = self.root / "INCAR_second"
        second_incar.write_text(TEMPLATE.replace("0 0 2", "0 2 0"))
        self.jdata["spin_incar"] = [str(self.incar), str(second_incar)]
        from dpgen.data.spin_perturb import build_spin_perturbation

        count, provider = build_spin_perturbation(
            [
                {"Rotation": {"angle": 90, "axis": [0, 1, 0]}},
                {"Scale": {"pert": 0.1, "pert_step": 0.1}},
            ]
        )
        self.jdata["spin_pert_numb"] = count

        paths = module.make_spin_tasks(self.jdata, {}, perturb=provider)

        self.assertIn("scale-1.000/000000/00/incar-000/origin/000000", paths)
        self.assertIn("scale-1.000/000000/00/incar-001/origin/000000", paths)
        self.assertIn("scale-1.000/000000/00/incar-001/Rotation-000/R1", paths)
        self.assertEqual(len(paths), 16)
        saved = json.loads((self.root / "03.spin/tasks.json").read_text())
        self.assertEqual(saved["tasks"], paths)
        self.assertEqual(module._saved_tasks(self.jdata)[1], paths)
        for task in paths:
            path = self.root / "03.spin" / task
            self.assertTrue((path / "POSCAR").is_symlink())
            self.assertTrue((path / "POTCAR").is_symlink())
            self.assertFalse((path / "INCAR").is_symlink())
            self.assertFalse(Path((path / "POSCAR").readlink()).is_absolute())
            self.assertFalse(Path((path / "POTCAR").readlink()).is_absolute())
            _, moments = module.read_spin_incar(path / "INCAR", 2)
            if task.endswith("incar-000/Rotation-000/R1"):
                np.testing.assert_allclose(moments, [[2, 0, 0], [0, 0, 0]], atol=1e-14)
            if task.endswith("incar-000/Scale-001/S1"):
                np.testing.assert_allclose(moments, [[0, 0, 1.8], [0, 0, 0]])

    def test_run_materializes_forward_files_missing_from_existing_tasks(self):
        module = self.module()
        path = self.prepare_run_inputs()
        (path / "KPOINTS").unlink()
        machine = {
            "fp_machine": {
                "batch_type": "Shell",
                "context_type": "LocalContext",
                "local_root": "./",
                "remote_root": str(self.root / "remote"),
            },
            "fp_resources": {
                "number_node": 1,
                "cpu_per_node": 1,
                "group_size": 1,
            },
            "fp_command": "vasp_ncl",
            "fp_group_size": 1,
            "fp_user_forward_files": [str(self.kpoints)],
        }
        with mock.patch.object(module, "_make_relative_symlink") as make_link:
            module.make_spin_submission(self.jdata, machine)
        make_link.assert_called_once_with(path / "KPOINTS", self.kpoints.resolve())

    def prepare_run_inputs(self, task="scale-1.000/000000/00/origin/000000"):
        stage = self.root / "03.spin"
        path = stage / task
        path.mkdir(parents=True)
        (path / "POSCAR").write_text(POSCAR_TEXT)
        (path / "INCAR").write_text(TEMPLATE)
        (path / "POTCAR").write_bytes(self.potcar.read_bytes())
        (path / "KPOINTS").write_text(self.kpoints.read_text())
        (stage / "tasks.json").write_text(json.dumps({"tasks": [task]}))
        return path

    def test_grouped_and_multi_incar_tasks_are_submitted_from_manifest(self):
        module = self.module()
        paths = [
            "scale-1.000/000000/00/Rotation-000/R1",
            "scale-1.000/000000/00/incar-001/Scale-001/S1",
            "scale-1.000/000000/01/incar-000/000000",
            "scale-1.000/000000/01/R1-S1",
            "scale-1.000/000000/00/origin/000000",
            "scale-1.000/000000/00/incar-001/origin/000000",
            "scale-1.000/000000/01/000000",
        ]
        for task in paths:
            path = self.prepare_run_inputs(task)
            (path / "OUTCAR").write_text("TOTAL-FORCE\nElapse\n")
            (path / "OSZICAR").write_text(" 1 F= -1 E0= -1\n")
        (self.root / "03.spin/tasks.json").write_text(json.dumps({"tasks": paths}))
        machine = {
            "fp_machine": {
                "batch_type": "Shell",
                "context_type": "LocalContext",
                "local_root": "./",
                "remote_root": str(self.root / "remote"),
            },
            "fp_resources": {"number_node": 1, "cpu_per_node": 1, "group_size": 1},
            "fp_command": "vasp_ncl",
            "fp_group_size": 1,
        }

        submission = module.make_spin_submission(self.jdata, machine)

        self.assertEqual(
            [task.task_work_path for task in submission.belonging_tasks], paths
        )
        self.assertTrue(
            all(
                task.backward_files == ["OUTCAR", "OSZICAR"]
                for task in submission.belonging_tasks
            )
        )
        self.assertEqual(module.check_spin_results(self.jdata), 7)

    def test_task_manifest_rejects_path_escape_and_unrecognized_group(self):
        stage = self.root / "03.spin"
        stage.mkdir()
        for task in (
            "scale-1.000/000000/00/incar-000/Rotation-000/../outside",
            "scale-1.000/000000/00/Unknown-000/R1",
            "scale-1.000/000000/00/Rotation-000/R1/extra",
            "scale-1.000/000000/00/origin/R1",
            "scale-1.000/000000/00/origin/000001",
            "scale-1.000/000000/00/origin/000000/extra",
            "scale-1.000/000000/00/incar-000/origin/../outside",
        ):
            (stage / "tasks.json").write_text(json.dumps({"tasks": [task]}))
            with (
                self.subTest(task=task),
                self.assertRaisesRegex(ValueError, "tasks.json"),
            ):
                self.module()._saved_tasks(self.jdata)

    def test_real_dispatcher_objects_use_spin_outputs_and_command(self):
        module = self.module()
        self.prepare_run_inputs()
        machine = {
            "fp_machine": {
                "batch_type": "Shell",
                "context_type": "LocalContext",
                "local_root": "./",
                "remote_root": str(self.root / "remote"),
            },
            "fp_resources": {
                "number_node": 1,
                "cpu_per_node": 1,
                "gpu_per_node": 0,
                "group_size": 1,
            },
            "fp_command": "vasp_ncl",
            "fp_group_size": 1,
            "fp_user_forward_files": [str(self.kpoints)],
            "fp_user_backward_files": ["vasprun.xml"],
        }
        submission = module.make_spin_submission(self.jdata, machine)
        task = submission.belonging_tasks[0]
        self.assertEqual(task.command, "vasp_ncl")
        self.assertEqual(task.task_work_path, "scale-1.000/000000/00/origin/000000")
        self.assertEqual(task.backward_files, ["OUTCAR", "OSZICAR", "vasprun.xml"])
        self.assertEqual(task.forward_files, ["POSCAR", "INCAR", "POTCAR", "KPOINTS"])

    def test_result_check_requires_both_outputs_and_normal_termination(self):
        module = self.module()
        path = self.prepare_run_inputs()
        with self.assertRaisesRegex(
            FileNotFoundError, "scale-1.000/000000/00/origin/000000"
        ):
            module.check_spin_results(self.jdata)
        (path / "OUTCAR").write_text("TOTAL-FORCE\nElapse\n")
        (path / "OSZICAR").write_text("")
        with self.assertRaises(RuntimeError):
            module.check_spin_results(self.jdata)
        (path / "OSZICAR").write_text(" 1 F= -1 E0= -1\n")
        self.assertEqual(module.check_spin_results(self.jdata), 1)

    def prepare_relax_results(self, converged=True):
        path = self.prepare_run_inputs()
        (path / "INCAR").write_text(
            TEMPLATE.replace("NSW = 0", "NSW = 50").replace("IBRION = -1", "IBRION = 2")
            + "ISIF = 3\n"
        )
        marker = (
            "reached required accuracy - stopping structural energy minimisation\n"
            if converged
            else ""
        )
        (path / "OUTCAR").write_text("TOTAL-FORCE\nTOTAL-FORCE\n" + marker + "Elapse\n")
        (path / "OSZICAR").write_text(" 1 F= -1\n 2 F= -2\n")
        (path / "CONTCAR").write_text(POSCAR_TEXT.replace("2.0", "2.1"))
        return path

    def test_relaxation_early_convergence_accepts_changed_cell(self):
        path = self.prepare_relax_results()
        with self.assertLogs("dpgen", level="INFO") as logs:
            self.assertEqual(self.module().check_spin_results(self.jdata), 1)
        self.assertTrue(any("convergence confirmed" in line for line in logs.output))
        self.assertEqual((path / "POSCAR").read_text(), POSCAR_TEXT)

    def test_normal_relaxation_exit_without_convergence_warns(self):
        path = self.prepare_relax_results(converged=False)
        with self.assertLogs("dpgen", level="WARNING") as logs:
            self.assertEqual(self.module().check_spin_results(self.jdata), 1)
        self.assertTrue(
            any(
                "not confirmed" in line
                and "OUTCAR" in line
                and "scale-1.000/000000/00/origin/000000" in line
                for line in logs.output
            )
        )
        self.assertEqual((path / "POSCAR").read_text(), POSCAR_TEXT)

    def test_relaxation_requires_valid_matching_contcar(self):
        path = self.prepare_relax_results()
        contcar = path / "CONTCAR"
        for text in (
            None,
            "",
            "not a POSCAR",
            POSCAR_TEXT.replace("Fe\n", "Ni\n"),
            POSCAR_TEXT.replace("Fe\n2\n", "Fe\n1\n"),
        ):
            with self.subTest(text=text):
                if text is None:
                    contcar.unlink()
                else:
                    contcar.write_text(text)
                with self.assertRaisesRegex(
                    (FileNotFoundError, ValueError),
                    "scale-1.000/000000/00/origin/000000.*CONTCAR",
                ):
                    self.module().check_spin_results(self.jdata)

    def test_relaxation_still_requires_normal_termination_and_force_output(self):
        path = self.prepare_relax_results()
        for text in ("TOTAL-FORCE\n", "Elapse\n", "TOTAL-FORCE\nElapse\nElapse\n"):
            with self.subTest(text=text):
                (path / "OUTCAR").write_text(text)
                with self.assertRaisesRegex(RuntimeError, "OUTCAR"):
                    self.module().check_spin_results(self.jdata)

    def test_mixed_submission_automatically_retrieves_contcar_once(self):
        relax = self.prepare_relax_results()
        static_task = "scale-1.000/000000/01/origin/000000"
        self.prepare_run_inputs(static_task)
        paths = [relax.relative_to(self.root / "03.spin").as_posix(), static_task]
        (self.root / "03.spin/tasks.json").write_text(json.dumps({"tasks": paths}))
        machine = {
            "fp_machine": {
                "batch_type": "Shell",
                "context_type": "LocalContext",
                "local_root": "./",
                "remote_root": str(self.root / "remote"),
            },
            "fp_resources": {"number_node": 1, "cpu_per_node": 1, "group_size": 1},
            "fp_command": "vasp_ncl",
            "fp_group_size": 1,
            "fp_user_backward_files": ["CONTCAR", "vasprun.xml"],
        }
        for extras in ([], ["CONTCAR", "vasprun.xml"]):
            machine["fp_user_backward_files"] = extras
            submission = self.module().make_spin_submission(self.jdata, machine)
            for task in submission.belonging_tasks:
                self.assertEqual(task.backward_files.count("CONTCAR"), 1)
                self.assertEqual(
                    task.backward_files[:3], ["OUTCAR", "OSZICAR", "CONTCAR"]
                )

    def test_relaxation_checks_do_not_replace_real_poscar_symlink(self):
        _skip_without_symlink_privilege()
        path = self.prepare_relax_results()
        source = self.root / "02.disp/scale-1.000/000000/00/POSCAR"
        (path / "POSCAR").unlink()
        self.module()._make_relative_symlink(path / "POSCAR", source)
        self.assertEqual(self.module().check_spin_results(self.jdata), 1)
        self.assertTrue((path / "POSCAR").is_symlink())
        self.assertEqual(source.read_text(), POSCAR_TEXT)

    def test_malformed_task_manifest_has_context(self):
        module = self.module()
        stage = self.root / "03.spin"
        stage.mkdir()
        manifest = stage / "tasks.json"
        for text in ("not json", "{}", '{"tasks": "not-a-list"}'):
            with self.subTest(text=text):
                manifest.write_text(text)
                with self.assertRaisesRegex(ValueError, "tasks.json"):
                    module.check_spin_results(self.jdata)

    def test_stage4_run_does_not_read_md_incar_or_recreate_tasks(self):
        module = self.module()
        path = self.prepare_run_inputs()
        self.jdata["spin_action"] = "run"
        param = self.root / "input.json"
        param.write_text(json.dumps(self.jdata))
        machine = self.root / "machine.json"
        machine.write_text(
            json.dumps({"fp": {"command": "vasp_ncl", "machine": {}, "resources": {}}})
        )
        with mock.patch.object(module, "run_spin_tasks") as run:
            spin_init.gen_spin_init(
                argparse.Namespace(PARAM=str(param), MACHINE=str(machine))
            )
        self.assertEqual(run.call_args.args[1]["fp_command"], "vasp_ncl")
        self.assertEqual((path / "INCAR").read_text(), TEMPLATE)

    def test_stage4_can_use_separate_spin_machine_entry(self):
        module = self.module()
        self.prepare_run_inputs()
        self.jdata["spin_action"] = "run"
        param = self.root / "input.json"
        param.write_text(json.dumps(self.jdata))
        machine = self.root / "machine.json"
        machine.write_text(
            json.dumps(
                {
                    "api_version": "1.0",
                    "fp": {
                        "command": "vasp_std",
                        "machine": {},
                        "resources": {},
                    },
                    "spin": {
                        "command": "vasp_ncl",
                        "machine": {},
                        "resources": {"group_size": 2},
                    },
                }
            )
        )
        with mock.patch.object(module, "run_spin_tasks") as run:
            spin_init.gen_spin_init(
                argparse.Namespace(PARAM=str(param), MACHINE=str(machine))
            )
        selected = run.call_args.args[1]
        self.assertEqual(selected["fp_command"], "vasp_ncl")
        self.assertEqual(selected["fp_group_size"], 2)

    def test_stage4_make_run_without_machine_only_prepares_tasks(self):
        module = self.module()
        self.jdata["spin_action"] = "make_run"
        param = self.root / "input.json"
        param.write_text(json.dumps(self.jdata))
        with (
            mock.patch.object(module, "make_spin_tasks") as make,
            mock.patch.object(module, "run_spin_tasks") as run,
        ):
            spin_init.gen_spin_init(argparse.Namespace(PARAM=str(param), MACHINE=None))
        make.assert_called_once()
        run.assert_not_called()

    def test_stage4_wires_independent_group_count_and_names_into_task_generation(self):
        module = self.module()
        self.jdata["spin_action"] = "make"
        self.jdata["pert_spin"] = [
            {"Rotation": {"angle": [90, 180], "axis": [0, 0, 1]}},
            {"Scale": {"pert": 0.1, "pert_step": 0.1}},
        ]
        param = self.root / "input.json"
        param.write_text(json.dumps(self.jdata))
        with mock.patch.object(module, "make_spin_tasks") as make:
            spin_init.gen_spin_init(argparse.Namespace(PARAM=str(param), MACHINE=None))

        generated_jdata, _ = make.call_args.args
        provider = make.call_args.kwargs["perturb"]
        self.assertEqual(generated_jdata["spin_pert_numb"], 4)
        self.assertEqual(
            list(provider(np.array([[1.0, 0, 0]]), 4)),
            ["Rotation-000/R1", "Rotation-000/R2", "Scale-001/S1", "Scale-001/S2"],
        )

    def test_invalid_spin_parameters_fail_before_an_earlier_stage_runs(self):
        self.jdata["stages"] = [1, 4]
        self.jdata["spin_action"] = "make"
        self.jdata["md_incar"] = str(self.incar)
        self.jdata["pert_spin"] = [{"Rotation": {"angle": 30, "axis": [0, 0, 0]}}]
        param = self.root / "input.json"
        param.write_text(json.dumps(self.jdata))

        with mock.patch.object(spin_init, "make_spin_init_structures") as make:
            with self.assertRaisesRegex(ValueError, r"pert_spin\[0\]\.Rotation\.axis"):
                spin_init.gen_spin_init(
                    argparse.Namespace(PARAM=str(param), MACHINE=None)
                )
        make.assert_not_called()

    def test_run_action_still_validates_spin_parameters(self):
        self.jdata["stages"] = [4]
        self.jdata["spin_action"] = "run"
        self.jdata["pert_spin"] = [{"Rotation": {"angle": 30, "axis": [0, 0, 0]}}]
        param = self.root / "input.json"
        machine = self.root / "machine.json"
        param.write_text(json.dumps(self.jdata))
        machine.write_text(json.dumps({}))

        with self.assertRaisesRegex(ValueError, r"pert_spin\[0\]\.Rotation\.axis"):
            spin_init.gen_spin_init(
                argparse.Namespace(PARAM=str(param), MACHINE=str(machine))
            )


if __name__ == "__main__":
    unittest.main()
