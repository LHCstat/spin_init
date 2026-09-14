"""Stage 4 contracts, without requiring a VASP executable."""

import argparse
import json
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
            TEMPLATE.replace("NSW = 0", "NSW = 5"),
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

    def test_all_snapshots_receive_named_configs_without_writing(self):
        module = self.module()
        self.jdata["spin_pert_numb"] = 1

        def provider(moments, count):
            return {"R1": [[2, 0, 0], [0, 0, 0]]}

        tasks = module.plan_spin_tasks(self.jdata, {}, perturb=provider)
        self.assertEqual(
            [t["task"] for t in tasks],
            [
                "scale-1.000/000000/00/000000",
                "scale-1.000/000000/00/R1",
                "scale-1.000/000000/01/000000",
                "scale-1.000/000000/01/R1",
            ],
        )
        self.assertFalse((self.root / "03.spin").exists())
        np.testing.assert_array_equal(tasks[1]["moments"], [[2, 0, 0], [0, 0, 0]])

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
                "scale-1.000/000000/00/000000",
                "scale-1.000/000000/01/000000",
                "scale-1.000/000000/11/000000",
                "scale-1.000/000000/100/000000",
                "scale-1.000/000001/00/000000",
            ],
        )

    def test_no_placeholder_perturbations_or_path_escape(self):
        module = self.module()
        self.jdata["spin_pert_numb"] = 1
        with self.assertRaises(NotImplementedError):
            module.plan_spin_tasks(self.jdata, {})
        for name in ("../outside", "000000", "R/1"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                module.plan_spin_tasks(self.jdata, {}, perturb=lambda m, n: {name: m})
        self.assertFalse((self.root / "03.spin").exists())

    def test_canting_failure_identifies_snapshot_config_and_atom(self):
        from dpgen.data.spin_perturb import build_spin_perturbation

        module = self.module()
        count, provider = build_spin_perturbation(
            [{"Canting": {"angle": 30, "Rcut": 3.0}}]
        )
        self.jdata["spin_pert_numb"] = count

        with self.assertRaisesRegex(ValueError, "POSCAR.*C1.*atom 0.*Rcut"):
            module.plan_spin_tasks(self.jdata, {}, perturb=provider)

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

    def prepare_run_inputs(self):
        stage = self.root / "03.spin"
        task = "scale-1.000/000000/00/000000"
        path = stage / task
        path.mkdir(parents=True)
        (path / "POSCAR").write_text(POSCAR_TEXT)
        (path / "INCAR").write_text(TEMPLATE)
        (path / "POTCAR").write_bytes(self.potcar.read_bytes())
        (path / "KPOINTS").write_text(self.kpoints.read_text())
        (stage / "tasks.json").write_text(json.dumps({"tasks": [task]}))
        return path

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
        self.assertEqual(task.task_work_path, "scale-1.000/000000/00/000000")
        self.assertEqual(task.backward_files, ["OUTCAR", "OSZICAR", "vasprun.xml"])
        self.assertEqual(task.forward_files, ["POSCAR", "INCAR", "POTCAR", "KPOINTS"])

    def test_result_check_requires_both_outputs_and_normal_termination(self):
        module = self.module()
        path = self.prepare_run_inputs()
        with self.assertRaisesRegex(FileNotFoundError, "scale-1.000/000000/00/000000"):
            module.check_spin_results(self.jdata)
        (path / "OUTCAR").write_text("TOTAL-FORCE\nElapse\n")
        (path / "OSZICAR").write_text("")
        with self.assertRaises(RuntimeError):
            module.check_spin_results(self.jdata)
        (path / "OSZICAR").write_text(" 1 F= -1 E0= -1\n")
        self.assertEqual(module.check_spin_results(self.jdata), 1)

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

    def test_stage4_wires_canting_combinations_into_task_generation(self):
        module = self.module()
        self.jdata["spin_action"] = "make"
        self.jdata["pert_spin"] = [
            {
                "Canting": {
                    "angle": [30, 60],
                    "Rcut": [0.4, 0.5],
                    "direction": [[1, 0, 0], [0, 1, 0]],
                }
            }
        ]
        param = self.root / "input.json"
        param.write_text(json.dumps(self.jdata))
        with mock.patch.object(module, "make_spin_tasks") as make:
            spin_init.gen_spin_init(argparse.Namespace(PARAM=str(param), MACHINE=None))

        generated_jdata, _ = make.call_args.args
        provider = make.call_args.kwargs["perturb"]
        self.assertEqual(generated_jdata["spin_pert_numb"], 8)
        self.assertEqual(
            list(provider(np.array([[0, 0, 1.0]]), 8)),
            [f"C{index}" for index in range(1, 9)],
        )

    def test_invalid_spin_mode_fails_before_an_earlier_stage_runs(self):
        self.jdata["stages"] = [1, 4]
        self.jdata["spin_action"] = "make"
        self.jdata["md_incar"] = str(self.incar)
        self.jdata["pert_spin"] = [{"Rotation": {"angle": 30, "axis": [0, 0, 1]}}]
        param = self.root / "input.json"
        param.write_text(json.dumps(self.jdata))

        with mock.patch.object(spin_init, "make_spin_init_structures") as make:
            with self.assertRaisesRegex(NotImplementedError, "Rotation"):
                spin_init.gen_spin_init(
                    argparse.Namespace(PARAM=str(param), MACHINE=None)
                )
        make.assert_not_called()

    def test_run_action_still_rejects_an_unsupported_spin_mode(self):
        self.jdata["stages"] = [4]
        self.jdata["spin_action"] = "run"
        self.jdata["pert_spin"] = [{"Rotation": {"angle": 30, "axis": [0, 0, 1]}}]
        param = self.root / "input.json"
        machine = self.root / "machine.json"
        param.write_text(json.dumps(self.jdata))
        machine.write_text(json.dumps({}))

        with self.assertRaisesRegex(NotImplementedError, "Rotation"):
            spin_init.gen_spin_init(
                argparse.Namespace(PARAM=str(param), MACHINE=str(machine))
            )


if __name__ == "__main__":
    unittest.main()
