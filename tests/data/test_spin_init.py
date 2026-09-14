import argparse
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from dpgen.data import arginfo as data_arginfo
from dpgen.data import spin_init
from dpgen.util import normalize


def _skip_without_symlink_privilege():
    """Skip tests that need a real symbolic link when the OS cannot make one.

    On Windows without SeCreateSymbolicLinkPrivilege ``os.symlink`` raises
    WinError 1314. That is an environment limitation, not a DP-GEN defect, so
    the calling test is reported as skipped rather than failing. The real
    symlink assertions below still run unchanged on symlink-capable systems
    (for example Linux). No copy fallback is introduced.
    """
    with tempfile.TemporaryDirectory() as temporary_directory:
        target = Path(temporary_directory) / "probe-target"
        target.write_text("x")
        link = Path(temporary_directory) / "probe-link"
        try:
            os.symlink(target, link)
        except OSError as error:
            if getattr(error, "winerror", None) == 1314:
                raise unittest.SkipTest(
                    "real symbolic link unsupported in this environment "
                    "(WinError 1314: privilege not held); run on a "
                    "symlink-capable system (e.g. Linux)"
                ) from error
            raise


POSCAR_TEXT = """Fe2
1.0
2.0 0.0 0.0
0.0 2.0 0.0
0.0 0.0 2.0
Fe
2
Direct
0.0 0.0 0.0
0.5 0.5 0.5
"""


class TestSpinInitArginfo(unittest.TestCase):
    def test_normalizes_minimal_parameters_without_output_suffix(self):
        self.assertTrue(hasattr(data_arginfo, "spin_init_jdata_arginfo"))

        parameters = {
            "stages": [1, 2, 3],
            "from_poscar_path": "POSCAR",
            "super_cell": [1, 1, 1],
            "scale": [1.0],
            "pert_numb": 1,
            "pert_box": 0.03,
            "pert_atom": 0.01,
            "md_incar": "INCAR.md",
            "md_nstep": 3,
            "potcars": ["POTCAR"],
        }

        normalized = normalize(data_arginfo.spin_init_jdata_arginfo(), parameters)

        self.assertEqual(normalized["out_dir"], ".")
        self.assertEqual(normalized["pert_spin"], [])
        self.assertNotIn("init_fp_style", normalized)
        self.assertNotIn("skip_relax", normalized)

    def test_normalizes_canting_parameter_blocks(self):
        parameters = {
            "stages": [4],
            "from_poscar_path": "POSCAR",
            "super_cell": [1, 1, 1],
            "scale": [1.0],
            "pert_numb": 0,
            "pert_box": 0.0,
            "pert_atom": 0.0,
            "md_incar": "INCAR.md",
            "md_nstep": 0,
            "potcars": ["POTCAR"],
            "spin_incar": "INCAR.spin",
            "pert_spin": [
                {
                    "Canting": {
                        "angle": [30, 60],
                        "Rcut": 0.4,
                    }
                }
            ],
        }

        normalized = normalize(data_arginfo.spin_init_jdata_arginfo(), parameters)

        self.assertEqual(normalized["pert_spin"][0]["Canting"]["Rcut"], 0.4)

    def test_normalizes_fp_machine_parameters(self):
        self.assertTrue(hasattr(data_arginfo, "spin_init_mdata_arginfo"))

        machine = {
            "fp": {
                "command": "vasp_std",
                "machine": {
                    "batch_type": "shell",
                    "context_type": "local",
                    "local_root": "./",
                },
                "resources": {"batch_type": "shell", "group_size": 1},
            }
        }

        normalized = normalize(data_arginfo.spin_init_mdata_arginfo(), machine)

        self.assertEqual(normalized["api_version"], "1.0")
        self.assertEqual(normalized["fp"]["command"], "vasp_std")

    def test_normalizes_optional_spin_machine_parameters(self):
        machine = {
            "fp": {
                "command": "vasp_std",
                "machine": {
                    "batch_type": "shell",
                    "context_type": "local",
                    "local_root": "./",
                },
                "resources": {"batch_type": "shell", "group_size": 1},
            },
            "spin": {
                "command": "vasp_ncl",
                "machine": {
                    "batch_type": "shell",
                    "context_type": "local",
                    "local_root": "./",
                },
                "resources": {"batch_type": "shell", "group_size": 1},
            },
        }

        normalized = normalize(data_arginfo.spin_init_mdata_arginfo(), machine)

        self.assertEqual(normalized["fp"]["command"], "vasp_std")
        self.assertEqual(normalized["spin"]["command"], "vasp_ncl")

    def test_rejects_invalid_structure_and_md_parameters(self):
        self.assertTrue(hasattr(spin_init, "validate_spin_init_parameters"))
        valid = {
            "stages": [1],
            "super_cell": [1, 1, 1],
            "scale": [1.0],
            "pert_numb": 0,
            "pert_box": 0.0,
            "pert_atom": 0.0,
            "md_nstep": 0,
            "potcars": ["POTCAR"],
        }
        invalid_values = {
            "stages": [],
            "super_cell": [1, 1],
            "scale": [],
            "pert_numb": -1,
            "pert_box": -0.1,
            "pert_atom": -0.1,
            "md_nstep": -1,
            "potcars": [],
        }

        for key, value in invalid_values.items():
            with self.subTest(key=key):
                parameters = dict(valid)
                parameters[key] = value
                with self.assertRaises(ValueError) as caught:
                    spin_init.validate_spin_init_parameters(parameters)
                self.assertIn(key, str(caught.exception))


class TestSpinInitStructures(unittest.TestCase):
    def test_generates_unperturbed_and_requested_perturbations_without_sys_layer(self):
        self.assertTrue(hasattr(spin_init, "make_spin_init_structures"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "POSCAR.source"
            source.write_text(POSCAR_TEXT)
            output = root / "output"
            parameters = {
                "from_poscar_path": str(source),
                "out_dir": str(output),
                "super_cell": [1, 1, 1],
                "scale": [1.0],
                "pert_numb": 2,
                "pert_box": 0.03,
                "pert_atom": 0.01,
            }

            tasks = spin_init.make_spin_init_structures(parameters)

            self.assertEqual(
                tasks,
                [
                    "scale-1.000/000000",
                    "scale-1.000/000001",
                    "scale-1.000/000002",
                ],
            )
            stage = output / "00.scale_pert"
            self.assertFalse(any(stage.glob("sys-*")))
            for task in tasks:
                structure_path = stage / task / "POSCAR"
                self.assertTrue(structure_path.is_file())

                from pymatgen.core import Structure

                self.assertEqual(len(Structure.from_file(structure_path)), 2)

    def test_rejects_existing_scale_perturbation_stage(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "POSCAR.source"
            source.write_text(POSCAR_TEXT)
            output = root / "output"
            parameters = {
                "from_poscar_path": str(source),
                "out_dir": str(output),
                "super_cell": [1, 1, 1],
                "scale": [1.0],
                "pert_numb": 0,
                "pert_box": 0.03,
                "pert_atom": 0.01,
            }
            stage = output / "00.scale_pert"
            stage.mkdir(parents=True)

            try:
                spin_init.make_spin_init_structures(parameters)
            except Exception as error:
                self.assertIs(type(error), RuntimeError)
                caught = error
            else:
                self.fail("existing stage was overwritten")
            self.assertIn(
                f"output stage already exists: {stage.resolve()}", str(caught)
            )

    def test_reports_missing_poscar_source_before_creating_stage(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            source = Path(temporary_directory) / "missing" / "POSCAR"
            parameters = {
                "from_poscar_path": str(source),
                "out_dir": str(output),
                "super_cell": [1, 1, 1],
                "scale": [1.0],
                "pert_numb": 0,
                "pert_box": 0.0,
                "pert_atom": 0.0,
            }

            with self.assertRaises(FileNotFoundError) as caught:
                spin_init.make_spin_init_structures(parameters)

            self.assertIn("POSCAR source does not exist", str(caught.exception))
            self.assertIn(str(source.resolve()), str(caught.exception))
            self.assertFalse((output / "00.scale_pert").exists())


class TestSpinInitMD(unittest.TestCase):
    def test_rejects_existing_symlink_path_without_overwriting_it(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            link = root / "POSCAR"
            link.write_text("keep me")
            target = root / "source"
            target.write_text(POSCAR_TEXT)

            with self.assertRaises(RuntimeError) as caught:
                spin_init._make_relative_symlink(link, target)

            self.assertIn(str(link), str(caught.exception))
            self.assertEqual(link.read_text(), "keep me")

    def test_creates_relative_poscar_incar_and_potcar_symlinks(self):
        self.assertTrue(hasattr(spin_init, "make_spin_init_md"))

        _skip_without_symlink_privilege()
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            structure_task = output / "00.scale_pert" / "scale-1.000" / "000000"
            structure_task.mkdir(parents=True)
            (structure_task / "POSCAR").write_text(POSCAR_TEXT)
            md_incar = output / "INCAR.source"
            md_incar.write_text("IBRION = 0\nNSW = 3\n")
            potcar = output / "POTCAR.source"
            potcar.write_text("POTCAR DATA\n")
            parameters = {
                "out_dir": str(output),
                "scale": [1.0],
                "pert_numb": 0,
                "md_incar": str(md_incar),
                "potcars": [str(potcar)],
            }

            tasks = spin_init.make_spin_init_md(parameters, {})

            self.assertEqual(tasks, ["scale-1.000/000000"])
            md_task = output / "01.md" / tasks[0]
            for name in ("POSCAR", "INCAR", "POTCAR"):
                self.assertTrue((md_task / name).is_symlink(), name)
                self.assertFalse(os.path.isabs(os.readlink(md_task / name)), name)
            self.assertEqual(
                os.readlink(md_task / "POSCAR"),
                os.path.join(
                    "..", "..", "..", "00.scale_pert", "scale-1.000", "000000", "POSCAR"
                ),
            )
            self.assertEqual(
                os.readlink(md_task / "INCAR"), os.path.join("..", "..", "INCAR")
            )
            self.assertEqual(
                os.readlink(md_task / "POTCAR"), os.path.join("..", "..", "POTCAR")
            )

    def test_missing_user_forward_file_does_not_leave_md_stage(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            structure_task = output / "00.scale_pert" / "scale-1.000" / "000000"
            structure_task.mkdir(parents=True)
            (structure_task / "POSCAR").write_text(POSCAR_TEXT)
            md_incar = output / "INCAR.source"
            md_incar.write_text("IBRION = 0\nNSW = 3\n")
            potcar = output / "POTCAR.source"
            potcar.write_text("POTCAR DATA\n")
            missing = output / "missing" / "vasp.slurm"
            parameters = {
                "out_dir": str(output),
                "scale": [1.0],
                "pert_numb": 0,
                "md_incar": str(md_incar),
                "potcars": [str(potcar)],
            }
            machine = {"fp_user_forward_files": [str(missing)]}

            with self.assertRaises(FileNotFoundError) as caught:
                spin_init.make_spin_init_md(parameters, machine)

            self.assertIn("user forward file does not exist", str(caught.exception))
            self.assertIn(str(missing.resolve()), str(caught.exception))
            self.assertFalse((output / "01.md").exists())


class TestSpinInitRunner(unittest.TestCase):
    def _make_md_task(self, root):
        task = root / "01.md" / "scale-1.000" / "000000"
        task.mkdir(parents=True)
        for name in ("POSCAR", "INCAR", "POTCAR"):
            (task / name).write_text(name)
        return task

    @mock.patch("dpgen.data.spin_init.make_submission", create=True)
    def test_always_backwards_xdatcar_and_does_not_infer_sbatch(self, make_submission):
        self.assertTrue(hasattr(spin_init, "run_spin_init_md"))

        submission = make_submission.return_value
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            self._make_md_task(output)
            slurm_file = output / "vasp.slurm"
            slurm_file.write_text("#!/bin/sh\nvasp_std\n")
            parameters = {"out_dir": str(output)}
            machine = {
                "fp_command": "vasp_std",
                "fp_group_size": 2,
                "fp_machine": {"batch_type": "shell"},
                "fp_resources": {},
                "fp_user_forward_files": [str(slurm_file)],
                "fp_user_backward_files": ["REPORT", "XDATCAR"],
            }

            tasks = spin_init.run_spin_init_md(parameters, machine)

        self.assertEqual(tasks, ["scale-1.000/000000"])
        kwargs = make_submission.call_args.kwargs
        self.assertEqual(kwargs["commands"], ["vasp_std"])
        self.assertEqual(
            kwargs["forward_files"], ["POSCAR", "INCAR", "POTCAR", "vasp.slurm"]
        )
        self.assertEqual(kwargs["backward_files"], ["OUTCAR", "XDATCAR", "REPORT"])
        submission.run_submission.assert_called_once_with()

    @mock.patch("dpgen.data.spin_init.make_submission", create=True)
    def test_runs_sbatch_only_when_fp_command_requests_it(self, make_submission):
        self.assertTrue(hasattr(spin_init, "run_spin_init_md"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            self._make_md_task(output)
            machine = {
                "fp_command": "sbatch vasp.slurm",
                "fp_group_size": 1,
                "fp_machine": {"batch_type": "shell"},
                "fp_resources": {},
            }

            spin_init.run_spin_init_md({"out_dir": str(output)}, machine)

        self.assertEqual(
            make_submission.call_args.kwargs["commands"], ["sbatch vasp.slurm"]
        )


class TestSpinInitOutputs(unittest.TestCase):
    fixture = Path(__file__).parent / "spin_init" / "XDATCAR"

    def test_parses_all_xdatcar_frames_without_outcar(self):
        self.assertTrue(hasattr(spin_init, "parse_xdatcar_snapshots"))

        structures = spin_init.parse_xdatcar_snapshots(self.fixture)

        self.assertEqual(len(structures), 3)
        self.assertEqual([site.specie.symbol for site in structures[0]], ["Fe", "O"])
        self.assertAlmostEqual(structures[0].lattice.a, 2.0)
        self.assertAlmostEqual(structures[2].frac_coords[0][0], 0.2)
        self.assertAlmostEqual(structures[2].frac_coords[1][1], 0.7)

    def test_rejects_missing_empty_and_malformed_xdatcar(self):
        self.assertTrue(hasattr(spin_init, "parse_xdatcar_snapshots"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            cases = {
                "missing": root / "missing" / "XDATCAR",
                "empty": root / "scale-1.000" / "000000" / "XDATCAR",
                "malformed": root / "scale-1.000" / "000001" / "XDATCAR",
            }
            cases["empty"].parent.mkdir(parents=True)
            cases["empty"].write_text("")
            cases["malformed"].parent.mkdir(parents=True)
            cases["malformed"].write_text("not an XDATCAR\n")

            for case, path in cases.items():
                with self.subTest(case=case):
                    with self.assertRaises((FileNotFoundError, RuntimeError)) as caught:
                        spin_init.parse_xdatcar_snapshots(path)
                    self.assertIn(str(path.resolve()), str(caught.exception))

    @mock.patch("pymatgen.io.vasp.outputs.Xdatcar")
    def test_rejects_inconsistent_frames_returned_by_pymatgen(self, xdatcar):
        from pymatgen.core import Lattice, Structure

        first = Structure(Lattice.cubic(2.0), ["Fe", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        second = Structure(Lattice.cubic(2.0), ["Fe"], [[0.1, 0, 0]])
        xdatcar.return_value.structures = [first, second]

        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "scale-1.000" / "000003" / "XDATCAR"
            path.parent.mkdir(parents=True)
            path.write_text("synthetic non-empty content\n")

            with self.assertRaises(RuntimeError) as caught:
                spin_init.parse_xdatcar_snapshots(path)

        self.assertIn("inconsistent frame 2", str(caught.exception))
        self.assertIn(str(path.resolve()), str(caught.exception))

    def test_checks_md_completion_without_requiring_xdatcar(self):
        self.assertTrue(hasattr(spin_init, "check_vasp_md_complete"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            task = Path(temporary_directory) / "scale-1.000" / "000000"
            task.mkdir(parents=True)
            (task / "OUTCAR").write_text(
                "TOTAL-FORCE\nTOTAL-FORCE\nTOTAL-FORCE\nElapsed time (sec)\n"
            )

            spin_init.check_vasp_md_complete(task, 3)

    def test_rejects_incomplete_md_independently_of_valid_xdatcar(self):
        self.assertTrue(hasattr(spin_init, "check_vasp_md_complete"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            task = Path(temporary_directory) / "scale-1.000" / "000002"
            task.mkdir(parents=True)
            (task / "OUTCAR").write_text("TOTAL-FORCE\nElapsed time (sec)\n")
            shutil.copy2(self.fixture, task / "XDATCAR")

            with self.assertRaises(RuntimeError) as caught:
                spin_init.check_vasp_md_complete(task, 3)
            self.assertIn(
                "scale-1.000/000002", str(caught.exception).replace("\\", "/")
            )
            self.assertIn(str((task / "OUTCAR").resolve()), str(caught.exception))


class TestSpinInitCollector(unittest.TestCase):
    fixture = Path(__file__).parent / "spin_init" / "XDATCAR"

    def _make_complete_task(self, output, task="scale-1.000/000000"):
        task_path = output / "01.md" / task
        task_path.mkdir(parents=True)
        (task_path / "OUTCAR").write_text(
            "TOTAL-FORCE\nTOTAL-FORCE\nTOTAL-FORCE\nElapsed time (sec)\n"
        )
        shutil.copy2(self.fixture, task_path / "XDATCAR")
        return task_path

    def test_writes_every_xdatcar_frame_as_an_independent_poscar(self):
        self.assertTrue(hasattr(spin_init, "collect_xdatcar_snapshots"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            self._make_complete_task(output)
            parameters = {
                "out_dir": str(output),
                "scale": [1.0],
                "pert_numb": 0,
                "md_nstep": 3,
            }

            frame_count = spin_init.collect_xdatcar_snapshots(parameters)

            self.assertEqual(frame_count, 3)
            frame_paths = [
                output
                / "02.disp"
                / "scale-1.000"
                / "000000"
                / f"{index:02d}"
                / "POSCAR"
                for index in range(3)
            ]
            self.assertTrue(all(path.is_file() for path in frame_paths))

            from pymatgen.io.vasp.inputs import Poscar

            structures = [Poscar.from_file(path).structure for path in frame_paths]
            self.assertEqual(
                [site.specie.symbol for site in structures[0]], ["Fe", "O"]
            )
            self.assertAlmostEqual(structures[0].lattice.a, 2.0)
            self.assertAlmostEqual(structures[2].frac_coords[0][0], 0.2)
            self.assertAlmostEqual(structures[2].frac_coords[1][1], 0.7)

    def test_preflight_failure_does_not_create_disp_stage(self):
        self.assertTrue(hasattr(spin_init, "collect_xdatcar_snapshots"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            task = self._make_complete_task(output)
            (task / "XDATCAR").unlink()
            parameters = {
                "out_dir": str(output),
                "scale": [1.0],
                "pert_numb": 0,
                "md_nstep": 3,
            }

            with self.assertRaises(FileNotFoundError) as caught:
                spin_init.collect_xdatcar_snapshots(parameters)

            self.assertIn("scale-1.000", str(caught.exception))
            self.assertIn("000000", str(caught.exception))
            self.assertFalse((output / "02.disp").exists())

    def test_rejects_existing_disp_stage(self):
        self.assertTrue(hasattr(spin_init, "collect_xdatcar_snapshots"))

        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "output"
            (output / "02.disp").mkdir(parents=True)
            parameters = {
                "out_dir": str(output),
                "scale": [1.0],
                "pert_numb": 0,
                "md_nstep": 3,
            }

            with self.assertRaises(RuntimeError) as caught:
                spin_init.collect_xdatcar_snapshots(parameters)

            self.assertIn(str((output / "02.disp").resolve()), str(caught.exception))


class TestSpinInitWorkflow(unittest.TestCase):
    @mock.patch("dpgen.data.spin_init.run_spin_init_md")
    @mock.patch("dpgen.data.spin_init.make_spin_init_md")
    def test_accepts_the_same_flat_machine_format_as_init_bulk(self, make_md, run_md):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            poscar = root / "POSCAR"
            poscar.write_text(POSCAR_TEXT)
            incar = root / "INCAR.md"
            incar.write_text("IBRION = 0\nNSW = 3\n")
            potcar = root / "POTCAR"
            potcar.write_text("POTCAR DATA\n")
            param_path = root / "param.json"
            param_path.write_text(
                json.dumps(
                    {
                        "stages": [2],
                        "from_poscar_path": str(poscar),
                        "out_dir": str(root / "output"),
                        "super_cell": [1, 1, 1],
                        "scale": [1.0],
                        "pert_numb": 0,
                        "pert_box": 0.0,
                        "pert_atom": 0.0,
                        "md_incar": str(incar),
                        "md_nstep": 3,
                        "potcars": [str(potcar)],
                    }
                )
            )
            flat_machine = {
                "api_version": "1.0",
                "fp_command": "vasp_std",
                "fp_machine": {
                    "batch_type": "shell",
                    "context_type": "local",
                    "local_root": "./",
                },
                "fp_resources": {"batch_type": "shell"},
                "fp_group_size": 1,
                "fp_user_forward_files": [],
                "fp_user_backward_files": [],
            }
            machine_path = root / "machine.json"
            machine_path.write_text(json.dumps(flat_machine))

            try:
                spin_init.gen_spin_init(
                    argparse.Namespace(PARAM=str(param_path), MACHINE=str(machine_path))
                )
            except Exception as error:
                self.fail(
                    "spin_init rejected a flat machine.json accepted by init_bulk: "
                    f"{error}"
                )

            make_md.assert_called_once()
            run_md.assert_called_once()
            self.assertEqual(make_md.call_args.args[1]["fp_command"], "vasp_std")

    @mock.patch("dpgen.data.spin_init.collect_xdatcar_snapshots")
    @mock.patch("dpgen.data.spin_init.run_spin_init_md")
    @mock.patch("dpgen.data.spin_init.make_spin_init_md")
    @mock.patch("dpgen.data.spin_init.make_spin_init_structures")
    def test_runs_three_stages_in_order_with_normalized_machine_data(
        self, make_structures, make_md, run_md, collect
    ):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            output = root / "output"
            poscar = root / "POSCAR"
            poscar.write_text(POSCAR_TEXT)
            incar = root / "INCAR.md"
            incar.write_text("IBRION = 0\nNSW = 3\n")
            potcar = root / "POTCAR"
            potcar.write_text("POTCAR DATA\n")
            param_path = root / "param.json"
            param_path.write_text(
                json.dumps(
                    {
                        "stages": [1, 2, 3],
                        "from_poscar_path": str(poscar),
                        "out_dir": str(output),
                        "super_cell": [1, 1, 1],
                        "scale": [1.0],
                        "pert_numb": 0,
                        "pert_box": 0.0,
                        "pert_atom": 0.0,
                        "md_incar": str(incar),
                        "md_nstep": 3,
                        "potcars": [str(potcar)],
                    }
                )
            )
            machine_path = root / "machine.json"
            machine_path.write_text(
                json.dumps(
                    {
                        "fp": {
                            "command": "vasp_std",
                            "machine": {
                                "batch_type": "shell",
                                "context_type": "local",
                                "local_root": "./",
                            },
                            "resources": {"batch_type": "shell", "group_size": 1},
                        }
                    }
                )
            )

            spin_init.gen_spin_init(
                argparse.Namespace(PARAM=str(param_path), MACHINE=str(machine_path))
            )

            self.assertTrue((output / "param.json").is_file())
            normalized_jdata = make_structures.call_args.args[0]
            normalized_mdata = make_md.call_args.args[1]
            self.assertEqual(normalized_jdata["out_dir"], str(output))
            self.assertEqual(normalized_mdata["fp_command"], "vasp_std")
            make_structures.assert_called_once_with(normalized_jdata)
            make_md.assert_called_once_with(normalized_jdata, normalized_mdata)
            run_md.assert_called_once_with(normalized_jdata, normalized_mdata)
            collect.assert_called_once_with(normalized_jdata)

    def test_rejects_unknown_stage(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            param_path = root / "param.json"
            param_path.write_text(
                json.dumps(
                    {
                        "stages": [5],
                        "from_poscar_path": "POSCAR",
                        "super_cell": [1, 1, 1],
                        "scale": [1.0],
                        "pert_numb": 0,
                        "pert_box": 0.0,
                        "pert_atom": 0.0,
                        "md_incar": "INCAR.md",
                        "md_nstep": 3,
                        "potcars": ["POTCAR"],
                    }
                )
            )

            with self.assertRaises(RuntimeError) as caught:
                spin_init.gen_spin_init(
                    argparse.Namespace(PARAM=str(param_path), MACHINE=None)
                )

            self.assertIn("unknown spin_init stage 5", str(caught.exception))

    @mock.patch("dpgen.data.spin_init.make_spin_init_structures")
    def test_uses_nsw_from_md_incar_as_upstream_init_bulk_does(self, make_structures):
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            poscar = root / "POSCAR"
            poscar.write_text(POSCAR_TEXT)
            incar = root / "INCAR.md"
            incar.write_text("IBRION = 0\nNSW = 3\n")
            potcar = root / "POTCAR"
            potcar.write_text("POTCAR DATA\n")
            param_path = root / "param.json"
            param_path.write_text(
                json.dumps(
                    {
                        "stages": [1],
                        "from_poscar_path": str(poscar),
                        "out_dir": str(root / "output"),
                        "super_cell": [1, 1, 1],
                        "scale": [1.0],
                        "pert_numb": 0,
                        "pert_box": 0.0,
                        "pert_atom": 0.0,
                        "md_incar": str(incar),
                        "md_nstep": 99,
                        "potcars": [str(potcar)],
                    }
                )
            )

            spin_init.gen_spin_init(
                argparse.Namespace(PARAM=str(param_path), MACHINE=None)
            )

            self.assertEqual(make_structures.call_args.args[0]["md_nstep"], 3)


if __name__ == "__main__":
    unittest.main()
