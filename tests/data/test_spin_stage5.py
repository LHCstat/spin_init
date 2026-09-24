"""Stage 5 contracts independent of VASP and Linux symlink privileges."""

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from dpgen.data import spin_init, spin_stage5
from tests.data.test_spin_data import INCAR_TEXT, _magnetization_outcar
from tests.data.test_spin_init import POSCAR_TEXT
from tests.data.test_spin_init import _skip_without_symlink_privilege
from tests.data.test_spin_raw_conversion import EXTXYZ_TWO_FRAMES


class TestSpinStage5(unittest.TestCase):
    def _make_tasks(self, root):
        paths = [
            "scale-1.000/000000/00/origin/000000",
            "scale-1.000/000000/00/000-rotation/R1",
            "scale-1.020/000000/00/origin/000000",
        ]
        stage = root / "03.spin"
        stage.mkdir()
        (stage / "tasks.json").write_text(json.dumps({"tasks": paths}))
        for task, magnitude in zip(paths, (1.003, 1.02, 0.997)):
            folder = stage / task
            folder.mkdir(parents=True)
            (folder / "POSCAR").write_text(POSCAR_TEXT)
            (folder / "INCAR").write_text(INCAR_TEXT)
            (folder / "OUTCAR").write_text(
                _magnetization_outcar((0.0, 0.0, magnitude))
            )
            (folder / "OSZICAR").write_text("magnetic force data\n")
        return stage, paths

    def test_all_completed_tasks_are_numbered_per_scale_without_rmse(self):
        """Every completed task must enter conversion regardless of moment drift."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, paths = self._make_tasks(root)
            with mock.patch.object(spin_stage5, "check_spin_results", return_value=3):
                selection = spin_stage5.select_spin_tasks({"out_dir": str(root)})

        self.assertEqual(
            [(row["scale"], row["index"], row["task"]) for row in selection["selected"]],
            [
                ("scale-1.000", 1, paths[0]),
                ("scale-1.000", 2, paths[1]),
                ("scale-1.020", 1, paths[2]),
            ],
        )
        self.assertEqual(selection["rejected"], [])
        self.assertNotIn("rmse_limit", selection)
        self.assertTrue(all("rmse" not in row for row in selection["selected"]))

    def test_converter_input_uses_get_outcar_flat_pair_names(self):
        """The converter needs data/OUTCAR-N and data/OSZICAR-N pairs."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage, _ = self._make_tasks(root)
            with mock.patch.object(spin_stage5, "check_spin_results", return_value=3):
                selection = spin_stage5.select_spin_tasks({"out_dir": str(root)})
            scratch = root / "scratch"
            scales = spin_stage5.prepare_converter_inputs(stage, scratch, selection)

            self.assertEqual(scales, ["scale-1.000", "scale-1.020"])
            for scale in scales:
                folder = scratch / scale / "data"
                count = 2 if scale == "scale-1.000" else 1
                self.assertEqual(
                    sorted(path.name for path in folder.iterdir()),
                    sorted(
                        f"{name}-{index}"
                        for index in range(1, count + 1)
                        for name in ("OUTCAR", "OSZICAR")
                    ),
                )
                for index in range(1, count + 1):
                    self.assertTrue((folder / f"OUTCAR-{index}").is_file())
                    self.assertTrue((folder / f"OSZICAR-{index}").is_file())
            self.assertNotIn("1.02", (scratch / "scale-1.000/data/OUTCAR-1").read_text())

    @mock.patch("dpgen.data.spin_stage5.make_submission")
    def test_convert_data_machine_entry_submits_one_task_per_scale(self, make_submission):
        """The machine command runs inside each scale and retrieves extxyz."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mdata = {
                "api_version": "1.0",
                "convert-data": [
                    {
                        "command": "source activate py39 && nequip-data -m -z 8",
                        "machine": {"batch_type": "shell", "local_root": "./"},
                        "resources": {"group_size": 1},
                    }
                ],
            }
            spin_stage5.run_convert_data(root, ["scale-1.000", "scale-1.020"], mdata)
            self.assertTrue((root / "scale-1.000/out").is_dir())
            self.assertTrue((root / "scale-1.020/out").is_dir())

        kwargs = make_submission.call_args.kwargs
        self.assertEqual(kwargs["run_tasks"], ["scale-1.000", "scale-1.020"])
        self.assertEqual(kwargs["forward_files"], ["data"])
        self.assertEqual(kwargs["backward_files"], ["out/data.extxyz"])
        self.assertEqual(
            kwargs["commands"],
            [
                "mkdir -p out && source activate py39 && "
                "nequip-data -m -z 8 -p data -o out/data.extxyz"
            ],
        )
        make_submission.return_value.run_submission.assert_called_once_with()

    @mock.patch("dpgen.data.spin_stage5.make_submission")
    def test_convert_data_does_not_duplicate_explicit_paths(self, make_submission):
        """Catches duplicate options when a machine already follows the contract."""
        with tempfile.TemporaryDirectory() as temporary:
            command = (
                "source activate py39 && nequip-data -m -z 8 "
                "-p data -o out/data.extxyz"
            )
            spin_stage5.run_convert_data(
                temporary,
                ["scale-1.000"],
                {
                    "convert-data": [
                        {
                            "command": command,
                            "machine": {"batch_type": "shell", "local_root": "./"},
                            "resources": {"group_size": 1},
                        }
                    ]
                },
            )

        actual = make_submission.call_args.kwargs["commands"][0]
        self.assertEqual(
            actual,
            "mkdir -p out && " + command,
        )
        self.assertEqual(actual.count("-p data"), 1)
        self.assertEqual(actual.count("-o out/data.extxyz"), 1)

    @mock.patch("dpgen.data.spin_stage5.make_submission")
    def test_convert_data_rejects_ambiguous_shell_commands(self, make_submission):
        """Catches path flags being attached to a command other than nequip-data."""
        machine = {
            "convert-data": [
                {
                    "command": "unused",
                    "machine": {"batch_type": "shell", "local_root": "./"},
                    "resources": {"group_size": 1},
                }
            ]
        }
        for command in (
            "echo nequip-data",
            "nequip-data -m -z 8 && echo done",
            "nequip-data -m -z 8 | tee conversion.log",
            "nequip-data -m -z 8\necho done",
            "nequip-data -m -z 8 # paths added too late",
        ):
            with self.subTest(command=command):
                machine["convert-data"][0]["command"] = command
                with self.assertRaisesRegex(ValueError, "nequip-data.*terminal"):
                    spin_stage5.run_convert_data(
                        "unused", ["scale-1.000"], machine
                    )
        make_submission.assert_not_called()

    def test_bohrium_convert_data_reports_missing_optional_oss_dependency(self):
        """Lebesgue upload must fail clearly before an opaque oss2 NameError."""
        machine = {
            "convert-data": [
                {
                    "command": "nequip-data -m -z 8",
                    "machine": {
                        "batch_type": "Lebesgue",
                        "context_type": "LebesgueContext",
                        "local_root": "./",
                    },
                    "resources": {"group_size": 1},
                }
            ]
        }
        with mock.patch.dict(sys.modules, {"oss2": None}):
            for context in ("LebesgueContext", "BohriumContext"):
                machine["convert-data"][0]["machine"]["context_type"] = context
                with self.subTest(context=context):
                    with self.assertRaisesRegex(RuntimeError, "dpdispatcher.*bohrium"):
                        spin_stage5.run_convert_data(
                            "unused", ["scale-1.000"], machine
                        )

    def test_invalid_convert_data_config_reports_required_fields(self):
        """A malformed MACHINE must fail clearly before contacting dpdispatcher."""
        cases = [
            ({"convert-data": []}, "convert-data.*nonempty"),
            (
                {"convert-data": {"machine": {}, "resources": {}}},
                "convert-data.*command",
            ),
            (
                {"convert-data": {"command": "nequip-data", "resources": {}}},
                "convert-data.*machine",
            ),
            (
                {"convert-data": {"command": "nequip-data", "machine": {}}},
                "convert-data.*resources",
            ),
        ]
        with mock.patch.object(spin_stage5, "make_submission") as submit:
            for machine, message in cases:
                with self.subTest(machine=machine):
                    with self.assertRaisesRegex(ValueError, message):
                        spin_stage5.run_convert_data(
                            "unused", ["scale-1.000"], machine
                        )
            submit.assert_not_called()

    def test_extxyz_raw_and_npy_share_pdf_out_directory(self):
        """The final scale/out folder retains extxyz and raw intermediates."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            out = root / "scale-1.000/out"
            out.mkdir(parents=True)
            (out / "data.extxyz").write_text(EXTXYZ_TWO_FRAMES)
            counts = spin_stage5.convert_scale_outputs(root, ["scale-1.000"])

            self.assertEqual(counts, {"scale-1.000": 2})
            self.assertTrue((out / "data.extxyz").is_file())
            self.assertTrue((out / "spin.raw").is_file())
            self.assertTrue((out / "force_mag.raw").is_file())
            energy = np.load(out / "set.000/energy.npy")
            self.assertEqual(energy.shape, (2,))
            np.testing.assert_allclose(energy, [-1, -2])
            self.assertFalse((out / "set").exists())

    def test_final_data_links_match_pdf_and_get_outcar_layout(self):
        """The published data dirs and numbered inputs must be real relative links."""
        _skip_without_symlink_privilege()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage, _ = self._make_tasks(root)
            with mock.patch.object(spin_stage5, "check_spin_results", return_value=3):
                selection = spin_stage5.select_spin_tasks({"out_dir": str(root)})
            scratch = root / "scratch"
            spin_stage5.prepare_converter_inputs(stage, scratch, selection)
            final = root / "04.data"
            spin_stage5.publish_scale_links(stage, scratch, final, selection)
            scratch.replace(final)

            for scale in ("scale-1.000", "scale-1.020"):
                source = stage / scale / "data"
                linked = final / scale / "data"
                self.assertTrue(linked.is_symlink())
                self.assertFalse(os.path.isabs(os.readlink(linked)))
                self.assertEqual(linked.resolve(), source.resolve())
                count = 2 if scale == "scale-1.000" else 1
                for name in (
                    f"{kind}-{index}"
                    for index in range(1, count + 1)
                    for kind in ("OUTCAR", "OSZICAR")
                ):
                    self.assertTrue((source / name).is_symlink())
                    self.assertFalse(os.path.isabs(os.readlink(source / name)))
                    self.assertTrue((linked / name).is_file())

    def test_existing_output_stage_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            final = root / "04.data"
            final.mkdir()
            (final / "keep.txt").write_text("user data")
            with self.assertRaisesRegex(RuntimeError, "already exists"):
                spin_stage5.collect_spin_data({"out_dir": str(root)}, {})
            self.assertEqual((final / "keep.txt").read_text(), "user data")

    def test_cli_rejects_existing_stage_before_replacing_saved_param(self):
        """A refused Stage 5 rerun must not change earlier run provenance."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "run_spin"
            (output / "04.data").mkdir(parents=True)
            saved = output / "param.json"
            saved.write_text("original run")
            param = root / "new_param.json"
            param.write_text(
                json.dumps(
                    {
                        "stages": [5],
                        "out_dir": str(output),
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

            with self.assertRaisesRegex(RuntimeError, "already exists"):
                spin_init.gen_spin_init(
                    __import__("argparse").Namespace(PARAM=str(param), MACHINE=None)
                )

            self.assertEqual(saved.read_text(), "original run")

    def test_stage_five_converts_returned_extxyz_and_publishes_raw_npy(self):
        """The CLI collector must use the returned extxyz, not legacy DPData."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._make_tasks(root)

            def fake_converter(work_path, scales, _mdata):
                for scale in scales:
                    out = Path(work_path) / scale / "out"
                    out.mkdir()
                    (out / "data.extxyz").write_text(EXTXYZ_TWO_FRAMES)

            with mock.patch.object(spin_stage5, "check_spin_results", return_value=3):
                with mock.patch.object(
                    spin_stage5, "run_convert_data", side_effect=fake_converter
                ):
                    # Windows cannot create real symlinks. Link publication has
                    # its own privilege-gated test above; keep the conversion real.
                    with mock.patch.object(spin_stage5, "publish_scale_links"):
                        count = spin_stage5.collect_spin_data(
                            {"out_dir": str(root)}, {"convert-data": [{}]}
                        )

            self.assertEqual(count, 4)
            for scale in ("scale-1.000", "scale-1.020"):
                out = root / "04.data" / scale / "out"
                self.assertTrue((out / "data.extxyz").is_file())
                self.assertTrue((out / "spin.raw").is_file())
                energy = np.load(out / "set.000/energy.npy")
                self.assertEqual(energy.shape, (2,))
                np.testing.assert_allclose(energy, [-1, -2])
            selection = json.loads((root / "04.data/selection.json").read_text())
            self.assertEqual(len(selection["selected"]), 3)
            self.assertEqual(selection["rejected"], [])
            self.assertNotIn("rmse_limit", selection)

    def test_missing_converter_output_does_not_publish_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stage, _ = self._make_tasks(root)
            with mock.patch.object(spin_stage5, "check_spin_results", return_value=3):
                with mock.patch.object(spin_stage5, "run_convert_data"):
                    with self.assertRaisesRegex(FileNotFoundError, "data.extxyz"):
                        spin_stage5.collect_spin_data(
                            {"out_dir": str(root)}, {"convert-data": [{}]}
                        )
            self.assertFalse((root / "04.data").exists())
            self.assertFalse((stage / "scale-1.000/data").exists())

    def test_stage_five_cli_passes_convert_data_machine_config(self):
        """The Stage 5 entrypoint must not invoke the legacy DPData collector."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            param = root / "param.json"
            param.write_text(
                json.dumps(
                    {
                        "stages": [5],
                        "out_dir": str(root / "output"),
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
            machine = root / "machine.json"
            machine.write_text(
                json.dumps(
                    {
                        "convert-data": [
                            {
                                "command": "nequip-data -m -z 8",
                                "machine": {"local_root": "./"},
                                "resources": {"group_size": 1},
                            }
                        ]
                    }
                )
            )
            with mock.patch.object(spin_stage5, "collect_spin_data") as collector:
                spin_init.gen_spin_init(
                    __import__("argparse").Namespace(PARAM=str(param), MACHINE=str(machine))
                )
            self.assertEqual(collector.call_count, 1)
            self.assertEqual(
                collector.call_args.args[1]["convert-data"][0]["command"],
                "nequip-data -m -z 8",
            )


if __name__ == "__main__":
    unittest.main()
