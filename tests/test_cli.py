import argparse
import subprocess as sp
import unittest

from dpgen.main import main_parser


class TestCLI(unittest.TestCase):
    def test_spin_init(self):
        parser = main_parser()
        subparsers = next(
            action
            for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        )
        self.assertIn("spin_init", subparsers.choices)

        args = parser.parse_args(["spin_init", "param.json"])
        self.assertEqual(args.func.__name__, "gen_spin_init")
        self.assertEqual(args.PARAM, "param.json")
        self.assertIsNone(args.MACHINE)

    def test_cli(self):
        sp.check_output(["dpgen", "-h"])
        for subcommand in (
            "run",
            "simplify",
            "init_surf",
            "init_bulk",
            "spin_init",
            "init_reaction",
            "autotest",
        ):
            sp.check_output(["dpgen", subcommand, "-h"])
