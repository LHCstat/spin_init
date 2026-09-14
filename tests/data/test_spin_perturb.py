"""Mathematical contracts for magnetic-moment perturbations."""

import unittest

import numpy as np


class TestCanting(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb

        return spin_perturb

    def test_matches_documented_z_axis_example(self):
        actual = self.module().cant_moments(
            [[0.0, 0.0, 1.0]], angle=30, rcut=0.5, direction=[1, 0, 0]
        )

        np.testing.assert_allclose(
            actual,
            [[np.sqrt(3) / 4, 0.25, np.sqrt(3) / 2]],
            atol=1e-14,
        )

    def test_zero_moments_stay_zero_and_all_nonzero_moments_are_canted(self):
        moments = np.array([[0, 0, 0], [0, 0, 1], [2, 0, 0]], dtype=float)

        actual = self.module().cant_moments(
            moments, angle=45, rcut=0.25, direction=[0, 1, 0]
        )

        np.testing.assert_array_equal(actual[0], [0, 0, 0])
        self.assertFalse(np.allclose(actual[1], moments[1]))
        self.assertFalse(np.allclose(actual[2], moments[2]))
        np.testing.assert_allclose(
            np.linalg.norm(actual, axis=1), np.linalg.norm(moments, axis=1)
        )

    def test_tiny_but_nonzero_moment_is_not_treated_as_zero(self):
        moments = np.array([[0.0, 0.0, 1e-13]])

        actual = self.module().cant_moments(
            moments, angle=30, rcut=0.5e-13, direction=[1, 0, 0]
        )

        self.assertFalse(np.array_equal(actual, moments))
        np.testing.assert_allclose(
            np.linalg.norm(actual, axis=1),
            np.linalg.norm(moments, axis=1),
            rtol=1e-14,
            atol=0,
        )

    def test_extreme_finite_moment_scales_do_not_underflow_or_overflow(self):
        normalized_expected = np.array([np.sqrt(3) / 4, 0.25, np.sqrt(3) / 2])
        for magnitude in (1e-200, 1e200):
            with self.subTest(magnitude=magnitude):
                actual = self.module().cant_moments(
                    [[0, 0, magnitude]],
                    angle=30,
                    rcut=0.5 * magnitude,
                    direction=[1, 0, 0],
                )
                np.testing.assert_allclose(
                    actual[0] / magnitude,
                    normalized_expected,
                    rtol=1e-14,
                    atol=0,
                )

    def test_missing_parallel_and_zero_directions_use_same_deterministic_axis(self):
        module = self.module()
        moments = [[0, 0, 1]]

        missing = module.cant_moments(moments, angle=30, rcut=0.5)
        parallel = module.cant_moments(moments, angle=30, rcut=0.5, direction=[0, 0, 2])
        zero = module.cant_moments(moments, angle=30, rcut=0.5, direction=[0, 0, 0])

        np.testing.assert_allclose(missing, parallel)
        np.testing.assert_allclose(missing, zero)
        np.testing.assert_allclose(missing, [[np.sqrt(3) / 4, 0.25, np.sqrt(3) / 2]])

    def test_extreme_direction_scales_preserve_direction_orientation(self):
        module = self.module()
        expected = module.cant_moments(
            [[0, 0, 1]], angle=30, rcut=0.5, direction=[0, 1, 0]
        )
        for scale in (1e-200, 1e200):
            with self.subTest(scale=scale):
                actual = module.cant_moments(
                    [[0, 0, 1]], angle=30, rcut=0.5, direction=[0, scale, 0]
                )
                np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0)

    def test_nearly_parallel_nonzero_projection_is_not_replaced_by_fallback(self):
        actual = self.module().cant_moments(
            [[0, 0, 1]], angle=0, rcut=0.5, direction=[0, 1e-13, 1]
        )

        np.testing.assert_allclose(actual, [[0, 0.5, np.sqrt(3) / 2]])

    def test_scalar_and_list_parameters_expand_as_a_cartesian_product(self):
        count, provider = self.module().build_spin_perturbation(
            [
                {
                    "Canting": {
                        "angle": [30, 60],
                        "Rcut": [0.4, 0.5],
                        "direction": [[1, 0, 0], [0, 1, 0]],
                    }
                }
            ]
        )

        configurations = provider(np.array([[0, 0, 1.0]]), count)

        self.assertEqual(count, 8)
        self.assertEqual(list(configurations), [f"C{index}" for index in range(1, 9)])
        np.testing.assert_allclose(
            configurations["C1"],
            [[0.4 * np.sqrt(3) / 2, 0.2, np.sqrt(1 - 0.4**2)]],
        )
        np.testing.assert_allclose(
            configurations["C8"],
            [[-0.5 * np.sqrt(3) / 2, 0.25, np.sqrt(1 - 0.5**2)]],
        )

    def test_multiple_canting_blocks_continue_configuration_numbering(self):
        count, provider = self.module().build_spin_perturbation(
            [
                {"Canting": {"angle": 10, "Rcut": 0.1}},
                {"Canting": {"angle": [20, 30], "Rcut": 0.2}},
            ]
        )

        configurations = provider(np.array([[0, 0, 1.0]]), count)

        self.assertEqual(list(configurations), ["C1", "C2", "C3"])

    def test_invalid_parameters_fail_before_generating_configurations(self):
        module = self.module()
        invalid = [
            {"Canting": {"Rcut": 0.1}},
            {"Canting": {"angle": 30}},
            {"Canting": {"angle": [], "Rcut": 0.1}},
            {"Canting": {"angle": 30, "Rcut": -0.1}},
            {"Canting": {"angle": True, "Rcut": 0.1}},
            {"Canting": {"angle": 30, "Rcut": False}},
            {"Canting": {"angle": 30, "Rcut": 0.1, "direction": [1, 0]}},
            {"Canting": {"angle": 30, "Rcut": 0.1, "direction": [True, 0, 0]}},
            {"Canting": {"angle": 30, "Rcut": 0.1, "unknown": 1}},
            {"Rotation": {"angle": 30, "axis": [0, 0, 1]}},
        ]
        for block in invalid:
            with (
                self.subTest(block=block),
                self.assertRaises((ValueError, NotImplementedError)),
            ):
                module.build_spin_perturbation([block])

    def test_rcut_larger_than_a_nonzero_moment_reports_atom_and_config(self):
        count, provider = self.module().build_spin_perturbation(
            [{"Canting": {"angle": 30, "Rcut": 1.1}}]
        )

        with self.assertRaisesRegex(ValueError, "C1.*atom 0.*Rcut"):
            provider(np.array([[0, 0, 1.0], [0, 0, 0]]), count)

    def test_rcut_above_magnitude_is_not_silently_clamped(self):
        with self.assertRaisesRegex(ValueError, "atom 0.*Rcut"):
            self.module().cant_moments([[0, 0, 1.0]], angle=30, rcut=1.0 + 5e-13)


if __name__ == "__main__":
    unittest.main()
