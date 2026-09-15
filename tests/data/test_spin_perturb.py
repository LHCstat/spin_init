"""Mathematical contracts for magnetic-moment perturbations."""

import unittest

import numpy as np


class SequenceRng:
    def __init__(self, values):
        self.values = iter(values)

    def uniform(self, low, high):
        value = next(self.values)
        if not low <= value < high:
            raise AssertionError(f"test azimuth {value} is outside [{low}, {high})")
        return value


class TestRotation(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb

        return spin_perturb

    def test_right_hand_rotation_uses_normalized_global_axis(self):
        moments = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        actual = self.module().rotate_moments(moments, 90, [0, 0, 2])

        np.testing.assert_allclose(actual, [[0, 1, 0], [0, 0, 0]], atol=1e-14)

    def test_rotation_endpoints_preserve_expected_vectors(self):
        module = self.module()
        moments = np.array([[1.0, 2.0, 3.0]])

        np.testing.assert_array_equal(
            module.rotate_moments(moments, 0, [0, 0, 1]), moments
        )
        np.testing.assert_allclose(
            module.rotate_moments(moments, 360, [0, 0, 1]),
            moments,
            atol=1e-14,
        )
        np.testing.assert_allclose(
            module.rotate_moments([[1, 0, 0]], 180, [0, 0, 1]),
            [[-1, 0, 0]],
            atol=1e-14,
        )

    def test_rotation_cartesian_product_is_angle_outer_axis_inner(self):
        count, provider = self.module().build_spin_perturbation(
            [
                {
                    "Rotation": {
                        "angle": [90, 180],
                        "axis": [[0, 0, 1], [0, 1, 0]],
                    }
                }
            ]
        )

        configurations = provider(np.array([[1.0, 0.0, 0.0]]), count)

        self.assertEqual(count, 4)
        self.assertEqual(list(configurations), ["R1", "R2", "R3", "R4"])
        np.testing.assert_allclose(configurations["R1"], [[0, 1, 0]], atol=1e-14)
        np.testing.assert_allclose(configurations["R2"], [[0, 0, -1]], atol=1e-14)

    def test_rotation_rejects_zero_axis_and_unknown_parameters(self):
        invalid = [
            ({"Rotation": {"angle": 30, "axis": [0, 0, 0]}}, "axis"),
            (
                {
                    "Rotation": {
                        "angle": 30,
                        "axis": [0, 0, 1],
                        "seed": 1,
                    }
                },
                "seed",
            ),
        ]
        for block, message in invalid:
            with (
                self.subTest(block=block),
                self.assertRaisesRegex(ValueError, message),
            ):
                self.module().build_spin_perturbation([block])


class TestCanting(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb

        return spin_perturb

    def test_angle_is_between_initial_and_perturbed_moments(self):
        moments = np.array([[0, 0, 0], [0, 0, 2], [1, 2, 3]], dtype=float)

        actual = self.module().cant_moments(
            moments, angle=30, rng=np.random.default_rng(12345)
        )

        np.testing.assert_array_equal(actual[0], [0, 0, 0])
        for atom_index in (1, 2):
            initial = moments[atom_index]
            perturbed = actual[atom_index]
            np.testing.assert_allclose(
                np.linalg.norm(perturbed), np.linalg.norm(initial), rtol=1e-14
            )
            cosine = np.dot(initial, perturbed) / (
                np.linalg.norm(initial) * np.linalg.norm(perturbed)
            )
            np.testing.assert_allclose(cosine, np.cos(np.deg2rad(30)), atol=1e-14)

    def test_each_nonzero_atom_gets_an_independent_random_azimuth(self):
        moments = np.tile([0.0, 0.0, 1.0], (3, 1))

        actual = self.module().cant_moments(
            moments, angle=60, rng=SequenceRng([0.0, np.pi / 2, np.pi])
        )

        np.testing.assert_allclose(
            actual,
            [
                [np.sqrt(3) / 2, 0.0, 0.5],
                [0.0, np.sqrt(3) / 2, 0.5],
                [-np.sqrt(3) / 2, 0.0, 0.5],
            ],
            atol=1e-14,
        )

    def test_zero_and_180_degree_endpoints_are_exact_and_do_not_consume_rng(self):
        moments = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])

        unchanged = self.module().cant_moments(moments, angle=0, rng=SequenceRng([]))
        reversed_moments = self.module().cant_moments(
            moments, angle=180, rng=SequenceRng([])
        )

        np.testing.assert_array_equal(unchanged, moments)
        np.testing.assert_array_equal(reversed_moments[0], moments[0])
        np.testing.assert_array_equal(reversed_moments[1], -moments[1])

    def test_low_level_canting_requires_one_scalar_angle(self):
        with self.assertRaisesRegex(ValueError, "scalar"):
            self.module().cant_moments(
                [[0, 0, 1]], angle=[30], rng=np.random.default_rng(1)
            )

    def test_seed_reproduces_the_complete_configuration_sequence(self):
        parameters = [{"Canting": {"angle": [30, 60], "seed": 9876}}]
        first_count, first_provider = self.module().build_spin_perturbation(parameters)
        second_count, second_provider = self.module().build_spin_perturbation(
            parameters
        )
        moments = np.tile([0.0, 0.0, 1.0], (3, 1))

        first_snapshot = first_provider(moments, first_count)
        first_next_snapshot = first_provider(moments, first_count)
        second_snapshot = second_provider(moments, second_count)
        second_next_snapshot = second_provider(moments, second_count)

        self.assertEqual(first_count, 2)
        for name in ("C1", "C2"):
            np.testing.assert_array_equal(first_snapshot[name], second_snapshot[name])
            np.testing.assert_array_equal(
                first_next_snapshot[name], second_next_snapshot[name]
            )
            self.assertFalse(
                np.array_equal(first_snapshot[name], first_next_snapshot[name])
            )

    def test_angle_list_creates_one_configuration_per_angle(self):
        count, provider = self.module().build_spin_perturbation(
            [{"Canting": {"angle": [0, 30, 180], "seed": 123}}]
        )

        configurations = provider(np.array([[0, 0, 2.0]]), count)

        self.assertEqual(count, 3)
        self.assertEqual(list(configurations), ["C1", "C2", "C3"])
        expected_cosines = [1.0, np.sqrt(3) / 2, -1.0]
        for name, expected_cosine in zip(configurations, expected_cosines):
            perturbed = configurations[name][0]
            np.testing.assert_allclose(np.linalg.norm(perturbed), 2.0)
            np.testing.assert_allclose(perturbed[2] / 2.0, expected_cosine, atol=1e-14)

    def test_multiple_canting_blocks_continue_configuration_numbering(self):
        count, provider = self.module().build_spin_perturbation(
            [
                {"Canting": {"angle": 10, "seed": 1}},
                {"Canting": {"angle": [20, 30], "seed": 2}},
            ]
        )

        configurations = provider(np.array([[0, 0, 1.0]]), count)

        self.assertEqual(count, 3)
        self.assertEqual(list(configurations), ["C1", "C2", "C3"])

    def test_seed_is_optional(self):
        count, provider = self.module().build_spin_perturbation(
            [{"Canting": {"angle": 45}}]
        )

        configurations = provider(np.array([[0, 0, 1.0]]), count)

        self.assertEqual(list(configurations), ["C1"])
        np.testing.assert_allclose(np.linalg.norm(configurations["C1"][0]), 1.0)

    def test_extreme_finite_moment_scales_do_not_underflow_or_overflow(self):
        for magnitude in (1e-200, 1e200):
            with self.subTest(magnitude=magnitude):
                actual = self.module().cant_moments(
                    [[0, 0, magnitude]],
                    angle=30,
                    rng=np.random.default_rng(123),
                )[0]
                normalized = actual / magnitude
                np.testing.assert_allclose(np.linalg.norm(normalized), 1.0, rtol=1e-14)
                np.testing.assert_allclose(normalized[2], np.sqrt(3) / 2, rtol=1e-14)

    def test_invalid_parameters_are_rejected(self):
        module = self.module()
        invalid = [
            {"Canting": {}},
            {"Canting": {"angle": []}},
            {"Canting": {"angle": -0.1}},
            {"Canting": {"angle": 180.1}},
            {"Canting": {"angle": True}},
            {"Canting": {"angle": 30, "seed": -1}},
            {"Canting": {"angle": 30, "seed": 1.5}},
            {"Canting": {"angle": 30, "seed": True}},
            {"Canting": {"angle": 30, "seed": None}},
            {"Canting": {"angle": 30, "Rcut": 0.1}},
            {"Canting": {"angle": 30, "direction": [1, 0, 0]}},
        ]
        for block in invalid:
            with (
                self.subTest(block=block),
                self.assertRaises((ValueError, NotImplementedError)),
            ):
                module.build_spin_perturbation([block])


if __name__ == "__main__":
    unittest.main()
