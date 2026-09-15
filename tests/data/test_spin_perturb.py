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

    def test_extreme_finite_axis_normalizes_without_overflow(self):
        actual = self.module().rotate_moments(
            [[1.0, 0.0, 0.0]], 120, [1.7e308, 1.7e308, 1.7e308]
        )

        np.testing.assert_allclose(actual, [[0.0, 1.0, 0.0]], atol=1e-14)

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


class TestRandom(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb

        return spin_perturb

    def test_uniform_sphere_formula_preserves_each_magnitude_and_zero(self):
        rng = SequenceRng([0.0, 0.0, 0.0, np.pi / 2])

        actual = self.module().randomize_moments(
            [[2, 0, 0], [0, 0, 0], [0, 0, 3]], rng=rng
        )

        np.testing.assert_allclose(
            actual,
            [[2, 0, 0], [0, 0, 0], [0, 3, 0]],
            atol=1e-14,
        )

    def test_random_num_creates_seeded_sequence_without_reseeding(self):
        parameters = [{"Random": {"num": 3, "seed": 789}}]
        count_a, provider_a = self.module().build_spin_perturbation(parameters)
        count_b, provider_b = self.module().build_spin_perturbation(parameters)
        moments = np.array([[0.0, 0.0, 2.0], [0.0, 0.0, 0.0]])

        first_a = provider_a(moments, count_a)
        next_a = provider_a(moments, count_a)
        first_b = provider_b(moments, count_b)

        self.assertEqual(list(first_a), ["Rand1", "Rand2", "Rand3"])
        for name in first_a:
            np.testing.assert_array_equal(first_a[name], first_b[name])
            np.testing.assert_allclose(np.linalg.norm(first_a[name][0]), 2.0)
            self.assertFalse(np.array_equal(first_a[name], next_a[name]))

    def test_random_requires_positive_integer_num_and_valid_seed(self):
        invalid = [
            {"num": 0},
            {"num": -1},
            {"num": 1.5},
            {"num": True},
            {"num": None},
            {"num": 1, "seed": None},
            {"num": 1, "extra": 1},
        ]
        for parameters in invalid:
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                self.module().build_spin_perturbation([{"Random": parameters}])


class TestScale(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb

        return spin_perturb

    def test_grid_excludes_zero_and_scale_is_relative(self):
        module = self.module()

        np.testing.assert_allclose(
            module._scale_deltas(0.25, 0.05),
            [
                -0.25,
                -0.20,
                -0.15,
                -0.10,
                -0.05,
                0.05,
                0.10,
                0.15,
                0.20,
                0.25,
            ],
            rtol=0,
            atol=1e-15,
        )
        np.testing.assert_allclose(
            module.scale_moments([[2, -4, 0]], -0.25), [[1.5, -3, 0]]
        )

    def test_scale_provider_names_follow_delta_order(self):
        count, provider = self.module().build_spin_perturbation(
            [{"Scale": {"pert": 0.10, "pert_step": 0.05}}]
        )

        configurations = provider(np.array([[0.0, 0.0, 2.0]]), count)

        self.assertEqual(count, 4)
        self.assertEqual(list(configurations), ["S1", "S2", "S3", "S4"])
        np.testing.assert_allclose(
            [configurations[name][0, 2] for name in configurations],
            [1.8, 1.9, 2.1, 2.2],
        )

    def test_scale_rejects_invalid_range_step_ratio_and_fields(self):
        invalid = [
            {"pert": 0, "pert_step": 0.1},
            {"pert": 1, "pert_step": 0.1},
            {"pert": -0.1, "pert_step": 0.1},
            {"pert": 0.25, "pert_step": 0},
            {"pert": 0.25, "pert_step": 0.06},
            {"pert": 0.25},
            {"pert_step": 0.05},
            {"pert": 0.25, "pert_step": 0.05, "seed": 1},
        ]
        for parameters in invalid:
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                self.module().build_spin_perturbation([{"Scale": parameters}])


class TestRotaCant(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb

        return spin_perturb

    def test_rotation_is_applied_before_canting(self):
        actual = self.module().rota_cant_moments(
            [[1.0, 0.0, 0.0]],
            rotation_angle=90,
            axis=[0, 0, 1],
            canting_angle=90,
            rng=SequenceRng([0.0]),
        )

        np.testing.assert_allclose(actual, [[0.0, 0.0, 1.0]], atol=1e-14)

    def test_rota_cant_cartesian_product_and_seed_are_stable(self):
        parameters = [
            {
                "Rota_Cant": {
                    "R_angle": [30, 60],
                    "axis": [[0, 1, 0], [0, 0, 1]],
                    "C_angle": [45, 90],
                    "seed": 12345,
                }
            }
        ]
        count_a, provider_a = self.module().build_spin_perturbation(parameters)
        count_b, provider_b = self.module().build_spin_perturbation(parameters)
        moments = np.array([[0.0, 0.0, 2.0]])

        configurations_a = provider_a(moments, count_a)
        configurations_b = provider_b(moments, count_b)

        self.assertEqual(count_a, 8)
        self.assertEqual(
            list(configurations_a), [f"RC{index}" for index in range(1, 9)]
        )
        for name in configurations_a:
            np.testing.assert_array_equal(
                configurations_a[name], configurations_b[name]
            )
            np.testing.assert_allclose(np.linalg.norm(configurations_a[name][0]), 2.0)

    def test_rota_cant_rejects_missing_invalid_and_unknown_fields(self):
        invalid = [
            {"axis": [0, 0, 1], "C_angle": 30},
            {"R_angle": 30, "C_angle": 30},
            {"R_angle": 30, "axis": [0, 0, 1]},
            {"R_angle": 361, "axis": [0, 0, 1], "C_angle": 30},
            {"R_angle": 30, "axis": [0, 0, 0], "C_angle": 30},
            {"R_angle": 30, "axis": [0, 0, 1], "C_angle": 181},
            {
                "R_angle": 30,
                "axis": [0, 0, 1],
                "C_angle": 30,
                "seed": None,
            },
            {
                "R_angle": 30,
                "axis": [0, 0, 1],
                "C_angle": 30,
                "extra": 1,
            },
        ]
        for parameters in invalid:
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                self.module().build_spin_perturbation([{"Rota_Cant": parameters}])


class TestOrderedPipeline(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb

        return spin_perturb

    def test_operations_expand_branches_and_emit_only_final_leaves(self):
        count, provider = self.module().build_spin_perturbation(
            [
                {"Rotation": {"angle": [90, 180], "axis": [0, 0, 1]}},
                {"Scale": {"pert": 0.1, "pert_step": 0.1}},
            ]
        )

        configurations = provider(np.array([[1.0, 0.0, 0.0]]), count)

        self.assertEqual(count, 4)
        self.assertEqual(
            list(configurations),
            ["R1-S1", "R1-S2", "R2-S1", "R2-S2"],
        )
        self.assertNotIn("R1", configurations)
        np.testing.assert_allclose(configurations["R1-S1"], [[0, 0.9, 0]], atol=1e-14)
        np.testing.assert_allclose(configurations["R1-S2"], [[0, 1.1, 0]], atol=1e-14)

    def test_repeated_modes_compose_in_list_order(self):
        count, provider = self.module().build_spin_perturbation(
            [
                {"Rotation": {"angle": 90, "axis": [0, 0, 1]}},
                {"Canting": {"angle": 0, "seed": 1}},
                {"Rotation": {"angle": 90, "axis": [0, 0, 1]}},
            ]
        )

        configurations = provider(np.array([[1.0, 0.0, 0.0]]), count)

        self.assertEqual(count, 1)
        self.assertEqual(list(configurations), ["R1-C1-R1"])
        np.testing.assert_allclose(configurations["R1-C1-R1"], [[-1, 0, 0]], atol=1e-14)

    def test_each_item_has_exactly_one_known_operation(self):
        invalid = [
            {},
            {
                "Rotation": {"angle": 30, "axis": [0, 0, 1]},
                "Scale": {"pert": 0.1, "pert_step": 0.1},
            },
            {"rotation": {"angle": 30, "axis": [0, 0, 1]}},
            "Canting",
        ]
        for block in invalid:
            with (
                self.subTest(block=block),
                self.assertRaises((ValueError, NotImplementedError)),
            ):
                self.module().build_spin_perturbation([block])

    def test_errors_include_operation_index_mode_and_field(self):
        with self.assertRaisesRegex(ValueError, r"pert_spin\[0\]\.Rotation\.axis"):
            self.module().build_spin_perturbation(
                [{"Rotation": {"angle": 30, "axis": [0, 0, 0]}}]
            )

    def test_stochastic_state_is_reproducible_across_branches_and_snapshots(self):
        parameters = [
            {"Rotation": {"angle": [0, 90], "axis": [0, 0, 1]}},
            {"Canting": {"angle": [30, 60], "seed": 2468}},
            {"Random": {"num": 2, "seed": 1357}},
        ]
        count_a, provider_a = self.module().build_spin_perturbation(parameters)
        count_b, provider_b = self.module().build_spin_perturbation(parameters)
        moments = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])

        first_a = provider_a(moments, count_a)
        next_a = provider_a(moments, count_a)
        first_b = provider_b(moments, count_b)
        next_b = provider_b(moments, count_b)

        self.assertEqual(count_a, 8)
        self.assertEqual(list(first_a), list(first_b))
        for name in first_a:
            np.testing.assert_array_equal(first_a[name], first_b[name])
            np.testing.assert_array_equal(next_a[name], next_b[name])
        self.assertTrue(
            any(not np.array_equal(first_a[name], next_a[name]) for name in first_a)
        )


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

    def test_multiple_canting_blocks_compose_in_pipeline_order(self):
        count, provider = self.module().build_spin_perturbation(
            [
                {"Canting": {"angle": 10, "seed": 1}},
                {"Canting": {"angle": [20, 30], "seed": 2}},
            ]
        )

        configurations = provider(np.array([[0, 0, 1.0]]), count)

        self.assertEqual(count, 2)
        self.assertEqual(list(configurations), ["C1-C1", "C1-C2"])

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
