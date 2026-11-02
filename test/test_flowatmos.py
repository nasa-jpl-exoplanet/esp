import unittest
from pathlib import Path

import numpy as np
import torch

from tools.flowatmos import (
    CONDITION_NAMES,
    MOLECULES,
    TARGET_NAMES,
    FlowAtmos,
    ForwardInputs,
    UnconditionalPriorFlow,
    _bh_mask,
    _global_ranking_metrics,
    _normalised_arrays,
    _prepare_inputs,
    _rank_molecule_samples,
    _sample_mixture_abundances,
    _sample_scenario,
    _stratified_ranking_metrics,
    _target_bounds,
    _uniform_prior_nll,
)


class FlowAtmosTests(unittest.TestCase):
    def test_spline_flow_and_parameter_transform(self):
        lower = np.asarray([-12.0] * 9 + [300.0, 0.1], dtype=np.float32)
        upper = np.asarray([-2.0] * 9 + [800.0, 0.2], dtype=np.float32)
        model = FlowAtmos(329, 6, lower, upper)
        spectra = torch.randn(2, 329)
        noise = torch.full((2,), -1.5)
        conditions = torch.randn(2, 6)
        targets = torch.tensor(
            np.tile(np.asarray([-7.0] * 9 + [500.0, 0.15], dtype=np.float32), (2, 1))
        )

        log_prob = model.log_prob(spectra, noise, conditions, targets)
        latent, _ = model.to_latent(targets)
        samples = model.flow(model.encode(spectra, noise, conditions)).sample((5,))

        self.assertEqual(tuple(log_prob.shape), (2,))
        self.assertTrue(torch.isfinite(log_prob).all())
        self.assertEqual(tuple(samples.shape), (5, 2, len(TARGET_NAMES)))
        torch.testing.assert_close(model.from_latent(latent), targets)

    def test_unconditional_prior_flow_scores_targets(self):
        lower = np.asarray([-12.0] * 9 + [300.0, 0.1], dtype=np.float32)
        upper = np.asarray([-2.0] * 9 + [800.0, 0.2], dtype=np.float32)
        model = UnconditionalPriorFlow(lower, upper)
        targets = torch.tensor(
            np.tile(np.asarray([-7.0] * 9 + [500.0, 0.15], dtype=np.float32), (3, 1))
        )

        log_prob = model.log_prob(targets)

        self.assertEqual(tuple(log_prob.shape), (3,))
        self.assertTrue(torch.isfinite(log_prob).all())

    def test_detection_mask_uses_bh_cutoff(self):
        selected = _bh_mask(np.asarray([0.001, 0.02, 0.8, 0.9]), q=0.05)
        np.testing.assert_array_equal(selected, [True, True, False, False])

    def test_global_ranking_counts_undetected_true_top_as_a_miss(self):
        truths = np.asarray(
            [
                [9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
                [9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
            ]
        )
        medians = truths.copy()
        selected = np.zeros((2, 9), dtype=bool)
        selected[0, 1] = True
        selected[1, 0] = True
        point_order = np.full((2, 9), -1, dtype=np.int16)
        point_order[:, 0] = [1, 0]

        metrics = _global_ranking_metrics(truths, medians, selected, point_order)

        self.assertEqual(metrics["global_top1_accuracy"], 0.5)
        self.assertEqual(metrics["global_top1_detection_rate"], 0.5)
        self.assertEqual(metrics["global_top1_accuracy_when_detected"], 1.0)
        self.assertEqual(metrics["top1_accuracy_among_detected_species"], 1.0)
        self.assertEqual(metrics["global_top5_overlap"], 0.2)

    def test_fast_rank_summary_has_valid_rank_probabilities(self):
        samples = np.tile(np.arange(len(MOLECULES), dtype=np.float32), (32, 1))
        medians, rank_probabilities, point_order = _rank_molecule_samples(samples)

        self.assertEqual(medians.shape, (len(MOLECULES),))
        self.assertEqual(rank_probabilities.shape, (len(MOLECULES), len(MOLECULES)))
        np.testing.assert_allclose(rank_probabilities.sum(axis=1), 1.0)
        self.assertEqual(MOLECULES[point_order[0]], MOLECULES[-1])

    def test_stratified_ranking_reports_each_scenario_group(self):
        truths = np.asarray([[9.0, 8.0, 7.0], [1.0, 2.0, 3.0]])
        medians = truths.copy()
        selected = np.ones((2, 3), dtype=bool)
        point_order = np.asarray([[0, 1, 2], [2, 1, 0]])
        counts = np.asarray([1, 3])

        result = _stratified_ranking_metrics(
            truths, medians, selected, point_order, counts
        )

        self.assertEqual(set(result), {"1", "3"})
        self.assertEqual(result["1"]["global_top1_accuracy"], 1.0)
        self.assertEqual(result["3"]["global_top1_accuracy"], 1.0)

    def test_uniform_prior_nll_is_positive_for_wide_priors(self):
        score = _uniform_prior_nll(np.asarray([-2.0, 0.0]), np.asarray([8.0, 4.0]))
        self.assertAlmostEqual(score, np.log(40.0))

    def test_spectrum_scale_is_preserved_in_conditioning_context(self):
        spectra = np.asarray(
            [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [5.0, 5.0, 5.0]],
            dtype=np.float32,
        )
        dataset = {
            "spectra": spectra,
            "conditions": np.zeros((3, 4), dtype=np.float32),
            "noise_token": np.zeros(3, dtype=np.float32),
        }

        normalized, _, conditions, preprocessing = _prepare_inputs(
            dataset, np.asarray([0, 1])
        )
        _, _, evaluation_conditions = _normalised_arrays(dataset, preprocessing)

        self.assertEqual(conditions.shape, (3, 6))
        np.testing.assert_allclose(normalized[:2].mean(axis=1), 0.0, atol=1e-6)
        np.testing.assert_allclose(normalized[:2].std(axis=1), 1.0, atol=1e-6)
        self.assertLess(conditions[0, -1], conditions[1, -1])
        np.testing.assert_allclose(conditions, evaluation_conditions)
        _, _, legacy_conditions = _normalised_arrays(
            dataset,
            {
                "condition_mean": np.zeros(4, dtype=np.float32),
                "condition_std": np.ones(4, dtype=np.float32),
            },
        )
        self.assertEqual(legacy_conditions.shape, (3, 4))

    def test_scenario_samples_stay_within_prior_bounds(self):
        priors = {
            "R*": 0.607,
            "R*_uperr": 0.014,
            "R*_lowerr": -0.014,
            "T*": 4026.0,
            "T*_uperr": 14.0,
            "T*_lowerr": -14.0,
            "b": {
                "rp": 0.1525,
                "rp_uperr": 0.007,
                "rp_lowerr": -0.007,
                "logg": 3.15,
                "logg_uperr": 0.14,
                "logg_lowerr": -0.21,
                "sma": 0.0943,
            },
        }
        inputs = ForwardInputs(
            priors=priors,
            wavelengths=np.arange(329, dtype=np.float64),
            relative_noise=np.full(329, 0.05),
            cross_sections={},
            qtgrid={},
            haze_library={},
            pickle_dir=Path("."),
            planet="b",
        )
        lower, upper = _target_bounds(priors, "b")
        targets, conditions = _sample_scenario(
            np.random.default_rng(42), inputs, lower, upper
        )
        mixture_targets, _ = _sample_scenario(
            np.random.default_rng(43),
            inputs,
            lower,
            upper,
            scenario_strategy="mixture_balanced_v2",
        )

        self.assertEqual(len(MOLECULES), 9)
        self.assertEqual(len(CONDITION_NAMES), 4)
        self.assertTrue(np.all(targets >= lower))
        self.assertTrue(np.all(targets <= upper))
        self.assertTrue(np.all(np.isfinite(conditions)))
        self.assertLessEqual(np.power(10.0, targets[:9]).sum(), 0.1)
        mixture_total = np.power(10.0, mixture_targets[:9]).sum()
        self.assertGreater(mixture_total, 1e-4)
        self.assertLessEqual(mixture_total, 0.1)

    def test_mixture_sampler_covers_sparse_and_crowded_compositions(self):
        lower = np.asarray([-12.0] * 9 + [300.0, 0.1], dtype=np.float32)
        upper = np.asarray([-2.0] * 9 + [800.0, 0.2], dtype=np.float32)
        rng = np.random.default_rng(42)
        samples = np.stack(
            [_sample_mixture_abundances(rng, lower, upper) for _ in range(5000)]
        )
        stronger_samples = np.stack(
            [
                _sample_mixture_abundances(
                    rng, lower, upper, total_log_vmr_bounds=(-3.5, -2.0)
                )
                for _ in range(1000)
            ]
        )
        active_counts = (samples > -8.0).sum(axis=1)
        total_vmr = np.power(10.0, samples.astype(np.float64)).sum(axis=1)
        stronger_total_vmr = np.power(
            10.0, stronger_samples.astype(np.float64)
        ).sum(axis=1)

        self.assertTrue(np.all(samples >= lower[:9]))
        self.assertTrue(np.all(samples <= upper[:9]))
        self.assertTrue(np.all(total_vmr <= 0.1))
        self.assertGreaterEqual(active_counts.min(), 1)
        self.assertLessEqual(active_counts.max(), 9)
        self.assertGreater(np.unique(active_counts).size, 5)
        self.assertTrue(np.all(np.isfinite(samples)))
        self.assertTrue(np.all(stronger_total_vmr <= 0.1))
        self.assertGreater(
            np.median(stronger_total_vmr), 5.0 * np.median(total_vmr)
        )


if __name__ == "__main__":
    unittest.main()
