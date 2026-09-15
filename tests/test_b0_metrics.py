import math
import unittest

from analysis.b0_metrics import (
    distribution_summary,
    evaluate_regression_policy,
    stage_efficiencies,
    wilson_efficiency,
)


class B0MetricsTest(unittest.TestCase):
    def test_wilson_efficiency_boundaries(self):
        eff = wilson_efficiency(10, 10)
        self.assertEqual(eff.value, 1.0)
        self.assertGreaterEqual(eff.low, 0.0)
        self.assertLessEqual(eff.high, 1.0)

    def test_empty_efficiency_is_nan(self):
        eff = wilson_efficiency(0, 0)
        self.assertTrue(math.isnan(eff.value))

    def test_distribution_ignores_nonfinite(self):
        summary = distribution_summary([1.0, 2.0, float("nan"), float("inf")])
        self.assertEqual(summary["n"], 2)
        self.assertAlmostEqual(summary["mean"], 1.5)

    def test_stage_efficiencies_use_eligible_denominator(self):
        result = stage_efficiencies(
            eligible=[True, True, False, True],
            seeded=[True, False, True, True],
            unfiltered=[True, False, True, False],
            filtered=[True, False, False, False],
            truth_matched=[True, False, False, False],
        )
        self.assertEqual(result["eligible"]["events"], 3)
        self.assertAlmostEqual(result["seeded_over_eligible"]["value"], 2 / 3)
        self.assertAlmostEqual(result["unfiltered_over_seeded"]["value"], 1 / 2)

    def test_regression_policy(self):
        baseline = {"metrics": {"eff": 0.90, "sigma": 0.10}}
        candidate = {"metrics": {"eff": 0.88, "sigma": 0.105}}
        checks = [
            {"path": "metrics.eff", "mode": "max_drop", "tolerance": 0.01},
            {"path": "metrics.sigma", "mode": "max_fractional_increase", "tolerance": 0.10},
        ]
        failures = evaluate_regression_policy(baseline, candidate, checks)
        self.assertEqual(len(failures), 1)
        self.assertEqual(failures[0]["path"], "metrics.eff")


if __name__ == "__main__":
    unittest.main()
