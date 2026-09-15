# Reproducible B0 validation workflow

A physics regression is meaningful only when the baseline and candidate use compatible inputs, geometry, material, and reconstruction settings. B0Trackers therefore treats provenance as part of the result rather than as external notebook context.

## 1. Create a run manifest

```bash
python3 analysis/make_run_manifest.py \
  --input sim.edm4hep.root \
  --dataset proton-gun-8-41GeV-4-22mrad \
  --material-map "$DETECTOR_PATH/calibrations/<validated-map>.cbor" \
  --eicrecon-repo ~/eic/dev/EICrecon \
  --epic-repo ~/eic/dev/epic \
  --events 10000 \
  --arg=-PB0TrackerStubSeeder:minStations=3 \
  --arg=-PB0TrackerCKFTracking:numMeasurementsMin=3
```

The manifest records:

- input path and SHA256;
- stable dataset/sample name;
- event count;
- selected-primary PDG/status;
- `DETECTOR_CONFIG` and `DETECTOR_PATH`;
- B0Trackers, EICrecon, and epic commit SHAs when the repositories are supplied;
- material-map path and SHA256;
- reconstruction arguments chosen for the study.

## 2. Produce the tracking report

```bash
python3 analysis/b0_report.py b0trackers.root \
  --manifest b0-run-manifest.json \
  --output-dir report
```

The report includes global and phase-space diagnostics:

- truth-matched efficiency vs truth momentum;
- truth-matched efficiency vs scattering angle;
- truth-matched efficiency vs azimuth `phi`;
- relative momentum bias/resolution vs momentum and scattering angle;
- selected-primary station/front-back crossing patterns;
- seed survival and ambiguity survival;
- mapping and residual diagnostics;
- association-defined fake and selected-primary duplicate metrics.

The fake/duplicate definition is explicit and configurable with `--min-truth-weight` (default `0.5`). Within truth-reachable events, a track whose dominant MC association is absent or below the threshold is an association-defined fake; more than one track above threshold associated to the selected primary makes that event a selected-primary duplicate.

## 3. Compare only compatible reports

```bash
python3 analysis/compare_b0_metrics.py \
  baseline/summary.json candidate/summary.json \
  --policy analysis/regression_policy.example.json
```

The comparator refuses incompatible schema/geometry/selection provenance by default. The policy may additionally require equality of manifest fields such as dataset ID, EICrecon/epic commit, and material-map SHA256. Use `--allow-incompatible` only for an intentional cross-configuration study, not for a regression gate.

## 4. Interpret binomial uncertainties

Efficiency bins use Wilson intervals rather than symmetric Gaussian errors. The JSON stores bin edges separately from interval bounds:

```text
bin_low / bin_high
value
interval_low / interval_high
numerator / denominator
```

This separation is important for small-occupancy edge bins and prevents plotting code from confusing physics-bin boundaries with statistical confidence bounds.
