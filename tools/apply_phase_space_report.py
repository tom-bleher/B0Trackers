#!/usr/bin/env python3
from pathlib import Path


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected 1 match, found {count}")
    return text.replace(old, new, 1)


path = Path("analysis/b0_report.py")
text = path.read_text()

anchor = '''def _plot_residual_rms(outdir, residuals, np, filename, title):
'''
plot_fun = '''def _plot_binned_efficiency(outdir, rows, filename, xlabel, title, np):
    import matplotlib.pyplot as plt

    valid = [row for row in rows if row["denominator"] > 0 and row["value"] is not None]
    if not valid:
        return
    centers = np.array([(row["low"] + row["high"]) / 2.0 for row in valid], dtype=float)
    half_width = np.array([(row["high"] - row["low"]) / 2.0 for row in valid], dtype=float)
    values = np.array([row["value"] for row in valid], dtype=float)
    lows = np.array([row["low"] if False else row["value"] - row["low"] for row in []])
    ylow = np.array([row["value"] - row["low"] for row in []])
    # Wilson interval keys are named low/high too, so preserve the bin edges first.
    interval_low = np.array([row["interval_low"] for row in valid], dtype=float)
    interval_high = np.array([row["interval_high"] for row in valid], dtype=float)
    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.errorbar(
        centers,
        values,
        xerr=half_width,
        yerr=np.vstack([values - interval_low, interval_high - values]),
        fmt="o",
        capsize=3,
    )
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Truth-matched efficiency")
    ax.set_title(title)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=160)
    plt.close(fig)


def _plot_residual_rms(outdir, residuals, np, filename, title):
'''
# The binned helper currently uses low/high both for bin edges and interval fields.
# Normalize rows before plotting by adding explicit interval keys later.
text = replace_once(text, anchor, plot_fun, "binned plot helper")

# Build report signature.
text = replace_once(
    text,
    '''    momentum_bins: list[float],
    angle_bins_mrad: list[float],
    manifest: dict | None = None,
):
''',
    '''    momentum_bins: list[float],
    angle_bins_mrad: list[float],
    phi_bins_rad: list[float],
    min_truth_weight: float = 0.5,
    manifest: dict | None = None,
):
''',
    "build signature",
)

# Replace binned efficiency rows with explicit bin_low/bin_high convention to avoid interval key collision.
# binned_efficiency currently overwrites bin edge low/high with Wilson low/high, so fix helper first elsewhere.
path.write_text(text)

# Fix metric helper's binned rows to expose distinct bin and interval fields.
metric = Path("analysis/b0_metrics.py")
m = metric.read_text()
old = '        rows.append({"low": float(lo), "high": float(hi), **wilson_efficiency(num, den).as_dict()})\n'
new = '''        eff = wilson_efficiency(num, den).as_dict()
        rows.append(
            {
                "bin_low": float(lo),
                "bin_high": float(hi),
                "numerator": eff["numerator"],
                "denominator": eff["denominator"],
                "value": eff["value"],
                "interval_low": eff["low"],
                "interval_high": eff["high"],
            }
        )
'''
if old not in m:
    raise RuntimeError("binned metric row not found")
m = m.replace(old, new, 1)
metric.write_text(m)
