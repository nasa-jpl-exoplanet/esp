"""Create presentation-ready diagnostic figures from saved FlowAtmos runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "blue": "#2878B5",
    "orange": "#E07A2D",
    "green": "#39966B",
    "purple": "#8264A2",
    "gray": "#64748B",
    "light": "#E8EEF4",
    "red": "#C75450",
}


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Required evaluation report is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _read_history(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f"Required training history is missing: {path}")
    saved = np.genfromtxt(path, delimiter=",", names=True)
    if saved.size == 0:
        raise ValueError(f"Training history is empty: {path}")
    return {name: np.atleast_1d(saved[name]) for name in saved.dtype.names or ()}


def _style_axis(axis: plt.Axes) -> None:
    axis.set_facecolor("white")
    axis.grid(axis="y", color=COLORS["light"], linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.tick_params(colors="#334155")
    axis.title.set_color("#102A43")


def _plot_architecture(output: Path) -> None:
    fig, axis = plt.subplots(figsize=(13, 5.2))
    axis.set_xlim(0, 13)
    axis.set_ylim(0, 5.2)
    axis.axis("off")

    boxes = [
        (0.25, 2.45, 1.65, 1.2, "Spectrum\n329 channels", "#E8F1FA"),
        (2.25, 2.45, 2.15, 1.2, "Residual CNN\n64 → 128 → 256", "#E6F4EF"),
        (4.75, 2.45, 1.7, 1.2, "Attention\npooling", "#F5EFFA"),
        (6.8, 2.45, 1.85, 1.2, "Context vector\n256 values", "#FFF2E5"),
        (9.0, 2.45, 1.8, 1.2, "Conditional\nspline flow", "#E8F1FA"),
        (11.15, 2.45, 1.6, 1.2, "Posterior\nsamples", "#E6F4EF"),
    ]
    for x, y, width, height, label, color in boxes:
        axis.add_patch(
            plt.Rectangle(
                (x, y),
                width,
                height,
                facecolor=color,
                edgecolor="#486581",
                linewidth=1.2,
                joinstyle="round",
            )
        )
        axis.text(
            x + width / 2,
            y + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=10,
            color="#102A43",
            weight="medium",
        )
    for start, end in ((1.9, 2.25), (4.4, 4.75), (6.45, 6.8), (8.65, 9.0), (10.8, 11.15)):
        axis.annotate(
            "",
            xy=(end, 3.05),
            xytext=(start, 3.05),
            arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#486581"},
        )

    axis.text(
        3.3,
        4.25,
        "16 positional channels + normalized spectrum",
        ha="center",
        fontsize=9,
        color=COLORS["gray"],
    )
    axis.text(
        7.72,
        1.75,
        "noise embedding + stellar / orbital context + original spectrum scale",
        ha="center",
        fontsize=9,
        color=COLORS["gray"],
    )
    axis.annotate(
        "",
        xy=(7.72, 2.4),
        xytext=(7.72, 2.0),
        arrowprops={"arrowstyle": "->", "lw": 1.2, "color": COLORS["gray"]},
    )
    axis.text(
        11.95,
        1.3,
        "Posterior medians → point ranking\nPosterior ranks → probability matrix",
        ha="center",
        va="center",
        fontsize=9,
        color="#102A43",
    )
    axis.annotate(
        "",
        xy=(11.95, 2.4),
        xytext=(11.95, 1.85),
        arrowprops={"arrowstyle": "->", "lw": 1.2, "color": COLORS["gray"]},
    )
    fig.suptitle(
        "FlowAtmos: amortized posterior inference for one spectrum",
        fontsize=16,
        weight="bold",
        color="#102A43",
        y=0.96,
    )
    fig.text(
        0.5,
        0.06,
        "The trained network is reused for each new spectrum; it does not run a fresh MCMC chain.",
        ha="center",
        color=COLORS["gray"],
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(output / "01_flowatmos_architecture.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_training(root: Path, output: Path) -> None:
    runs = (
        ("Uniform-log baseline", root / "flowatmos", COLORS["blue"]),
        ("Scale-context model", root / "flowatmos_context", COLORS["green"]),
        ("Mixture v1", root / "flowatmos_mixture_v1", COLORS["orange"]),
        ("Mixture v2", root / "flowatmos_mixture_v2", COLORS["purple"]),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False)
    for axis, (label, run_dir, color) in zip(axes.flat, runs):
        history = _read_history(run_dir / "training_history.csv")
        axis.plot(
            history["epoch"],
            history["train_nll"],
            color=color,
            linewidth=1.8,
            label="Training",
        )
        axis.plot(
            history["epoch"],
            history["validation_nll"],
            color="#334155",
            linewidth=1.8,
            label="Validation",
        )
        best_epoch = int(np.argmin(history["validation_nll"]))
        axis.scatter(
            history["epoch"][best_epoch],
            history["validation_nll"][best_epoch],
            color=COLORS["red"],
            zorder=3,
            label="Selected checkpoint",
        )
        _style_axis(axis)
        axis.set_title(label)
        axis.set_xlabel("Epoch")
        axis.set_ylabel("Negative log likelihood (lower is better)")
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Training diagnostics: validation loss chooses the checkpoint",
        fontsize=15,
        weight="bold",
        color="#102A43",
    )
    fig.text(
        0.5,
        0.01,
        "A widening train/validation gap is a warning sign for overfitting. "
        "Mixture-v2 held-out metrics are not included in this training-curve comparison.",
        ha="center",
        color=COLORS["gray"],
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(output / "02_training_diagnostics.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_ranking(root: Path, output: Path) -> None:
    reports = {
        "Uniform-log model": _read_json(
            root / "flowatmos" / "evaluation_metrics.json"
        ),
        "Scale-context model": _read_json(
            root / "flowatmos_context" / "evaluation_metrics.json"
        ),
        "Uniform model on mixture-v1 test": _read_json(
            root
            / "flowatmos_mixture_v1"
            / "uniform_context_comparison"
            / "evaluation_metrics.json"
        ),
        "Mixture-v1 model": _read_json(
            root / "flowatmos_mixture_v1" / "evaluation_metrics.json"
        ),
    }
    cohorts = (
        (
            "Original uniform-log test set",
            ("Uniform-log model", "Scale-context model"),
        ),
        (
            "Mixture-v1 test set",
            ("Uniform model on mixture-v1 test", "Mixture-v1 model"),
        ),
    )
    metric_names = (
        ("global_top1_accuracy", "Global top-1 accuracy"),
        ("global_top5_overlap", "Global top-5 overlap"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    colors = (COLORS["blue"], COLORS["orange"])
    for axis, (cohort, model_names) in zip(axes, cohorts):
        x = np.arange(len(metric_names))
        width = 0.34
        for offset, model_name, color in zip(
            (-width / 2, width / 2), model_names, colors
        ):
            metrics = reports[model_name]["ranking_quality_over_all_spectra"]
            values = [metrics[key] for key, _ in metric_names]
            bars = axis.bar(x + offset, values, width, label=model_name, color=color)
            axis.bar_label(bars, fmt="%.2f", padding=3, fontsize=9)
        _style_axis(axis)
        axis.set_title(cohort)
        axis.set_xticks(x, [label for _, label in metric_names])
        axis.set_ylim(0, 0.85)
        axis.set_ylabel("Fraction")
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Detection-aware ranking on held-out synthetic spectra",
        fontsize=15,
        weight="bold",
        color="#102A43",
    )
    fig.text(
        0.5,
        0.015,
        "Left: identical original test split. Right: identical mixture-v1 test split. "
        "Undetected true leaders count as misses; mixture-v2 is omitted from this comparison.",
        ha="center",
        color=COLORS["gray"],
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.92))
    fig.savefig(output / "03_ranking_comparison.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_calibration(root: Path, output: Path) -> None:
    report_specs = (
        ("Uniform-log model", root / "flowatmos"),
        ("Scale-context model", root / "flowatmos_context"),
        (
            "Mixture-v1 model",
            root / "flowatmos_mixture_v1",
        ),
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    colors = (COLORS["blue"], COLORS["green"], COLORS["orange"])
    for axis, coverage_key, target in zip(axes, ("68%", "95%"), (0.68, 0.95)):
        axis.axhline(target, color="#1F2937", linestyle="--", linewidth=1.2, label="Ideal")
        for index, ((label, run_dir), color) in enumerate(zip(report_specs, colors)):
            report = _read_json(run_dir / "evaluation_metrics.json")
            value = report["credible_interval_coverage"][coverage_key]["mean_coverage"]
            axis.scatter([index], [value], s=75, color=color, label=f"{label}: {value:.3f}")
            axis.annotate(
                f"{value:.3f}",
                (index, value),
                xytext=(0, 8),
                textcoords="offset points",
                ha="center",
                fontsize=8,
            )
        _style_axis(axis)
        axis.set_title(f"Nominal {coverage_key} intervals")
        axis.set_ylim(0.45, 1.0)
        axis.set_ylabel("Observed coverage on synthetic test set")
        axis.set_xticks(np.arange(len(report_specs)), [label for label, _ in report_specs])
        axis.tick_params(axis="x", labelrotation=18)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Posterior interval calibration (closer to the dashed line is better)",
        fontsize=14,
        weight="bold",
        color="#102A43",
    )
    fig.text(
        0.5,
        0.015,
        "Mixture-v1 uses a different test distribution. Mixture-v2 is omitted from this figure.",
        ha="center",
        color=COLORS["gray"],
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))
    fig.savefig(output / "04_uncertainty_calibration.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_scenario_counts(root: Path, output: Path) -> None:
    datasets = (
        ("Uniform-log baseline", root / "flowatmos"),
        ("Mixture v1", root / "flowatmos_mixture_v1"),
        ("Mixture v2", root / "flowatmos_mixture_v2"),
    )
    counts_by_run: list[np.ndarray] = []
    for _, run_dir in datasets:
        with np.load(run_dir / "synthetic_dataset.npz") as saved:
            if "active_species_count" in saved:
                counts_by_run.append(saved["active_species_count"])
            else:
                counts_by_run.append(
                    (saved["targets"][:, :9] > -8.0).sum(axis=1)
                )
    fig, axis = plt.subplots(figsize=(10, 5))

    for (label, _), counts, color in zip(
        datasets,
        counts_by_run,
        (COLORS["blue"], COLORS["orange"], COLORS["purple"]),
    ):
        fractions = np.asarray(
            [np.mean(counts == value) for value in range(10)], dtype=np.float64
        )
        axis.plot(
            np.arange(10),
            fractions,
            marker="o",
            linewidth=2,
            color=color,
            label=label,
        )
    _style_axis(axis)
    axis.set_xticks(np.arange(10))
    axis.set_xlabel("Number of molecules above the -8 log10 VMR training floor")
    axis.set_ylabel("Fraction of generated simulations")
    axis.set_title("Training-prior design changes which ranking scenarios are learned")
    axis.legend(frameon=False)
    fig.text(
        0.5,
        0.015,
        "This describes generated labels, not observational detectability. "
        "Mixture-v1 is the failed lower-concentration experiment; v2 results are not plotted here.",
        ha="center",
        color=COLORS["gray"],
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(output / "05_scenario_composition.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_inference_latency(output: Path) -> None:
    labels = ["64 draws\nCPU", "64 draws\nMPS", "256 draws\nMPS", "Load checkpoint"]
    latency_ms = np.asarray([55, 52, 82, 218])
    colors = [COLORS["blue"], COLORS["green"], COLORS["green"], COLORS["gray"]]
    fig, axis = plt.subplots(figsize=(8, 4.8))
    bars = axis.bar(labels, latency_ms, color=colors, width=0.62)
    axis.bar_label(bars, fmt="%.0f ms", padding=3, fontsize=9)
    _style_axis(axis)
    axis.set_ylabel("Measured time (milliseconds)")
    axis.set_title("Example amortized inference latency by backend")
    fig.text(
        0.5,
        0.015,
        "Example warm timings; results vary by hardware and software. "
        "Checkpoint loading is a one-time cost. No MCMC baseline was benchmarked.",
        ha="center",
        color=COLORS["gray"],
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(output / "06_inference_latency.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def create_professor_figures(
    project_root: str | Path,
    output_dir: str | Path | None = None,
) -> Path:
    root = Path(project_root).expanduser().resolve() / "datasets"
    destination = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else root / "flowatmos_professor_figures"
    )
    destination.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    _plot_architecture(destination)
    _plot_training(root, destination)
    _plot_ranking(root, destination)
    _plot_calibration(root, destination)
    _plot_scenario_counts(root, destination)
    _plot_inference_latency(destination)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build professor-ready FlowAtmos diagnostic figures."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--output-dir", type=Path)
    arguments = parser.parse_args()
    output = create_professor_figures(arguments.project_root, arguments.output_dir)
    print(f"Saved presentation figures to {output}")
    for figure in sorted(output.glob("*.png")):
        print(f"  {figure.name}")


if __name__ == "__main__":
    main()
