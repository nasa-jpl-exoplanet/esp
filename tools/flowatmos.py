"""Synthetic CERBERUS data and conditional FlowAtmos training helpers."""

from __future__ import annotations

import csv
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import torch
from scipy.stats import chisquare, kendalltau, spearmanr
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from zuko.flows import NSF


MOLECULES = ("H2O", "CO2", "CH4", "CO", "HCN", "C2H2", "NH3", "SO2", "H2S")
CONDITION_NAMES = (
    "stellar_radius_rsun",
    "stellar_temperature_K",
    "planet_logg_cgs",
    "semi_major_axis_au",
)
MODEL_CONDITION_NAMES = (
    *CONDITION_NAMES,
    "spectrum_mean_depth",
    "spectrum_standard_deviation",
)
DETECTION_FLOOR_LOG_VMR = -8.0
TARGET_NAMES = (*MOLECULES, "temperature_K", "planet_radius_Rjup")


@dataclass(frozen=True)
class ForwardInputs:
    priors: dict[str, Any]
    wavelengths: np.ndarray
    relative_noise: np.ndarray
    cross_sections: dict[str, Any]
    qtgrid: dict[str, Any]
    haze_library: dict[str, Any]
    pickle_dir: Path
    planet: str


def forward_inputs_from_products(
    system: dict[str, Any],
    observed: dict[str, Any],
    xslib: Any,
    pickle_dir: str | Path,
    planet: str = "b",
) -> ForwardInputs:
    """Build the simulation inputs from the products already loaded in the notebook."""
    root = Path(pickle_dir).expanduser().resolve()
    haze_path = root / "Jup-ISS-aerosol.dat"
    if planet not in system["priors"] or planet not in observed["data"]:
        raise KeyError(f"Planet {planet!r} is missing from the local pickle products.")
    if planet not in xslib["data"]:
        raise KeyError(f"Planet {planet!r} is missing from the opacity library.")

    detector_data = [
        observed["data"][planet][detector]["1"]
        for detector in ("NRS1", "NRS2")
    ]
    wavelengths = np.concatenate([np.asarray(row["WB"], dtype=np.float64) for row in detector_data])
    observed_depth = np.concatenate([np.asarray(row["ES"], dtype=np.float64) for row in detector_data])
    observed_error = np.concatenate([np.asarray(row["ESerr"], dtype=np.float64) for row in detector_data])
    relative_noise = observed_error / np.maximum(np.abs(observed_depth), 1e-12)
    if not (
        wavelengths.ndim == 1
        and wavelengths.size == observed_depth.size == observed_error.size
        and np.all(np.isfinite(wavelengths))
        and np.all(np.diff(wavelengths) > 0)
        and np.all(np.isfinite(relative_noise))
        and np.all(relative_noise > 0)
    ):
        raise ValueError("The selected NRS1/NRS2 order-1 channels or errors are invalid.")

    xsec_data = xslib["data"][planet]
    missing = set(MOLECULES) - set(xsec_data["XSECS"])
    if missing:
        raise KeyError(f"The cached opacity library is missing species: {sorted(missing)}")

    import excalibur.cerberus.core as crbcore

    haze_library: dict[str, Any] = {"PROFILE": []}
    crbcore.hazelib(
        haze_library,
        hazedir=str(root),
        datafile=haze_path.name,
        verbose=False,
    )
    return ForwardInputs(
        priors=system["priors"],
        wavelengths=wavelengths,
        relative_noise=relative_noise,
        cross_sections=xsec_data["XSECS"],
        qtgrid=xsec_data["QTGRID"],
        haze_library=haze_library,
        pickle_dir=root,
        planet=planet,
    )


def _target_bounds(priors: dict[str, Any], planet: str) -> tuple[np.ndarray, np.ndarray]:
    star = priors
    planet_prior = priors[planet]
    solar_to_au = 0.004650467260962157
    eq_temperature = float(star["T*"]) * math.sqrt(
        float(star["R*"]) * solar_to_au / (2.0 * float(planet_prior["sma"]))
    )
    temp_bounds = (max(300.0, 0.75 * eq_temperature), 1.5 * eq_temperature)
    radius = float(planet_prior["rp"])
    radius_sigma = max(
        abs(float(planet_prior.get("rp_uperr", 0.0))),
        abs(float(planet_prior.get("rp_lowerr", 0.0))),
    )
    radius_bounds = (max(radius * 0.5, radius - 3 * radius_sigma), radius + 3 * radius_sigma)
    lower = np.asarray([-12.0] * len(MOLECULES) + [temp_bounds[0], radius_bounds[0]], dtype=np.float32)
    upper = np.asarray([-2.0] * len(MOLECULES) + [temp_bounds[1], radius_bounds[1]], dtype=np.float32)
    if not np.all(np.isfinite(lower)) or not np.all(upper > lower):
        raise ValueError("Could not construct finite, non-empty parameter bounds.")
    return lower, upper


def _sample_mixture_abundances(
    rng: np.random.Generator,
    lower: np.ndarray,
    upper: np.ndarray,
    total_log_vmr_bounds: tuple[float, float] = (-5.0, -2.2),
) -> np.ndarray:
    """Sample sparse-to-crowded molecular mixtures with varied abundance ratios."""
    abundance_lower = lower[: len(MOLECULES)]
    abundance_upper = upper[: len(MOLECULES)]
    # Balance the number of active labels instead of making crowded mixtures dominate.
    n_active = int(rng.integers(1, len(MOLECULES) + 1))
    active_ids = rng.choice(len(MOLECULES), size=n_active, replace=False)
    # The Dirichlet concentration controls mixture shape: small values are spiky,
    # while large values distribute the total VMR more evenly.
    dirichlet_concentration = float(rng.choice((0.25, 1.0, 5.0)))
    fractions = rng.dirichlet(
        np.full(n_active, dirichlet_concentration, dtype=np.float64)
    )
    total_vmr = 10.0 ** rng.uniform(*total_log_vmr_bounds)
    minimum_fraction = 10.0 ** (
        float(np.max(abundance_lower)) - math.log10(total_vmr)
    )
    fractions = minimum_fraction + (1.0 - n_active * minimum_fraction) * fractions
    abundances = rng.uniform(
        abundance_lower,
        np.minimum(abundance_upper, DETECTION_FLOOR_LOG_VMR - 0.5),
    ).astype(np.float32)
    abundances[active_ids] = np.log10(total_vmr * fractions).astype(np.float32)
    if np.any(abundances < abundance_lower) or np.any(abundances > abundance_upper):
        raise ValueError("Mixture abundances exceed the configured molecular bounds.")
    if np.power(10.0, abundances.astype(np.float64)).sum() > 0.1:
        raise ValueError("Sampled molecular mixture exceeds the 0.1 total-VMR limit.")
    return abundances


def _sample_scenario(
    rng: np.random.Generator,
    inputs: ForwardInputs,
    lower: np.ndarray,
    upper: np.ndarray,
    scenario_strategy: str = "uniform_log",
) -> tuple[np.ndarray, np.ndarray]:
    priors = inputs.priors
    planet_prior = priors[inputs.planet]
    radius_star_sigma = max(
        abs(float(priors.get("R*_uperr", 0.0))),
        abs(float(priors.get("R*_lowerr", 0.0))),
    )
    temperature_star_sigma = max(
        abs(float(priors.get("T*_uperr", 0.0))),
        abs(float(priors.get("T*_lowerr", 0.0))),
    )
    logg_sigma = max(
        abs(float(planet_prior.get("logg_uperr", 0.0))),
        abs(float(planet_prior.get("logg_lowerr", 0.0))),
    )

    conditions = np.asarray(
        [
            np.clip(
                rng.normal(float(priors["R*"]), radius_star_sigma),
                float(priors["R*"]) - 3 * radius_star_sigma,
                float(priors["R*"]) + 3 * radius_star_sigma,
            ),
            np.clip(
                rng.normal(float(priors["T*"]), temperature_star_sigma),
                float(priors["T*"]) - 3 * temperature_star_sigma,
                float(priors["T*"]) + 3 * temperature_star_sigma,
            ),
            np.clip(
                rng.normal(float(planet_prior["logg"]), logg_sigma),
                float(planet_prior["logg"]) - 3 * logg_sigma,
                float(planet_prior["logg"]) + 3 * logg_sigma,
            ),
            float(planet_prior["sma"]),
        ],
        dtype=np.float32,
    )

    targets = np.empty(len(TARGET_NAMES), dtype=np.float32)
    if scenario_strategy == "uniform_log":
        targets[: len(MOLECULES)] = rng.uniform(
            lower[: len(MOLECULES)], upper[: len(MOLECULES)]
        )
    elif scenario_strategy == "mixture_balanced_v1":
        targets[: len(MOLECULES)] = _sample_mixture_abundances(rng, lower, upper)
    elif scenario_strategy == "mixture_balanced_v2":
        targets[: len(MOLECULES)] = _sample_mixture_abundances(
            rng,
            lower,
            upper,
            total_log_vmr_bounds=(-3.5, -2.0),
        )
    else:
        raise ValueError(
            "scenario_strategy must be 'uniform_log', 'mixture_balanced_v1', "
            "or 'mixture_balanced_v2'."
        )
    targets[-2:] = rng.uniform(lower[-2:], upper[-2:])
    return targets, conditions


def _simulate_spectrum(
    targets: np.ndarray,
    conditions: np.ndarray,
    inputs: ForwardInputs,
) -> np.ndarray:
    import excalibur.system.core as syscore
    import excalibur.cerberus.forward_model as crbfm

    temperature = float(targets[-2])
    radius = float(targets[-1])
    sample_priors = dict(inputs.priors)
    sample_priors["R*"] = float(conditions[0])
    sample_priors["T*"] = float(conditions[1])
    sample_planet = dict(inputs.priors[inputs.planet])
    sample_planet["logg"] = float(conditions[2])
    sample_priors[inputs.planet] = sample_planet

    mixratio = {
        molecule: np.full(100, float(abundance) + 6.0)
        for molecule, abundance in zip(MOLECULES, targets[: len(MOLECULES)])
    }
    forward_model = crbfm.crbFM().crbmodel(
        np.full(100, temperature, dtype=np.float64),
        1,
        cheq=None,
        mixratio=mixratio,
        hazescale=-10,
        hazethick=3.1,
        hazeloc=-3.43,
        hazeprof="AVERAGE",
        hzlib=inputs.haze_library,
        planet=inputs.planet,
        rp0=radius * syscore.ssconstants(mks=True)["Rjup"],
        orbp=sample_priors,
        wgrid=inputs.wavelengths,
        xsecs=inputs.cross_sections,
        qtgrid=inputs.qtgrid,
        cialist=["H2-H", "H2-H2", "H2-He", "He-H"],
        xmollist=[
            "TIO", "H2O", "H2CO", "HCN", "CO", "CO2", "NH3",
            "CH4", "C2H2", "C2H6", "C3H8", "CH3CHO", "SO2", "H2S",
        ],
        nlevels=100,
        Hsmax=25,
        solrad=1,
        break_down_by_molecule=False,
        logx=False,
        verbose=False,
        debug=False,
    )
    spectrum = np.asarray(forward_model.spectrum, dtype=np.float32)
    if spectrum.shape != inputs.wavelengths.shape or not np.all(np.isfinite(spectrum)):
        raise FloatingPointError("CERBERUS returned an invalid synthetic spectrum.")
    return spectrum


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    def clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        if isinstance(value, list):
            return [clean(item) for item in value]
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            return None
        if isinstance(value, np.integer):
            return int(value)
        return value

    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(clean(payload), indent=2), encoding="utf-8")
    temporary.replace(path)


def generate_dataset(
    inputs: ForwardInputs,
    output_dir: str | Path,
    n_samples: int = 10240,
    seed: int = 42,
    checkpoint_every: int = 128,
    scenario_strategy: str = "uniform_log",
) -> Path:
    """Generate resumable, noisy CERBERUS spectra and disjoint split indices."""
    if n_samples < 100:
        raise ValueError("Use at least 100 simulations for a meaningful split.")
    if checkpoint_every < 1:
        raise ValueError("checkpoint_every must be positive.")
    if scenario_strategy not in {
        "uniform_log",
        "mixture_balanced_v1",
        "mixture_balanced_v2",
    }:
        raise ValueError(
            "scenario_strategy must be 'uniform_log', 'mixture_balanced_v1', "
            "or 'mixture_balanced_v2'."
        )
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    dataset_path = output / "synthetic_dataset.npz"
    split_path = output / "splits.npz"
    manifest_path = output / "manifest.json"
    work = output / ".generation"
    work.mkdir(exist_ok=True)

    lower, upper = _target_bounds(inputs.priors, inputs.planet)
    config = {
        "n_samples": int(n_samples),
        "seed": int(seed),
        "target": inputs.planet,
        "n_channels": int(inputs.wavelengths.size),
        "target_names": list(TARGET_NAMES),
        "scenario_strategy": scenario_strategy,
    }
    if dataset_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        existing_config = dict(manifest.get("generation", {}))
        existing_config.setdefault("scenario_strategy", "uniform_log")
        if existing_config != config:
            raise ValueError(
                f"{dataset_path} already exists with different settings; move it or choose a new output directory."
            )
        return dataset_path

    state_path = work / "progress.json"
    array_paths = {
        "spectra": work / "spectra.npy",
        "targets": work / "targets.npy",
        "noise_sigma": work / "noise_sigma.npy",
        "noise_token": work / "noise_token.npy",
        "conditions": work / "conditions.npy",
    }
    resume = False
    start = 0
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("generation") != config:
            raise ValueError(
                "A partial generation used different settings. Remove only the .generation folder to restart."
            )
        start = int(state["next_index"])
        resume = all(path.exists() for path in array_paths.values())
        if not resume:
            start = 0

    shapes = {
        "spectra": (n_samples, inputs.wavelengths.size),
        "targets": (n_samples, len(TARGET_NAMES)),
        "noise_sigma": (n_samples, inputs.wavelengths.size),
        "noise_token": (n_samples,),
        "conditions": (n_samples, len(CONDITION_NAMES)),
    }
    arrays: dict[str, np.memmap] = {}
    for name, path in array_paths.items():
        if not resume:
            path.unlink(missing_ok=True)
        arrays[name] = np.lib.format.open_memmap(
            path,
            mode="r+" if resume else "w+",
            dtype=np.float32,
            shape=shapes[name],
        )

    for index in range(start, n_samples):
        rng = np.random.default_rng(np.random.SeedSequence([seed, index]))
        targets, conditions = _sample_scenario(
            rng, inputs, lower, upper, scenario_strategy=scenario_strategy
        )
        noiseless = _simulate_spectrum(targets, conditions, inputs)
        noise_multiplier = 10.0 ** rng.uniform(-0.3, 0.3)
        relative_sigma = inputs.relative_noise * noise_multiplier
        sigma = np.maximum(np.abs(noiseless) * relative_sigma, 1e-12)
        spectra = noiseless + rng.normal(0.0, sigma).astype(np.float32)

        arrays["spectra"][index] = spectra
        arrays["targets"][index] = targets
        arrays["noise_sigma"][index] = sigma
        arrays["noise_token"][index] = np.log10(np.median(relative_sigma))
        arrays["conditions"][index] = conditions

        if (index + 1) % checkpoint_every == 0 or index + 1 == n_samples:
            for array in arrays.values():
                array.flush()
            _write_json(
                state_path,
                {"generation": config, "next_index": index + 1},
            )
            print(f"Generated {index + 1:,}/{n_samples:,} spectra")

    indices = np.random.default_rng(seed).permutation(n_samples)
    n_train = int(0.60 * n_samples)
    n_validation = int(0.10 * n_samples)
    n_calibration = int(0.15 * n_samples)
    splits = {
        "train": indices[:n_train],
        "validation": indices[n_train : n_train + n_validation],
        "calibration": indices[n_train + n_validation : n_train + n_validation + n_calibration],
        "test": indices[n_train + n_validation + n_calibration :],
    }
    if sum(len(part) for part in splits.values()) != n_samples:
        raise RuntimeError("The generated data splits do not cover the dataset.")
    flattened = np.concatenate(list(splits.values()))
    if np.unique(flattened).size != n_samples:
        raise RuntimeError("The generated data splits overlap.")

    np.savez_compressed(
        dataset_path,
        spectra=arrays["spectra"],
        targets=arrays["targets"],
        active_species_count=np.sum(
            arrays["targets"][:, : len(MOLECULES)] > DETECTION_FLOOR_LOG_VMR,
            axis=1,
        ).astype(np.int8),
        noise_sigma=arrays["noise_sigma"],
        noise_token=arrays["noise_token"],
        conditions=arrays["conditions"],
        wavelengths=inputs.wavelengths.astype(np.float32),
        target_names=np.asarray(TARGET_NAMES),
        condition_names=np.asarray(CONDITION_NAMES),
        lower_bounds=lower,
        upper_bounds=upper,
    )
    np.savez_compressed(split_path, **splits)
    target_values = np.asarray(arrays["targets"])
    active_species_counts = (
        target_values[:, : len(MOLECULES)] > DETECTION_FLOOR_LOG_VMR
    ).sum(axis=1)
    total_molecular_vmr = np.power(
        10.0, target_values[:, : len(MOLECULES)].astype(np.float64)
    ).sum(axis=1)
    _write_json(
        manifest_path,
        {
            "generation": config,
            "split_sizes": {name: int(len(part)) for name, part in splits.items()},
            "noise": "per-channel fractional JWST ESerr/ES, scaled by a sampled factor",
            "abundance_units": "log10 VMR",
            "CERBERUS_units": "log10 ppm; log10 VMR + 6 passed to mixratio",
            "scenario_strategy": scenario_strategy,
            "mixture_sampling": (
                "Uniform log10 VMR for each molecule."
                if scenario_strategy == "uniform_log"
                else (
                    "Uniform active-molecule count 1-9; symmetric Dirichlet "
                    "mixture fractions with concentration 0.25, 1, or 5; "
                    "log10 total molecular VMR uniform from -5 to -2.2; "
                    "inactive labels sampled below the -8 detection floor."
                    if scenario_strategy == "mixture_balanced_v1"
                    else (
                        "Uniform active-molecule count 1-9; symmetric Dirichlet "
                        "mixture fractions with concentration 0.25, 1, or 5; "
                        "log10 total molecular VMR uniform from -3.5 to -2; "
                        "inactive labels sampled below the -8 detection floor."
                    )
                )
            ),
            "scenario_summary": {
                "active_species_count_histogram": {
                    str(count): int(np.count_nonzero(active_species_counts == count))
                    for count in range(len(MOLECULES) + 1)
                },
                "total_molecular_vmr_quantiles_0_50_95_100": [
                    float(value)
                    for value in np.quantile(
                        total_molecular_vmr, [0.0, 0.5, 0.95, 1.0]
                    )
                ],
            },
            "outputs_are_synthetic": True,
        },
    )
    arrays.clear()
    for path in array_paths.values():
        path.unlink(missing_ok=True)
    state_path.unlink(missing_ok=True)
    work.rmdir()
    return dataset_path


def sinusoidal_position_encoding(
    length: int, dimensions: int, device: torch.device, dtype: torch.dtype
) -> torch.Tensor:
    if dimensions % 2:
        raise ValueError("The sinusoidal encoding dimension must be even.")
    position = torch.arange(length, device=device, dtype=dtype).unsqueeze(1)
    frequency = torch.exp(
        torch.arange(0, dimensions, 2, device=device, dtype=dtype)
        * (-math.log(10000.0) / dimensions)
    )
    encoding = torch.zeros(length, dimensions, device=device, dtype=dtype)
    encoding[:, 0::2] = torch.sin(position * frequency)
    encoding[:, 1::2] = torch.cos(position * frequency)
    return encoding


def _transform_to_latent(
    targets: torch.Tensor,
    lower_bounds: torch.Tensor,
    upper_bounds: torch.Tensor,
    logistic_scale: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    width = upper_bounds - lower_bounds
    unit = (targets - lower_bounds) / width
    epsilon = torch.finfo(targets.dtype).eps
    unit = unit.clamp(epsilon, 1.0 - epsilon)
    latent = torch.logit(unit) / logistic_scale
    log_jacobian = (
        -math.log(logistic_scale)
        - torch.log(width)
        - torch.log(unit)
        - torch.log1p(-unit)
    ).sum(dim=-1)
    return latent, log_jacobian


class ResidualBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int, stride: int = 1):
        super().__init__()
        self.conv1 = nn.Conv1d(
            in_channels, out_channels, kernel_size=7, stride=stride, padding=3, bias=False
        )
        self.norm1 = nn.BatchNorm1d(out_channels)
        self.conv2 = nn.Conv1d(out_channels, out_channels, kernel_size=3, padding=1, bias=False)
        self.norm2 = nn.BatchNorm1d(out_channels)
        self.activation = nn.GELU()
        self.skip = (
            nn.Identity()
            if in_channels == out_channels and stride == 1
            else nn.Sequential(
                nn.Conv1d(in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm1d(out_channels),
            )
        )

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        residual = self.skip(inputs)
        hidden = self.activation(self.norm1(self.conv1(inputs)))
        hidden = self.norm2(self.conv2(hidden))
        return self.activation(hidden + residual)


class FlowAtmos(nn.Module):
    """Encode a spectrum and model its physical-parameter posterior.

    The CNN turns wavelength-dependent spectral features into a compact context.
    A conditional spline flow then maps that context to a flexible joint posterior
    over molecular abundances and planet parameters.
    """

    def __init__(
        self,
        n_channels: int,
        n_conditions: int,
        lower_bounds: np.ndarray,
        upper_bounds: np.ndarray,
    ):
        super().__init__()
        self.n_channels = int(n_channels)
        self.encoder_stem = nn.Sequential(
            nn.Conv1d(17, 64, kernel_size=7, padding=3, bias=False),
            nn.BatchNorm1d(64),
            nn.GELU(),
        )
        self.stage1 = nn.Sequential(ResidualBlock(64, 64), ResidualBlock(64, 64))
        self.stage2 = nn.Sequential(ResidualBlock(64, 128, stride=2), ResidualBlock(128, 128))
        self.stage3 = nn.Sequential(ResidualBlock(128, 256, stride=2), ResidualBlock(256, 256))
        self.attention = nn.Conv1d(256, 1, kernel_size=1)
        self.noise_embedding = nn.Sequential(
            nn.Linear(1, 32),
            nn.GELU(),
            nn.Linear(32, 64),
            nn.GELU(),
        )
        self.context_head = nn.Sequential(
            nn.Linear(256 + 64 + n_conditions, 384),
            nn.GELU(),
            nn.Linear(384, 256),
            nn.GELU(),
        )
        self.flow = NSF(
            features=len(lower_bounds),
            context=256,
            transforms=10,
            bins=16,
            hidden_features=[384, 384, 384],
            activation=nn.GELU,
        )
        self.register_buffer("lower_bounds", torch.as_tensor(lower_bounds, dtype=torch.float32))
        self.register_buffer("upper_bounds", torch.as_tensor(upper_bounds, dtype=torch.float32))
        self.logistic_scale = math.pi / math.sqrt(3.0)

    def encode(
        self,
        spectra: torch.Tensor,
        noise_tokens: torch.Tensor,
        conditions: torch.Tensor,
    ) -> torch.Tensor:
        if spectra.ndim != 2 or spectra.shape[1] != self.n_channels:
            raise ValueError(f"Expected spectra with shape (batch, {self.n_channels}).")
        if noise_tokens.ndim == 1:
            noise_tokens = noise_tokens.unsqueeze(1)
        # Position features let the shared convolution filters distinguish wavelengths.
        position = sinusoidal_position_encoding(
            self.n_channels, 16, spectra.device, spectra.dtype
        ).transpose(0, 1)
        position = position.unsqueeze(0).expand(spectra.shape[0], -1, -1)
        hidden = torch.cat((spectra.unsqueeze(1), position), dim=1)
        hidden = self.encoder_stem(hidden)
        hidden = self.stage1(hidden)
        hidden = self.stage2(hidden)
        hidden = self.stage3(hidden)
        weights = torch.softmax(self.attention(hidden), dim=-1)
        pooled = torch.sum(hidden * weights, dim=-1)
        noise = self.noise_embedding(noise_tokens)
        return self.context_head(torch.cat((pooled, noise, conditions), dim=-1))

    def to_latent(self, targets: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return _transform_to_latent(
            targets,
            self.lower_bounds,
            self.upper_bounds,
            self.logistic_scale,
        )

    def from_latent(self, latent: torch.Tensor) -> torch.Tensor:
        unit = torch.sigmoid(self.logistic_scale * latent)
        return self.lower_bounds + (self.upper_bounds - self.lower_bounds) * unit

    def log_prob(
        self,
        spectra: torch.Tensor,
        noise_tokens: torch.Tensor,
        conditions: torch.Tensor,
        targets: torch.Tensor,
    ) -> torch.Tensor:
        context = self.encode(spectra, noise_tokens, conditions)
        latent, log_jacobian = self.to_latent(targets)
        # Add the change-of-variables term so density is reported in physical units.
        return self.flow(context).log_prob(latent) + log_jacobian


class UnconditionalPriorFlow(nn.Module):
    """Fit the target prior without spectrum context as a matched NLL baseline."""

    def __init__(self, lower_bounds: np.ndarray, upper_bounds: np.ndarray):
        super().__init__()
        self.register_buffer(
            "lower_bounds", torch.as_tensor(lower_bounds, dtype=torch.float32)
        )
        self.register_buffer(
            "upper_bounds", torch.as_tensor(upper_bounds, dtype=torch.float32)
        )
        self.logistic_scale = math.pi / math.sqrt(3.0)
        self.flow = NSF(
            features=len(lower_bounds),
            context=0,
            transforms=10,
            bins=16,
            hidden_features=[384, 384, 384],
            activation=nn.GELU,
        )

    def log_prob(self, targets: torch.Tensor) -> torch.Tensor:
        latent, log_jacobian = _transform_to_latent(
            targets,
            self.lower_bounds,
            self.upper_bounds,
            self.logistic_scale,
        )
        return self.flow().log_prob(latent) + log_jacobian


def _prepare_inputs(
    dataset: dict[str, np.ndarray],
    train_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """Normalize spectra and retain their original scale in the context."""
    spectra = dataset["spectra"].astype(np.float32)
    spectrum_mean = spectra.mean(axis=1, keepdims=True)
    spectrum_std = spectra.std(axis=1, keepdims=True)
    spectra = (spectra - spectrum_mean) / np.maximum(spectrum_std, 1e-8)
    conditions = np.concatenate(
        (
            dataset["conditions"].astype(np.float32),
            spectrum_mean.astype(np.float32),
            spectrum_std.astype(np.float32),
        ),
        axis=1,
    )
    condition_mean = conditions[train_indices].mean(axis=0)
    condition_std = conditions[train_indices].std(axis=0)
    condition_std = np.maximum(condition_std, 1e-8)
    conditions = (conditions - condition_mean) / condition_std
    noise_tokens = dataset["noise_token"].astype(np.float32)
    stats = {
        "condition_mean": condition_mean.astype(np.float32),
        "condition_std": condition_std.astype(np.float32),
        "spectrum_mean": np.asarray([0.0], dtype=np.float32),
        "spectrum_std": np.asarray([1.0], dtype=np.float32),
    }
    return spectra, noise_tokens, conditions, stats


def _choose_device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def _seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def _epoch_nll(
    model: FlowAtmos,
    loader: DataLoader,
    device: torch.device,
    optimizer: torch.optim.Optimizer | None = None,
) -> float:
    training = optimizer is not None
    model.train(training)
    total = 0.0
    count = 0
    for spectra, noise, conditions, targets in loader:
        spectra = spectra.to(device=device, dtype=torch.float32)
        noise = noise.to(device=device, dtype=torch.float32)
        conditions = conditions.to(device=device, dtype=torch.float32)
        targets = targets.to(device=device, dtype=torch.float32)
        if optimizer is not None:
            optimizer.zero_grad(set_to_none=True)
        with torch.set_grad_enabled(training):
            loss = -model.log_prob(spectra, noise, conditions, targets).mean()
            if not torch.isfinite(loss):
                raise FloatingPointError("FlowAtmos produced a non-finite NLL.")
            if optimizer is not None:
                loss.backward()
                nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
        total += float(loss.detach().cpu()) * len(spectra)
        count += len(spectra)
    if count == 0:
        raise ValueError("An NLL data loader is empty.")
    return total / count


def _prior_epoch_nll(
    model: UnconditionalPriorFlow,
    loader: DataLoader,
    device: torch.device,
    optimizer: torch.optim.Optimizer | None = None,
) -> float:
    training = optimizer is not None
    model.train(training)
    total = 0.0
    count = 0
    for (targets,) in loader:
        targets = targets.to(device=device, dtype=torch.float32)
        if optimizer is not None:
            optimizer.zero_grad(set_to_none=True)
        with torch.set_grad_enabled(training):
            loss = -model.log_prob(targets).mean()
            if not torch.isfinite(loss):
                raise FloatingPointError("The unconditional prior flow produced non-finite NLL.")
            if optimizer is not None:
                loss.backward()
                nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()
        total += float(loss.detach().cpu()) * len(targets)
        count += len(targets)
    if count == 0:
        raise ValueError("The prior-flow data loader is empty.")
    return total / count


def _fit_prior_baseline(
    targets: np.ndarray,
    train_indices: np.ndarray,
    validation_indices: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    output: Path,
    seed: int,
    max_epochs: int = 300,
    patience: int = 25,
) -> dict[str, Any]:
    """Fit an unconditional flow to estimate the dataset's sampling prior."""
    _seed_everything(seed)
    train_loader = DataLoader(
        TensorDataset(torch.from_numpy(targets[train_indices])),
        batch_size=512,
        shuffle=True,
        generator=torch.Generator().manual_seed(seed),
    )
    validation_loader = DataLoader(
        TensorDataset(torch.from_numpy(targets[validation_indices])),
        batch_size=512,
        shuffle=False,
    )
    device = _choose_device()
    model = UnconditionalPriorFlow(lower, upper).to(
        device=device, dtype=torch.float32
    )
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-5)
    best_loss = math.inf
    best_epoch = 0
    stale_epochs = 0
    history: list[dict[str, float]] = []
    checkpoint_path = output / "flowatmos_prior_baseline.pt"

    for epoch in range(max_epochs):
        if epoch < 5:
            learning_rate = 3e-4 * (epoch + 1) / 5
        else:
            progress = (epoch - 5) / max(1, max_epochs - 6)
            learning_rate = 5e-6 + (3e-4 - 5e-6) * 0.5 * (
                1.0 + math.cos(math.pi * progress)
            )
        for group in optimizer.param_groups:
            group["lr"] = learning_rate
        train_nll = _prior_epoch_nll(model, train_loader, device, optimizer)
        validation_nll = _prior_epoch_nll(model, validation_loader, device)
        history.append(
            {
                "epoch": float(epoch + 1),
                "learning_rate": float(learning_rate),
                "train_nll": train_nll,
                "validation_nll": validation_nll,
            }
        )
        print(
            f"prior epoch {epoch + 1:03d} | train NLL {train_nll:.4f} | "
            f"validation NLL {validation_nll:.4f}",
            flush=True,
        )
        if validation_nll < best_loss:
            best_loss = validation_nll
            best_epoch = epoch + 1
            stale_epochs = 0
            torch.save(
                {
                    "state_dict": {
                        name: value.detach().cpu()
                        for name, value in model.state_dict().items()
                    },
                    "epoch": best_epoch,
                    "validation_nll": best_loss,
                    "seed": seed,
                    "lower_bounds": torch.as_tensor(lower),
                    "upper_bounds": torch.as_tensor(upper),
                },
                checkpoint_path,
            )
        else:
            stale_epochs += 1
        if stale_epochs >= patience:
            print(f"Prior-flow early stopping after {epoch + 1} epochs.")
            break

    with (output / "prior_baseline_history.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=history[0].keys())
        writer.writeheader()
        writer.writerows(history)
    summary = {
        "best_epoch": best_epoch,
        "best_validation_nll": best_loss,
        "epochs_run": len(history),
        "device": str(device),
        "precision": "float32",
        "seed": seed,
        "train_count": int(len(train_indices)),
        "validation_count": int(len(validation_indices)),
        "checkpoint": str(checkpoint_path),
    }
    _write_json(output / "prior_baseline_summary.json", summary)
    return summary


def train_flowatmos(
    dataset_path: str | Path,
    output_dir: str | Path,
    seed: int = 42,
    max_epochs: int = 300,
    patience: int = 25,
    progress_callback: Callable[[dict[str, float]], None] | None = None,
) -> dict[str, Any]:
    """Train FlowAtmos in fp32 with a validation-based early-stopping checkpoint."""
    dataset_path = Path(dataset_path).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    if max_epochs < 6 or patience < 1:
        raise ValueError("Use at least six epochs and a positive early-stopping patience.")
    _seed_everything(seed)

    with np.load(dataset_path) as saved:
        dataset = {name: saved[name] for name in saved.files}
    with np.load(output / "splits.npz") as saved_splits:
        splits = {name: saved_splits[name] for name in saved_splits.files}
    train_indices = splits["train"]
    validation_indices = splits["validation"]
    lower = dataset["lower_bounds"].astype(np.float32)
    upper = dataset["upper_bounds"].astype(np.float32)
    if len(train_indices) < 512 or len(validation_indices) == 0:
        raise ValueError("The training split must contain at least one full batch.")

    spectra, noise_tokens, conditions, preprocessing = _prepare_inputs(dataset, train_indices)
    targets = dataset["targets"].astype(np.float32)
    def tensor_data(indices: np.ndarray) -> TensorDataset:
        return TensorDataset(
            torch.from_numpy(spectra[indices]),
            torch.from_numpy(noise_tokens[indices]),
            torch.from_numpy(conditions[indices]),
            torch.from_numpy(targets[indices]),
        )

    generator = torch.Generator().manual_seed(seed)
    train_loader = DataLoader(
        tensor_data(train_indices), batch_size=512, shuffle=True, generator=generator
    )
    validation_loader = DataLoader(
        tensor_data(validation_indices), batch_size=512, shuffle=False
    )
    device = _choose_device()
    model = FlowAtmos(
        n_channels=spectra.shape[1],
        n_conditions=conditions.shape[1],
        lower_bounds=lower,
        upper_bounds=upper,
    ).to(device=device, dtype=torch.float32)
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-5)
    best_loss = math.inf
    best_epoch = 0
    stale_epochs = 0
    history: list[dict[str, float]] = []
    checkpoint_path = output / "flowatmos_best.pt"

    for epoch in range(max_epochs):
        if epoch < 5:
            learning_rate = 3e-4 * (epoch + 1) / 5
        else:
            progress = (epoch - 5) / max(1, max_epochs - 6)
            learning_rate = 5e-6 + (3e-4 - 5e-6) * 0.5 * (
                1.0 + math.cos(math.pi * progress)
            )
        for group in optimizer.param_groups:
            group["lr"] = learning_rate

        train_nll = _epoch_nll(model, train_loader, device, optimizer)
        validation_nll = _epoch_nll(model, validation_loader, device)
        record = {
            "epoch": float(epoch + 1),
            "learning_rate": float(learning_rate),
            "train_nll": float(train_nll),
            "validation_nll": float(validation_nll),
        }
        history.append(record)
        if progress_callback is not None:
            progress_callback(record)
        print(
            f"epoch {epoch + 1:03d} | train NLL {train_nll:.4f} | "
            f"validation NLL {validation_nll:.4f}"
        )

        if validation_nll < best_loss:
            best_loss = validation_nll
            best_epoch = epoch + 1
            stale_epochs = 0
            torch.save(
                {
                    "state_dict": {
                        name: value.detach().cpu()
                        for name, value in model.state_dict().items()
                    },
                    "epoch": best_epoch,
                    "validation_nll": best_loss,
                    "seed": seed,
                    "lower_bounds": torch.as_tensor(lower),
                    "upper_bounds": torch.as_tensor(upper),
                    "conditioning_features": list(MODEL_CONDITION_NAMES),
                    "preprocessing": {
                        name: torch.as_tensor(value)
                        for name, value in preprocessing.items()
                    },
                },
                checkpoint_path,
            )
        else:
            stale_epochs += 1
        if stale_epochs >= patience:
            print(f"Early stopping after {epoch + 1} epochs.")
            break

    with (output / "training_history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=history[0].keys())
        writer.writeheader()
        writer.writerows(history)
    summary = {
        "best_epoch": best_epoch,
        "best_validation_nll": best_loss,
        "epochs_run": len(history),
        "device": str(device),
        "precision": "float32",
        "seed": seed,
        "max_epochs": max_epochs,
        "patience": patience,
        "training_sample_count": int(len(train_indices)),
        "conditioning_features": list(MODEL_CONDITION_NAMES),
    }
    manifest_path = output / "manifest.json"
    scenario_strategy = "uniform_log"
    if manifest_path.exists():
        generation = json.loads(manifest_path.read_text(encoding="utf-8")).get(
            "generation", {}
        )
        scenario_strategy = generation.get("scenario_strategy", "uniform_log")
    summary["scenario_strategy"] = scenario_strategy
    if scenario_strategy != "uniform_log":
        summary["matched_prior_baseline"] = _fit_prior_baseline(
            targets,
            train_indices,
            validation_indices,
            lower,
            upper,
            output,
            seed + 1,
            max_epochs=max_epochs,
            patience=patience,
        )
    _write_json(output / "training_summary.json", summary)
    return summary


def _load_model(
    dataset: dict[str, np.ndarray],
    checkpoint: dict[str, Any],
    device: torch.device,
) -> tuple[FlowAtmos, dict[str, np.ndarray]]:
    model = FlowAtmos(
        n_channels=dataset["spectra"].shape[1],
        n_conditions=int(checkpoint["preprocessing"]["condition_mean"].numel()),
        lower_bounds=checkpoint["lower_bounds"].cpu().numpy(),
        upper_bounds=checkpoint["upper_bounds"].cpu().numpy(),
    ).to(device=device, dtype=torch.float32)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    preprocessing = {
        name: value.cpu().numpy()
        for name, value in checkpoint["preprocessing"].items()
    }
    return model, preprocessing


def _normalised_arrays(
    dataset: dict[str, np.ndarray], preprocessing: dict[str, np.ndarray]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    spectra = dataset["spectra"].astype(np.float32)
    spectrum_mean = spectra.mean(axis=1, keepdims=True)
    spectrum_std = spectra.std(axis=1, keepdims=True)
    spectra = (spectra - spectrum_mean) / np.maximum(spectrum_std, 1e-8)
    conditions = dataset["conditions"].astype(np.float32)
    conditioning_width = len(preprocessing["condition_mean"])
    if conditioning_width == conditions.shape[1] + 2:
        conditions = np.concatenate(
            (conditions, spectrum_mean.astype(np.float32), spectrum_std.astype(np.float32)),
            axis=1,
        )
    elif conditioning_width != conditions.shape[1]:
        raise ValueError("Checkpoint conditioning features do not match the dataset.")
    conditions = (
        conditions - preprocessing["condition_mean"]
    ) / preprocessing["condition_std"]
    return spectra, dataset["noise_token"].astype(np.float32), conditions


@torch.no_grad()
def _posterior_draws(
    model: FlowAtmos,
    spectra: np.ndarray,
    noise: np.ndarray,
    conditions: np.ndarray,
    indices: np.ndarray,
    n_draws: int,
    device: torch.device,
    batch_size: int = 8,
) -> np.ndarray:
    draws = np.empty(
        (len(indices), n_draws, len(TARGET_NAMES)),
        dtype=np.float32,
    )
    for start in range(0, len(indices), batch_size):
        batch_indices = indices[start : start + batch_size]
        x = torch.as_tensor(spectra[batch_indices], device=device, dtype=torch.float32)
        n = torch.as_tensor(noise[batch_indices], device=device, dtype=torch.float32)
        c = torch.as_tensor(conditions[batch_indices], device=device, dtype=torch.float32)
        context = model.encode(x, n, c)
        latent = model.flow(context).sample((n_draws,))
        samples = model.from_latent(latent).permute(1, 0, 2).contiguous()
        draws[start : start + len(batch_indices)] = samples.cpu().numpy()
        completed = start + len(batch_indices)
        if completed == len(indices) or completed % 100 == 0:
            print(
                f"Posterior sampling: {completed:,}/{len(indices):,} spectra "
                f"at {n_draws:,} draws each",
                flush=True,
            )
    return draws


def _rank_molecule_samples(
    samples: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if samples.ndim != 2 or samples.shape[1] != len(MOLECULES):
        raise ValueError(
            f"Expected molecule samples with shape (draws, {len(MOLECULES)})."
        )
    if samples.shape[0] < 1 or not np.all(np.isfinite(samples)):
        raise ValueError("Molecule samples must be non-empty and finite.")
    medians = np.median(samples, axis=0)
    ranks = np.argsort(np.argsort(-samples, axis=1), axis=1)
    rank_probabilities = np.stack(
        [
            np.bincount(ranks[:, molecule], minlength=len(MOLECULES))
            / samples.shape[0]
            for molecule in range(len(MOLECULES))
        ]
    )
    point_order = np.argsort(medians)[::-1]
    return medians, rank_probabilities, point_order


class FlowAtmosPredictor:
    """Load a trained model once and reuse it for low-latency spectrum ranking.

    Keep one predictor alive for a batch of observed spectra: checkpoint loading
    and model setup are fixed overheads that should not be repeated per target.
    """

    def __init__(
        self,
        checkpoint_path: str | Path,
        n_channels: int,
        device: str = "auto",
    ):
        if n_channels < 1:
            raise ValueError("n_channels must be positive.")
        if device == "auto":
            chosen_device = _choose_device()
        elif device in {"cpu", "mps", "cuda"}:
            chosen_device = torch.device(device)
        else:
            raise ValueError("device must be 'auto', 'cpu', 'mps', or 'cuda'.")
        if chosen_device.type == "mps" and not (
            hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
        ):
            raise RuntimeError("MPS inference was requested but is unavailable.")
        if chosen_device.type == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA inference was requested but is unavailable.")

        checkpoint = torch.load(
            Path(checkpoint_path).expanduser().resolve(),
            map_location="cpu",
            weights_only=True,
        )
        model_data = {
            "spectra": np.empty((1, n_channels), dtype=np.float32),
            "conditions": np.empty((1, len(CONDITION_NAMES)), dtype=np.float32),
        }
        self.model, self.preprocessing = _load_model(
            model_data, checkpoint, chosen_device
        )
        self.device = chosen_device
        self.n_channels = n_channels

    @torch.no_grad()
    def predict_molecule_ranks(
        self,
        spectrum: np.ndarray,
        noise_token: float,
        conditions: np.ndarray,
        posterior_draws: int = 256,
    ) -> dict[str, Any]:
        """Return a fast top-five posterior ranking and rank-probability matrix.

        Detection gating and formal interval coverage remain part of the
        separately calibrated offline evaluation.
        """
        spectrum = np.asarray(spectrum, dtype=np.float32)
        conditions = np.asarray(conditions, dtype=np.float32)
        if spectrum.shape != (self.n_channels,) or not np.all(np.isfinite(spectrum)):
            raise ValueError(
                f"spectrum must contain {self.n_channels} finite channel values."
            )
        if conditions.shape != (len(CONDITION_NAMES),) or not np.all(
            np.isfinite(conditions)
        ):
            raise ValueError(
                f"conditions must contain {len(CONDITION_NAMES)} finite values in "
                f"the order {CONDITION_NAMES}."
            )
        if not np.isfinite(noise_token):
            raise ValueError("noise_token must be finite.")
        if posterior_draws < 1:
            raise ValueError("posterior_draws must be positive.")

        spectrum_mean = float(spectrum.mean())
        spectrum_std = float(spectrum.std())
        normalized_spectrum = (spectrum - spectrum_mean) / max(spectrum_std, 1e-8)
        model_conditions = np.concatenate(
            (
                conditions,
                np.asarray([spectrum_mean, spectrum_std], dtype=np.float32),
            )
        )
        normalized_conditions = (
            model_conditions - self.preprocessing["condition_mean"]
        ) / self.preprocessing["condition_std"]
        x = torch.as_tensor(
            normalized_spectrum[None, :], device=self.device, dtype=torch.float32
        )
        n = torch.as_tensor([noise_token], device=self.device, dtype=torch.float32)
        c = torch.as_tensor(
            normalized_conditions[None, :], device=self.device, dtype=torch.float32
        )
        context = self.model.encode(x, n, c)
        samples = self.model.from_latent(
            self.model.flow(context).sample((posterior_draws,))
        )
        molecule_samples = samples[:, 0, : len(MOLECULES)].cpu().numpy()
        medians, rank_probabilities, point_order = _rank_molecule_samples(
            molecule_samples
        )
        probability_above_detection_floor = np.mean(
            molecule_samples > DETECTION_FLOOR_LOG_VMR, axis=0
        )
        intervals = np.quantile(
            molecule_samples, [0.025, 0.16, 0.84, 0.975], axis=0
        )
        return {
            "posterior_draws": posterior_draws,
            "device": str(self.device),
            "detection_gating": "not applied; use offline calibrated evaluation",
            "point_ranking": [
                {
                    "rank": rank + 1,
                    "molecule": MOLECULES[molecule_index],
                    "posterior_median_log10_vmr": float(medians[molecule_index]),
                    "posterior_probability_above_detection_floor": float(
                        probability_above_detection_floor[molecule_index]
                    ),
                    "credible_interval_68_log10_vmr": [
                        float(intervals[1, molecule_index]),
                        float(intervals[2, molecule_index]),
                    ],
                    "credible_interval_95_log10_vmr": [
                        float(intervals[0, molecule_index]),
                        float(intervals[3, molecule_index]),
                    ],
                    "rank_probabilities": rank_probabilities[
                        molecule_index
                    ].tolist(),
                }
                for rank, molecule_index in enumerate(point_order[:5])
            ],
            "rank_probability_matrix": {
                molecule: rank_probabilities[molecule_index].tolist()
                for molecule_index, molecule in enumerate(MOLECULES)
            },
        }


def _physical_nll(
    model: FlowAtmos,
    spectra: np.ndarray,
    noise: np.ndarray,
    conditions: np.ndarray,
    targets: np.ndarray,
    indices: np.ndarray,
    device: torch.device,
    batch_size: int = 512,
) -> np.ndarray:
    results = np.empty(len(indices), dtype=np.float64)
    model.eval()
    with torch.no_grad():
        for start in range(0, len(indices), batch_size):
            chosen = indices[start : start + batch_size]
            x = torch.as_tensor(spectra[chosen], device=device, dtype=torch.float32)
            n = torch.as_tensor(noise[chosen], device=device, dtype=torch.float32)
            c = torch.as_tensor(conditions[chosen], device=device, dtype=torch.float32)
            y = torch.as_tensor(targets[chosen], device=device, dtype=torch.float32)
            results[start : start + len(chosen)] = (
                -model.log_prob(x, n, c, y).cpu().numpy()
            )
    return results


def _finite_sample_quantile(scores: np.ndarray, coverage: float) -> float:
    ordered = np.sort(np.asarray(scores, dtype=np.float64))
    rank = int(math.ceil((len(ordered) + 1) * coverage))
    if rank > len(ordered):
        return float("inf")
    return float(ordered[rank - 1])


def _uniform_prior_nll(lower: np.ndarray, upper: np.ndarray) -> float:
    widths = np.asarray(upper, dtype=np.float64) - np.asarray(lower, dtype=np.float64)
    if widths.size == 0 or np.any(~np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("Uniform-prior bounds must have finite, positive widths.")
    return float(np.log(widths).sum())


def _prior_physical_nll(
    model: UnconditionalPriorFlow,
    targets: np.ndarray,
    indices: np.ndarray,
    device: torch.device,
    batch_size: int = 512,
) -> np.ndarray:
    results = np.empty(len(indices), dtype=np.float64)
    model.eval()
    with torch.no_grad():
        for start in range(0, len(indices), batch_size):
            chosen = indices[start : start + batch_size]
            y = torch.as_tensor(targets[chosen], device=device, dtype=torch.float32)
            results[start : start + len(chosen)] = (
                -model.log_prob(y).cpu().numpy()
            )
    return results


def _bh_mask(p_values: np.ndarray, q: float = 0.05) -> np.ndarray:
    order = np.argsort(p_values)
    thresholds = q * np.arange(1, len(p_values) + 1) / len(p_values)
    passing = np.flatnonzero(p_values[order] <= thresholds)
    selected = np.zeros(len(p_values), dtype=bool)
    if passing.size:
        selected[order[: passing[-1] + 1]] = True
    return selected


def _global_ranking_metrics(
    truths: np.ndarray,
    medians: np.ndarray,
    selected: np.ndarray,
    point_order: np.ndarray,
) -> dict[str, float | int]:
    """Score ranks across all test spectra, counting undetected leaders as misses."""
    global_top1_correct: list[float] = []
    global_top1_detected: list[float] = []
    global_top5_overlap: list[float] = []
    global_top5_detection_recall: list[float] = []
    top1_correct_when_detected: list[float] = []
    top1_correct_given_any_detection: list[float] = []
    top5_overlap_given_any_detection: list[float] = []

    for row in range(len(truths)):
        true_top1 = int(np.argmax(truths[row]))
        true_top5 = set(np.argsort(truths[row])[::-1][:5])
        selected_ids = np.flatnonzero(selected[row])
        ranked_ids = point_order[row][point_order[row] >= 0]
        top1_detected = true_top1 in selected_ids
        top1_predicted_correctly = bool(ranked_ids.size and ranked_ids[0] == true_top1)
        global_top1_detected.append(float(top1_detected))
        global_top1_correct.append(float(top1_predicted_correctly))

        predicted_top5 = set(ranked_ids[:5])
        global_top5_overlap.append(len(true_top5 & predicted_top5) / 5.0)
        global_top5_detection_recall.append(len(true_top5 & set(selected_ids)) / 5.0)

        if selected_ids.size:
            true_top_among_detected = selected_ids[
                np.argmax(truths[row, selected_ids])
            ]
            predicted_top_among_detected = selected_ids[
                np.argmax(medians[row, selected_ids])
            ]
            correct_among_detected = (
                true_top_among_detected == predicted_top_among_detected
            )
            top1_correct_given_any_detection.append(float(correct_among_detected))
            top5_overlap_given_any_detection.append(
                len(true_top5 & predicted_top5) / 5.0
            )
            if top1_detected:
                top1_correct_when_detected.append(
                    float(top1_predicted_correctly)
                )

    return {
        "global_top1_accuracy": float(np.mean(global_top1_correct)),
        "global_top1_accuracy_when_detected": (
            float(np.mean(top1_correct_when_detected))
            if top1_correct_when_detected
            else float("nan")
        ),
        "global_top5_overlap": float(np.mean(global_top5_overlap)),
        "global_top1_detection_rate": float(np.mean(global_top1_detected)),
        "mean_global_top5_detection_recall": float(
            np.mean(global_top5_detection_recall)
        ),
        "mean_top5_overlap_given_any_detection": (
            float(np.mean(top5_overlap_given_any_detection))
            if top5_overlap_given_any_detection
            else float("nan")
        ),
        "top1_accuracy_among_detected_species": (
            float(np.mean(top1_correct_given_any_detection))
            if top1_correct_given_any_detection
            else float("nan")
        ),
        "spectra_with_at_least_one_detection": int(
            len(top1_correct_given_any_detection)
        ),
        "test_spectra": int(len(truths)),
    }


def _stratified_ranking_metrics(
    truths: np.ndarray,
    medians: np.ndarray,
    selected: np.ndarray,
    point_order: np.ndarray,
    active_species_count: np.ndarray,
) -> dict[str, dict[str, float | int]]:
    """Report top-k ranking quality separately for each active-species count."""
    if not (
        len(truths)
        == len(medians)
        == len(selected)
        == len(point_order)
        == len(active_species_count)
    ):
        raise ValueError("Scenario labels and ranking arrays must have matching lengths.")
    summaries: dict[str, dict[str, float | int]] = {}
    for count in sorted(np.unique(active_species_count).tolist()):
        mask = active_species_count == count
        summaries[str(int(count))] = _global_ranking_metrics(
            truths[mask], medians[mask], selected[mask], point_order[mask]
        )
    return summaries


def _reliability_summary(
    predicted: np.ndarray, observed: np.ndarray, n_bins: int = 10
) -> dict[str, Any]:
    finite = np.isfinite(predicted) & np.isfinite(observed)
    predicted = predicted[finite]
    observed = observed[finite]
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bins: list[dict[str, float | int]] = []
    ece = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (predicted >= left) & (
            (predicted < right) if right < 1.0 else (predicted <= right)
        )
        count = int(mask.sum())
        mean_probability = float(predicted[mask].mean()) if count else float("nan")
        observed_frequency = float(observed[mask].mean()) if count else float("nan")
        if count:
            ece += count / max(1, len(predicted)) * abs(
                mean_probability - observed_frequency
            )
        bins.append(
            {
                "lower": float(left),
                "upper": float(right),
                "count": count,
                "mean_probability": mean_probability,
                "observed_frequency": observed_frequency,
            }
        )
    brier = float(np.mean((predicted - observed) ** 2)) if len(predicted) else float("nan")
    return {"brier_score": brier, "expected_calibration_error": float(ece), "bins": bins}


def evaluate_flowatmos(
    dataset_path: str | Path,
    output_dir: str | Path,
    posterior_draws: int = 5000,
    seed: int = 42,
    posterior_batch_size: int = 8,
    sampling_device: str = "cpu",
    cpu_threads: int = 4,
) -> dict[str, Any]:
    """Evaluate held-out posterior, conformal, detection, and ranking behavior."""
    if posterior_draws < 5000:
        raise ValueError("At least 5,000 posterior draws are required for test ranking.")
    if posterior_batch_size < 1:
        raise ValueError("posterior_batch_size must be positive.")
    if sampling_device not in {"cpu", "mps", "cuda"}:
        raise ValueError("sampling_device must be one of: cpu, mps, cuda.")
    if cpu_threads < 1:
        raise ValueError("cpu_threads must be positive.")
    if sampling_device == "mps" and not (
        hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
    ):
        raise RuntimeError("MPS posterior sampling was requested but is unavailable.")
    if sampling_device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA posterior sampling was requested but is unavailable.")
    output = Path(output_dir).expanduser().resolve()
    dataset_path = Path(dataset_path).expanduser().resolve()
    with np.load(dataset_path) as saved:
        dataset = {name: saved[name] for name in saved.files}
    with np.load(output / "splits.npz") as saved_splits:
        splits = {name: saved_splits[name] for name in saved_splits.files}
    checkpoint = torch.load(output / "flowatmos_best.pt", map_location="cpu", weights_only=True)
    device = torch.device(sampling_device)
    if device.type == "cpu":
        torch.set_num_threads(cpu_threads)
    model, preprocessing = _load_model(dataset, checkpoint, device)
    spectra, noise, conditions = _normalised_arrays(dataset, preprocessing)
    targets = dataset["targets"].astype(np.float32)
    validation_indices = splits["validation"]
    calibration_indices = splits["calibration"]
    test_indices = splits["test"]

    print("Scoring validation, calibration, and test NLL.")
    validation_nll = _physical_nll(
        model, spectra, noise, conditions, targets, validation_indices, device
    )
    calibration_nll = _physical_nll(
        model, spectra, noise, conditions, targets, calibration_indices, device
    )
    test_nll = _physical_nll(
        model, spectra, noise, conditions, targets, test_indices, device
    )
    lower = dataset["lower_bounds"].astype(np.float64)
    upper = dataset["upper_bounds"].astype(np.float64)
    prior_nll = _uniform_prior_nll(lower, upper)
    manifest_path = dataset_path.parent / "manifest.json"
    generation = (
        json.loads(manifest_path.read_text(encoding="utf-8")).get("generation", {})
        if manifest_path.exists()
        else {}
    )
    scenario_strategy = generation.get("scenario_strategy", "uniform_log")
    prior_baseline_comparable = scenario_strategy == "uniform_log"
    prior_checkpoint_path = dataset_path.parent / "flowatmos_prior_baseline.pt"
    prior_validation_nll: float | None = None
    prior_test_nll: float | None = None
    prior_baseline_method = "unavailable"
    if prior_checkpoint_path.is_file():
        prior_checkpoint = torch.load(
            prior_checkpoint_path, map_location="cpu", weights_only=True
        )
        prior_model = UnconditionalPriorFlow(
            prior_checkpoint["lower_bounds"].cpu().numpy(),
            prior_checkpoint["upper_bounds"].cpu().numpy(),
        ).to(device=device, dtype=torch.float32)
        prior_model.load_state_dict(prior_checkpoint["state_dict"])
        prior_validation_scores = _prior_physical_nll(
            prior_model,
            targets,
            splits["validation"],
            device,
        )
        prior_test_scores = _prior_physical_nll(
            prior_model,
            targets,
            test_indices,
            device,
        )
        prior_validation_nll = float(prior_validation_scores.mean())
        prior_test_nll = float(prior_test_scores.mean())
        prior_baseline_method = "unconditional target-only spline flow trained on train split"
    elif prior_baseline_comparable:
        prior_validation_nll = prior_nll
        prior_test_nll = prior_nll
        prior_baseline_method = "analytic uniform box prior"
    conformal: dict[str, Any] = {}
    for coverage in (0.68, 0.95):
        quantile = _finite_sample_quantile(calibration_nll, coverage)
        conformal[f"{int(coverage * 100)}%"] = {
            "calibration_quantile_nll": quantile,
            "test_coverage": float(np.mean(test_nll <= quantile)),
            "calibration_count": int(len(calibration_nll)),
        }

    print("Drawing calibration posteriors for null-score thresholds.")
    cal_draws = _posterior_draws(
        model, spectra, noise, conditions, calibration_indices,
        n_draws=min(512, posterior_draws), device=device,
        batch_size=posterior_batch_size,
    )
    print("Drawing test posteriors for intervals, SBC, and ranking.")
    test_draws = _posterior_draws(
        model, spectra, noise, conditions, test_indices,
        n_draws=posterior_draws, device=device,
        batch_size=posterior_batch_size,
    )
    molecule_count = len(MOLECULES)
    cal_abundances = cal_draws[:, :, :molecule_count]
    test_abundances = test_draws[:, :, :molecule_count]
    calibration_truth = targets[calibration_indices, :molecule_count]
    test_truth = targets[test_indices]

    null_scores: list[np.ndarray] = []
    for molecule_index in range(molecule_count):
        is_null = calibration_truth[:, molecule_index] <= DETECTION_FLOOR_LOG_VMR
        if not np.any(is_null):
            raise ValueError(f"No null calibration spectra for {MOLECULES[molecule_index]}.")
        null_scores.append(np.median(cal_abundances[is_null, :, molecule_index], axis=1))

    median = np.median(test_draws, axis=1)
    intervals: dict[str, Any] = {}
    for coverage in (0.68, 0.95):
        tail = (1.0 - coverage) / 2.0
        quantiles = np.quantile(test_draws, [tail, 0.5, 1.0 - tail], axis=1)
        covered = (test_truth >= quantiles[0]) & (test_truth <= quantiles[2])
        intervals[f"{int(coverage * 100)}%"] = {
            "coverage_by_parameter": {
                name: float(value)
                for name, value in zip(TARGET_NAMES, covered.mean(axis=0))
            },
            "mean_coverage": float(covered.mean()),
        }

    ranking_probabilities = np.full(
        (len(test_indices), molecule_count, molecule_count), np.nan, dtype=np.float32
    )
    selected_mask = np.zeros((len(test_indices), molecule_count), dtype=bool)
    detection_p = np.ones((len(test_indices), molecule_count), dtype=np.float64)
    point_order = np.full((len(test_indices), molecule_count), -1, dtype=np.int16)
    taus: list[float] = []
    spearmans: list[float] = []
    reliability_probabilities: list[np.ndarray] = []
    reliability_outcomes: list[np.ndarray] = []

    for row in range(len(test_indices)):
        scores = np.median(test_abundances[row], axis=0)
        p_values = np.asarray(
            [
                (1 + np.count_nonzero(null_scores[molecule] >= scores[molecule]))
                / (len(null_scores[molecule]) + 1)
                for molecule in range(molecule_count)
            ],
            dtype=np.float64,
        )
        selected = _bh_mask(p_values, q=0.05)
        selected_mask[row] = selected
        detection_p[row] = p_values
        selected_ids = np.flatnonzero(selected)
        if not selected_ids.size:
            continue

        point_order[row, : len(selected_ids)] = selected_ids[
            np.argsort(scores[selected_ids])[::-1]
        ]
        posterior = test_abundances[row][:, selected_ids]
        ranks = np.argsort(np.argsort(-posterior, axis=1), axis=1)
        probabilities = np.zeros((len(selected_ids), molecule_count), dtype=np.float32)
        for molecule_position in range(len(selected_ids)):
            probabilities[molecule_position] = np.bincount(
                ranks[:, molecule_position], minlength=molecule_count
            ) / posterior_draws
            ranking_probabilities[row, selected_ids[molecule_position]] = probabilities[
                molecule_position
            ]

        truth_order = np.argsort(np.argsort(-test_truth[row, selected_ids]))
        predicted_order = np.argsort(np.argsort(-scores[selected_ids]))
        if len(selected_ids) >= 2:
            tau = float(kendalltau(truth_order, predicted_order).statistic)
            rho = float(spearmanr(truth_order, predicted_order).statistic)
            if np.isfinite(tau):
                taus.append(tau)
            if np.isfinite(rho):
                spearmans.append(rho)
        actual_rank = np.argsort(np.argsort(-test_truth[row, selected_ids]))
        reliability_probabilities.append(probabilities.reshape(-1))
        reliability_outcomes.append(
            np.eye(molecule_count, dtype=np.float32)[actual_rank].reshape(-1)
        )

    sbc: dict[str, Any] = {}
    rank_histograms = np.zeros((len(TARGET_NAMES), 30), dtype=np.int64)
    for parameter, name in enumerate(TARGET_NAMES):
        ranks = np.mean(
            test_draws[:, :, parameter] < test_truth[:, None, parameter], axis=1
        )
        counts, _ = np.histogram(ranks, bins=np.linspace(0.0, 1.0, 31))
        rank_histograms[parameter] = counts
        statistic, p_value = chisquare(counts)
        sbc[name] = {
            "bin_counts": counts.tolist(),
            "chi_square_statistic": float(statistic),
            "chi_square_p_value": float(p_value),
        }

    if reliability_probabilities:
        rank_reliability = _reliability_summary(
            np.concatenate(reliability_probabilities),
            np.concatenate(reliability_outcomes),
        )
    else:
        rank_reliability = _reliability_summary(
            np.asarray([], dtype=np.float32), np.asarray([], dtype=np.float32)
        )

    residual = test_truth - median
    target_variance = np.sum((test_truth - test_truth.mean(axis=0)) ** 2, axis=0)
    r2 = 1.0 - np.sum(residual**2, axis=0) / np.maximum(target_variance, 1e-12)
    rmse = np.sqrt(np.mean(residual**2, axis=0))
    parameter_metrics = {
        name: {"r2": float(r2[i]), "rmse": float(rmse[i])}
        for i, name in enumerate(TARGET_NAMES)
    }
    all_spectrum_ranking = _global_ranking_metrics(
        test_truth[:, :molecule_count],
        median[:, :molecule_count],
        selected_mask,
        point_order,
    )
    active_species_count = dataset.get(
        "active_species_count",
        (targets[:, :molecule_count] > DETECTION_FLOOR_LOG_VMR).sum(axis=1),
    )[test_indices]
    scenario_ranking = _stratified_ranking_metrics(
        test_truth[:, :molecule_count],
        median[:, :molecule_count],
        selected_mask,
        point_order,
        active_species_count,
    )
    mean_r2 = float(np.mean(r2))
    validation_nll_better_than_prior = (
        bool(validation_nll.mean() < prior_validation_nll)
        if prior_validation_nll is not None
        else None
    )
    results = {
        "evaluation_type": "synthetic held-out CERBERUS simulations; not independent observed-planet validation",
        "dataset_scenario_strategy": scenario_strategy,
        "seed": seed,
        "device": str(device),
        "posterior_draws_per_spectrum": posterior_draws,
        "test_count": int(len(test_indices)),
        "calibration_count": int(len(calibration_indices)),
        "prior_nll_baseline": prior_validation_nll,
        "prior_test_nll_baseline": prior_test_nll,
        "prior_baseline_method": prior_baseline_method,
        "uniform_box_reference_nll": prior_nll,
        "uniform_box_reference_comparable": prior_baseline_comparable,
        "mean_test_nll": float(test_nll.mean()),
        "mean_validation_nll": float(validation_nll.mean()),
        "checkpoint_validation_nll": float(checkpoint["validation_nll"]),
        "validation_nll_gain_over_prior": (
            float(prior_validation_nll - validation_nll.mean())
            if prior_validation_nll is not None
            else None
        ),
        "test_nll_gain_over_prior": (
            float(prior_test_nll - test_nll.mean())
            if prior_test_nll is not None
            else None
        ),
        "validation_nll_better_than_prior": validation_nll_better_than_prior,
        "mean_parameter_r2": mean_r2,
        "mean_parameter_r2_above_zero": bool(mean_r2 > 0.0),
        "sanity_checks_passed": bool(
            validation_nll_better_than_prior is True and mean_r2 > 0.0
        ),
        "sanity_check_status": (
            "passed"
            if validation_nll_better_than_prior is True and mean_r2 > 0.0
            else (
                "matched prior baseline unavailable"
                if prior_validation_nll is None
                else "failed"
            )
        ),
        "parameter_metrics": parameter_metrics,
        "credible_interval_coverage": intervals,
        "split_conformal_nll": conformal,
        "simulation_based_calibration": {
            "bins": 30,
            "parameters": sbc,
            "chi_square_tests_with_p_below_0.05": int(
                sum(item["chi_square_p_value"] < 0.05 for item in sbc.values())
            ),
        },
        "detection": {
            "method": "empirical low-abundance null calibration and per-spectrum Benjamini-Hochberg FDR",
            "null_threshold_log10_vmr": DETECTION_FLOOR_LOG_VMR,
            "fdr_q": 0.05,
            "mean_detected_species": float(selected_mask.sum(axis=1).mean()),
            "fraction_with_no_detections": float(np.mean(selected_mask.sum(axis=1) == 0)),
            "true_global_top1_detection_rate": all_spectrum_ranking[
                "global_top1_detection_rate"
            ],
            "mean_true_top5_detection_recall": all_spectrum_ranking[
                "mean_global_top5_detection_recall"
            ],
            "null_calibration_counts": {
                molecule: int(len(null_scores[i]))
                for i, molecule in enumerate(MOLECULES)
            },
        },
        "ranking_quality_among_detected_species": {
            "mean_kendall_tau": float(np.mean(taus)) if taus else float("nan"),
            "mean_spearman_correlation": float(np.mean(spearmans)) if spearmans else float("nan"),
            "top1_accuracy_among_detected_species": all_spectrum_ranking[
                "top1_accuracy_among_detected_species"
            ],
            "mean_top5_overlap_given_any_detection": all_spectrum_ranking[
                "mean_top5_overlap_given_any_detection"
            ],
            "spectra_with_at_least_one_detection": all_spectrum_ranking[
                "spectra_with_at_least_one_detection"
            ],
            "rank_probability_reliability": rank_reliability,
        },
        "ranking_quality_over_all_spectra": all_spectrum_ranking,
        "ranking_quality_by_active_species_count": scenario_ranking,
    }
    _write_json(output / "evaluation_metrics.json", results)
    np.savez_compressed(
        output / "test_posterior_summary.npz",
        test_indices=test_indices,
        true_targets=test_truth,
        posterior_medians=median,
        selected_mask=selected_mask,
        detection_p_values=detection_p,
        point_order=point_order,
        rank_probabilities=ranking_probabilities,
        sbc_rank_histograms=rank_histograms,
    )
    _write_ranking_csvs(
        output, test_indices, test_truth, median, selected_mask,
        detection_p, point_order, ranking_probabilities,
    )
    return results


def _write_ranking_csvs(
    output: Path,
    test_indices: np.ndarray,
    truths: np.ndarray,
    medians: np.ndarray,
    selected: np.ndarray,
    p_values: np.ndarray,
    point_order: np.ndarray,
    rank_probabilities: np.ndarray,
) -> None:
    with (output / "test_rankings.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["spectrum_id", "molecule", "posterior_median_log10_vmr", "detection_p_value", "selected", "point_rank", "in_reported_top5", "true_log10_vmr"]
        )
        for row, sample_index in enumerate(test_indices):
            ranks = {
                molecule_index: rank + 1
                for rank, molecule_index in enumerate(point_order[row])
                if molecule_index >= 0
            }
            for molecule_index, molecule in enumerate(MOLECULES):
                writer.writerow(
                    [
                        int(sample_index),
                        molecule,
                        float(medians[row, molecule_index]),
                        float(p_values[row, molecule_index]),
                        bool(selected[row, molecule_index]),
                        ranks.get(molecule_index, "below detection"),
                        bool(ranks.get(molecule_index, 6) <= 5),
                        float(truths[row, molecule_index]),
                    ]
                )

    with (output / "rank_probability.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["spectrum_id", "molecule", *[f"rank_{rank}_probability" for rank in range(1, len(MOLECULES) + 1)]]
        )
        for row, sample_index in enumerate(test_indices):
            for molecule_index, molecule in enumerate(MOLECULES):
                probabilities = rank_probabilities[row, molecule_index]
                writer.writerow(
                    [
                        int(sample_index),
                        molecule,
                        *[
                            float(probability) if np.isfinite(probability) else "below detection"
                            for probability in probabilities
                        ],
                    ]
                )
