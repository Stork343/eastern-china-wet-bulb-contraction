#!/usr/bin/env python3
"""Size stress test for cross-record dependence in the product-shift test.

This script targets a departure that is not covered by the product-invariance
proposition.  It uses the 121 fixed ERA5-Land coordinates and the five fixed
Gaussian graph operators.  Each simulated panel contains June--August records
from two summers.  Within a summer, all three records share the phase of a
latent Gaussian process on the unit circle.  Each record remains marginally
cyclic stationary, but independently rotating its outcome sequence changes
the cross-record covariance whenever the shared loading is positive.

The regional-mean series used for type-7 quartile labels are generated from
separate random-number streams and are independent of the spatial fields.
Only lower-tail null rejection is studied.  The defaults are the manuscript
design: 2,000 replications per loading, 999 product shifts and alpha=0.05.
Environment variables ESH_CROSS_RECORD_N, ESH_CROSS_RECORD_B and
ESH_CROSS_RECORD_BATCH may be used for development runs.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_DIR = SCRIPT_PATH.parent.parent
OUTPUT_DIR = PROJECT_DIR / "output_extension_methods"
REFERENCE_PATH = (
    PROJECT_DIR / "data" / "era5_consistent" / "era5_201506_consistent.csv"
)

SEED = 20260812
REPLICATIONS = int(os.environ.get("ESH_CROSS_RECORD_N", "2000"))
PRODUCT_SHIFTS = int(os.environ.get("ESH_CROSS_RECORD_B", "999"))
BATCH_SIZE = int(os.environ.get("ESH_CROSS_RECORD_BATCH", "100"))
ALPHA = 0.05
LAMBDAS = np.array([0.0, 0.3, 0.6, 0.8], dtype=float)
H_FACTORS = np.array([0.125, 0.25, 0.5, 1.0, 2.0], dtype=float)
RECORD_LENGTHS = np.array([30, 31, 31, 30, 31, 31], dtype=int)
RECORD_MONTHS = np.array([6, 7, 8, 6, 7, 8], dtype=int)
RECORD_SUMMERS = np.array([1, 1, 1, 2, 2, 2], dtype=int)
SPATIAL_NOISE_RANK = 12
SPATIAL_RANGE_KM = 450.0
SPATIAL_NUGGET = 0.05

# The weights below define a unit-variance Gaussian process on normalized
# circular phase.  Sharing its sine/cosine coefficients, rather than a raw
# day index, retains exact marginal cyclic stationarity for both 30- and
# 31-day records.
FOURIER_VARIANCE_WEIGHTS = np.array([1.0, 0.45, 0.20], dtype=float)
FOURIER_VARIANCE_WEIGHTS /= FOURIER_VARIANCE_WEIGHTS.sum()

# Circular moving-average filters for the independent regional-mean and
# record-specific spatial-score series.  Dividing by the Euclidean norm makes
# the output variance one at every day.
MEAN_FILTER = np.array([0.25, 0.60, 1.00, 0.60, 0.25], dtype=float)
MEAN_FILTER /= np.linalg.norm(MEAN_FILTER)
NOISE_FILTER = np.array([0.20, 0.50, 1.00, 0.50, 0.20], dtype=float)
NOISE_FILTER /= np.linalg.norm(NOISE_FILTER)
FILTER_OFFSETS = np.arange(-2, 3, dtype=int)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def projected_coordinates_km(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lat0 = np.deg2rad(lat.mean())
    return np.column_stack((lon * 111.32 * np.cos(lat0), lat * 110.57))


@dataclass(frozen=True)
class GraphSystem:
    bandwidths_km: np.ndarray
    laplacians: np.ndarray
    divisors: np.ndarray
    gradient: np.ndarray
    noise_basis: np.ndarray
    noise_quadratics: np.ndarray
    noise_gradient_cross: np.ndarray
    gradient_quadratic: np.ndarray
    distances_km: np.ndarray


def build_graph_system() -> tuple[GraphSystem, pd.DataFrame]:
    reference = pd.read_csv(
        REFERENCE_PATH, usecols=["site_id", "lon", "lat"]
    ).drop_duplicates().sort_values("site_id").reset_index(drop=True)
    if len(reference) != 121 or reference.site_id.nunique() != 121:
        raise RuntimeError("Expected 121 unique fixed reference sites")

    coordinates = projected_coordinates_km(
        reference.lon.to_numpy(), reference.lat.to_numpy()
    )
    delta = coordinates[:, None, :] - coordinates[None, :, :]
    distance = np.sqrt(np.sum(delta * delta, axis=2))
    median_distance = float(np.median(distance[np.tril_indices(121, -1)]))
    bandwidths = H_FACTORS * median_distance

    laplacians = []
    divisors = []
    for bandwidth in bandwidths:
        weights = np.exp(-(distance * distance) / (2 * bandwidth * bandwidth))
        np.fill_diagonal(weights, 0.0)
        laplacians.append(np.diag(weights.sum(axis=1)) - weights)
        divisors.append(2.0 * weights[np.triu_indices(121, 1)].sum())
    laplacians_array = np.stack(laplacians)
    divisors_array = np.asarray(divisors)

    raw_gradient = 0.8 * coordinates[:, 0] + coordinates[:, 1]
    gradient = raw_gradient - raw_gradient.mean()
    gradient /= math.sqrt(float(np.mean(gradient * gradient)))

    # A deliberately specified low-rank spatial Gaussian noise process keeps
    # the full 8,000-cell experiment tractable.  The basis comprises the first
    # 12 centred eigenmodes of the same exponential covariance used in the
    # existing field simulations.  It is rescaled so the mean site variance
    # is one; it is part of this stress-test DGP, not a numerical approximation
    # applied after simulation.
    covariance = (
        (1.0 - SPATIAL_NUGGET) * np.exp(-distance / SPATIAL_RANGE_KM)
        + SPATIAL_NUGGET * np.eye(121)
    )
    centering = np.eye(121) - np.ones((121, 121)) / 121.0
    covariance = centering @ covariance @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1][:SPATIAL_NOISE_RANK]
    if np.min(eigenvalues[order]) <= 0:
        raise RuntimeError("Centred spatial covariance has insufficient rank")
    basis = eigenvectors[:, order] * np.sqrt(eigenvalues[order])[None, :]
    basis *= math.sqrt(121.0 / float(np.sum(basis * basis)))

    noise_quadratics = np.empty(
        (len(H_FACTORS), SPATIAL_NOISE_RANK, SPATIAL_NOISE_RANK), dtype=float
    )
    noise_gradient_cross = np.empty(
        (len(H_FACTORS), SPATIAL_NOISE_RANK), dtype=float
    )
    gradient_quadratic = np.empty(len(H_FACTORS), dtype=float)
    for h, (laplacian, divisor) in enumerate(
        zip(laplacians_array, divisors_array, strict=True)
    ):
        noise_quadratics[h] = basis.T @ laplacian @ basis / divisor
        noise_gradient_cross[h] = basis.T @ laplacian @ gradient / divisor
        gradient_quadratic[h] = gradient @ laplacian @ gradient / divisor

    graph = GraphSystem(
        bandwidths_km=bandwidths,
        laplacians=laplacians_array,
        divisors=divisors_array,
        gradient=gradient,
        noise_basis=basis,
        noise_quadratics=noise_quadratics,
        noise_gradient_cross=noise_gradient_cross,
        gradient_quadratic=gradient_quadratic,
        distances_km=distance,
    )
    return graph, reference


def circular_filter(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Apply the fixed unit-variance moving average along axis 1."""
    result = np.zeros_like(values)
    for offset, weight in zip(FILTER_OFFSETS, weights, strict=True):
        result += weight * np.roll(values, shift=offset, axis=1)
    return result


def fourier_design(n_days: int) -> tuple[np.ndarray, np.ndarray]:
    phase = np.arange(n_days, dtype=float) / n_days
    harmonics = np.arange(1, len(FOURIER_VARIANCE_WEIGHTS) + 1, dtype=float)
    angle = 2.0 * np.pi * np.outer(phase, harmonics)
    scale = np.sqrt(FOURIER_VARIANCE_WEIGHTS)[None, :]
    return np.cos(angle) * scale, np.sin(angle) * scale


FOURIER_DESIGNS = {
    n: fourier_design(int(n)) for n in np.unique(RECORD_LENGTHS)
}


def shared_process(
    cosine_coefficients: np.ndarray,
    sine_coefficients: np.ndarray,
    n_days: int,
) -> np.ndarray:
    cosine, sine = FOURIER_DESIGNS[n_days]
    return cosine_coefficients @ cosine.T + sine_coefficients @ sine.T


def classify_type7(regional_mean: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    quartiles = np.quantile(
        regional_mean, [0.25, 0.75], axis=1, method="linear"
    )
    lower = quartiles[0, :, None]
    upper = quartiles[1, :, None]
    high = regional_mean >= upper
    middle = (regional_mean >= lower) & (regional_mean < upper)
    if np.any(high.sum(axis=1) < 2) or np.any(middle.sum(axis=1) < 2):
        raise RuntimeError("Unexpected empty type-7 quartile group")
    return high, middle


def component_metrics(
    spatial_scores: np.ndarray, shared: np.ndarray, graph: GraphSystem
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return epsilon energy, epsilon-gradient cross term and U^2 energy."""
    noise_energy = np.einsum(
        "bti,hij,btj->bth",
        spatial_scores,
        graph.noise_quadratics,
        spatial_scores,
        optimize=True,
    )
    noise_gradient = np.einsum(
        "bti,hi->bth",
        spatial_scores,
        graph.noise_gradient_cross,
        optimize=True,
    )
    shared_energy = shared[:, :, None] ** 2 * graph.gradient_quadratic
    return noise_energy, noise_gradient, shared_energy


def effect_lookup(
    metrics: np.ndarray, high: np.ndarray, middle: np.ndarray
) -> np.ndarray:
    """Record-level five-scale effects for every circular offset."""
    batch, n_days, n_scales = metrics.shape
    lookup = np.empty((batch, n_days, n_scales), dtype=float)
    high_count = high.sum(axis=1)[:, None]
    middle_count = middle.sum(axis=1)[:, None]
    for offset in range(n_days):
        shifted = np.roll(metrics, shift=-offset, axis=1)
        high_mean = (shifted * high[:, :, None]).sum(axis=1) / high_count
        middle_mean = (
            (shifted * middle[:, :, None]).sum(axis=1) / middle_count
        )
        lookup[:, offset, :] = high_mean / middle_mean - 1.0
    return lookup


class OnlineCorrelation:
    def __init__(self) -> None:
        self.n = 0
        self.sum_x = 0.0
        self.sum_y = 0.0
        self.sum_xx = 0.0
        self.sum_yy = 0.0
        self.sum_xy = 0.0

    def update(self, x: np.ndarray, y: np.ndarray) -> None:
        x = np.asarray(x, dtype=float).ravel()
        y = np.asarray(y, dtype=float).ravel()
        if x.size != y.size:
            raise RuntimeError("Correlation inputs have different lengths")
        self.n += int(x.size)
        self.sum_x += float(x.sum())
        self.sum_y += float(y.sum())
        self.sum_xx += float(x @ x)
        self.sum_yy += float(y @ y)
        self.sum_xy += float(x @ y)

    def correlation(self) -> float:
        covariance = self.sum_xy - self.sum_x * self.sum_y / self.n
        variance_x = self.sum_xx - self.sum_x * self.sum_x / self.n
        variance_y = self.sum_yy - self.sum_y * self.sum_y / self.n
        return covariance / math.sqrt(variance_x * variance_y)


def covariance_checks(graph: GraphSystem) -> dict[str, object]:
    checks: dict[str, object] = {}
    maximum_marginal_error = 0.0
    marginal_details = {}
    for n_days in np.unique(RECORD_LENGTHS):
        cosine, sine = FOURIER_DESIGNS[int(n_days)]
        covariance = cosine @ cosine.T + sine @ sine.T
        shifted = np.roll(np.roll(covariance, 1, axis=0), 1, axis=1)
        error = float(np.max(np.abs(shifted - covariance)))
        maximum_marginal_error = max(maximum_marginal_error, error)
        marginal_details[str(int(n_days))] = error

    cosine_30, sine_30 = FOURIER_DESIGNS[30]
    cosine_31, sine_31 = FOURIER_DESIGNS[31]
    cross_covariance = cosine_30 @ cosine_31.T + sine_30 @ sine_31.T
    independently_shifted = np.roll(cross_covariance, 1, axis=0)
    difference = independently_shifted - cross_covariance
    cross_relative_change = float(
        np.linalg.norm(difference) / np.linalg.norm(cross_covariance)
    )
    # For a centred Gaussian shared process, the cross-covariance between
    # graph-dispersion contributions is proportional to
    # Cov(U_r(t)^2, U_s(u)^2) = 2 Cov(U_r(t), U_s(u))^2.  This check therefore
    # targets the dispersion collection D in Proposition 1 directly.
    energy_cross_covariance = 2.0 * cross_covariance * cross_covariance
    shifted_energy_cross_covariance = np.roll(
        energy_cross_covariance, 1, axis=0
    )
    energy_cross_relative_change = float(
        np.linalg.norm(
            shifted_energy_cross_covariance - energy_cross_covariance
        )
        / np.linalg.norm(energy_cross_covariance)
    )

    lambda_checks = {}
    for loading in LAMBDAS:
        lambda_checks[f"{loading:.1f}"] = {
            "shared_variance_fraction_spatial_average": float(loading**2),
            "max_abs_cross_covariance_change_after_one_independent_shift": (
                float(loading**2 * np.max(np.abs(difference)))
            ),
            "frobenius_cross_covariance_change": float(
                loading**2 * np.linalg.norm(difference)
            ),
        }

    checks["marginal_shared_process_one_step_rotation_max_abs_error"] = (
        maximum_marginal_error
    )
    checks["marginal_error_by_record_length"] = marginal_details
    checks["unscaled_30_by_31_cross_covariance_relative_change_after_one_30_day_shift"] = (
        cross_relative_change
    )
    checks[
        "unscaled_30_by_31_graph_dispersion_cross_covariance_relative_change_after_one_30_day_shift"
    ] = energy_cross_relative_change
    checks["joint_product_invariance_departure_by_lambda"] = lambda_checks
    checks["laplacian_constant_annihilation_max_abs"] = float(
        np.max(np.abs(graph.laplacians @ np.ones(121)))
    )
    checks["gradient_spatial_mean"] = float(graph.gradient.mean())
    checks["gradient_spatial_mean_square"] = float(
        np.mean(graph.gradient * graph.gradient)
    )
    checks["noise_basis_spatial_mean_max_abs"] = float(
        np.max(np.abs(graph.noise_basis.mean(axis=0)))
    )
    checks["noise_process_spatial_average_variance"] = float(
        np.sum(graph.noise_basis * graph.noise_basis) / 121.0
    )
    return checks


def run_simulation(graph: GraphSystem) -> tuple[list[dict[str, object]], dict]:
    if REPLICATIONS <= 0 or PRODUCT_SHIFTS <= 0 or BATCH_SIZE <= 0:
        raise ValueError("Replication, shift and batch counts must be positive")

    field_seed, label_seed, shift_seed = np.random.SeedSequence(SEED).spawn(3)
    field_rng = np.random.default_rng(field_seed)
    label_rng = np.random.default_rng(label_seed)
    shift_rng = np.random.default_rng(shift_seed)

    reject_counts = np.zeros(len(LAMBDAS), dtype=int)
    p_value_sums = np.zeros(len(LAMBDAS), dtype=float)
    estimate_sums = np.zeros(len(LAMBDAS), dtype=float)
    dependence_correlations = [OnlineCorrelation() for _ in LAMBDAS]
    label_counts: dict[str, set[tuple[int, int]]] = {"30": set(), "31": set()}
    completed = 0

    while completed < REPLICATIONS:
        batch = min(BATCH_SIZE, REPLICATIONS - completed)
        # One isotropic sine/cosine coefficient pair per summer is shared by
        # its June, July and August records.
        shared_cosine = field_rng.standard_normal(
            (batch, 2, len(FOURIER_VARIANCE_WEIGHTS))
        )
        shared_sine = field_rng.standard_normal(
            (batch, 2, len(FOURIER_VARIANCE_WEIGHTS))
        )

        records = []
        offsets = []
        for record, n_days_value in enumerate(RECORD_LENGTHS):
            n_days = int(n_days_value)
            summer = int(RECORD_SUMMERS[record] - 1)
            shared = shared_process(
                shared_cosine[:, summer, :], shared_sine[:, summer, :], n_days
            )

            raw_scores = field_rng.standard_normal(
                (batch, n_days, SPATIAL_NOISE_RANK)
            )
            spatial_scores = circular_filter(raw_scores, NOISE_FILTER)
            noise_energy, noise_gradient, shared_energy = component_metrics(
                spatial_scores, shared, graph
            )

            raw_mean = label_rng.standard_normal((batch, n_days))
            regional_mean = circular_filter(raw_mean, MEAN_FILTER)
            high, middle = classify_type7(regional_mean)
            for high_n, middle_n in zip(
                high.sum(axis=1), middle.sum(axis=1), strict=True
            ):
                label_counts[str(n_days)].add((int(high_n), int(middle_n)))

            records.append({
                "shared": shared,
                "noise_energy": noise_energy,
                "noise_gradient": noise_gradient,
                "shared_energy": shared_energy,
                "high": high,
                "middle": middle,
            })
            offsets.append(
                shift_rng.integers(
                    0, n_days, size=(batch, PRODUCT_SHIFTS), endpoint=False
                )
            )

        for lambda_index, loading in enumerate(LAMBDAS):
            noise_weight = math.sqrt(1.0 - loading * loading)
            observed = np.zeros((batch, len(H_FACTORS)), dtype=float)
            shifted = np.zeros(
                (batch, PRODUCT_SHIFTS, len(H_FACTORS)), dtype=float
            )
            record_metrics = []
            for record_index, record in enumerate(records):
                exact_cross = (
                    record["shared"][:, :, None]
                    * record["noise_gradient"]
                )
                metrics = (
                    (1.0 - loading * loading) * record["noise_energy"]
                    + loading * loading * record["shared_energy"]
                    + 2.0 * loading * noise_weight * exact_cross
                )
                if float(metrics.min()) < -1e-10:
                    raise RuntimeError("Negative graph energy beyond tolerance")
                metrics = np.maximum(metrics, 0.0)
                lookup = effect_lookup(
                    metrics, record["high"], record["middle"]
                )
                observed += lookup[:, 0, :]
                selected = np.take_along_axis(
                    lookup,
                    offsets[record_index][:, :, None],
                    axis=1,
                )
                shifted += selected
                record_metrics.append(metrics)

            observed /= len(records)
            shifted /= len(records)
            observed_profile = observed.mean(axis=1)
            shifted_profile = shifted.mean(axis=2)
            p_values = (
                1
                + np.sum(
                    shifted_profile <= observed_profile[:, None], axis=1
                )
            ) / (PRODUCT_SHIFTS + 1)
            reject_counts[lambda_index] += int(np.sum(p_values <= ALPHA))
            p_value_sums[lambda_index] += float(p_values.sum())
            estimate_sums[lambda_index] += float(observed_profile.sum())

            # July and August use the same 31-point normalized phase grid.
            # Their daily graph-profile correlation is a direct empirical
            # check that the shared phase produces cross-record dependence.
            for left, right in ((1, 2), (4, 5)):
                dependence_correlations[lambda_index].update(
                    record_metrics[left].mean(axis=2),
                    record_metrics[right].mean(axis=2),
                )

        completed += batch
        print(
            f"completed {completed}/{REPLICATIONS} paired replications",
            flush=True,
        )

    summary = []
    for index, loading in enumerate(LAMBDAS):
        rejection = reject_counts[index] / REPLICATIONS
        summary.append({
            "shared_loading_lambda": float(loading),
            "shared_variance_fraction_spatial_average": float(loading**2),
            "replications": REPLICATIONS,
            "product_shifts_per_replication": PRODUCT_SHIFTS,
            "nominal_alpha": ALPHA,
            "rejections": int(reject_counts[index]),
            "empirical_rejection_rate": rejection,
            "monte_carlo_se": math.sqrt(
                rejection * (1.0 - rejection) / REPLICATIONS
            ),
            "mean_p_value": p_value_sums[index] / REPLICATIONS,
            "mean_observed_statistic": estimate_sums[index] / REPLICATIONS,
            "same_summer_july_august_daily_profile_correlation": (
                dependence_correlations[index].correlation()
            ),
        })

    diagnostics = {
        "type7_group_sizes_by_record_length": {
            length: [list(pair) for pair in sorted(pairs)]
            for length, pairs in label_counts.items()
        },
        "spawned_rng_states": {
            "field_entropy": field_seed.entropy,
            "field_spawn_key": list(field_seed.spawn_key),
            "label_entropy": label_seed.entropy,
            "label_spawn_key": list(label_seed.spawn_key),
            "shift_entropy": shift_seed.entropy,
            "shift_spawn_key": list(shift_seed.spawn_key),
        },
    }
    return summary, diagnostics


def write_outputs(
    summary: list[dict[str, object]],
    diagnostics: dict,
    graph: GraphSystem,
    reference: pd.DataFrame,
) -> tuple[Path, Path]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_path = OUTPUT_DIR / "extension_cross_record_stress_summary.csv"
    audit_path = OUTPUT_DIR / "extension_cross_record_stress_audit.json"

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    checks = covariance_checks(graph)
    checks.update(diagnostics)
    if checks["marginal_shared_process_one_step_rotation_max_abs_error"] > 1e-12:
        raise RuntimeError("Marginal Fourier covariance is not cyclic invariant")
    if (
        checks[
            "unscaled_30_by_31_cross_covariance_relative_change_after_one_30_day_shift"
        ]
        <= 0.01
    ):
        raise RuntimeError("Independent shift does not alter cross covariance")
    if checks["laplacian_constant_annihilation_max_abs"] > 1e-10:
        raise RuntimeError("Graph Laplacian does not annihilate constants")

    audit = {
        "script": str(SCRIPT_PATH.relative_to(PROJECT_DIR)),
        "analysis_role": (
            "targeted size stress test; not an extension of the product-"
            "invariance proposition"
        ),
        "seed": SEED,
        "replications_per_lambda": REPLICATIONS,
        "product_shift_draws_per_replication": PRODUCT_SHIFTS,
        "nominal_alpha": ALPHA,
        "shared_loadings_lambda": LAMBDAS.tolist(),
        "paired_base_fields_labels_and_offsets_across_lambda": True,
        "site_count": int(len(reference)),
        "record_structure": {
            "summers": 2,
            "records_per_summer": 3,
            "months": RECORD_MONTHS.tolist(),
            "record_lengths": RECORD_LENGTHS.tolist(),
        },
        "statistic": (
            "equal-record, equal-five-bandwidth mean of record-specific "
            "high/middle graph-dispersion ratios minus one"
        ),
        "test": {
            "tail": "lower",
            "p_value": (
                "(1 + number shifted statistics <= observed)/"
                f"({PRODUCT_SHIFTS} + 1)"
            ),
            "product_action": (
                "one independently sampled uniform cyclic offset per monthly "
                "record; each five-scale daily profile shifted intact"
            ),
        },
        "dgp": {
            "field": "y_ymt = X_ymt*1 + z_ymt",
            "anomaly": (
                "z_ymt = lambda*U_y(t/n_m)*g + "
                "sqrt(1-lambda^2)*epsilon_ymt"
            ),
            "shared_process": {
                "definition": (
                    "three-harmonic isotropic Gaussian Fourier process on "
                    "normalized circular phase; sine/cosine coefficients "
                    "shared by the three records within each summer"
                ),
                "harmonics": [1, 2, 3],
                "variance_weights": FOURIER_VARIANCE_WEIGHTS.tolist(),
                "pointwise_variance": 1.0,
            },
            "spatial_pattern_g": (
                "centred 0.8*x+y coordinate gradient, rescaled to spatial "
                "mean square one"
            ),
            "month_specific_noise": {
                "spatial_basis": (
                    "12 leading eigenmodes of the centred 0.95*exp(-d/450 "
                    "km)+0.05*I covariance, rescaled to mean site variance one"
                ),
                "temporal_filter_offsets": FILTER_OFFSETS.tolist(),
                "temporal_filter_weights": NOISE_FILTER.tolist(),
                "independence": "independent across month-year records",
            },
            "regional_mean_X": {
                "temporal_filter_offsets": FILTER_OFFSETS.tolist(),
                "temporal_filter_weights": MEAN_FILTER.tolist(),
                "independence": (
                    "separate SeedSequence child stream; independent of U "
                    "and epsilon"
                ),
                "labels": (
                    "within-record NumPy method='linear' quartiles, equivalent "
                    "to Hyndman-Fan/R type 7; high >= q75, middle q25 <= X < q75"
                ),
            },
        },
        "graphs": {
            "h_factors": H_FACTORS.tolist(),
            "bandwidths_km": graph.bandwidths_km.tolist(),
            "weights": "exp(-distance_km^2/(2*h^2)), diagonal zero",
            "dispersion": "z' L_h z/(2*sum_{i<j} w_ij)",
        },
        "numerical_checks": checks,
        "summary": summary,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "platform": platform.platform(),
        },
        "inputs": [{
            "path": str(REFERENCE_PATH.relative_to(PROJECT_DIR)),
            "bytes": REFERENCE_PATH.stat().st_size,
            "sha256": sha256(REFERENCE_PATH),
        }],
        "outputs": [{
            "path": str(summary_path.relative_to(PROJECT_DIR)),
            "bytes": summary_path.stat().st_size,
            "sha256": sha256(summary_path),
        }],
    }
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return summary_path, audit_path


def main() -> None:
    graph, reference = build_graph_system()
    summary, diagnostics = run_simulation(graph)
    summary_path, audit_path = write_outputs(
        summary, diagnostics, graph, reference
    )
    print(pd.DataFrame(summary).to_string(index=False))
    print(f"wrote {summary_path.relative_to(PROJECT_DIR)}")
    print(f"wrote {audit_path.relative_to(PROJECT_DIR)}")


if __name__ == "__main__":
    main()
