"""Utilities for comparing APWP construction and rate-estimation methods.

The functions in this module support the revision analysis that compares
traditional running means, site-level bootstrap APWPs, and Bayesian Euler-pole
paths.  Plotting and project-specific file loading belong in the accompanying
notebook; this module contains only reusable numerical operations.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from vgptools.utils import GCD_cartesian, spherical2cartesian


D2R = np.pi / 180.0
R2D = 180.0 / np.pi


def circular_mean(longitudes: Iterable[float]) -> float:
    """Return the circular mean of longitudes in degrees on ``[0, 360)``.

    Parameters
    ----------
    longitudes
        Longitude values in degrees. Missing values are ignored.

    Returns
    -------
    float
        Circular mean longitude in degrees. This avoids the erroneous result
        produced by an arithmetic mean when values straddle 0/360 degrees.
    """

    lon_rad = np.deg2rad(np.asarray(longitudes, dtype=float))
    mean_lon = np.arctan2(np.nanmean(np.sin(lon_rad)), np.nanmean(np.cos(lon_rad)))
    return float(np.rad2deg(mean_lon) % 360.0)


def align_longitudes(longitudes: Iterable[float], center: float | None = None) -> np.ndarray:
    """Unwrap longitudes onto the branch centered on a reference longitude.

    Parameters
    ----------
    longitudes
        Longitude values in degrees.
    center
        Center of the desired 360-degree branch. If omitted, the circular mean
        of ``longitudes`` is used.

    Returns
    -------
    numpy.ndarray
        Longitudes expressed between ``center - 180`` and ``center + 180``.
        The aligned values are useful for scalar quantiles around a pole cloud
        that crosses the conventional longitude discontinuity.
    """

    lon = np.asarray(longitudes, dtype=float)
    if center is None:
        center = circular_mean(lon)
    return center + ((lon - center + 180.0) % 360.0 - 180.0)


def cartesian_mean_path(
    paths: pd.DataFrame,
    age_col: str = "age",
    longitude_col: str = "plon",
    latitude_col: str = "plat",
) -> pd.DataFrame:
    """Summarize an ensemble of paleopoles at each age on the unit sphere.

    Each pole is converted to a Cartesian unit vector before averaging. This
    avoids coordinate-dependent bias from separately averaging latitude and
    longitude. Empirical 2.5, 50, and 97.5 percent quantiles are also returned
    for plotting; longitude quantiles are calculated after branch alignment.

    Parameters
    ----------
    paths
        Long-form table containing one or more pole positions at each age.
    age_col, longitude_col, latitude_col
        Column names for age, longitude, and latitude.

    Returns
    -------
    pandas.DataFrame
        One row per age with the spherical Cartesian mean (``plon``, ``plat``)
        and marginal empirical quantiles for longitude and latitude.
    """

    rows: list[dict[str, float]] = []
    for age, group in paths.groupby(age_col):
        xyz = np.array(
            [
                spherical2cartesian(
                    [np.radians(row[latitude_col]), np.radians(row[longitude_col])]
                )
                for _, row in group.iterrows()
            ]
        )
        mean_xyz = xyz.mean(axis=0)
        mean_xyz /= np.linalg.norm(mean_xyz)
        latitude = np.rad2deg(np.arcsin(mean_xyz[2]))
        longitude = np.rad2deg(np.arctan2(mean_xyz[1], mean_xyz[0])) % 360.0
        aligned_lon = align_longitudes(group[longitude_col].to_numpy(), longitude)
        rows.append(
            {
                "age": float(age),
                "plon": longitude,
                "plat": latitude,
                "plon_q025": np.nanquantile(aligned_lon, 0.025),
                "plon_q50": np.nanquantile(aligned_lon, 0.5),
                "plon_q975": np.nanquantile(aligned_lon, 0.975),
                "plat_q025": np.nanquantile(group[latitude_col], 0.025),
                "plat_q50": np.nanquantile(group[latitude_col], 0.5),
                "plat_q975": np.nanquantile(group[latitude_col], 0.975),
            }
        )
    return pd.DataFrame(rows).sort_values("age").reset_index(drop=True)


def summarize_quantiles(
    values: pd.DataFrame,
    age_col: str = "age",
    value_col: str = "APW_rate",
) -> pd.DataFrame:
    """Calculate median and empirical 95 percent limits by age.

    Parameters
    ----------
    values
        Long-form ensemble table.
    age_col
        Column used to group values in time.
    value_col
        Numeric column to summarize.

    Returns
    -------
    pandas.DataFrame
        Columns ``age``, ``q025``, ``q50``, and ``q975``, sorted by age.
    """

    return (
        values.groupby(age_col)[value_col]
        .quantile([0.025, 0.5, 0.975])
        .unstack()
        .rename(columns={0.025: "q025", 0.5: "q50", 0.975: "q975"})
        .reset_index()
        .rename(columns={age_col: "age"})
        .sort_values("age")
        .reset_index(drop=True)
    )


def calculate_coherent_path_rates(
    paths: pd.DataFrame,
    path_col: str = "run",
    age_col: str = "age",
    longitude_col: str = "plon",
    latitude_col: str = "plat",
    output_col: str = "APW_rate",
) -> pd.DataFrame:
    """Calculate segment rates while preserving each sampled path trajectory.

    This is the trajectory-coherent alternative to the pooled effective-age
    resampling used by Gallo et al. (2023). Adjacent poles always come from the
    same bootstrap or posterior realization. It is retained because it is
    useful for sensitivity tests, even when figures reproduce the published
    Gallo convention.

    Parameters
    ----------
    paths
        Long-form table with one pole per age and path realization.
    path_col, age_col, longitude_col, latitude_col
        Column names identifying paths, ages, and pole coordinates.
    output_col
        Name assigned to the calculated rate column.

    Returns
    -------
    pandas.DataFrame
        Copy of ``paths`` sorted by path and age, with great-circle segment
        rates in degrees per Myr. The first pole in each path has a missing
        rate because no preceding segment exists.
    """

    frames: list[pd.DataFrame] = []
    for _, group in paths.sort_values([path_col, age_col]).groupby(path_col):
        group = group.copy()
        vectors = [
            spherical2cartesian(
                [np.radians(row[latitude_col]), np.radians(row[longitude_col])]
            )
            for _, row in group.iterrows()
        ]
        distances = [np.nan]
        for index in range(1, len(vectors)):
            distances.append(np.degrees(GCD_cartesian(vectors[index - 1], vectors[index])))
        group[output_col] = np.asarray(distances) / group[age_col].diff().to_numpy()
        frames.append(group)
    return pd.concat(frames, ignore_index=True)


def gallo_effective_age_resample(
    ensemble: pd.DataFrame,
    ages: Sequence[float] | np.ndarray,
    n_paths: int,
    random_seed: int | None = None,
    effective_age_col: str = "effective_age",
    longitude_col: str = "plon",
    latitude_col: str = "plat",
) -> pd.DataFrame:
    """Reproduce the APW-rate resampling convention of Gallo et al. (2023).

    For every requested rounded effective age, one pole is sampled from all
    ensemble rows assigned that effective age. Sampling is repeated
    independently at the next age, so a synthetic path generally combines
    poles from different original bootstrap runs. Great-circle rates are then
    calculated between adjacent draws and divided by their effective-age
    difference.

    This behavior intentionally matches cell 19 of
    ``2_MC_uncertainty_propagation.ipynb`` in the Gallo et al. (2023) archive.
    It should not be confused with preserving the covariance of one bootstrap
    trajectory; use :func:`calculate_coherent_path_rates` for that alternative.

    Parameters
    ----------
    ensemble
        Long-form moving-window ensemble containing rounded effective ages and
        paleopole coordinates.
    ages
        Effective ages to sample, in the order used to construct each path.
        Every requested age must occur in ``effective_age_col``.
    n_paths
        Number of synthetic paths to generate.
    random_seed
        Optional seed for reproducible row selection.
    effective_age_col, longitude_col, latitude_col
        Column names for rounded effective age and pole coordinates.

    Returns
    -------
    pandas.DataFrame
        Sampled source rows plus ``source_run`` (when the input has ``run``),
        synthetic ``run``, ``GCD``, and ``APW_rate`` columns.

    Raises
    ------
    ValueError
        If ``n_paths`` is not positive, fewer than two ages are requested, ages
        are duplicated or unordered, or an age has no candidate pole.
    """

    requested_ages = np.asarray(ages, dtype=float)
    if n_paths <= 0:
        raise ValueError("n_paths must be positive")
    if requested_ages.size < 2:
        raise ValueError("At least two effective ages are required")
    if np.any(np.diff(requested_ages) <= 0):
        raise ValueError("Effective ages must be unique and strictly increasing")

    pools: dict[float, pd.DataFrame] = {}
    for age in requested_ages:
        pool = ensemble.loc[np.isclose(ensemble[effective_age_col], age)].reset_index(drop=True)
        if pool.empty:
            raise ValueError(f"No ensemble rows have {effective_age_col}={age:g}")
        pools[float(age)] = pool

    rng = np.random.default_rng(random_seed)
    frames: list[pd.DataFrame] = []
    for synthetic_run in range(n_paths):
        rows = []
        for age in requested_ages:
            pool = pools[float(age)]
            rows.append(pool.iloc[int(rng.integers(len(pool)))].copy())
        path = pd.DataFrame(rows).reset_index(drop=True)
        if "run" in path:
            path = path.rename(columns={"run": "source_run"})
        path[effective_age_col] = requested_ages
        path["run"] = synthetic_run

        vectors = [
            spherical2cartesian(
                [np.radians(row[latitude_col]), np.radians(row[longitude_col])]
            )
            for _, row in path.iterrows()
        ]
        distances = [np.nan]
        for index in range(1, len(vectors)):
            distances.append(np.degrees(GCD_cartesian(vectors[index - 1], vectors[index])))
        path["GCD"] = distances
        path["APW_rate"] = path["GCD"] / path[effective_age_col].diff()
        frames.append(path)

    return pd.concat(frames, ignore_index=True)


def rotation_matrix(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Construct the Euler rotation matrix used by the Bayesian PEP model.

    Angles are supplied in radians. The matrix follows the same ordered
    z-y-z rotations implemented in ``Bayesian_PEP_inversion.py`` so archived
    posterior parameters can be evaluated without importing legacy PyMC3.
    """

    rot_alpha = np.array(
        [[np.cos(alpha), -np.sin(alpha), 0.0], [np.sin(alpha), np.cos(alpha), 0.0], [0.0, 0.0, 1.0]]
    )
    rot_beta = np.array(
        [[np.cos(beta), 0.0, np.sin(beta)], [0.0, 1.0, 0.0], [-np.sin(beta), 0.0, np.cos(beta)]]
    )
    rot_gamma = np.array(
        [[np.cos(gamma), -np.sin(gamma), 0.0], [np.sin(gamma), np.cos(gamma), 0.0], [0.0, 0.0, 1.0]]
    )
    return rot_gamma @ rot_beta @ rot_alpha


def spherical_to_cartesian(longitude: float, latitude: float) -> np.ndarray:
    """Convert longitude and latitude in degrees to a Cartesian unit vector."""

    colatitude = 90.0 - latitude
    return np.array(
        [
            np.sin(colatitude * D2R) * np.cos(longitude * D2R),
            np.sin(colatitude * D2R) * np.sin(longitude * D2R),
            np.cos(colatitude * D2R),
        ]
    )


def cartesian_to_spherical(vector: np.ndarray) -> tuple[float, float]:
    """Convert a Cartesian vector to longitude and latitude in degrees."""

    vector = np.asarray(vector, dtype=float)
    vector /= np.linalg.norm(vector)
    longitude = np.arctan2(vector[1], vector[0]) * R2D
    latitude = 90.0 - np.arccos(vector[2]) * R2D
    return float(longitude % 360.0), float(latitude)


def rotate_about_euler_pole(
    pole_vector: np.ndarray,
    euler_longitude: float,
    euler_latitude: float,
    angle: float,
) -> np.ndarray:
    """Rotate a paleopole vector about an Euler pole by an angle in degrees."""

    euler_vector = spherical_to_cartesian(euler_longitude, euler_latitude)
    normalized_lon, normalized_lat = cartesian_to_spherical(euler_vector)
    colatitude = 90.0 - normalized_lat
    first = rotation_matrix(-normalized_lon * D2R, -colatitude * D2R, angle * D2R)
    second = rotation_matrix(0.0, colatitude * D2R, normalized_lon * D2R)
    return second @ (first @ pole_vector)


def evaluate_two_euler_position(parameters: pd.Series, age: float) -> tuple[float, float]:
    """Evaluate one two-stage Euler posterior sample at a requested age.

    ``parameters`` must contain the transformed posterior variables
    ``start_lon``, ``start_lat``, ``start_pole_age``, ``switchpoint``,
    ``euler1_lon``, ``euler1_lat``, ``rate_1``, ``euler2_lon``,
    ``euler2_lat``, and ``rate_2``. The first stage applies before the inferred
    switchpoint and the second stage continues from that position afterward.
    """

    pole = spherical_to_cartesian(parameters.start_lon, parameters.start_lat)
    if age > parameters.switchpoint:
        pole = rotate_about_euler_pole(
            pole,
            parameters.euler1_lon,
            parameters.euler1_lat,
            parameters.rate_1 * (parameters.start_pole_age - age),
        )
    else:
        pole = rotate_about_euler_pole(
            pole,
            parameters.euler1_lon,
            parameters.euler1_lat,
            parameters.rate_1 * (parameters.start_pole_age - parameters.switchpoint),
        )
        pole = rotate_about_euler_pole(
            pole,
            parameters.euler2_lon,
            parameters.euler2_lat,
            parameters.rate_2 * (parameters.switchpoint - age),
        )
    return cartesian_to_spherical(pole)


def load_two_euler_npz_paths(
    chain_paths: Sequence[str | Path],
    ages: Sequence[float] | np.ndarray,
    thin: int = 100,
) -> pd.DataFrame:
    """Evaluate archived two-Euler Bayesian PEP chains as paleopole paths.

    Parameters
    ----------
    chain_paths
        Paths to ``samples.npz`` files containing transformed posterior draws.
    ages
        Ages in Ma at which every retained posterior draw is evaluated.
    thin
        Retain every ``thin``-th draw from each chain. A value of 100 converts
        four 200,000-draw chains into 8,000 path realizations.

    Returns
    -------
    pandas.DataFrame
        Long-form table with ``sample``, ``age``, ``plon``, ``plat``, and
        coherent ``APW_rate`` columns. Rates connect adjacent ages from the
        same posterior draw because each draw is a complete kinematic model.
    """

    if thin <= 0:
        raise ValueError("thin must be positive")
    requested_ages = np.asarray(ages, dtype=float)
    rows: list[dict[str, float | int]] = []
    sample_id = 0
    for chain_path in sorted(Path(path) for path in chain_paths):
        with np.load(chain_path) as samples:
            take = np.arange(0, len(samples["rate_1"]), thin)
            trace = pd.DataFrame(
                {
                    "euler1_lon": samples["euler_1"][take, 0],
                    "euler1_lat": samples["euler_1"][take, 1],
                    "rate_1": samples["rate_1"][take],
                    "euler2_lon": samples["euler_2"][take, 0],
                    "euler2_lat": samples["euler_2"][take, 1],
                    "rate_2": samples["rate_2"][take],
                    "start_pole_age": samples["start_pole_age"][take],
                    "start_lon": samples["start_pole"][take, 0],
                    "start_lat": samples["start_pole"][take, 1],
                    "switchpoint": samples["switchpoint"][take],
                }
            )
        for _, parameters in trace.iterrows():
            previous_vector = None
            previous_age = None
            for age in requested_ages:
                longitude, latitude = evaluate_two_euler_position(parameters, age)
                current_vector = spherical2cartesian(
                    [np.radians(latitude), np.radians(longitude)]
                )
                rate = np.nan
                if previous_vector is not None:
                    distance = np.degrees(GCD_cartesian(previous_vector, current_vector))
                    rate = distance / (age - previous_age)
                rows.append(
                    {
                        "sample": sample_id,
                        "age": age,
                        "plon": longitude,
                        "plat": latitude,
                        "APW_rate": rate,
                    }
                )
                previous_vector = current_vector
                previous_age = age
            sample_id += 1
    return pd.DataFrame(rows)
