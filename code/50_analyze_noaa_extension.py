#!/usr/bin/env python3
"""Validate ERA5-Land in ten frozen non-development NOAA station years."""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT = Path(__file__).resolve().parent.parent
PROTOCOL = PROJECT / "EXTENSION_ANALYSIS_PROTOCOL.md"
BASE = PROJECT / "data" / "noaa_isd_extension"
ERA_DIR = BASE / "era5_land_points" / "trimmed_points"
ERA_PROVENANCE = BASE / "era5_land_points" / "era5_station_point_provenance.json"
MANIFEST = BASE / "noaa_extension_station_manifest.csv"
QUALIFICATION = BASE / "noaa_extension_station_year_qualification.csv"
EVENT_FILE = PROJECT / "output_confirmatory" / "sensitivity_event_manifest.csv"
OUTPUT = PROJECT / "output_noaa_extension"
MATCHED_FILE = OUTPUT / "noaa_extension_era5_matched.csv.gz"
MEASUREMENT_YEAR_FILE = OUTPUT / "noaa_extension_measurement_by_year.csv"
MEASUREMENT_OVERALL_FILE = OUTPUT / "noaa_extension_measurement_overall.csv"
FIELD_FILE = OUTPUT / "noaa_extension_graph_fields.csv.gz"
RECORD_FILE = OUTPUT / "noaa_extension_record_effects.csv"
YEAR_EFFECT_FILE = OUTPUT / "noaa_extension_year_scale_effects.csv"
SCALE_FILE = OUTPUT / "noaa_extension_scale_summary.csv"
PROFILE_FILE = OUTPUT / "noaa_extension_broad_profile_summary.csv"
TEX_FILE = OUTPUT / "noaa_extension_tables.tex"
AUDIT_FILE = OUTPUT / "noaa_extension_analysis_audit.json"
YEARS = (1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020, 2023, 2025)
BANDWIDTHS = np.array([125.799765, 251.599530, 503.199060, 1006.398120, 2012.796241])
BROAD = {503.199060, 1006.398120, 2012.796241}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def saturation_vapor_pressure_pa(temperature_k: np.ndarray) -> np.ndarray:
    temperature_c = temperature_k - 273.15
    return 611.2 * np.exp(17.67 * temperature_c / (temperature_c + 243.5))


def theta_e(temperature_k: np.ndarray, dewpoint_k: np.ndarray,
            pressure_pa: np.ndarray) -> np.ndarray:
    vapor = saturation_vapor_pressure_pa(dewpoint_k)
    mixing = 0.622 * vapor / (pressure_pa - vapor)
    t_lcl = 1 / (1 / (dewpoint_k - 56) + np.log(temperature_k / dewpoint_k) / 800) + 56
    theta_l = (temperature_k * (100000 / (pressure_pa - vapor)) ** 0.2854 *
               (temperature_k / t_lcl) ** (0.28 * mixing))
    return theta_l * np.exp((3036 / t_lcl - 1.78) * mixing * (1 + 0.448 * mixing))


def wet_bulb_c(temperature_k: np.ndarray, dewpoint_k: np.ndarray,
               pressure_pa: np.ndarray) -> np.ndarray:
    dewpoint_k = np.minimum(dewpoint_k, temperature_k)
    target = theta_e(temperature_k, dewpoint_k, pressure_pa)
    lower, upper = dewpoint_k.copy(), temperature_k.copy()
    for _ in range(40):
        middle = (lower + upper) / 2
        move = theta_e(middle, middle, pressure_pa) < target
        lower = np.where(move, middle, lower)
        upper = np.where(move, upper, middle)
    return (lower + upper) / 2 - 273.15


def load_noaa() -> tuple[pd.DataFrame, list[Path]]:
    paths = [BASE / f"noaa_isd_extension_{year}_jja_exact_hours.csv.gz" for year in YEARS]
    if not all(path.exists() for path in paths):
        raise RuntimeError("NOAA extension download is incomplete")
    data = pd.concat([pd.read_csv(path, dtype={"station": str}) for path in paths],
                     ignore_index=True)
    data["time_utc"] = pd.to_datetime(data.time_utc, utc=True).dt.tz_localize(None)
    qualification = pd.read_csv(QUALIFICATION, dtype={"station": str})
    qualified = qualification.loc[qualification.qualified, ["station", "year"]]
    data = data.merge(qualified, on=["station", "year"], how="inner", validate="many_to_one")
    data["observed_wbt_c"] = wet_bulb_c(
        data.temperature_c.to_numpy() + 273.15,
        data.dewpoint_c.to_numpy() + 273.15,
        data.station_pressure_pa.to_numpy())
    return data, paths


def load_era() -> tuple[pd.DataFrame, list[Path]]:
    import xarray as xr

    provenance = json.loads(ERA_PROVENANCE.read_text(encoding="utf-8"))
    available = {
        str(item["station"])
        for item in provenance["stations"]
        if item["status"] == "available"
    }
    paths = sorted(ERA_DIR.glob("station_*_evaluation_jja.nc"))
    path_stations = {path.name.split("_")[1] for path in paths}
    if len(paths) != provenance["available_station_series"] or path_stations != available:
        raise RuntimeError(
            "Finite ERA5-Land station files do not match their provenance"
        )
    pieces = []
    for path in paths:
        station = path.name.split("_")[1]
        with xr.open_dataset(path, engine="h5netcdf") as dataset:
            frame = dataset[["t2m", "d2m", "sp"]].to_dataframe().reset_index()
        frame["station"] = station
        frame["time_utc"] = pd.to_datetime(frame.time, utc=True).dt.tz_localize(None)
        frame["era_temperature_c"] = frame.t2m - 273.15
        frame["era_dewpoint_c"] = np.minimum(frame.d2m, frame.t2m) - 273.15
        frame["era_wbt_c"] = wet_bulb_c(
            frame.t2m.to_numpy(), np.minimum(frame.d2m, frame.t2m).to_numpy(),
            frame.sp.to_numpy())
        pieces.append(frame[["station", "time_utc", "era_temperature_c",
                             "era_dewpoint_c", "sp", "era_wbt_c"]])
    return pd.concat(pieces, ignore_index=True), paths


def correlation(x: pd.Series, y: pd.Series) -> float:
    return float(x.corr(y)) if len(x) > 1 else math.nan


def measurement_summary(data: pd.DataFrame) -> dict:
    error = data.era_wbt_c - data.observed_wbt_c
    station_correlations = data.groupby("station").apply(
        lambda group: correlation(group.observed_wbt_c, group.era_wbt_c),
        include_groups=False)
    centered = data.copy()
    centered["observed_centered"] = centered.observed_wbt_c - centered.groupby(
        "station").observed_wbt_c.transform("mean")
    centered["era_centered"] = centered.era_wbt_c - centered.groupby(
        "station").era_wbt_c.transform("mean")
    return {
        "matched_hours": len(data), "stations": data.station.nunique(),
        "bias_c": float(error.mean()), "mae_c": float(error.abs().mean()),
        "rmse_c": float(np.sqrt(np.square(error).mean())),
        "pooled_correlation": correlation(data.observed_wbt_c, data.era_wbt_c),
        "within_station_centered_correlation": correlation(
            centered.observed_centered, centered.era_centered),
        "equal_station_mean_correlation": float(station_correlations.mean()),
    }


def haversine(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    longitude, latitude = np.deg2rad(lon), np.deg2rad(lat)
    dlat = latitude[:, None] - latitude[None, :]
    dlon = longitude[:, None] - longitude[None, :]
    a = (np.sin(dlat / 2) ** 2 + np.cos(latitude[:, None]) *
         np.cos(latitude[None, :]) * np.sin(dlon / 2) ** 2)
    return 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def graph_fields(matched: pd.DataFrame, event: pd.DataFrame) -> pd.DataFrame:
    event = event.copy()
    event["peak_time"] = pd.to_datetime(event.peak_time)
    event = event.loc[event.year.isin(YEARS), ["peak_time", "record_id", "year", "month", "regime"]]
    event = event.loc[event.regime.isin(["high", "middle"])]
    selected = matched.merge(event, left_on="time_utc", right_on="peak_time",
                             how="inner", validate="many_to_one")
    rows = []
    for time, group in selected.groupby("time_utc", sort=True):
        group = group.sort_values("station").drop_duplicates("station")
        if len(group) < 10:
            continue
        distance = haversine(group.LONGITUDE.to_numpy(), group.LATITUDE.to_numpy())
        upper = np.triu_indices(len(group), 1)
        observed, era = group.observed_wbt_c.to_numpy(), group.era_wbt_c.to_numpy()
        for bandwidth in BANDWIDTHS:
            weights = np.exp(-np.square(distance[upper]) / (2 * bandwidth**2))
            rows.append({
                "time_utc": time, "record_id": int(group.record_id.iloc[0]),
                "year": int(group.year_y.iloc[0]), "month": int(group.month.iloc[0]),
                "regime": group.regime.iloc[0], "bandwidth_km": bandwidth,
                "stations": len(group),
                "observed_q": float(np.sum(weights * np.square(
                    observed[upper[0]] - observed[upper[1]])) / (2 * weights.sum())),
                "era_q": float(np.sum(weights * np.square(
                    era[upper[0]] - era[upper[1]])) / (2 * weights.sum())),
            })
    return pd.DataFrame(rows)


def summarize_effects(fields: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    records = []
    for (record_id, bandwidth), group in fields.groupby(["record_id", "bandwidth_km"], sort=True):
        high, middle = group.regime.eq("high"), group.regime.eq("middle")
        if high.sum() < 1 or middle.sum() < 1:
            continue
        records.append({
            "record_id": int(record_id), "year": int(record_id) // 100,
            "month": int(record_id) % 100, "bandwidth_km": bandwidth,
            "high_fields": int(high.sum()), "middle_fields": int(middle.sum()),
            "mean_stations": float(group.stations.mean()),
            "station_effect": group.loc[high, "observed_q"].mean() /
                              group.loc[middle, "observed_q"].mean() - 1,
            "era_effect": group.loc[high, "era_q"].mean() /
                          group.loc[middle, "era_q"].mean() - 1,
        })
    records = pd.DataFrame(records)
    years = records.groupby(["year", "bandwidth_km"], as_index=False).agg(
        station_effect=("station_effect", "mean"), era_effect=("era_effect", "mean"),
        records=("record_id", "nunique"), high_fields=("high_fields", "sum"),
        middle_fields=("middle_fields", "sum"), mean_stations=("mean_stations", "mean"))
    scale_rows = []
    for bandwidth, group in years.groupby("bandwidth_km", sort=True):
        row = {"bandwidth_km": bandwidth, "years": group.year.nunique(),
               "records": int(group.records.sum()),
               "mean_stations": float(group.mean_stations.mean())}
        for name in ("station_effect", "era_effect"):
            values = group[name].to_numpy()
            n = len(values); estimate = values.mean(); se = values.std(ddof=1) / math.sqrt(n)
            critical = stats.t.ppf(0.975, n - 1)
            row.update({f"{name}_estimate": estimate, f"{name}_se": se,
                        f"{name}_ci_lower": estimate - critical * se,
                        f"{name}_ci_upper": estimate + critical * se,
                        f"{name}_negative_years": int((values < 0).sum())})
        scale_rows.append(row)
    scales = pd.DataFrame(scale_rows)
    broad_year = years.loc[years.bandwidth_km.isin(BROAD)].groupby("year", as_index=False).agg(
        station_effect=("station_effect", "mean"), era_effect=("era_effect", "mean"),
        scales=("bandwidth_km", "nunique"))
    profile_rows = []
    for name in ("station_effect", "era_effect"):
        values = broad_year[name].to_numpy(); n = len(values)
        estimate = values.mean(); se = values.std(ddof=1) / math.sqrt(n)
        critical = stats.t.ppf(0.975, n - 1)
        profile_rows.append({
            "field": "NOAA station" if name == "station_effect" else "ERA5-Land",
            "years": n, "estimate": estimate, "standard_error": se,
            "ci_lower": estimate - critical * se, "ci_upper": estimate + critical * se,
            "negative_years": int((values < 0).sum()),
        })
    return records, years, scales, pd.DataFrame(profile_rows)


def write_tex(by_year: pd.DataFrame, overall: pd.DataFrame,
              scales: pd.DataFrame, profile: pd.DataFrame) -> None:
    lines = [
        "% Generated by code/50_analyze_noaa_extension.py",
        r"\begin{table}[H]\centering",
        r"\caption{Measurement agreement in the ten non-development NOAA years. Bias is ERA5-Land minus station WBT. The centred correlation removes each station's mean.}",
        r"\label{tab:noaa-extension-measurement}\footnotesize",
        r"\begin{tabular}{rrrrrrr}\toprule",
        r"Year & Stations & Hours & Bias & MAE & RMSE & Centred correlation \\",
        r"\midrule",
    ]
    for row in by_year.itertuples(index=False):
        lines.append(f"{row.year} & {row.stations} & {row.matched_hours:,} & "
                     f"{row.bias_c:.2f} & {row.mae_c:.2f} & {row.rmse_c:.2f} & "
                     f"{row.within_station_centered_correlation:.3f} " + r"\\")
    row = overall.iloc[0]
    lines += [r"\midrule", f"All & {int(row.stations)} & {int(row.matched_hours):,} & "
              f"{row.bias_c:.2f} & {row.mae_c:.2f} & {row.rmse_c:.2f} & "
              f"{row.within_station_centered_correlation:.3f} " + r"\\",
              r"\bottomrule\end{tabular}\end{table}", ""]
    lines += [
        r"\begin{table}[H]\centering",
        r"\caption{High--middle graph contrasts at frozen ERA5-Land peak times and labels in ten non-development station years. Values are percentage changes. The three broadest scales were fixed as station-supported before the extension was evaluated.}",
        r"\label{tab:noaa-extension-effect}\footnotesize",
        r"\begin{tabular}{rrrrrrr}\toprule",
        r"Bandwidth & Years & Records & Mean stations & NOAA effect & ERA effect & NOAA negative \\",
        r"(km) & & & & (\%) & (\%) & summers \\", r"\midrule",
    ]
    for row in scales.itertuples(index=False):
        lines.append(f"{row.bandwidth_km:,.0f} & {row.years} & {row.records} & "
                     f"{row.mean_stations:.1f} & {100 * row.station_effect_estimate:.2f} & "
                     f"{100 * row.era_effect_estimate:.2f} & "
                     f"{row.station_effect_negative_years}/{row.years} " + r"\\")
    lines += [r"\midrule"]
    station = profile.loc[profile.field.eq("NOAA station")].iloc[0]
    era = profile.loc[profile.field.eq("ERA5-Land")].iloc[0]
    lines += [f"Broad-profile mean & {int(station.years)} & -- & -- & "
              f"{100 * station.estimate:.2f} & {100 * era.estimate:.2f} & "
              f"{int(station.negative_years)}/{int(station.years)} " + r"\\",
              r"\bottomrule\end{tabular}\end{table}", ""]
    TEX_FILE.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    observed, noaa_paths = load_noaa()
    era, era_paths = load_era()
    matched = observed.merge(era, on=["station", "time_utc"], how="inner", validate="one_to_one")
    if matched.empty or not np.isfinite(matched[["observed_wbt_c", "era_wbt_c"]]).all().all():
        raise RuntimeError("NOAA--ERA matched panel is empty or nonfinite")
    matched.to_csv(MATCHED_FILE, index=False, compression="gzip")
    year_rows = []
    for year, group in matched.groupby("year", sort=True):
        year_rows.append({"year": year, **measurement_summary(group)})
    by_year = pd.DataFrame(year_rows)
    overall = pd.DataFrame([measurement_summary(matched)])
    by_year.to_csv(MEASUREMENT_YEAR_FILE, index=False)
    overall.to_csv(MEASUREMENT_OVERALL_FILE, index=False)
    event = pd.read_csv(EVENT_FILE)
    fields = graph_fields(matched, event)
    records, years, scales, profile = summarize_effects(fields)
    fields.to_csv(FIELD_FILE, index=False, compression="gzip")
    records.to_csv(RECORD_FILE, index=False)
    years.to_csv(YEAR_EFFECT_FILE, index=False)
    scales.to_csv(SCALE_FILE, index=False)
    profile.to_csv(PROFILE_FILE, index=False)
    write_tex(by_year, overall, scales, profile)
    outputs = (MATCHED_FILE, MEASUREMENT_YEAR_FILE, MEASUREMENT_OVERALL_FILE,
               FIELD_FILE, RECORD_FILE, YEAR_EFFECT_FILE, SCALE_FILE, PROFILE_FILE, TEX_FILE)
    audit = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "post-analysis non-development external measurement and effect check",
        "years": list(YEARS), "qualified_station_year_threshold": 400,
        "effect_labels": "frozen 121-site ERA5-Land peak times and quartile labels",
        "minimum_stations_per_effect_field": 10,
        "minimum_available_fields_per_regime_record": 1,
        "station_supported_bandwidths_km": sorted(BROAD),
        "matched_hours": len(matched), "matched_stations": matched.station.nunique(),
        "selected_stations": 30,
        "era_available_stations": len(era_paths),
        "era_unavailable_stations": [
            item for item in json.loads(
                ERA_PROVENANCE.read_text(encoding="utf-8")
            )["stations"] if item["status"] != "available"
        ],
        "effect_fields": fields.time_utc.nunique(), "effect_records": records.record_id.nunique(),
        "protocol_sha256": sha256(PROTOCOL),
        "inputs": [{"file": str(path.relative_to(PROJECT)), "sha256": sha256(path)}
                   for path in (*noaa_paths, *era_paths, QUALIFICATION, MANIFEST,
                                ERA_PROVENANCE, EVENT_FILE)],
        "outputs": [{"file": str(path.relative_to(PROJECT)), "sha256": sha256(path)}
                    for path in outputs],
    }
    AUDIT_FILE.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(overall.to_string(index=False))
    print(scales.to_string(index=False))
    print(profile.to_string(index=False))


if __name__ == "__main__":
    main()
