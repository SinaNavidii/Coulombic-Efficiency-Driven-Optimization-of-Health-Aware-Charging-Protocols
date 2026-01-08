from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, List, Tuple
import numpy as np
import pandas as pd

from utils.soc_sweep_analysis import get_coulombic_efficiency, shift_ce_to_100


def build_ce_dataframe(
    groups: List[Dict[str, Any]],
    processed_root: Path,
) -> pd.DataFrame:
    """
    df_ce construction:
    - Loads CE for each cell
    - Converts to %
    - Shifts to start at 100%
    - Stores long-form rows for each SOC point
    """
    rows = []

    for grp in groups:
        group_name = grp["group"]
        soh = grp["SOH"]
        legend_labels = grp["legend_labels"]
        channels_cells = grp["channels_cells"]
        cell_to_legend_label = grp["cell_to_legend_label"]

        for channel_id, cell_ids in channels_cells.items():
            for cell_id in cell_ids:
                dir_processed = Path(processed_root) / f"SOH{soh}"

                ret = get_coulombic_efficiency(
                    dir_processed=dir_processed,
                    channel_id=int(channel_id),
                    cell_id=int(cell_id),
                    soc=None,
                )

                ce_percent = ret["CE"] * 100.0
                u_ce_percent = ret["u_CE"] * 100.0
                ce_percent, u_ce_percent = shift_ce_to_100(ce_percent, u_ce_percent)

                legend_label = cell_to_legend_label[(channel_id, cell_id)]
                c_rate = legend_labels[legend_label]  # "2C", "2.5C", ...

                for soc_val, ce_val, ce_err in zip(ret["SOC"], ce_percent, u_ce_percent):
                    rows.append(
                        {
                            "Group": group_name,
                            "Channel": int(channel_id),
                            "Cell": int(cell_id),
                            "SOC (%)": float(soc_val),
                            "CE (%)": float(ce_val),
                            "CE Uncertainty (%)": float(ce_err),
                            "C-rate": c_rate,
                        }
                    )

    df_ce = pd.DataFrame(rows)
    df_ce.sort_values(by=["Group", "Channel", "Cell", "SOC (%)"], inplace=True)
    df_ce.reset_index(drop=True, inplace=True)
    return df_ce


def find_crossing_soc(socs, ces, threshold=99.0) -> float:
    """
    Finds SOC where CE drops from >= threshold to < threshold (linear interpolation).
    Returns np.nan if no crossing.
    """
    socs = np.array(socs, dtype=float)
    ces = np.array(ces, dtype=float)

    for i in range(len(ces) - 1):
        if ces[i] >= threshold and ces[i + 1] < threshold:
            frac = (ces[i] - threshold) / (ces[i] - ces[i + 1])
            return socs[i] + frac * (socs[i + 1] - socs[i])

    return np.nan


def build_summary_dataframe(
    df_ce: pd.DataFrame,
    *,
    threshold: float = 99.3,
    soh_mapping: Dict[str, float] | None = None,
) -> pd.DataFrame:
    """
    Produces df_summary:
    group by (Group, Channel, Cell, C-rate) and compute Crossing_SOC using find_crossing_soc.
    """
    if soh_mapping is None:
        soh_mapping = {"G1": 100, "G2": 90, "G3": 80, "G4": 70}

    rows = []
    for (group, channel, cell, c_rate_str), df_cell in df_ce.groupby(["Group", "Channel", "Cell", "C-rate"]):
        df_cell = df_cell.sort_values("SOC (%)")
        socs = df_cell["SOC (%)"].values
        ces = df_cell["CE (%)"].values

        crossing_soc = find_crossing_soc(socs, ces, threshold=threshold)
        numeric_c_rate = float(str(c_rate_str).replace("C", ""))

        rows.append(
            {
                "Group": group,
                "Channel": int(channel),
                "Cell": int(cell),
                "C-rate": c_rate_str,
                "Numeric C-rate": numeric_c_rate,
                "Crossing_SOC": crossing_soc,
                "SOH": soh_mapping.get(group, np.nan),
            }
        )

    df_summary = pd.DataFrame(rows)
    df_summary.sort_values(by=["Numeric C-rate", "Group", "Channel", "Cell"], inplace=True)
    df_summary.reset_index(drop=True, inplace=True)
    return df_summary


def build_soc_table(df_summary: pd.DataFrame) -> pd.DataFrame:
    """
    soc_table:
      df_summary.groupby(["Numeric C-rate","Discrete SOH"])["Crossing_SOC"].mean().unstack()
    """
    discrete_soh = {"G1": 100, "G2": 90, "G3": 80, "G4": 70}
    df = df_summary.copy()
    df["Discrete SOH"] = df["Group"].map(discrete_soh)

    soc_table = df.groupby(["Numeric C-rate", "Discrete SOH"])["Crossing_SOC"].mean().unstack()
    return soc_table


def fit_quadratic_boundaries(
    soc_table: pd.DataFrame,
    *,
    t_fit: np.ndarray | None = None,
) -> Tuple[np.ndarray, np.ndarray, Dict[float, np.ndarray]]:
    """
    Quadratic parametric fit:
      t_points = [0,1,2,3] maps to SOH = [100,90,80,70]
      fit SOC(t) quadratic per C-rate
      soh_fit = 100 - 10*t_fit
    Returns: (t_fit, soh_fit, curve_dict[c_rate] = soc_fit)
    """
    soc_table = soc_table.dropna()

    if t_fit is None:
        t_fit = np.linspace(-0.5, 5.0, 200)

    t_points = np.array([0, 1, 2, 3], dtype=float)
    soh_points = np.array([100, 90, 80, 70], dtype=float)
    soh_fit = 100.0 - 10.0 * t_fit

    curve_dict: Dict[float, np.ndarray] = {}
    for cr in soc_table.index.values.astype(float):
        soc_points = np.array(
            [soc_table.loc[cr, 100], soc_table.loc[cr, 90], soc_table.loc[cr, 80], soc_table.loc[cr, 70]],
            dtype=float,
        )
        coeffs = np.polyfit(t_points, soc_points, 2)
        curve_dict[float(cr)] = np.polyval(coeffs, t_fit)

    return t_fit, soh_fit, curve_dict


def find_boundary_crossings(soh_target, soh_fit, curve_dict):
    idx = np.abs(soh_fit - soh_target).argmin()
    return {c_rate: float(curve_dict[c_rate][idx]) for c_rate in sorted(curve_dict.keys(), reverse=True)}


def build_adaptive_profile(crossings, final_soc=90.0, final_c_rate=1.5):
    soc_steps = [0]
    c_rate_steps = [max(crossings.keys())]

    for c_rate in sorted(crossings.keys(), reverse=True):
        soc = crossings[c_rate]
        if np.isnan(soc):
            continue

        soc_steps.append(soc)
        c_rate_steps.append(c_rate)

        if c_rate != min(crossings.keys()):
            soc_steps.append(soc)
            next_c_rate = sorted([c for c in crossings.keys() if c < c_rate])[-1]
            c_rate_steps.append(next_c_rate)

    last_soc = soc_steps[-1]
    soc_steps += [last_soc, final_soc]
    c_rate_steps += [final_c_rate, final_c_rate]

    return soc_steps, c_rate_steps


def build_profiles_over_soh(soh_levels, soh_fit, curve_dict, final_soc=90.0, final_c_rate=1.5):
    profiles = {}
    for soh in soh_levels:
        crossings = find_boundary_crossings(soh, soh_fit, curve_dict)
        soc, c_rate = build_adaptive_profile(crossings, final_soc=final_soc, final_c_rate=final_c_rate)
        profiles[float(soh)] = {"SOC": soc, "C-rate": c_rate}
    return profiles


def profiles_to_df(profiles):
    rows = []
    for soh, p in profiles.items():
        rows += [{"SOH": soh, "SOC": s, "C-rate": c} for s, c in zip(p["SOC"], p["C-rate"])]
    return pd.DataFrame(rows)


def profile_at_soh(soh, soh_fit, curve_dict, final_soc=90.0, final_c_rate=1.5):
    crossings = find_boundary_crossings(soh, soh_fit, curve_dict)
    soc, c_rate = build_adaptive_profile(crossings, final_soc=final_soc, final_c_rate=final_c_rate)
    return {"SOH": float(soh), "SOC": soc, "C-rate": c_rate}
