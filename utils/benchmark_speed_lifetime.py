"""
Benchmark charging speed + lifetime metrics.

- Extracts average charging time from "Cycle N" folders.
- Extracts RPT capacity from RPT folders.
- Computes:
  (1) per-cell lifetime to 70% (with Constant 1C extrapolation if needed),
  (2) protocol lifetime from mean trajectory,
  (3) per-batch metrics, and
  (4) protocol-level metrics (batches combined).

"""

import os
import re
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


# -----------------------------
# Helpers
# -----------------------------
def _list_dirs(path, pattern=None, contains=None):
    out = []
    for d in os.listdir(path):
        p = os.path.join(path, d)
        if not os.path.isdir(p):
            continue
        if contains and contains not in d:
            continue
        if pattern and not re.match(pattern, d):
            continue
        out.append(d)
    return out


def _first_int(s):
    m = re.search(r"\d+", s)
    return int(m.group()) if m else None


def _protocol_label(profile: str) -> str:
    # keep identical to your current "Protocol" cleaning intent
    return (
        profile
        .replace(" (20% SOC)", "")
        .replace(" (10% SOC)", "")
    )


# -----------------------------
# Core extraction
# -----------------------------
def extract_charge_times(base_dir: str, batches: dict) -> pd.DataFrame:
    """
    Returns charge_df with columns:
      Batch, Profile, Cell, CycleNumber, AvgChargeTime_s
    """
    charge_records = []

    for batch_name, cfg in batches.items():
        cell_groups = cfg["cell_groups"]
        batch_path = os.path.join(base_dir, batch_name)

        cycle_folders = _list_dirs(batch_path, pattern=r"Cycle \d+")
        for folder in cycle_folders:
            cycle_num = _first_int(folder)
            folder_path = os.path.join(batch_path, folder)

            for profile, cells in cell_groups.items():
                for ch, cid in cells:
                    prefix = f"250046-{ch}-{cid}-"
                    for fn in os.listdir(folder_path):
                        if not (fn.startswith(prefix) and fn.endswith(".xlsx")):
                            continue

                        fp = os.path.join(folder_path, fn)
                        stats = pd.read_excel(fp, sheet_name=2)

                        # faithful to your logic: parse stats['Date'] as datetime
                        stats["DateTime"] = pd.to_datetime(stats["Date"], errors="coerce")

                        # charging steps: "chg" in State and Original step >= 5
                        chg = stats[
                            stats["State"].str.contains("chg", case=False, na=False)
                            & (stats["Original step"] >= 5)
                        ]
                        if chg.empty:
                            continue

                        durs = chg.groupby("Cycle")["DateTime"].agg(
                            lambda x: (x.max() - x.min()).total_seconds()
                        )
                        avg_time = durs.mean()

                        charge_records.append(
                            {
                                "Batch": batch_name,
                                "Profile": profile,
                                "Cell": f"{ch}-{cid}",
                                "CycleNumber": cycle_num,
                                "AvgChargeTime_s": avg_time,
                            }
                        )

    return pd.DataFrame(charge_records)


def extract_rpt_capacities(base_dir: str, batches: dict, nominal_capacity_ah: float) -> pd.DataFrame:
    """
    Returns cap_df with columns:
      Batch, Profile, Cell, RPT, Capacity%
    """
    cap_records = []

    for batch_name, cfg in batches.items():
        cell_groups = cfg["cell_groups"]
        use13 = cfg.get("use_step13_for_rpt0", False)
        batch_path = os.path.join(base_dir, batch_name)

        rpt_folders = _list_dirs(batch_path, contains="RPT")
        for rpt_folder in rpt_folders:
            rpt = _first_int(rpt_folder)
            if rpt is None:
                continue

            step = 13 if (rpt == 0 and use13) else 5
            folder_path = os.path.join(batch_path, rpt_folder)

            for profile, cells in cell_groups.items():
                for ch, cid in cells:
                    prefix = f"250046-{ch}-{cid}-"
                    for fn in os.listdir(folder_path):
                        if not (fn.startswith(prefix) and fn.endswith(".xlsx")):
                            continue

                        fp = os.path.join(folder_path, fn)
                        df = pd.read_excel(fp, sheet_name=3)

                        df_f = df[df["Steps"] == step].copy()
                        if df_f.empty:
                            continue

                        # faithful to your logic: pick Date(h:min:s.ms) else Date
                        if "Date(h:min:s.ms)" in df_f.columns:
                            date_col = "Date(h:min:s.ms)"
                        elif "Date" in df_f.columns:
                            date_col = "Date"
                        else:
                            continue

                        t_as_dt = pd.to_datetime(df_f[date_col], errors="coerce")
                        if t_as_dt.isna().all():
                            t_as_td = pd.to_timedelta(df_f[date_col], errors="coerce")
                            if t_as_td.isna().all():
                                continue
                            times_sec = (t_as_td - t_as_td.iloc[0]).dt.total_seconds().values
                        else:
                            times_sec = (t_as_dt - t_as_dt.iloc[0]).dt.total_seconds().values

                        df_f["Time(s)"] = times_sec

                        cur = df_f["Current(A)"].to_numpy()
                        tm = df_f["Time(s)"].to_numpy()

                        # NOTE: your newer version uses trapz; we keep that
                        cap_Ah = abs(np.trapz(cur, tm)) / 3600.0
                        cap_pct = cap_Ah / nominal_capacity_ah * 100.0

                        cap_records.append(
                            {
                                "Batch": batch_name,
                                "Profile": profile,
                                "Cell": f"{ch}-{cid}",
                                "RPT": rpt,
                                "Capacity%": cap_pct,
                            }
                        )

    return pd.DataFrame(cap_records)


# -----------------------------
# Metric computation
# -----------------------------
def compute_benchmark_speed_metrics(
    base_dir: str,
    batches: dict,
    nominal_capacity_ah: float = 47 / 1000,
    eol_pct: float = 70.0,
) -> dict:
    """
    Returns a dict with:
      charge_df, cap_df, life_df,
      metrics_batch, metrics_combined,
      lifetimes_profile, life_std_profile
    """
    # 1) charging time stats per batch/profile
    charge_df = extract_charge_times(base_dir, batches)

    charge_stats = (
        charge_df
        .groupby(["Batch", "Profile"])["AvgChargeTime_s"]
        .agg(["mean", "std"])
        .rename(columns={"mean": "ChargeTimeMean_s", "std": "ChargeTimeStd_s"})
        .reset_index()
    )

    # 2) RPT capacities
    cap_df = extract_rpt_capacities(base_dir, batches, nominal_capacity_ah)

    # 3) per-cell lifetimes (for std), with Constant 1C extrapolation
    lifetime_records = []
    for (batch, profile, cell), group in cap_df.groupby(["Batch", "Profile", "Cell"]):
        group = group.sort_values("RPT")
        rpts = group["RPT"].to_numpy(dtype=float)
        caps = group["Capacity%"].to_numpy(dtype=float)

        below = rpts[caps <= eol_pct]

        if profile == "Constant 1C" and below.size == 0:
            n_fit = min(len(rpts), 10)
            x_fit = rpts[-n_fit:]
            y_fit = caps[-n_fit:]

            if len(x_fit) >= 2:
                a, b = np.polyfit(x_fit, y_fit, 1)
                if a < 0:
                    rpt_eol = (eol_pct - b) / a
                    rpt_eol = max(rpt_eol, rpts[-1])
                else:
                    rpt_eol = rpts[-1]
            else:
                rpt_eol = rpts[-1]
        else:
            rpt_eol = float(below[0]) if below.size > 0 else float(rpts[-1])

        lifetime_records.append(
            {
                "Batch": batch,
                "Profile": profile,
                "Cell": cell,
                "Lifetime_cycles": rpt_eol * 4.0,
            }
        )

    life_df = pd.DataFrame(lifetime_records)

    life_std_profile = (
        life_df
        .groupby("Profile")["Lifetime_cycles"]
        .std()
        .rename("LifeStd_cycles")
        .reset_index()
    )

    # 4) canonical lifetime from mean trajectory
    mean_caps_profile = (
        cap_df
        .groupby(["Profile", "RPT"])["Capacity%"]
        .mean()
        .reset_index()
    )

    lifetimes_profile = {}
    for profile, grp in mean_caps_profile.groupby("Profile"):
        grp = grp.sort_values("RPT")
        rpts = grp["RPT"].to_numpy(dtype=float)
        caps = grp["Capacity%"].to_numpy(dtype=float)

        below = rpts[caps <= eol_pct]

        if profile == "Constant 1C" and below.size == 0:
            n_fit = min(len(rpts), 10)
            x_fit = rpts[-n_fit:]
            y_fit = caps[-n_fit:]

            if len(x_fit) >= 2:
                a, b = np.polyfit(x_fit, y_fit, 1)
                if a < 0:
                    rpt_eol = (eol_pct - b) / a
                    rpt_eol = max(rpt_eol, rpts[-1])
                else:
                    rpt_eol = rpts[-1]
            else:
                rpt_eol = rpts[-1]
        else:
            rpt_eol = float(below[0]) if below.size > 0 else float(rpts[-1])

        lifetimes_profile[profile] = rpt_eol * 4.0

    # 5) per-batch metrics 
    life_stats_batch = (
        life_df
        .groupby(["Batch", "Profile"])
        .agg(
            LifeMean_cycles=("Lifetime_cycles", "mean"),
            LifeStd_cycles=("Lifetime_cycles", "std"),
        )
        .reset_index()
    )

    metrics_batch = pd.merge(
        charge_stats,
        life_stats_batch,
        on=["Batch", "Profile"],
        how="inner",
    )

    # 6) protocol-level metrics (batches combined + canonical lifetime)
    metrics_batch["Protocol"] = metrics_batch["Profile"].apply(_protocol_label)

    charge_protocol = (
        metrics_batch
        .groupby("Protocol", as_index=False)
        .agg(
            ChargeTimeMean_s=("ChargeTimeMean_s", "mean"),
            ChargeTimeStd_s=("ChargeTimeStd_s", "mean"),
        )
    )

    life_mean_protocol = (
        pd.Series(lifetimes_profile, name="LifeMean_cycles")
        .reset_index()
        .rename(columns={"index": "Protocol"})
    )
    life_std_protocol = life_std_profile.rename(columns={"Profile": "Protocol"})

    metrics_combined = (
        charge_protocol
        .merge(life_mean_protocol, on="Protocol", how="left")
        .merge(life_std_protocol, on="Protocol", how="left")
    )

    return {
        "charge_df": charge_df,
        "cap_df": cap_df,
        "life_df": life_df,
        "metrics_batch": metrics_batch,
        "metrics_combined": metrics_combined,
        "lifetimes_profile": lifetimes_profile,
        "life_std_profile": life_std_profile,
    }
