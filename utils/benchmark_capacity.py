import os
import re
import warnings
from typing import Dict, Any, Tuple, Optional, List

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


def _find_rpt_folders(batch_path: str) -> List[str]:
    return [
        d for d in os.listdir(batch_path)
        if os.path.isdir(os.path.join(batch_path, d)) and ("RPT" in d)
    ]


def _extract_rpt_num(folder_name: str) -> int:
    m = re.search(r"\d+", folder_name)
    if m is None:
        raise ValueError(f"Could not parse RPT number from folder name: {folder_name}")
    return int(m.group())


def _time_to_seconds(df_f: pd.DataFrame) -> np.ndarray:
    if "Date(h:min:s.ms)" in df_f.columns:
        date_col = "Date(h:min:s.ms)"
    elif "Date" in df_f.columns:
        date_col = "Date"
    else:
        raise KeyError("Neither 'Date(h:min:s.ms)' nor 'Date' column found")

    t_series = pd.to_datetime(df_f[date_col], errors="coerce")
    if t_series.isna().all():
        td = pd.to_timedelta(df_f[date_col], errors="coerce")
        if td.isna().all():
            raise ValueError(f"Could not parse time from column '{date_col}'")
        return (td - td.iloc[0]).dt.total_seconds().values

    return (t_series - t_series.iloc[0]).dt.total_seconds().values


def collect_benchmark_rpt_capacities(
    *,
    base_path: str,
    batches: Dict[str, Dict[str, Any]],
    nominal_capacity_ah: float = 47 / 1000,
) -> pd.DataFrame:
    """
    Collects capacity data from benchmark RPT folders.
    Returns a tidy DataFrame with:
      ['Batch','Profile','RPT','Cell','Capacity (Ah)','Capacity (%)']
    """
    results = []

    for batch_dir, cfg in batches.items():
        batch_path = os.path.join(base_path, batch_dir)
        cell_groups = cfg["cell_groups"]

        is_first = (batch_dir == "Benchmark_Tests_Batch1")

        rpt_folders = _find_rpt_folders(batch_path)

        for rpt_folder in rpt_folders:
            rpt_num = _extract_rpt_num(rpt_folder)

            # preserved exactly
            if is_first:
                step_to_check = 13 if rpt_num == 0 else 5
            else:
                step_to_check = 5

            rpt_path = os.path.join(batch_path, rpt_folder)

            for profile_name, cells in cell_groups.items():
                for ch, cid in cells:
                    prefix = f"250046-{ch}-{cid}-"

                    for fn in os.listdir(rpt_path):
                        if (not fn.startswith(prefix)) or (not fn.endswith(".xlsx")):
                            continue

                        fp = os.path.join(rpt_path, fn)

                        try:
                            df = pd.read_excel(fp, sheet_name=3)
                            df_f = df[df["Steps"] == step_to_check].copy()
                            if df_f.empty:
                                continue

                            times_sec = _time_to_seconds(df_f)
                            df_f["Time(s)"] = times_sec

                            current = df_f["Current(A)"].values
                            times = df_f["Time(s)"].values

                            cap_ah = abs((current * np.diff(np.append(times, times[-1]))).sum()) / 3600.0
                            cap_pct = (cap_ah / nominal_capacity_ah) * 100.0

                            results.append(
                                {
                                    "Batch": batch_dir,
                                    "Profile": profile_name,
                                    "RPT": rpt_num,
                                    "Cell": f"{ch}-{cid}",
                                    "Capacity (Ah)": cap_ah,
                                    "Capacity (%)": cap_pct,
                                }
                            )

                        except Exception as e:
                            print(f"Error processing {fp}: {e}")

    return pd.DataFrame(results)
