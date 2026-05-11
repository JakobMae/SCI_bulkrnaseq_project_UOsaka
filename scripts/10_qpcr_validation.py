#!/usr/bin/env python3

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


warnings.filterwarnings(
    "ignore",
    message="Precision loss occurred in moment calculation"
)


INPUT_DIR = Path("data/qpcr")
OUTPUT_DIR = Path("results/qpcr")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


SAMPLE_ORDER = [
    "KOsham_dpi7_1", "KOsham_dpi7_2", "KOsham_dpi7_3", "KOsham_dpi7_4", "KOsham_dpi7_5",
    "KOSCI_dpi7_1", "KOSCI_dpi7_2", "KOSCI_dpi7_3", "KOSCI_dpi7_4", "KOSCI_dpi7_5",
    "DTR_Ns_sham_dpi14_1", "DTR_Ns_sham_dpi14_2", "DTR_Ns_sham_dpi14_3", "DTR_Ns_sham_dpi14_4", "DTR_Ns_sham_dpi14_5",
    "DTR_Ns_SCI_dpi14_1", "DTR_Ns_SCI_dpi14_2", "DTR_Ns_SCI_dpi14_3", "DTR_Ns_SCI_dpi14_4", "DTR_Ns_SCI_dpi14_5",
    "DTR_DT_sham_dpi14_1", "DTR_DT_sham_dpi14_2", "DTR_DT_sham_dpi14_3", "DTR_DT_sham_dpi14_4", "DTR_DT_sham_dpi14_5",
    "DTR_DT_SCI_dpi14_1", "DTR_DT_SCI_dpi14_2", "DTR_DT_SCI_dpi14_3", "DTR_DT_SCI_dpi14_4", "DTR_DT_SCI_dpi14_5",
    "WS_1", "WS_2", "WS_3", "WS_4", "WS_5", "WS_6",
    "WI_1", "WI_2", "WI_3", "WI_4", "WI_5", "WI_6",
]

SAMPLE_MAP = {f"Sample {i + 1}": sample for i, sample in enumerate(SAMPLE_ORDER)}

TARGET_RENAME = {
    "Slcla2": "Slc1a2",
    "NxF1": "Nxf1",
    "TFNAIP": "Tnfaip3",
    "CALN": "Caln1",
    "CCNB": "Ccnb1",
}

HOUSEKEEPERS = ["GAPDH", "RPS18"]

PAIRWISE_COMPARISONS = [
    ("WT_sham_dpi7", "WT_SCI_dpi7"),
    ("KO_sham_dpi7", "KO_SCI_dpi7"),
    ("DTR_Ns_SCI_dpi14", "DTR_DT_SCI_dpi14"),
]


def parse_group(sample_id: str) -> str:
    if sample_id.startswith("KOsham"):
        return "KO_sham_dpi7"
    if sample_id.startswith("KOSCI"):
        return "KO_SCI_dpi7"
    if sample_id.startswith("DTR_Ns_sham"):
        return "DTR_Ns_sham_dpi14"
    if sample_id.startswith("DTR_Ns_SCI"):
        return "DTR_Ns_SCI_dpi14"
    if sample_id.startswith("DTR_DT_sham"):
        return "DTR_DT_sham_dpi14"
    if sample_id.startswith("DTR_DT_SCI"):
        return "DTR_DT_SCI_dpi14"
    if sample_id.startswith("WS_"):
        return "WT_sham_dpi7"
    if sample_id.startswith("WI_"):
        return "WT_SCI_dpi7"
    raise ValueError(f"Could not parse group for sample: {sample_id}")


def load_results_table(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Results", header=None)

    header_row = None
    for i in range(min(100, len(raw))):
        vals = raw.iloc[i].astype(str).tolist()
        if "Well" in vals and "Sample Name" in vals and "Target Name" in vals:
            header_row = i
            break

    if header_row is None:
        raise ValueError(f"Could not find Results header row in {path}")

    return pd.read_excel(path, sheet_name="Results", header=header_row)


files = sorted([p for p in INPUT_DIR.glob("*.xls*") if "Primer" not in p.name])

rows = []

for file in files:
    df = load_results_table(file)

    df = df[df["Task"].astype(str).str.upper() == "UNKNOWN"].copy()
    df = df[df["Sample Name"].isin(SAMPLE_MAP)].copy()

    df["source_file"] = file.name
    df["sample_id"] = df["Sample Name"].map(SAMPLE_MAP)
    df["group"] = df["sample_id"].map(parse_group)
    df["target"] = df["Target Name"].replace(TARGET_RENAME).astype(str)
    df["ct"] = pd.to_numeric(df["CT"], errors="coerce")
    df["tm1"] = pd.to_numeric(df["Tm1"], errors="coerce")
    df["efficiency"] = pd.to_numeric(df["Efficiency"], errors="coerce")

    if "NOAMP" in df.columns:
        df = df[df["NOAMP"].fillna("N").astype(str) != "Y"]
    if "EXPFAIL" in df.columns:
        df = df[df["EXPFAIL"].fillna("N").astype(str) != "Y"]

    rows.append(
        df[
            [
                "source_file",
                "sample_id",
                "group",
                "target",
                "ct",
                "tm1",
                "efficiency",
            ]
        ]
    )


raw_df = pd.concat(rows, ignore_index=True)
raw_df = raw_df[raw_df["ct"].notna()].copy()

processed = (
    raw_df
    .groupby(["sample_id", "group", "target"], as_index=False)
    .agg(
        ct_mean=("ct", "mean"),
        ct_sd=("ct", "std"),
        n_reps=("ct", "size"),
        tm1_mean=("tm1", "mean"),
        efficiency_mean=("efficiency", "mean"),
    )
)

hk = (
    processed[processed["target"].isin(HOUSEKEEPERS)]
    .groupby("sample_id", as_index=False)
    .agg(hk_ct=("ct_mean", "mean"))
)

processed = processed.merge(hk, on="sample_id", how="left")
processed["delta_ct"] = processed["ct_mean"] - processed["hk_ct"]

control_group_by_group = {
    "WT_sham_dpi7": "WT_sham_dpi7",
    "WT_SCI_dpi7": "WT_sham_dpi7",
    "KO_sham_dpi7": "KO_sham_dpi7",
    "KO_SCI_dpi7": "KO_sham_dpi7",
    "DTR_Ns_sham_dpi14": "DTR_Ns_sham_dpi14",
    "DTR_Ns_SCI_dpi14": "DTR_Ns_sham_dpi14",
    "DTR_DT_sham_dpi14": "DTR_DT_sham_dpi14",
    "DTR_DT_SCI_dpi14": "DTR_DT_sham_dpi14",
}

control_means = (
    processed
    .assign(control_group=lambda x: x["group"].map(control_group_by_group))
    .groupby(["target", "control_group"], as_index=False)
    .agg(control_delta_ct=("delta_ct", "mean"))
)

processed["control_group"] = processed["group"].map(control_group_by_group)

processed = processed.merge(
    control_means,
    on=["target", "control_group"],
    how="left",
)

processed["delta_delta_ct"] = processed["delta_ct"] - processed["control_delta_ct"]
processed["fold_change"] = 2 ** (-processed["delta_delta_ct"])

stats_rows = []

for target in sorted(processed["target"].unique()):
    for group_1, group_2 in PAIRWISE_COMPARISONS:
        vals_1 = processed[
            (processed["target"] == target) &
            (processed["group"] == group_1)
        ]["delta_ct"]

        vals_2 = processed[
            (processed["target"] == target) &
            (processed["group"] == group_2)
        ]["delta_ct"]

        if len(vals_1) < 2 or len(vals_2) < 2:
            continue

        _, p_value = ttest_ind(vals_1, vals_2, equal_var=False)

        stats_rows.append({
            "target": target,
            "group_1": group_1,
            "group_2": group_2,
            "n_1": len(vals_1),
            "n_2": len(vals_2),
            "mean_delta_ct_1": vals_1.mean(),
            "mean_delta_ct_2": vals_2.mean(),
            "log2_fc_qpcr": -(vals_2.mean() - vals_1.mean()),
            "p_value": p_value,
        })

stats = pd.DataFrame(stats_rows)

raw_df.to_csv(OUTPUT_DIR / "qpcr_long_raw.tsv", sep="\t", index=False)
processed.to_csv(OUTPUT_DIR / "qpcr_processed.tsv", sep="\t", index=False)
stats.to_csv(OUTPUT_DIR / "qpcr_stats.tsv", sep="\t", index=False)

print("qPCR validation analysis complete.")
