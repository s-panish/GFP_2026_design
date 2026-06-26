# ============================================================
# GeneMeow v3: process ColabFold results archive for WT + final 6
# Memory-safe Colab script
# ============================================================
# Upload the ColabFold results zip generated from:
# outputs/colabfold_input_WT_plus_6designs.fasta
# Optional but recommended: upload design_report_brightness_gate.csv too.
# This script creates:
#   stage4_colabfold_metrics_final6.csv
#   final6_structural_summary.csv
#   submission_GeneMeow.csv
#   figures zip
#   minimal package zip
# ============================================================

import os
import re
import glob
import json
import shutil
import zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from google.colab import files
    IN_COLAB = True
except Exception:
    IN_COLAB = False

TEAM_NAME = "GeneMeow"
OUTPUT_DIR = "/content/GeneMeow_v3_colabfold_outputs"
EXTRACT_DIR = os.path.join(OUTPUT_DIR, "colabfold_results_extracted")
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")
RAW_FIGURES_DIR = os.path.join(FIGURES_DIR, "colabfold_raw_plots_final6")
MINIMAL_PACKAGE_DIR = "/content/GeneMeow_v3_minimal_submission_package"

CHROMOPHORE_POSITIONS = [65, 66, 67]
CORE_START = 1
CORE_END = 220
ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")

FINAL_ORDER = [
    "Seq1_conservative_surface_D19E_R73L_H231Y",
    "Seq2_balanced_core_D19E_R73L_H231Y_K166E_N212D",
    "Seq3_MAIN_balanced_D19E_R73L_H231Y_N198D_N212D_K166E",
    "Seq4_BEST_thermal_D19E_R73L_H231Y_N198D_N212D_K166E_K156E_K101E",
    "Seq5_BEST_brightness_D19E_R73L_H231Y_Y237N",
    "Seq6_INSURANCE_supercharge_K166E_N212D_K101E_K156E_N198D",
]

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(EXTRACT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RAW_FIGURES_DIR, exist_ok=True)


def normalize_sequence(seq):
    return re.sub(r"\s+", "", str(seq)).upper()


def validate_submission_sequence(seq_id, sequence):
    seq = normalize_sequence(sequence)
    errors = []
    if not seq.startswith("M"):
        errors.append("does not start with M")
    if not 220 <= len(seq) <= 250:
        errors.append(f"length is {len(seq)}, expected 220-250 aa")
    invalid = sorted(set(seq) - ALLOWED_AA)
    if invalid:
        errors.append("invalid characters: " + "".join(invalid))
    if seq[64:67] != "TYG":
        errors.append(f"chromophore 65-67 is {seq[64:67]}, expected TYG")
    if errors:
        raise ValueError(f"Seq_ID {seq_id}: " + "; ".join(errors))
    return seq


def recursively_extract_zip(zip_path, extract_to):
    if os.path.exists(extract_to):
        shutil.rmtree(extract_to)
    os.makedirs(extract_to, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(extract_to)
    nested = glob.glob(os.path.join(extract_to, "**", "*.zip"), recursive=True)
    for zp in nested:
        nested_dir = os.path.join(os.path.dirname(zp), os.path.splitext(os.path.basename(zp))[0] + "_extracted")
        os.makedirs(nested_dir, exist_ok=True)
        try:
            with zipfile.ZipFile(zp, "r") as z:
                z.extractall(nested_dir)
        except zipfile.BadZipFile:
            pass


def extract_query_id_from_filename(path):
    name = os.path.basename(path)
    name = re.sub(r"_scores_rank_001.*$", "", name)
    name = re.sub(r"_unrelaxed_rank_001.*$", "", name)
    name = re.sub(r"_relaxed_rank_001.*$", "", name)
    name = re.sub(r"_predicted_aligned_error.*$", "", name)
    name = re.sub(r"_pae.*$", "", name)
    name = re.sub(r"_plddt.*$", "", name)
    name = re.sub(r"_coverage.*$", "", name)
    name = re.sub(r"\.a3m$", "", name)
    return name


def find_files(result_dir, patterns):
    out = []
    for p in patterns:
        out.extend(glob.glob(os.path.join(result_dir, "**", p), recursive=True))
    return sorted(set(out))


def read_scores_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    plddt = np.array(data.get("plddt"), dtype=float)
    if plddt.size == 0:
        raise ValueError(f"No pLDDT in {path}")
    ptm = data.get("ptm", np.nan)
    chrom_idx = [p - 1 for p in CHROMOPHORE_POSITIONS]
    return {
        "Mean_pLDDT": float(np.mean(plddt)),
        "pTM": float(ptm) if ptm is not None else np.nan,
        "Chromophore_pLDDT_Mean_65_67": float(np.mean(plddt[chrom_idx])),
        "Chromophore_pLDDT_Min_65_67": float(np.min(plddt[chrom_idx])),
        "Core_pLDDT_Mean_1_220": float(np.mean(plddt[CORE_START - 1:CORE_END])),
        "CTerm_pLDDT_Mean_221_End": float(np.mean(plddt[220:])),
        "plddt_array": plddt,
    }


def parse_mutation_positions(mutation_string):
    return [int(x) for x in re.findall(r"[A-Z](\d+)[A-Z]", str(mutation_string))]


def read_ca_coordinates_from_pdb(path, start=1, end=220):
    coords = []
    residues = []
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                resseq = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            if start <= resseq <= end:
                residues.append(resseq)
                coords.append([x, y, z])
    return np.array(coords, dtype=float), residues


def kabsch_rmsd(P, Q):
    if P.shape != Q.shape:
        raise ValueError(f"Shape mismatch: {P.shape} vs {Q.shape}")
    P_cent = P - P.mean(axis=0)
    Q_cent = Q - Q.mean(axis=0)
    C = np.dot(P_cent.T, Q_cent)
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(np.dot(V, Wt)))
    D = np.diag([1.0, 1.0, d])
    U = np.dot(np.dot(V, D), Wt)
    P_rot = np.dot(P_cent, U)
    diff = P_rot - Q_cent
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}


def read_sequence_from_pdb(path):
    residues = []
    seen = set()
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            resname = line[17:20].strip()
            chain_id = line[21].strip()
            resseq = line[22:26].strip()
            key = (chain_id, resseq)
            if key in seen:
                continue
            seen.add(key)
            residues.append(THREE_TO_ONE.get(resname, "X"))
    seq = "".join(residues)
    return None if "X" in seq else seq


def find_matching_key(query_id, key_collection):
    if query_id in key_collection:
        return query_id
    contains = [key for key in key_collection if query_id in key or key in query_id]
    return contains[0] if contains else None


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = normalize_sequence("".join(buf))
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = normalize_sequence("".join(buf))
    return seqs


def save_barplot(df, x_col, y_col, title, ylabel, output_path, rotate=30):
    plt.figure(figsize=(9, 5))
    plt.bar(df[x_col].astype(str), df[y_col])
    plt.title(title)
    plt.xlabel(x_col)
    plt.ylabel(ylabel)
    plt.xticks(rotation=rotate, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
    plt.close()


print("Upload the ColabFold results zip. Optional: upload design_report_brightness_gate.csv and submission_GeneMeow.csv.")
if not IN_COLAB:
    raise RuntimeError("This script is intended for Google Colab.")

uploaded = files.upload()
uploaded_paths = [f"/content/{name}" for name in uploaded.keys()]
zip_files = [p for p in uploaded_paths if p.lower().endswith(".zip")]
csv_files = [p for p in uploaded_paths if p.lower().endswith(".csv")]
fasta_files = [p for p in uploaded_paths if p.lower().endswith((".fasta", ".fa"))]

if not zip_files:
    raise FileNotFoundError("No ColabFold results zip was uploaded.")

colabfold_zip = sorted(zip_files, key=lambda p: os.path.getsize(p), reverse=True)[0]
print("Selected ColabFold results archive:", colabfold_zip)

report_path = None
submission_path_uploaded = None
for p in csv_files:
    b = os.path.basename(p)
    if "design_report" in b:
        report_path = p
    if "submission" in b:
        submission_path_uploaded = p

fasta_seqs = {}
if fasta_files:
    fasta_seqs = read_fasta(fasta_files[0])

recursively_extract_zip(colabfold_zip, EXTRACT_DIR)
json_files = find_files(EXTRACT_DIR, ["*scores_rank_001*.json", "*_scores_rank_001_*.json"])
pdb_files = find_files(EXTRACT_DIR, ["*unrelaxed_rank_001*.pdb", "*_unrelaxed_rank_001_*.pdb", "*relaxed_rank_001*.pdb", "*_relaxed_rank_001_*.pdb"])
a3m_files = find_files(EXTRACT_DIR, ["*.a3m"])
print(f"Rank 001 JSON files found: {len(json_files)}")
print(f"Rank 001 PDB files found: {len(pdb_files)}")
print(f"A3M files found: {len(a3m_files)}")

json_by_key = {extract_query_id_from_filename(p): p for p in json_files}
pdb_by_key = {extract_query_id_from_filename(p): p for p in pdb_files}
wt_keys = [k for k in pdb_by_key if "WT_sfGFP" in k or "sfGFP_WT" in k]
if not wt_keys:
    raise ValueError("WT sfGFP PDB not found. Make sure WT_sfGFP was included in the ColabFold input FASTA.")
wt_key = wt_keys[0]
wt_coords, _ = read_ca_coordinates_from_pdb(pdb_by_key[wt_key], CORE_START, CORE_END)
if wt_coords.size == 0:
    raise ValueError("Could not read WT CA coordinates.")

report_df = pd.read_csv(report_path) if report_path else pd.DataFrame()
report_by_seq = {}
if not report_df.empty and "Seq_ID" in report_df.columns:
    for _, row in report_df.iterrows():
        report_by_seq[int(row["Seq_ID"])] = row.to_dict()

records = []
for seq_i, query_id in enumerate(FINAL_ORDER, start=1):
    json_key = find_matching_key(query_id, json_by_key.keys())
    pdb_key = find_matching_key(query_id, pdb_by_key.keys())
    if json_key is None or pdb_key is None:
        raise ValueError(f"Missing ColabFold rank_001 files for {query_id}")
    score = read_scores_json(json_by_key[json_key])
    plddt = score["plddt_array"]
    coords, _ = read_ca_coordinates_from_pdb(pdb_by_key[pdb_key], CORE_START, CORE_END)
    core_rmsd = kabsch_rmsd(coords, wt_coords) if coords.shape == wt_coords.shape else np.nan
    meta = report_by_seq.get(seq_i, {})
    muts = meta.get("Mutations", "")
    mut_pos = parse_mutation_positions(muts)
    mut_values = [plddt[p - 1] for p in mut_pos if 1 <= p <= len(plddt)]
    sequence = meta.get("Sequence", "")
    if not sequence and query_id in fasta_seqs:
        sequence = fasta_seqs[query_id]
    if not sequence:
        sequence = read_sequence_from_pdb(pdb_by_key[pdb_key])
    sequence = normalize_sequence(sequence)
    record = {
        "Seq_ID": seq_i,
        "ColabFold_Query_ID": query_id,
        "Design_Label": meta.get("Design_Label", query_id),
        "Mutations": muts,
        "Sequence": sequence,
        "Length": len(sequence),
        "Brightness_Delta_vs_sfGFP": meta.get("Brightness_Delta_vs_sfGFP", np.nan),
        "Brightness_Verdict": meta.get("Brightness_Verdict", ""),
        "Measured": meta.get("Measured", np.nan),
        "Estimated": meta.get("Estimated", np.nan),
        "Ref_Mismatch": meta.get("Ref_Mismatch", np.nan),
        "Near_Chromophore": meta.get("Near_Chromophore", np.nan),
        "Net_Charge": meta.get("Net_Charge", np.nan),
        "Delta_Charge": meta.get("Delta_Charge", np.nan),
        "Mean_pLDDT": score["Mean_pLDDT"],
        "pTM": score["pTM"],
        "Chromophore_pLDDT_Mean_65_67": score["Chromophore_pLDDT_Mean_65_67"],
        "Chromophore_pLDDT_Min_65_67": score["Chromophore_pLDDT_Min_65_67"],
        "Core_pLDDT_Mean_1_220": score["Core_pLDDT_Mean_1_220"],
        "CTerm_pLDDT_Mean_221_End": score["CTerm_pLDDT_Mean_221_End"],
        "Mutation_Site_pLDDT_Mean": float(np.mean(mut_values)) if mut_values else np.nan,
        "Mutation_Site_pLDDT_Min": float(np.min(mut_values)) if mut_values else np.nan,
        "Core_RMSD_to_WT_A": core_rmsd,
        "ColabFold_JSON": json_by_key[json_key],
        "ColabFold_PDB": pdb_by_key[pdb_key],
    }
    records.append(record)

metrics_df = pd.DataFrame(records)
metrics_df["Structural_Pass"] = (
    (metrics_df["Mean_pLDDT"] >= 90) &
    (metrics_df["Chromophore_pLDDT_Min_65_67"] >= 85) &
    (metrics_df["Core_RMSD_to_WT_A"] <= 1.0)
)

stage4_metrics_path = os.path.join(OUTPUT_DIR, "stage4_colabfold_metrics_final6.csv")
metrics_df.to_csv(stage4_metrics_path, index=False)
print("ColabFold metrics:")
display(metrics_df[["Seq_ID", "Design_Label", "Mutations", "Mean_pLDDT", "pTM", "Chromophore_pLDDT_Min_65_67", "Mutation_Site_pLDDT_Min", "Core_RMSD_to_WT_A", "Structural_Pass"]])

submission_records = []
for _, row in metrics_df.iterrows():
    seq = validate_submission_sequence(int(row["Seq_ID"]), row["Sequence"])
    submission_records.append({"Team_Name": TEAM_NAME, "Seq_ID": int(row["Seq_ID"]), "Sequence": seq})
submission_df = pd.DataFrame(submission_records)
submission_path = os.path.join(OUTPUT_DIR, "submission_GeneMeow.csv")
submission_df.to_csv(submission_path, index=False)
print("Submission:")
display(submission_df)

summary_path = os.path.join(OUTPUT_DIR, "final6_structural_summary.csv")
metrics_df.drop(columns=["ColabFold_JSON", "ColabFold_PDB"]).to_csv(summary_path, index=False)

plot_df = metrics_df.copy()
plot_df["Seq_Label"] = plot_df["Seq_ID"].apply(lambda x: f"Seq {x}")
save_barplot(plot_df, "Seq_Label", "Mean_pLDDT", "Mean pLDDT for final six", "Mean pLDDT", os.path.join(FIGURES_DIR, "final6_mean_plddt.png"))
save_barplot(plot_df, "Seq_Label", "Chromophore_pLDDT_Min_65_67", "Minimum pLDDT at chromophore positions 65-67", "Min pLDDT 65-67", os.path.join(FIGURES_DIR, "final6_chromophore_plddt_min.png"))
save_barplot(plot_df, "Seq_Label", "Core_RMSD_to_WT_A", "Core RMSD to WT sfGFP", "Core RMSD to WT, A", os.path.join(FIGURES_DIR, "final6_core_rmsd_to_wt.png"))
save_barplot(plot_df, "Seq_Label", "Brightness_Delta_vs_sfGFP", "Corrected brightness delta vs sfGFP", "Brightness delta", os.path.join(FIGURES_DIR, "final6_brightness_delta.png"))

for query_id in FINAL_ORDER:
    for pattern in [
        os.path.join(EXTRACT_DIR, "**", f"{query_id}*plddt*.png"),
        os.path.join(EXTRACT_DIR, "**", f"{query_id}*pae*.png"),
        os.path.join(EXTRACT_DIR, "**", f"{query_id}*coverage*.png"),
    ]:
        for src in glob.glob(pattern, recursive=True):
            dst = os.path.join(RAW_FIGURES_DIR, f"{query_id}_{os.path.basename(src)}")
            shutil.copy2(src, dst)

figures_zip_path = shutil.make_archive("/content/GeneMeow_v3_figures_for_PDF", "zip", FIGURES_DIR)
if os.path.exists(MINIMAL_PACKAGE_DIR):
    shutil.rmtree(MINIMAL_PACKAGE_DIR)
os.makedirs(MINIMAL_PACKAGE_DIR, exist_ok=True)
for src in [submission_path, summary_path, stage4_metrics_path, figures_zip_path]:
    shutil.copy2(src, os.path.join(MINIMAL_PACKAGE_DIR, os.path.basename(src)))
minimal_zip_path = shutil.make_archive("/content/GeneMeow_v3_minimal_submission_package", "zip", MINIMAL_PACKAGE_DIR)

print("Saved:")
for p in [submission_path, summary_path, stage4_metrics_path, figures_zip_path, minimal_zip_path]:
    print(p)

files.download(submission_path)
files.download(summary_path)
files.download(stage4_metrics_path)
files.download(figures_zip_path)
files.download(minimal_zip_path)
