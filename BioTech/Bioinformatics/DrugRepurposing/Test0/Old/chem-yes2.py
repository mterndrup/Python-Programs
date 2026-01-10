import argparse
import sqlite3
from Bio import SeqIO
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
from datetime import datetime
from rdkit import RDLogger
import itertools

# ======================
# Logging
# ======================
def create_logger(script_name):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"log_{timestamp}.txt"

    def log(msg):
        time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        entry = f"[{time_str}] {msg}"
        print(entry)
        with open(log_file, "a") as f:
            f.write(entry + "\n")
    log(f"Script: {script_name}")
    log(f"Start time: {datetime.now()}")
    return log, log_file

# ======================
# Fingerprint helper
# ======================
RDLogger.DisableLog('rdApp.*')  # suppress RDKit warnings

def mol_to_fp(mol):
    """Convert RDKit Mol object to Morgan fingerprint (bit vector)."""
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

# ======================
# Fetch TB drug fingerprints
# ======================
def get_tb_fps(tb_drugs, conn, log):
    cursor = conn.cursor()
    placeholders = ",".join("?" for _ in tb_drugs)
    query = f"""
        SELECT md.pref_name, cs.canonical_smiles, md.molregno
        FROM molecule_dictionary md
        JOIN compound_structures cs ON md.molregno = cs.molregno
        WHERE LOWER(md.pref_name) IN ({placeholders})
    """
    cursor.execute(query, [d.lower() for d in tb_drugs])
    tb_fps = {}
    log("\n=== TB DRUG FINGERPRINTS ===\n============================")
    for pref_name, smiles, molregno in cursor.fetchall():
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = mol_to_fp(mol)
            tb_fps[pref_name] = {"fp": fp, "molregno": molregno}
            log(f"[INFO] TB Drug: {pref_name}, Molregno: {molregno}, FP length: {fp.GetNumBits()}")
        else:
            log(f"[WARN] Invalid SMILES for {pref_name}")
    if not tb_fps:
        log("[WARN] No TB drug fingerprints found! Check drug names.")
    return tb_fps

# ======================
# Fetch molecules in batches
# ======================
def fetch_all_molecules(fasta_path, min_phase, conn, log, batch_size=1000):
    target_ids = [r.id for r in SeqIO.parse(fasta_path, "fasta") if r.id.startswith("CHEMBL")]
    molecules = []
    log(f"[INFO] Fetching molecules for {len(target_ids)} targets in batches of {batch_size}...")

    cursor = conn.cursor()
    for i in range(0, len(target_ids), batch_size):
        batch_targets = target_ids[i:i+batch_size]
        placeholders = ",".join("?" for _ in batch_targets)
        query = f"""
            SELECT td.chembl_id AS target_chembl_id,
                   md.molregno,
                   md.pref_name,
                   md.max_phase,
                   cs.canonical_smiles
            FROM activities a
            JOIN molecule_dictionary md ON a.molregno = md.molregno
            JOIN compound_structures cs ON md.molregno = cs.molregno
            JOIN assays ass ON a.assay_id = ass.assay_id
            JOIN target_dictionary td ON ass.tid = td.tid
            WHERE td.chembl_id IN ({placeholders})
              AND md.max_phase >= ?
        """
        cursor.execute(query, batch_targets + [min_phase])
        rows = cursor.fetchall()
        log(f"[INFO] Retrieved {len(rows)} molecules for targets {i+1}-{i+len(batch_targets)}")

        for row in rows:
            target_chembl_id, molregno, pref_name, max_phase, smiles = row
            mol_obj = Chem.MolFromSmiles(smiles)
            if mol_obj:
                fp = mol_to_fp(mol_obj)
                molecules.append({
                    "target_chembl_id": target_chembl_id,
                    "molregno": molregno,
                    "pref_name": pref_name,
                    "max_phase": max_phase,
                    "fp": fp
                })
    log(f"[INFO] Total molecules fetched: {len(molecules)}")
    return molecules

# ======================
# Compute dynamic threshold per TB drug
# ======================
def suggest_threshold(tb_fps, log):
    suggested = {}
    for tb_name, data in tb_fps.items():
        others = [d["fp"] for n,d in tb_fps.items() if n != tb_name]
        if others:
            sims = DataStructs.BulkTanimotoSimilarity(data["fp"], others)
            mean_sim = sum(sims) / len(sims)
            threshold = min(max(mean_sim + 0.474, 0.5), 0.95)  # ensure reasonable range
        else:
            threshold = 0.7
            mean_sim = 0.0
        suggested[tb_name] = threshold
        log(f"[INFO] {tb_name} - Mean similarity to other TB drugs: {mean_sim:.3f}, Suggested threshold: {threshold:.3f}")
    return suggested

# ======================
# Compute similarities
# ======================
def screen_molecules(tb_fps, molecules, thresholds, log):
    results = []
    for tb_name, data in tb_fps.items():
        tb_fp = data["fp"]
        tb_molregno = data["molregno"]
        fps_list = [m["fp"] for m in molecules]
        molregnos = [m["molregno"] for m in molecules]
        sims = DataStructs.BulkTanimotoSimilarity(tb_fp, fps_list)
        threshold = thresholds[tb_name]
        for m, sim in zip(molecules, sims):
            # skip TB drug itself
            if m["molregno"] == tb_molregno:
                continue
            if sim >= threshold:
                results.append({
                    "TB_Drug": tb_name,
                    "Target_CHEMBL_ID": m["target_chembl_id"],
                    "Drug_Name": m["pref_name"],
                    "Drug_molregno": m["molregno"],
                    "Max_Phase": m["max_phase"],
                    "Tanimoto": round(sim, 3)
                })
    log(f"[INFO] Total hits found above thresholds: {len(results)}")
    return results

# ======================
# Main
# ======================
def main():
    parser = argparse.ArgumentParser(description="Fast TB drug repurposing with ChEMBL SQLite")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--tb-drug", nargs="+", required=True)
    parser.add_argument("--min-phase", type=int, default=3)
    parser.add_argument("--sqlite-db", required=True)
    args = parser.parse_args()

    log, log_file = create_logger("chem.py")
    conn = sqlite3.connect(args.sqlite_db)

    # Fetch TB drug fingerprints
    tb_fps = get_tb_fps(args.tb_drug, conn, log)

    # Fetch molecules
    molecules = fetch_all_molecules(args.fasta, args.min_phase, conn, log)

    # Suggest thresholds per TB drug
    thresholds = suggest_threshold(tb_fps, log)

    # Screen molecules
    results = screen_molecules(tb_fps, molecules, thresholds, log)

    # Save results with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_csv = f"{args.out}_{timestamp}.csv"
    df = pd.DataFrame(results)
    df.to_csv(out_csv, index=False)
    log(f"[INFO] Results written to {out_csv}")
    conn.close()

if __name__ == "__main__":
    main()
