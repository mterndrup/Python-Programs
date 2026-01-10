import argparse
import sqlite3
import os
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from datetime import datetime
from rdkit import RDLogger
from multiprocessing import Pool, cpu_count
import numpy as np

RDLogger.DisableLog('rdApp.*')


# ======================
# Logging
# ======================
def create_logger(script_name):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"log_{timestamp}.txt"

    def log(msg):
        time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        entry = f"[{time_str}] [INFO] {msg}"
        print(entry)
        with open(log_file, "a") as f:
            f.write(entry + "\n")
    log(f"Script: {script_name}")
    log(f"Start time: {datetime.now()}")
    return log, log_file


# ======================
# Fingerprint helpers
# ======================
def mol_to_fp(mol, radius=2, nBits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)


def fp_to_numpy(fp):
    arr = np.zeros((1,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def compute_similarity_batch(args):
    batch, tb_fps, threshold, start_idx, log = args
    results = []
    for idx, (mol_id, fp) in enumerate(batch, start=start_idx):
        sims = DataStructs.BulkTanimotoSimilarity(fp, tb_fps)
        for i, sim in enumerate(sims):
            if sim >= threshold:
                results.append((mol_id, i, sim))
    log(f"Processed candidates {start_idx+1} to {start_idx+len(batch)} - Matches: {len(results)}")
    return results


# ======================
# Main
# ======================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sqlite-db", default="chembl_35.db", help="Path to ChEMBL SQLite DB")
    parser.add_argument("--tb-file", dest="tb_file", default=None, help="TB molecules input file")
    parser.add_argument("--organism", default="Mycobacterium tuberculosis")
    parser.add_argument("--threshold", type=float, default=0.8)
    parser.add_argument("--batch_size", type=int, default=5000)
    args = parser.parse_args()

    log, _ = create_logger("chem_organism_search.py")

    # ======================
    # Load TB molecules
    # ======================
    if args.tb_file and os.path.exists(args.tb_file):
        tb_smis = [line.strip() for line in open(args.tb_file)]
        tb_mols = [Chem.MolFromSmiles(smi) for smi in tb_smis]
        tb_mols = [m for m in tb_mols if m is not None]
        tb_fps = [mol_to_fp(m) for m in tb_mols]
        log(f"Loaded {len(tb_mols)} valid TB molecules from {args.tb_file}")
    else:
        tb_fps = []
        log("No TB molecule file provided or file not found. Skipping TB filtering.")

    # ======================
    # Connect to DB and fetch candidates
    # ======================
    conn = sqlite3.connect(args.sqlite_db)
    cur = conn.cursor()
    log(f"Connected to database: {args.sqlite_db}")

    cur.execute("SELECT molregno, canonical_smiles FROM compound_structures")
    candidates = cur.fetchall()
    conn.close()
    log(f"Fetched {len(candidates)} candidate molecules from DB")

    # Filter out TB molecules (if any)
    tb_smiles_set = set(tb_smis) if tb_fps else set()
    candidates = [(mid, smi) for mid, smi in candidates if smi not in tb_smiles_set]
    log(f"{len(candidates)} non-TB candidates remaining after filtering")

    # ======================
    # Precompute candidate fingerprints in parallel
    # ======================
    def candidate_fp_worker(args):
        mid, smi = args
        mol = Chem.MolFromSmiles(smi)
        if mol:
            return mid, mol_to_fp(mol)
        return None

    n_cpu = cpu_count()
    pool = Pool(n_cpu)
    candidate_fps = pool.map(candidate_fp_worker, candidates)
    pool.close()
    pool.join()
    candidate_fps = [c for c in candidate_fps if c is not None]
    log(f"Computed fingerprints for {len(candidate_fps)} candidates")

    # ======================
    # Parallel similarity search in batches
    # ======================
    batch_size = args.batch_size
    batches = [candidate_fps[i:i + batch_size] for i in range(0, len(candidate_fps), batch_size)]
    all_matches = []

    for i, batch in enumerate(batches):
        start_idx = i * batch_size
        args_list = [(batch, tb_fps, args.threshold, start_idx, log)]
        with Pool(n_cpu) as pool:
            results = pool.map(compute_similarity_batch, args_list)
        # Flatten and append to CSV incrementally
        flat_results = [r for batch_res in results for r in batch_res]
        all_matches.extend(flat_results)
        df = pd.DataFrame(flat_results, columns=["candidate_id", "tb_index", "similarity"])
        mode = "w" if i == 0 else "a"
        header = i == 0
        df.to_csv("similarity_results.csv", mode=mode, index=False, header=header)

    log(f"Similarity search finished. Total matches: {len(all_matches)}")


if __name__ == "__main__":
    main()
