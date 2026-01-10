import argparse
import sqlite3
import os
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from datetime import datetime
from rdkit import RDLogger
from rdkit.DataStructs import ExplicitBitVect

#Program created by Matthew Terndrup (Aug 2025)
#Pasadena City College

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


# ======================
# Fetch targets by organism
# ======================
def get_targets_by_organism(conn, organism, log, debug=False):
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='target_dictionary'")
    if not cursor.fetchone():
        log("[ERROR] Table 'target_dictionary' does not exist! Cannot fetch targets.")
        return []

    query = """
        SELECT td.chembl_id, td.pref_name, cs.accession AS uniprot_accession
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE td.target_type='SINGLE PROTEIN' AND td.organism=?
    """
    cursor.execute(query, (organism,))
    targets = cursor.fetchall()
    log(f"[INFO] Found {len(targets)} targets for organism '{organism}'")
    if debug:
        log(f"[DEBUG] Targets fetched: {targets}")
    return targets

# ======================
# Fetch molecules for targets
# ======================
def fetch_molecules_for_targets(conn, targets, log, debug=False, batch_size=1000, out_txt="tb_molecules.txt"):
    molecules = []
    if not targets:
        log("[WARN] No targets provided. Skipping molecule fetch.")
        return molecules

    cursor = conn.cursor()
    log(f"[INFO] Fetching molecules for {len(targets)} targets in batches of {batch_size}...")

    # Check required tables exist
    for table in ["activities", "assays", "target_dictionary", "compound_structures"]:
        cursor.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table}'")
        if not cursor.fetchone():
            log(f"[ERROR] Table '{table}' does not exist! Cannot fetch molecules.")
            return molecules

    # Open output file once
    with open(out_txt, "w") as f_out:

        for i in range(0, len(targets), batch_size):
            batch_targets = targets[i:i + batch_size]
            placeholders = ",".join("?" for _ in batch_targets)
            query = f"""
                SELECT td.chembl_id AS target_chembl_id,
                       td.pref_name,
                       csq.accession AS uniprot_accession,
                       a.molregno,
                       cs.canonical_smiles,
                       a.standard_type,
                       a.standard_value,
                       a.standard_units
                FROM activities a
                JOIN assays ass ON a.assay_id = ass.assay_id
                JOIN target_dictionary td ON ass.tid = td.tid
                JOIN target_components tc ON td.tid = tc.tid
                JOIN component_sequences csq ON tc.component_id = csq.component_id
                JOIN compound_structures cs ON a.molregno = cs.molregno
                WHERE td.chembl_id IN ({placeholders})
                  AND a.standard_type IN ('Inhibition', 'IC50', 'Ki', 'AC50')
            """
            target_ids = [t[0] for t in batch_targets]

            # Show SQL query and wait for confirmation
            print("\n====================================")
            print("[SQL QUERY TO BE EXECUTED]:")
            print(query)
            print("With target_ids:", target_ids)
            print("====================================")
            proceed = input("Do you want to execute this query? (y/n): ").strip().lower()
            if proceed != "y":
                log("[INFO] User chose not to execute query. Skipping this batch.")
                continue

            cursor.execute(query, target_ids)
            rows = cursor.fetchall()
            log(f"[INFO] Retrieved {len(rows)} molecules for targets {i+1}-{i+len(batch_targets)}")
            if debug:
                log(f"[DEBUG] Molecules fetched: {rows[:5]} ...")

            # Process molecules and write to file immediately
            for row in rows:
                target_chembl_id, pref_name, uniprot_accession, molregno, smiles, standard_type, standard_value, standard_units = row
                mol_obj = Chem.MolFromSmiles(smiles)
                if mol_obj:
                    fp = mol_to_fp(mol_obj)
                    bitstring = fp.ToBitString()
                    molecules.append({
                        "target_chembl_id": target_chembl_id,
                        "pref_name": pref_name,
                        "uniprot_accession": uniprot_accession,
                        "molregno": molregno,
                        "smiles": smiles,
                        "fp_bitstring": bitstring,
                        "standard_type": standard_type,
                        "standard_value": standard_value,
                        "standard_units": standard_units
                    })
                    # Write SMILES + optional info to txt file
                    f_out.write(f"{smiles}\t{pref_name}\t{target_chembl_id}\n")

    log(f"[INFO] Total molecules fetched: {len(molecules)}")
    log(f"[INFO] All molecules saved to {out_txt}")
    return molecules


def mol_to_fp(mol):
    """Generate Morgan fingerprint for a molecule."""
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def bitstring_to_fp(bitstring: str):
    """Convert a 0/1 bitstring back into RDKit ExplicitBitVect."""
    fp = ExplicitBitVect(len(bitstring))
    for i, bit in enumerate(bitstring):
        if bit == "1":
            fp.SetBit(i)
    return fp


def similarity_search(tb_drugs, candidates, threshold=0.7, batch_size=5000, log=None, out_file="similarity_results.csv"):
    """
    Compare TB reference drugs to candidate molecules using Tanimoto similarity in batches.
    Writes results incrementally to a CSV file.

    Args:
        tb_drugs: list of dicts, each containing 'smiles' and 'name'
        candidates: pandas DataFrame with 'fp_bitstring', 'smiles', 'pref_name', 'target_chembl_id', etc.
        threshold: similarity threshold (0.0-1.0)
        batch_size: number of candidate molecules to process at a time
        log: optional logging function
        out_file: CSV file to write incremental results

    Returns:
        pandas DataFrame of all matches
    """
    import os

    total_candidates = len(candidates)
    log(f"[INFO] Total candidate molecules: {total_candidates}") if log else None

    # Precompute fingerprints for TB reference drugs
    tb_with_fp = []
    for tb in tb_drugs:
        mol = Chem.MolFromSmiles(tb["smiles"])
        if mol:
            tb_with_fp.append((tb, mol_to_fp(mol)))

    # Initialize output file
    if os.path.exists(out_file):
        os.remove(out_file)
    first_batch = True

    all_results = []

    for start in range(0, total_candidates, batch_size):
        batch = candidates.iloc[start:start + batch_size]
        log(f"[INFO] Processing candidates {start + 1} to {start + len(batch)}") if log else None

        batch_results = []
        for _, cand in batch.iterrows():
            cand_fp = bitstring_to_fp(cand["fp_bitstring"])
            for tb, tb_fp in tb_with_fp:
                sim = DataStructs.TanimotoSimilarity(tb_fp, cand_fp)
                if sim >= threshold:
                    batch_results.append({
                        "ref_chembl_id": tb.get("target_chembl_id", ""),
                        "ref_smiles": tb.get("smiles", ""),
                        "candidate_target": cand.get("pref_name", ""),
                        "candidate_chembl_id": cand.get("target_chembl_id", ""),
                        "candidate_smiles": cand.get("smiles", ""),
                        "candidate_max_phase": cand.get("max_phase", ""),
                        # make sure candidates DataFrame includes max_phase
                        "similarity": sim,
                        "standard_type": cand.get("standard_type", ""),
                        "standard_value": cand.get("standard_value", ""),
                        "standard_units": cand.get("standard_units", "")
                    })

        if batch_results:
            df_batch = pd.DataFrame(batch_results)
            df_batch.to_csv(out_file, mode="a", header=first_batch, index=False)
            first_batch = False
            log(f"[INFO] Batch {start + 1}-{start + len(batch)}: {len(batch_results)} matches written to {out_file}") if log else None
            all_results.extend(batch_results)

    log(f"[INFO] Similarity search finished. Total matches: {len(all_results)}") if log else None
    return pd.DataFrame(all_results)




# ======================
# Main
# ======================
def main():
    parser = argparse.ArgumentParser(description="Drug search by organism (safe DB handling)")
    parser.add_argument("--out", required=False, help="Optional output CSV file for molecules")
    parser.add_argument("--in", dest="infile", required=False, help="Optional input CSV file with TB molecules")
    parser.add_argument("--targets-out", required=False, help="Output CSV file for proteins/targets")
    parser.add_argument("--sqlite-db", required=True)
    parser.add_argument("--organism", default="Mycobacterium tuberculosis")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode to stop after each step")
    args = parser.parse_args()

    log, log_file = create_logger("chem_organism_search.py")

    # ===== Check DB exists =====
    if not os.path.isfile(args.sqlite_db):
        log(f"[ERROR] Database file not found: {args.sqlite_db}")
        return
    conn = sqlite3.connect(args.sqlite_db)
    log(f"[INFO] Connected to database: {args.sqlite_db}")

    # ===== Load TB molecules from file or fetch from DB =====
    if args.infile and os.path.isfile(args.infile):
        log(f"[INFO] Loading TB molecules from file: {args.infile}")
        molecules = pd.read_csv(args.infile).to_dict(orient="records")
    else:
        log("[INFO] Fetching targets and molecules from database...")
        # Fetch targets
        targets = get_targets_by_organism(conn, args.organism, log, debug=args.debug)
        if not targets:
            log(f"[WARN] No targets found for organism '{args.organism}'")
            conn.close()
            return

        # Fetch molecules
        molecules = fetch_molecules_for_targets(conn, targets, log, debug=args.debug)
        if args.out:
            df = pd.DataFrame(molecules)
            df.to_csv(args.out, index=False)
            log(f"[INFO] Molecules output saved to {args.out}")

    if not molecules:
        log("[WARN] No TB molecules available for similarity search.")
        conn.close()
        return

    # ===== Fetch candidates for similarity search =====
    log("[INFO] Fetching non-TB candidate drugs...")
    cursor = conn.cursor()
    query = """
           SELECT m.molregno,
                  m.chembl_id,
                  cs.canonical_smiles,
                  td.pref_name,
                  td.organism,
                  a.standard_type,
                  a.standard_value,
                  a.standard_units
           FROM molecule_dictionary m
           JOIN compound_structures cs ON m.molregno = cs.molregno
           JOIN activities a ON m.molregno = a.molregno
           JOIN assays ass ON a.assay_id = ass.assay_id
           JOIN target_dictionary td ON ass.tid = td.tid
           WHERE m.max_phase >= 3
             AND td.organism != ?
             AND cs.canonical_smiles IS NOT NULL
       """
    cursor.execute(query, (args.organism,))
    rows = cursor.fetchall()
    log(f"[INFO] Retrieved {len(rows)} non-TB candidate drugs")

    candidates = []
    for molregno, chembl_id, smiles, target_name, organism, std_type, std_val, std_units in rows:
        mol_obj = Chem.MolFromSmiles(smiles)
        if mol_obj:
            fp = mol_to_fp(mol_obj)
            candidates.append({
                "molregno": molregno,
                "target_chembl_id": chembl_id,
                "smiles": smiles,
                "fp_bitstring": fp.ToBitString(),
                "pref_name": target_name,
                "organism": organism,
                "standard_type": std_type,
                "standard_value": std_val,
                "standard_units": std_units
            })
    df_candidates = pd.DataFrame(candidates)

    # ===== Run similarity search =====
    log("[INFO] Starting similarity search...")
    sim_results = similarity_search(molecules, df_candidates, threshold=0.7, log=log)
    sim_results.to_csv("similarity_results.csv", index=False)
    log(f"[INFO] Similarity results saved to similarity_results.csv")

    conn.close()
    log("[INFO] Database connection closed.")

if __name__ == "__main__":
    main()
