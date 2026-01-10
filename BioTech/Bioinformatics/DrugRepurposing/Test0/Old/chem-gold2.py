import argparse
import sqlite3
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
from datetime import datetime
from rdkit import RDLogger

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
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

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
def fetch_molecules_for_targets(conn, targets, log, debug=False, batch_size=1000):
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

        # ===== NEW: Show SQL query and wait for confirmation =====
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
    log(f"[INFO] Total molecules fetched: {len(molecules)}")
    return molecules

# ======================
# Main
# ======================
def main():
    parser = argparse.ArgumentParser(description="Drug search by organism (safe DB handling)")
    parser.add_argument("--out", required=False, help="Optional output CSV file for molecules")
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

    # ===== Fetch targets first =====
    log("[INFO] Fetching targets for organism...")
    targets = get_targets_by_organism(conn, args.organism, log, debug=args.debug)

    if targets:
        log(f"[INFO] Proteins for organism '{args.organism}':")
        for t in targets:
            chembl_id, pref_name, uniprot_accession = t
            log(f" - {chembl_id}: {pref_name} (UniProt: {uniprot_accession})")

        # ===== Save targets to CSV if requested =====
        if args.targets_out:
            df_targets = pd.DataFrame(targets, columns=["chembl_id", "pref_name", "uniprot_accession"])
            df_targets.to_csv(args.targets_out, index=False)
            log(f"[INFO] Proteins/targets saved to {args.targets_out}")
    else:
        log(f"[WARN] No targets found for organism '{args.organism}'")
        conn.close()
        return

    if args.debug:
        log("[DEBUG] Stopping after printing targets in debug mode")
        conn.close()
        return

    # ===== Fetch molecules for targets =====
    log("[INFO] Fetching molecules for targets...")
    molecules = fetch_molecules_for_targets(conn, targets, log, debug=args.debug)
    if args.debug:
        log("[DEBUG] Stopping after fetching molecules in debug mode")
        conn.close()
        return

    if molecules:
        log(f"[INFO] Total molecules fetched: {len(molecules)}")
        # Print first 10 molecules for immediate inspection
        log("[INFO] Sample of fetched molecules:")
        for m in molecules[:10]:
            log(m)

        # Optional: write CSV if --out provided
        if args.out:
            df = pd.DataFrame(molecules)
            df.to_csv(args.out, index=False)
            log(f"[INFO] Molecules output saved to {args.out}")
    else:
        log("[WARN] No molecules found. Output file not created.")

    conn.close()
    log("[INFO] Database connection closed.")

if __name__ == "__main__":
    main()
