import argparse
from Bio import SeqIO
from chembl_webresource_client.new_client import new_client
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
import os

# ======================
# Helper Functions
# ======================
def get_mol_smiles(tb_drugs):
    """Fetch canonical SMILES for known TB drugs by name"""
    molecule = new_client.molecule
    tb_fps = {}
    for drug in tb_drugs:
        try:
            res = molecule.filter(pref_name__iexact=drug).only(
                ["molecule_chembl_id", "pref_name", "molecule_structures"]
            )
            res = list(res)
            if res and res[0]["molecule_structures"]:
                smiles = res[0]["molecule_structures"]["canonical_smiles"]
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                    tb_fps[drug] = fp
        except Exception as e:
            print(f"[WARN] Could not fetch SMILES for {drug}: {e}")
    return tb_fps


def screen_targets(fasta_path, tb_fps, min_phase, threshold, out_csv):
    """Screen FASTA protein targets against ChEMBL for repurposing hits"""
    activity = new_client.activity
    molecule = new_client.molecule
    results = []

    # Check if output CSV exists and load existing results
    if os.path.exists(out_csv):
        try:
            df_existing = pd.read_csv(out_csv)
            processed_targets = set(df_existing["Target_CHEMBL_ID"].unique())
            results = df_existing.to_dict(orient="records")
        except Exception:
            processed_targets = set()
    else:
        processed_targets = set()

    # Extract ChEMBL target IDs from FASTA headers
    fasta_records = list(SeqIO.parse(fasta_path, "fasta"))
    for idx, record in enumerate(fasta_records, start=1):
        target_id = record.id
        if target_id in processed_targets:
            print(f"[INFO] Skipping already processed target {idx}: {target_id}")
            continue

        print(f"[INFO] Processing target {idx}: {target_id}")

        # Query bioactivities linked to this target
        try:
            acts = activity.filter(target_chembl_id=target_id)
            for act in acts:
                mol_id = act.get("molecule_chembl_id")
                if not mol_id:
                    continue

                mol_data = molecule.get(mol_id)
                if not mol_data:
                    continue

                # Only keep FDA/EMA approved or Phase III+
                phase = mol_data.get("max_phase")
                try:
                    phase = int(phase)
                except (TypeError, ValueError):
                    continue  # skip if max_phase is invalid
                if phase < min_phase:
                    continue

                if not mol_data.get("molecule_structures"):
                    continue

                smiles = mol_data["molecule_structures"]["canonical_smiles"]
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    continue

                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

                # Compare with TB reference drugs
                for tb_name, tb_fp in tb_fps.items():
                    sim = DataStructs.TanimotoSimilarity(fp, tb_fp)
                    if sim >= threshold:
                        results.append({
                            "TB_Drug": tb_name,
                            "Target_CHEMBL_ID": target_id,
                            "Drug_Name": mol_data.get("pref_name"),
                            "Drug_CHEMBL_ID": mol_id,
                            "Max_Phase": phase,
                            "Tanimoto": round(sim, 3),
                            "ATC_Codes": mol_data.get("atc_classifications"),
                            "Indications": mol_data.get("therapeutic_flag")
                        })

            # Save CSV incrementally after each target
            df = pd.DataFrame(results)
            df.to_csv(out_csv, index=False)

        except Exception as e:
            print(f"[WARN] Skipping target {target_id}: {e}")
            continue

    return results


# ======================
# Main
# ======================
def main():
    parser = argparse.ArgumentParser(description="TB drug repurposing pipeline using ChEMBL.")
    parser.add_argument("--fasta", required=True, help="Path to ChEMBL FASTA file (targets).")
    parser.add_argument("--out", required=True, help="Output CSV file.")
    parser.add_argument("--threshold", type=float, default=0.7, help="Tanimoto similarity cutoff.")
    parser.add_argument("--tb-drug", nargs="+", required=True, help="List of TB reference drugs.")
    parser.add_argument("--min-phase", type=int, default=3,
                        help="Minimum clinical phase (3=late stage, 4=approved).")
    args = parser.parse_args()

    print("[INFO] Fetching SMILES for TB reference drugs...")
    tb_fps = get_mol_smiles(args.tb_drug)

    print(f"[INFO] Screening targets from {args.fasta} ...")
    results = screen_targets(args.fasta, tb_fps, args.min_phase, args.threshold, args.out)

    print(f"[INFO] Found {len(results)} repurposing candidates.")
    print(f"[INFO] Results written to {args.out}")


if __name__ == "__main__":
    main()
