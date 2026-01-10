import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from datetime import datetime
import os
import re


# ======================
# Logging helper
# ======================
def create_logger(script_name="process_repurposing.py"):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"log_{script_name}_{timestamp}.txt"

    def log(msg):
        time_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        entry = f"[{time_str}] {msg}"
        print(entry)
        with open(log_file, "a") as f:
            f.write(entry + "\n")

    log(f"Script started: {script_name}")
    return log, log_file


# ======================
# Sanitize filenames
# ======================
def sanitize_filename(name):
    return re.sub(r'[<>:"/\\|?*]', '_', name)


# ======================
# Load similarity results
# ======================
def load_similarity_results(log, file_path="similarity_results_Pseudomonas_aeruginosa.csv"):
    df = pd.read_csv(file_path)
    log(f"Loaded {len(df)} similarity results from {file_path}.")
    return df


# ======================
# Deduplicate rows
# ======================
def deduplicate_rows(df, log):
    before = len(df)
    df_clean = df.drop_duplicates(subset=['ref_smiles', 'candidate_smiles', 'candidate_target'])
    removed = before - len(df_clean)
    log(f"Deduplicated rows: removed {removed} duplicates.")
    return df_clean


# ======================
# Merge metadata (optional)
# ======================
def merge_metadata(df, log, metadata_file="candidate_metadata.csv"):
    try:
        df_meta = pd.read_csv(metadata_file)
        df_merged = df.merge(df_meta, on="candidate_chembl_id", how="left")
        log(f"Metadata merged from {metadata_file}. Total rows after merge: {df_merged.shape[0]}")
        return df_merged
    except FileNotFoundError:
        log(f"Metadata file '{metadata_file}' not found. Skipping metadata merge.")
        return df


# ======================
# Summarize top hits per target
# ======================
def summarize_results(df, log, top_n_mols=5, top_n_targets=25):
    # Select top N targets by maximum similarity
    top_targets = df.groupby('candidate_target')['similarity'].max().sort_values(ascending=False).head(
        top_n_targets).index
    df_top_targets = df[df['candidate_target'].isin(top_targets)]

    summary = {}
    for target, group in df_top_targets.groupby('candidate_target'):
        top_hits = group.sort_values('similarity', ascending=False).head(top_n_mols)
        summary[target] = top_hits
        log(f"Top {top_n_mols} hits selected for target: {target}")
    return summary


# ======================
# Visualize top hits with labels and summary
# ======================
def visualize_top_hits(summary, log, output_dir="visualizations"):
    os.makedirs(output_dir, exist_ok=True)
    summary_path = os.path.join(output_dir, "visualization_summary_H.txt")

    with open(summary_path, "w") as f_summary:
        f_summary.write("Visualization Summary\n")
        f_summary.write("====================\n\n")

        for target, group in summary.items():
            mols = [Chem.MolFromSmiles(s) for s in group['candidate_smiles'] if Chem.MolFromSmiles(s)]
            if mols:
                # Create labels for each molecule
                legends = [
                    f"{row['candidate_chembl_id']}\nSim: {row['similarity']:.2f}\n{row['standard_type']}: {row['standard_value']}"
                    for idx, row in group.iterrows()
                ]
                img = Draw.MolsToGridImage(mols, molsPerRow=len(mols), legends=legends)
                safe_name = sanitize_filename(target)
                img_path = os.path.join(output_dir, f"{safe_name}_top_hits.png")
                img.save(img_path)
                log(f"Visualization for {target} saved to {img_path}")

                # Write summary info
                f_summary.write(f"Target: {target}\n")
                f_summary.write(f"Reason: Top similarity hits (included in top N targets)\n")
                f_summary.write(f"Image file: {img_path}\n")
                f_summary.write(f"Molecules included:\n")
                for idx, row in group.iterrows():
                    f_summary.write(
                        f"  - ChEMBL ID: {row['candidate_chembl_id']}, "
                        f"Similarity: {row['similarity']:.2f}, "
                        f"{row['standard_type']}: {row['standard_value']}, "
                        f"Max Phase: {row.get('candidate_max_phase', '')}\n"
                    )
                f_summary.write("\n")

    log(f"Visualization summary saved to {summary_path}")


# ======================
# Save final report
# ======================
def save_final_report(df, log, file_name="repurposing_final_report_PS.csv"):
    df.to_csv(file_name, index=False)
    log(f"Final repurposing report saved to {file_name}")


# ======================
# Main program
# ======================
def main():
    log, log_file = create_logger("process_repurposing.py")

    # Load similarity results
    df = load_similarity_results(log, "similarity_results_Pseudomonas_aeruginosa.csv")

    # Deduplicate rows
    df = deduplicate_rows(df, log)

    # Merge metadata (optional)
    df = merge_metadata(df, log)

    # Summarize top hits per target (limits top N targets and top M molecules)
    summary = summarize_results(df, log, top_n_mols=5, top_n_targets=25)
    for target, group in summary.items():
        print(f"\nTop hits for {target}:")
        print(group[['ref_smiles', 'candidate_smiles', 'similarity', 'candidate_max_phase']])

    # Visualize top hits with informative labels and summary
    visualize_top_hits(summary, log)

    # Save final CSV
    save_final_report(df, log)

    log("[INFO] Script completed successfully.")
    print(f"Logs saved to {log_file}")


if __name__ == "__main__":
    main()
