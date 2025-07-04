import subprocess
import os
from collections import Counter

# Paths to executables - update if needed
VSEARCH_EXE = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tools\VSEARCH\vsearch-2.30.0-win-x86_64\bin\vsearch.exe"
MAFFT_EXE = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Tools\MAFFT\mafft-7.526-win64-signed (1)\mafft-win\mafft.bat"

# Input ITS sequences fasta file
input_fasta = r"C:\Users\ketgl\OneDrive\Documents\GitHub\Python-Programs\BioTech\Bioinformatics\Plants\Aya\DNA_Sequences\barcode32\ITSx\ITSx_out_PITS.full.fasta"

# Output files
base_dir = os.path.dirname(input_fasta)
output_rep_seqs = os.path.join(base_dir, "rep_seqs.fasta")
output_clusters = os.path.join(base_dir, "clusters.uc")
output_aligned = os.path.join(base_dir, "rep_seqs_aligned.fasta")
output_consensus = os.path.join(base_dir, "consensus_sequence.fasta")

identity_threshold = 0.97  # clustering identity threshold

def run_vsearch_clustering(vsearch_path, input_file, output_fasta, output_uc, identity):
    command = [
        vsearch_path,
        "--cluster_fast", input_file,
        "--id", str(identity),
        "--centroids", output_fasta,
        "--uc", output_uc,
        "--strand", "plus",
        "--threads", "4"
    ]
    print(f"Running vsearch clustering:\n{' '.join(command)}")
    try:
        subprocess.run(command, check=True)
        print(f"✅ Clustering complete.")
        print(f"Representative sequences saved to: {output_fasta}")
        print(f"Cluster info saved to: {output_uc}")
    except subprocess.CalledProcessError as e:
        print(f"Error running vsearch: {e}")
        exit(1)

def run_mafft_alignment(mafft_path, input_fasta, output_fasta):
    command = [mafft_path, "--auto", input_fasta]
    print(f"Running MAFFT alignment:\n{' '.join(command)}")
    try:
        with open(output_fasta, "w") as outfile:
            subprocess.run(command, stdout=outfile, check=True)
        print(f"✅ MAFFT alignment complete.")
        print(f"Aligned sequences saved to: {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT: {e}")
        exit(1)

def read_fasta(fasta_path):
    sequences = {}
    with open(fasta_path) as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq_lines)
                seq_id = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line.lower())
        if seq_id:
            sequences[seq_id] = "".join(seq_lines)
    return sequences

def write_fasta(seq_id, sequence, output_path):
    with open(output_path, "w") as f:
        f.write(f">{seq_id}\n")
        # wrap lines at 60 chars per fasta convention
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")

def consensus_from_alignment(sequences_dict):
    # sequences_dict: {id: aligned_seq}
    # all sequences should be same length
    seqs = list(sequences_dict.values())
    length = len(seqs[0])
    consensus_chars = []

    for i in range(length):
        column = [seq[i] for seq in seqs if seq[i] != '-']  # ignore gaps
        if not column:
            consensus_chars.append('N')  # no info at this position
            continue
        count = Counter(column)
        most_common_base, count_most = count.most_common(1)[0]
        consensus_chars.append(most_common_base.upper())

    return "".join(consensus_chars)

if __name__ == "__main__":
    # Check executables and input file
    if not os.path.isfile(VSEARCH_EXE):
        print(f"❌ VSEARCH executable not found:\n{VSEARCH_EXE}")
        exit(1)
    if not os.path.isfile(MAFFT_EXE):
        print(f"❌ MAFFT executable not found:\n{MAFFT_EXE}")
        exit(1)
    if not os.path.isfile(input_fasta):
        print(f"❌ Input FASTA file not found:\n{input_fasta}")
        exit(1)

    # Step 1: Cluster sequences with vsearch
    run_vsearch_clustering(VSEARCH_EXE, input_fasta, output_rep_seqs, output_clusters, identity_threshold)

    # Step 2: Align representative sequences with mafft
    run_mafft_alignment(MAFFT_EXE, output_rep_seqs, output_aligned)

    # Step 3: Generate consensus sequence from alignment
    aligned_seqs = read_fasta(output_aligned)
    consensus_seq = consensus_from_alignment(aligned_seqs)
    write_fasta("consensus_sequence", consensus_seq, output_consensus)
    print(f"✅ Consensus sequence generated and saved to: {output_consensus}")
