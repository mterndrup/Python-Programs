import os
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def find_orfs(seq, min_pro_len=100):
    """Find ORFs in all 6 frames, return nucleotide and protein sequences as Seq objects."""
    if not isinstance(seq, Seq):
        seq = Seq(seq)  # convert to Seq if raw string passed
    orfs = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            # Translate from frame, don't stop at stop codons
            trans = nuc[frame:].translate(to_stop=False)
            trans_len = len(trans)
            aa_start = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len  # no stop codon found, go to end
                # Check if ORF length meets minimum protein length
                if aa_end - aa_start >= min_pro_len:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = frame + aa_end * 3
                        orf_nuc_seq = seq[start:end]
                    else:
                        # Reverse strand: coordinates reversed relative to original seq
                        start = seq_len - (frame + aa_end * 3)
                        end = seq_len - (frame + aa_start * 3)
                        orf_nuc_seq = seq[start:end].reverse_complement()
                    orf_prot_seq = trans[aa_start:aa_end]
                    orfs.append((orf_nuc_seq, orf_prot_seq))
                aa_start = aa_end + 1
    return orfs

def extract_orfs_from_assembly(assembly_file, min_pro_len=100):
    records = list(SeqIO.parse(assembly_file, "fasta"))
    all_orfs = []
    for record in records:
        orfs = find_orfs(record.seq, min_pro_len)
        for i, (nuc_seq, prot_seq) in enumerate(orfs):
            rec = SeqRecord(
                nuc_seq,
                id=f"{record.id}_orf{i}",
                description=""
            )
            rec.annotations['protein_seq'] = prot_seq  # store protein seq for scoring
            all_orfs.append(rec)
    return all_orfs

def score_orfs_against_target(orfs, target_seq):
    scores = []
    for rec in orfs:
        prot_seq = rec.annotations['protein_seq']
        # localxx is faster for alignment and good enough here
        alignments = pairwise2.align.localxx(target_seq, prot_seq, one_alignment_only=True)
        score = alignments[0].score if alignments else 0
        scores.append((rec.id, score))
    return sorted(scores, key=lambda x: x[1], reverse=True)

def main():
    base_path = r"C:\Users\ketgl\OneDrive\Documents\Glowing-Fungi"

    assembly_fasta = os.path.join(base_path, "assembly.fasta")
    target_protein_fasta = os.path.join(base_path, "Mycenaal.fasta")
    output_top_protein_fasta = os.path.join(base_path, "top_candidate_proteins.fasta")
    output_top_nuc_fasta = os.path.join(base_path, "top_candidate_nucleotides.fasta")
    min_orf_length = 100  # Minimum amino acid length

    print("Finding ORFs...")
    orfs = extract_orfs_from_assembly(assembly_fasta, min_orf_length)
    print(f"Total ORFs found: {len(orfs)}")

    print("Loading target protein...")
    target_record = SeqIO.read(target_protein_fasta, "fasta")
    target_seq = target_record.seq

    print("Scoring alignments with localxx...")
    scored_orfs = score_orfs_against_target(orfs, target_seq)

    print("\nTop matches:")
    for i, (orf_id, score) in enumerate(scored_orfs[:10], 1):
        print(f"{i:2d}. {orf_id:30} Score: {score}")

    top_ids = {orf_id for orf_id, _ in scored_orfs[:10]}
    top_orfs_nuc = [rec for rec in orfs if rec.id in top_ids]

    # Write nucleotide sequences of top hits
    SeqIO.write(top_orfs_nuc, output_top_nuc_fasta, "fasta")

    # Write protein sequences of top hits
    top_orfs_prot = [
        SeqRecord(
            rec.annotations['protein_seq'],
            id=rec.id,
            description=""
        )
        for rec in top_orfs_nuc
    ]
    SeqIO.write(top_orfs_prot, output_top_protein_fasta, "fasta")

    print(f"\nTop 10 matching nucleotide sequences saved to: {output_top_nuc_fasta}")
    print(f"Top 10 matching protein sequences saved to: {output_top_protein_fasta}")

if __name__ == "__main__":
    main()
