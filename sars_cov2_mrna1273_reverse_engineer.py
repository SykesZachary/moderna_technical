# Analysis of SARS CoV 2 Spike protein in an attempt to replicate mrna-1273 sequence
# 08.05.2021 - Zachary Sykes

import os
from Bio.Seq import Seq
from Bio.SeqUtils import GC123

from file_handling import data_extraction
from file_handling import file_output

# PATHs set
proj_path = os.path.abspath(os.getcwd())
sars_cov_2_aa = os.path.join(proj_path, 'SARS_CoV_2_aaSeq/')  # SARS CoV 2 S protein sequences
mrna_1273_sequence = os.path.join(proj_path, 'mrna1273_fastaSeq/')  # mRNA-1273 fasta file

# Human codon table
hsap_codon_table = {'F': ['TTT', 'TTC'],
                    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                    'Y': ['TAT', 'TAC'], '*': ['TAA', 'TAG', 'TGA'],
                    'C': ['TGT', 'TGC'], 'W': ['TGG'],
                    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                    'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
                    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'],
                    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
                    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
                    'G': ['GGT', 'GGC', 'GGA', 'GGG']}


# Builds a dynamic codon table based on codons present in NT sequence
# Takes a NT seq str and an AA seq str (or Seq obj from Biopython)
def codon_counter(nt, codons, nt_type='dna'):

    # Stores codons used for each amino acid and frequency used for said amino acid
    codon_table = dict()

    # Grabs the key (aa) for the given value (codon)
    def get_key(val):
        for key, value in codons.items():
            if val in value:
                return key

    # Handles a RNA string passed to the codon counter
    if nt_type == 'rna' and type(nt) is not Seq:
        nt = Seq(nt)
        nt = nt.back_transcribe()
    elif nt_type == 'rna' and type(nt) is Seq:
        nt = nt.back_transcribe()

    start = None
    stop = None

    # Start and stop codons identified for the sequence
    for frame in range(0, len(nt), 3):
        if nt[frame: frame + 3] == 'ATG' and not start:
            print(f'Start codon {nt[frame: frame + 3]} identified at position {frame}')
            start = frame
        # mRNA-1273 contains all three stop codons at the end of the sequence
        # TAG was the last one before the 3' UTR so all stop codons included in the codon table
        if nt[frame: frame + 3] == 'TAG' and not stop:
            print(f'Stop codon {nt[frame: frame + 3]} identified at position {frame}')
            stop = frame + 3

    # Trimmed nt sequence starting at ATG and ending at TAG
    nt_cds = nt[start:stop]
    prev_codon = ''
    # Counting codons used per amino acid
    for frame in range(3, (len(nt_cds) + 3), 3):
        aa = get_key(nt_cds[frame - 3: frame])
        codon_table.setdefault(aa, []).append(str(nt_cds[frame - 3: frame]))

    # Returns a list of tuples (codon, num times used to translate aa in nt seq provided / total codons for aa)
    for aa in codon_table.keys():
        codon_counts = {aa: [(codon, round(codon_table[aa].count(codon) / len(codon_table[aa]), 3))
                             for codon in set(codon_table[aa])]}
        codon_table.update(codon_counts)
    print(GC123(nt_cds))

    return codon_table


# Storing AA and NT sequence data for SARS CoV 2 Spike protein
aa_seq = data_extraction(sars_cov_2_aa, 'sars_cov_2_aa.fasta')
mrna1273_dna_seq = data_extraction(mrna_1273_sequence, 'mrna1273_spike_encoding.fasta')

# Substituting the 2 proline at residues 986 and 987 to spike protein sequence
aa_s2p_seq = aa_seq[0][:985] + 'PP' + aa_seq[0][987:]
file_output(
    sars_cov_2_aa,
    'sars_cov_2_s2p_aa.fasta',
    f'{aa_seq[1]}_S2P_986|987',
    aa_s2p_seq
)  # Copy of S2P variant saved to disk


# Generating a codon table with frequencies of use from an NT sequence
# Codon table passed is the human codon/amino acid associations table
# NT sequence is mRNA-1273 encoding the S2P perfusion stabilized variant of the SARS-CoV-2 spike protein
for k, v in codon_counter(mrna1273_dna_seq[0], hsap_codon_table).items():
    print(k, v)

# NOTE - There is a preferred codon for each amino acid in the encoded in mRNA-1273
# NOTE - In very low frequency a synonymous codon is selected.
# NOTE - Multiple readings on degenerate codons usage across genomes leads me to believe that the reason
# NOTE - for the low frequency synonymous codon usage is due to regional GC content
