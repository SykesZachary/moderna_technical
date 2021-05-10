# Analysis of SARS CoV 2 Spike protein in an attempt to replicate mrna-1273 sequence
# 08.05.2021 -- Zachary Sykes

import os

from file_handling import data_extraction
from file_handling import file_output

# PATH set
proj_path = os.path.abspath(os.getcwd())
sars_cov_2_aa = os.path.join(proj_path, 'SARS_CoV_2_aaSeq/')  # SARS CoV 2 S2P
sars_cov_2_s_gene = os.path.join(proj_path, 'SARS_CoV_2_S_Gene/')  # S gene sequences
possible_mrna_out = os.path.join(proj_path, 'SARS_CoV_2_bt_output/')  # Back translated mRNA

# Human codon table for amino acids -- Mostly for reference
hsap_codon_table = {'F': ['UUU', 'UUC'],
                    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
                    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                    'Y': ['UAU', 'UAC'], 'Stop': ['UAA', 'UAG', 'UGA'],
                    'C': ['UGU', 'UGC'], 'W': ['UGG'],
                    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                    'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'],
                    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                    'I': ['AUU', 'AUC', 'AUA'], 'M': ['AUG'],
                    'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'N': ['AAU', 'AAC'],
                    'K': ['AAA', 'AAG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                    'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'D': ['GAU', 'GAC'],
                    'E': ['GAA', 'GAG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG']}

# Storing aa sequence data for SARS CoV 2 Spike protein
aa_seq = data_extraction(sars_cov_2_aa, 'sars_cov_2_aa.fasta')
s_dna_seq = data_extraction(sars_cov_2_s_gene, 'sars_cov_2_S_Gene_dna.fasta')

# Substituting the 2 proline at residues 986 and 987 to spike protein sequence
aa_s2p_seq = aa_seq[0][:985] + 'PP' + aa_seq[0][987:]
file_output(
    sars_cov_2_aa,
    'sars_cov_2_s2p_aa.fasta',
    f'{aa_seq[1]}_S2P_986|987',
    aa_s2p_seq
)  # Copy of S2P variant saved

# Generating wild type spike protein mRNA sequence
wt_s_mrna = s_dna_seq.transcribe()

# Generating potential mRNA-1273 sequences with all possible proline codon combinations at residue 986/3 and 987/3
i = 1  # Counter for file naming
for codon_1 in hsap_codon_table['P']:
    for codon_2 in hsap_codon_table['P']:
        mrna_1273_candidate = wt_s_mrna[:985 * 3] + codon_1 + codon_2 + wt_s_mrna[987 * 3:]
        file_output(
            possible_mrna_out,
            f'mRNA1273_Candidate_{i}.fasta',
            f'{s_dna_seq[1]} | mRNA-1273 Candidate | Proline codon substitution start at {985 * 3} and {986 * 3}',
            mrna_1273_candidate,
        )
        i += 1
