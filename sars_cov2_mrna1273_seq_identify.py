# Analysis of SARS CoV 2 Spike protein in an attempt to replicate mrna-1273 sequence
# 08.05.2021 - Zachary Sykes

import os
from Bio.Seq import Seq

from file_handling import data_extraction
from file_handling import file_output
from seq_manipulation import SeqManipulation

# PATH set
proj_path = os.path.abspath(os.getcwd())
sars_cov_2_aa = os.path.join(proj_path, 'SARS_CoV_2_aaSeq/')  # SARS CoV 2 S2P
sars_cov_2_s_gene = os.path.join(proj_path, 'SARS_CoV_2_S_Gene/')  # S gene sequences
possible_mrna_out = os.path.join(proj_path, 'SARS_CoV_2_bt_output/')  # Back translated mRNA

# Human codon table for amino acids with usage frequency
hsap_codon_table = {'F': [('UUU', 0.45), ('UUC', 0.55)],
                    'L': [('UUA', 0.07), ('UUG', 0.13), ('CUU', 0.13), ('CUC', 0.20), ('CUA', 0.07), ('CUG', 0.41)],
                    'S': [('UCU', 0.18), ('UCC', 0.22), ('UCA', 0.15), ('UCG', 0.06), ('AGU', 0.15), ('AGC', 0.24)],
                    'Y': [('UAU', 0.43), ('UAC', 0.57)], 'Stop': [('UAA', 0.28), ('UAG', 0.20), ('UGA', 0.52)],
                    'C': [('UGU', 0.45), ('UGC', 0.55)], 'W': [('UGG', 1.00)],
                    'P': [('CCU', 0.28), ('CCC', 0.33), ('CCA', 0.27), ('CCG', 0.11)],
                    'H': [('CAU', 0.41), ('CAC', 0.59)], 'Q': [('CAA', 0.25), ('CAG', 0.75)],
                    'R': [('CGU', 0.08), ('CGC', 0.19), ('CGA', 0.11), ('CGG', 0.21), ('AGA', 0.20), ('AGG', 0.20)],
                    'I': [('AUU', 0.36), ('AUC', 0.48), ('AUA', 0.16)], 'M': [('AUG', 1.00)],
                    'T': [('ACU', 0.24), ('ACC', 0.36), ('ACA', 0.28), ('ACG', 0.12)],
                    'N': [('AAU', 0.46), ('AAC', 0.54)], 'K': [('AAA', 0.42), ('AAG', 0.58)],
                    'V': [('GUU', 0.18), ('GUC', 0.24), ('GUA', 0.11), ('GUG', 0.47)],
                    'A': [('GCU', 0.26), ('GCC', 0.40), ('GCA', 0.23), ('GCG', 0.11)],
                    'D': [('GAU', 0.46), ('GAC', 0.54)], 'E': [('GAA', 0.42), ('GAG', 0.58)],
                    'G': [('GGU', 0.16), ('GGC', 0.34), ('GGA', 0.25), ('GGG', 0.25)]}


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
wt_s_mrna = s_dna_seq[0].transcribe()
wt_s_mrna = str(wt_s_mrna)

# Generating potential mRNA-1273 sequences with all possible proline codon combinations at residue 986/3 and 987/3
i = 1  # Counter for file naming
for codon_1 in hsap_codon_table['P']:
    for codon_2 in hsap_codon_table['P']:
        mrna_1273_candidate = wt_s_mrna[:985 * 3] + codon_1[0] + codon_2[0] + wt_s_mrna[987 * 3:]
        file_output(
            possible_mrna_out,
            f'mRNA1273_Candidate_{i}.fasta',
            f'{s_dna_seq[1]} | mRNA-1273 Candidate | Proline codon substitution start at {985 * 3} and {986 * 3}',
            mrna_1273_candidate,
        )
        i += 1

# Back translating SARS-CoV-2 S2P based on codon usage frequencies
aa_s2p_sm = SeqManipulation(aa_s2p_seq)
potential_mrna_1273 = aa_s2p_sm.reverse_translate(hsap_codon_table)
file_output(
    possible_mrna_out,
    f'mRNA1273_Candidate_ReverseTranslate_CodonFrequency.fasta',
    f'{s_dna_seq[1]} | mRNA-1273 Candidate | Reverse translated based on Codon Usage Frequencies',
    potential_mrna_1273,
)