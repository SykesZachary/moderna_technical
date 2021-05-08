# Passes sequence files to MEGA for alignment and further analysis
# 08.05.2021 -- Zachary Sykes

import os

from file_handling import data_extraction
from file_handling import file_output
from seq_manipulation import SeqManipulation

# PATH set
proj_path = os.path.abspath(os.getcwd())
sars_cov_2_aa = os.path.join(proj_path, 'SARS_CoV_2_aaSeq/')
sars_cov_2_dna = os.path.join(proj_path, 'SARS_CoV_2_dnaSeq/')
possible_mrna_out = os.path.join(proj_path, 'SARS_CoV_2_bt_output/')

# Human codon table for amino acids
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

# Storing sequence data
seq = data_extraction(sars_cov_2_aa, 'sars_cov_2_s2p_aa.fasta')

# Class initialization
sm = SeqManipulation(str(seq[0]))

i = 1  # Counter to stop generator after 10 iterations
# Pulls potential mRNA sequences from back translation generator
for rna_seq in sm.reverse_translate(hsap_codon_table):
    if i <= 10:
        file_output(
            possible_mrna_out,
            f'>{seq[1]}_mRNA_Candidate_{i}.fasta',
            rna_seq,
        )
    else:
        break
    i += 1
