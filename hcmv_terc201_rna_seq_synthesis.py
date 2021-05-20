# Analysis of C5MKY7 HCMV EGF Protein Seq and TERC-201 cDNA to Generate RNA Sequences
# 09.05.2021 - Zachary Sykes

import os

from file_handling import data_extraction
from file_handling import file_output
from seq_manipulation import SeqManipulation

# PATHs set
proj_path = os.path.abspath(os.getcwd())
hcmv_egfp = os.path.join(proj_path, 'C5MKY7_HCMV_egfp_aaSeq/')
terc_201 = os.path.join(proj_path, 'TERC201_HSap')

# Human codon table for amino acids with usage frequency -- Hard coded table
# Frequencies found on https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N
# The link above sourced it's frequencies from NCBI's Genbank
hsap_codon_table_rna_freq = {
    'F': [('UUU', 0.45), ('UUC', 0.55)],
    'L': [('UUA', 0.07), ('UUG', 0.13), ('CUU', 0.13), ('CUC', 0.20), ('CUA', 0.07), ('CUG', 0.41)],
    'S': [('UCU', 0.18), ('UCC', 0.22), ('UCA', 0.15), ('UCG', 0.06), ('AGU', 0.15), ('AGC', 0.24)],
    'Y': [('UAU', 0.43), ('UAC', 0.57)], '*': [('UAA', 0.28), ('UAG', 0.20), ('UGA', 0.52)],
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
    'G': [('GGU', 0.16), ('GGC', 0.34), ('GGA', 0.25), ('GGG', 0.25)]
}

# Storing seq objects and fasta id's
hcmv_egfp_aaseq = data_extraction(hcmv_egfp, 'c5mky7_hcmv_egfp_aa.fasta')  # Amino acid sequence for hcmv
terc_201_seq = data_extraction(terc_201, 'Homo_sapiens_TERC_201_sequence.fasta')

# Initializing sequence manipulation class with HCMV egfp aa seq
sm_hcmv_aa = SeqManipulation(hcmv_egfp_aaseq[0])
file_output(
    hcmv_egfp,
    'C5MKY7_HCMV_EGFP_ReverseTranslate_CodonFrequency.fasta',
    f'{hcmv_egfp_aaseq[1]} | C5MKY7 HCMV egfp | Reverse Translated based on Codon Usage Bias',
    sm_hcmv_aa.reverse_translate(hsap_codon_table_rna_freq),  # Calling reverse_translate method
)  # RNA sequence output for HCMV egfp

# RNA synthesis of TERC-201 cDNA
terc_201_comp = terc_201_seq[0].reverse_complement()  # Taking the reverse complement of the cDNA
terc_201_rna = terc_201_comp.transcribe()  # Transcribing the complemented strand
file_output(
    terc_201,
    'Homo_sapiens_TERC_201_rna_transcript.fasta',
    f'{terc_201_seq[1]} | TERC-201 | RNA sequence for cDNA:lncRNA',
    terc_201_rna,
)  # RNA sequence output for TERC-201 cDNA:lncRNA
