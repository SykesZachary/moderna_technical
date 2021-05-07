# Handles file input and output for genetic sequences
# 07.05.2021 - Zachary Sykes

from Bio import SeqIO
import os

from seq_manipulation import SeqManipulation

# PATH set
proj_path = os.path.abspath(os.getcwd())
sars_cov_2_aa = os.path.join(proj_path, 'SARS_CoV_2_aaSeq/')
sars_cov_2_dna = os.path.join(proj_path, 'SARS_CoV_2_dnaSeq/')

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


# Pulls sequences from files, (.fasta, or .txt)
def data_extraction(path, file):
    # CWD to sequence file
    os.chdir(path)

    # Importing fasta file
    if file[-6:] == '.fasta' or file[-3:] == '.fa':
        fasta = SeqIO.read(file, 'fasta')
        return [fasta.seq, fasta.id]

    # Importing txt file
    if file[-4:] == '.txt':
        data = []
        string = ''
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    data.append(line)
                else:
                    string += line
        data.insert(0, string)
        return data


def file_output(path, filename):
    pass


# Storing sequence data
seq = data_extraction(sars_cov_2_dna, 'sars_cov_2_s2p_dna.fasta')

# Class initialization
sm = SeqManipulation(seq[0])

print(sm.translate(codons=hsap_codon_table, seq_type='dna'))
print()
print(sm.transcribe())
