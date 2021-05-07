# Class containing sequence manipulation methods
# 07.05.2021 - Zachary Sykes

from Bio import SeqIO
# import itertools
import os


class SeqManipulation:

    def __init__(self, seq):
        self.seq = seq

    # Transcribes DNA into RNA
    def transcribe(self):
        # Storing seq to be returned
        transcribed_seq = self.seq.replace('T', 'U')
        return transcribed_seq

    # Translates seq into protein
    def translate(self, codons, seq_type='rna'):
        # Storing aa seq to be returned
        prot = ''
        i = 0  # Counter for while loop

        # Ensures codon table provided is a dictionary
        if type(codons) is not dict:
            raise TypeError('Please provide a dictionary for the codons argument')

        # Transcribes DNA to RNA prior to translation
        if seq_type == 'dna':
            sm = SeqManipulation(self.seq)  # Initializing the class
            self.seq = sm.transcribe()  # Calling the previous method

        while i < len(self.seq):
            for k, v in codons.items():
                if self.seq[i:(i + 3)] in v:
                    if k == 'Stop':
                        break
                    else:
                        prot += k
            i += 1

        return prot

    # Reverse translates RNA into protein based on codon set provided
    def reverse_translate(self, codons):
        pass

