# Class containing sequence manipulation methods
# 07.05.2021 - Zachary Sykes

import itertools
from Bio.Seq import Seq


class SeqManipulation:

    def __init__(self, seq):
        self.seq = seq

    # Transcribes DNA into RNA
    def transcribe(self, rna_poly='t7'):
        self.seq = str(self.seq)

        # TODO - Add other rna polymerase promoter regions
        # TODO - Add dna manipulation via promotor
        if rna_poly == 't7':
            promoter = Seq('TAATACGACTCACTATAGGGAGA')  # Promoter sequence for T7 RNA polymerase

        # Storing seq to be returned
        transcribed_seq = self.seq.replace('T', 'U')
        return transcribed_seq

    # Translates seq into protein
    def translate(self, codons, seq_type='rna'):
        self.seq = str(self.seq)
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

    # Generator to reverse translate RNA into protein based on codon set provided
    def reverse_translate(self, codons):
        # identifies possible seq and stores it -- May be far too intense
        rna_seq_lst = [codons[aa] for aa in self.seq]
        for comb in itertools.product(*rna_seq_lst):
            yield ''.join(comb)

    # Returns the reverse complement of DNA for cDNA handling
    def reverse_complement(self):
        self.seq = str(self.seq)
        # Base complements
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

        bases = list(self.seq)
        bases = reversed([complement.get(base, base) for base in bases])
        bases = ''.join(bases)

        return bases

    # Identifies motif in a seq
    def motif_identify(self, substring):
        i = 0
        while i < (len(self.seq) - len(substring)):
            if self.seq[i:i + len(substring)] == substring:
                print('Substring found')
                return len(self.seq[:i + 1])  # Starting index for specified motif returned
            i += 1

    # Returns information about the sequence provided
    def seq_stats(self, seq_type):
        # Helps convert fasta files fed as Seq type

        if seq_type == 'aa':
            print(len(self.seq))  # Seq length
        elif seq_type == 'dna':
            seq_len = len(self.seq)  # Seq length

            # Calculating seq gc content
            g_content = self.seq.count('G')
            c_content = self.seq.count('C')
            total_gc = g_content + c_content
            gc_content = (total_gc/seq_len) * 100

            # Diagnostics returned
            print(f'Sequence length: {seq_len}')
            print(f'Sequence type: {seq_type.upper()}')
            print(f'GC content: {gc_content}%')

        elif seq_type == 'rna':
            seq_len = len(self.seq)  # Seq length

            # Calculating seq gc content
            g_content = self.seq.count('G')
            c_content = self.seq.count('C')
            total_gc = g_content + c_content
            gc_content = (total_gc / seq_len) * 100

            # Diagnostics returned
            print(f'Sequence length: {seq_len}')
            print(f'Sequence type: {seq_type.upper()}')
            print(f'GC content: {gc_content}%')

        else:
            raise TypeError('Please enter "aa", "dna", or "rna" for seq_type arg')
