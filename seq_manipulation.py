# Class containing sequence manipulation methods
# 07.05.2021 - Zachary Sykes

from Bio import SeqIO
import itertools


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

    # Generator to reverse translate RNA into protein based on codon set provided
    def reverse_translate(self, codons):
        # identifies possible seq and stores it -- May be far too intense
        rna_seq_lst = [codons[aa] for aa in self.seq]
        for comb in itertools.product(*rna_seq_lst):
            yield ''.join(comb)

    # Returns the reverse complement of DNA for cDNA handling
    def reverse_complement(self):
        # Base complements
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

        bases = list(self.seq)
        bases = reversed([complement.get(base, base) for base in bases])
        bases = ''.join(bases)

        return bases

    # Returns information about the sequence provided
    def seq_stats(self, seq_type):
        # Helps convert fasta files fed as Seq type
        self.seq = str(self.seq)

        # TODO - Need to return base densities for nuc seq
        if seq_type == 'aa':
            return len(self.seq)
        elif seq_type == 'dna':
            return len(self.seq)
        elif seq_type == 'rna':
            return len(self.seq)
        else:
            raise TypeError('Please enter "aa", "dna", or "rna" for seq_type arg')
