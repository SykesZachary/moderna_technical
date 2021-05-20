# Class containing sequence manipulation methods
# 07.05.2021 - Zachary Sykes

from Bio.Seq import Seq


class SeqManipulation:

    def __init__(self, seq):
        self.seq = seq

    # Transcribes DNA into RNA -- Work in Progress
    def transcribe(self, rna_poly=None):
        self.seq = str(self.seq)

        # TODO - Add other rna polymerase promoter regions
        if rna_poly == 't7':
            # TODO - Add dna manipulation via promoter
            promoter = Seq('TAATACGACTCACTATAGGGAGA')  # Promoter sequence for T7 RNA polymerase

            # Storing seq to be returned
            transcribed_seq = self.seq.replace('T', 'U')
            return transcribed_seq
        else:
            # Storing seq to be returned
            transcribed_seq = self.seq.replace('T', 'U')
            return transcribed_seq

    # Translates seq into protein
    def translate(self, seq_type='rna'):
        self.seq = str(self.seq)

        # Using Biopython for simplicity
        prot = self.seq.translate()

        return prot

    # Generator to reverse translate RNA into protein based on codon set provided
    def reverse_translate(self, codons):
        # Ensures codon table provided is a dictionary
        if type(codons) is not dict:
            raise TypeError('Please provide a dictionary for the codons argument')

        # Identifies possible seq and stores it based on codon usage frequency provided
        rna_seq_lst = list()
        for aa in self.seq:
            if aa in codons.keys():
                # TODO - Identify how GC content effects the codon selected
                highest_freq = 0.0
                likely_codon = None
                for codon in codons[aa]:
                    if codon[1] > highest_freq:
                        highest_freq = codon[1]
                        likely_codon = codon
                rna_seq_lst.append(likely_codon[0])

        # Adding all three stop codons to the final sequence
        # In order from highest to lowest frequency in Human codon usage bias
        for stop_codon in ['UGA', 'UAA', 'UAG']:
            rna_seq_lst.append(stop_codon)

        return ''.join(rna_seq_lst)

    # Returns the reverse complement of DNA
    def reverse_complement(self, reverse=False):
        self.seq = str(self.seq)
        # Base complements
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

        if reverse:
            bases = list(self.seq)
            bases = reversed([complement.get(base, base) for base in bases])
            bases = ''.join(bases)
        else:
            bases = list(self.seq)
            bases = [complement.get(base, base) for base in bases]
            bases = ''.join(bases)

        return bases

    # Identifies motif in a seq
    def motif_identify(self, substring):
        i = 0
        while i < (len(self.seq) - len(substring)):
            if self.seq[i:i + len(substring)] == substring:
                print('Substring found:', end=' ')
                return len(self.seq[:i])  # Starting index for specified motif returned
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
