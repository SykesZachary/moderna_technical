# Handles file input and output for genetic sequences
# 07.05.2021 - Zachary Sykes

from Bio import SeqIO
import os


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


def file_output(path, filename, data):
    # CWD to output directory
    os.chdir(path)

    # Opening output file and writing to it
    out_file = open(filename, 'w+')
    out_file.write(filename)
    out_file.write('\n')

    # Splitting the lines into lengths of 80 aa or bases
    for idx, item in enumerate(data):
        out_file.write(item)
        if (idx + 1) / 80 in range(1, 101):
            out_file.write('\n')

