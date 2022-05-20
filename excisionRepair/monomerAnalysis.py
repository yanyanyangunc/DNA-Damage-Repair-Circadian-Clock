#module load python/3.6.6
from sys import version_info
import argparse
from os import listdir
from statistics import mean
from Bio import SeqIO
import requests
import string

parser = argparse.ArgumentParser(description='fastqPut')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
fasta_in = open(inputFile, 'r')
bed_out = open(outputFile, "w")
    
desired_lengths = list(range(26, 28))
dict_of_monomers = {
	length: {
		i + 1: {
			'A': 0,
			'T': 0,
			'G': 0,
			'C': 0
		} for i in range(length)
	} for length in desired_lengths
}

with open(inputFile) as fasta:

    # Count the number of monomers
    for entry in SeqIO.parse(fasta, 'fasta'):
        seq = entry.seq.upper()
        length = len(seq)
        if length in desired_lengths:
            for pos, nt in enumerate(seq):
                pos += 1  # position is 1 based, not 0 based
                dict_of_monomers[length][pos][nt] += 1

    # Convert numbers to ratios
    for length in dict_of_monomers:
        for pos in dict_of_monomers[length]:
            count_sum = sum(dict_of_monomers[length][pos].values())
            for nt, count in dict_of_monomers[length][pos].items():
                try:
                    dict_of_monomers[length][pos][nt] = count / count_sum
                except ZeroDivisionError:
                    dict_of_monomers[length][pos][nt] = 0

    # Convert to R data frame and write
    #outputFile = 'results/' + fasta[:-3] + '_monomer_R_df.txt'
    with open(outputFile, 'a') as f:
        for length in dict_of_monomers:
            for pos in dict_of_monomers[length]:
                for nt, ratio in dict_of_monomers[length][pos].items():
                    #f.write(f'{length}\t{pos}\t{nt}\t{ratio}\n')
                    f.write('{0}\t{1}\t{2}\t{3}\n'.format(length, pos, nt, ratio))
