# -*- coding: utf-8 -*-
"""
Created on Thu OCT 12 2020

@author: Yanyan Yang UNC-CH
"""


#continue with remove50mer fasta file
import argparse
parser = argparse.ArgumentParser(description='nucleotideDistribution')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
input = open(inputFile, 'r')
out = open(outputFile, "w")

count = dict()

for line in input:
    if line.startswith('>') : continue
    length = len(line.strip())
    if length != 13 : continue
    if length == 13 :
       nuc = line.strip()
       
       #print(nuc)
       
       for i in range(0, len(nuc)):
           #print(nuc[i])
           n = nuc[i].upper() + str(i)
           #print(n)
          
           count[n] = count.get(n,0) + 1
#print(count.items()) 
for k, v in count.items():
    #print(k[0], k[1:3], v)
    out.write('{0}\t{1}\t{2}\n'.format(k[0].upper(), k[1:3], v))
out.close()

# bsub python nucleotideDistributionfor_damageSeq.py -i input.fa -o myOutput.txt
