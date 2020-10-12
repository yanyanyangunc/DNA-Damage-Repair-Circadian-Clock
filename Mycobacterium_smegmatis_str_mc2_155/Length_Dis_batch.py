# -*- coding: utf-8 -*-
"""
Created on Thu OCT 12 2020

@author: Yanyan Yang UNC-CH
"""

#!/usr/bin/python

#module needed python

#To analysis the length and the nucleotide distribution
#Input files: fasta

import os
from glob import glob

for file in glob("./*.fa"):
    file = file.strip()
    (x, name) = file.split('/', 1)
    (name, x) = name.split('.', 1)
    F_IN = open(file, "r")
    Len_DIS = {}
    for LINES in F_IN:
        if LINES[0] == ">":
            continue
        else:
            LEN = len(LINES)- 1
            try: Len_DIS[LEN] = Len_DIS[LEN] + 1
            except: Len_DIS[LEN] = 1
    F_OUT = open('{0}_lent_dis.txt'.format(name), 'w')
    F_OUT.write('Length\tReads\n')
    for X, V in Len_DIS.items():
        F_OUT.write('{0}\t{1}\n'.format(X, V))
    F_IN.close()
    F_OUT.close()
