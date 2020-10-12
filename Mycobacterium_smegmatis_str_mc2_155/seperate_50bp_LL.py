# -*- coding: utf-8 -*-
"""
Created on Thu OCT 12 2020

@author: Yanyan Yang UNC-CH
"""

#Input files: fasta
import os
from glob import glob
for LINES in glob("./*.fa"):
    FA = LINES.strip()
    (temp, Name) = FA.split('/',1)
    (Name, temp) = Name.split('_',1)
    F_IN = open(FA, "r")
    Len_DIS = {}

    F_IN = open(FA, "r")
    F_OUT = open('{0}_50bp.fa'.format(Name), 'w')
    F_OUT2 = open('{0}_No50.fa'.format(Name), 'w')
    for LINES in F_IN:
        if LINES[0] == ">":
            HEAD = LINES
        else:
            if len(LINES)- 1 == 50:
                F_OUT.write('{0}{1}'.format(HEAD, LINES)) 
            else:
                F_OUT2.write('{0}{1}'.format(HEAD, LINES))
        
    F_IN.close()
    F_OUT.close()
    F_OUT2.close()