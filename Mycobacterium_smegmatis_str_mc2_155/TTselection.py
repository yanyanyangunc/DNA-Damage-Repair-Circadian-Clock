# -*- coding: utf-8 -*-
"""
Created on Thu OCT 12 2020

@author: Yanyan Yang UNC-CH
"""

import argparse
parser = argparse.ArgumentParser(description='fastqPut')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
input = open(inputFile, 'r')
out = open(outputFile, "w")

input_mem=input.read()
input_mem=input_mem.split('\n')
input_len=len(input_mem)/2
scan_point=1

key_list = ['TT']


while scan_point<=input_len:
    judge=input_mem[scan_point*2-1][-6:-4]
    for key in key_list:
        if key in judge:
            out.writelines(input_mem[scan_point*2-2]+'\n'+input_mem[scan_point*2-1]+'\n')
            break
    scan_point=scan_point+1
out.close()
