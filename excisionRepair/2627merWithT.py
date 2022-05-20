import argparse
parser = argparse.ArgumentParser(description='fastqPut')
parser.add_argument('-i', required= True, help='input')
parser.add_argument('-o', required= True, help='output')
args = parser.parse_args()

inputFile = args.i
outputFile = args.o
input = open(inputFile, 'r')
out = open(outputFile, "w")

while 1:
    a=input.readline()
    if a=="":
        break
    if a[0]=='>':
        b=input.readline()
        if len(b)==27 and b[19].upper()=='T' or len(b)==28 and b[20].upper()=='T' or len(b)==27 and b[19].upper()!='T' and b[20].upper()=='T' or len(b)==27 and b[19].upper()!='T' and b[20].upper()!='T' and b[21].upper()=='T' or len(b)==28 and b[20].upper()!='T' and b[21].upper()=='T' or len(b)==28 and b[20].upper()!='T' and b[21].upper()!='T' and b[22].upper()=='T':
            out.writelines(a+b)
out.close()