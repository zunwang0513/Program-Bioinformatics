import os
import sys
import multiprocessing
import subprocess
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i1")
    parser.add_argument("-i2")
    parser.add_argument("-t")
    parser.add_argument("-o")
    args = parser.parse_args()

    if args.i1:
        input1 = args.i1

    if args.i2:
        input2 = args.i2

    if args.o:
        outputname = args.o

    if args.t:
        seqtype = args.t

    file1 = open(input1,'r')
    file2 = open(input2,'r')
    file1content = file1.read()

    if seqtype == "n":
        makedatabase1 = "makeblastdb -in " + file1 + " -dbtype nucl -parse_seqids -out database1 -title 'db1'"
        makedatabase2 = "makeblastdb -in " + file2 + " -dbtype nucl -parse_seqids -out database2 -title 'db2'"
    elif seqtype == "p":
        makedatabase1 = "makeblastdb -in " + file1 + " -dbtype prot -parse_seqids -out database1 -title 'db1'"
        makedatabase2 = "makeblastdb -in " + file2 + " -dbtype prot -parse_seqids -out database2 -title 'db2'"

    mkdb1 = subprocess.call(makedatabase1)
    mkdb2 = subprocess.call(makedatabase2)

    if seqtype == "n":
        blast1 = "blastn -query " + file1 + " -db db2 -out output1.txt"
        blast2 = "blastn -query " + file2 + " -db db1 -out output2.txt"
    elif seqtype == "p":
        blast1 = "blastp -query " + file1 + " -db db2 -out output1.txt"
        blast2 = "blastp -query " + file2 + " -db db1 -out output2.txt"

    output1 = open("output1.txt",'r')
    output2 = open("output2.txt",'r')
    genedict1 = {}
    genedict2 = {}
    reciprocal = {}
    for row in output1:
        genedict1[row.split()[0]] = row.split()[2]
    for row in output2:
        genedict2[row.split()[2]] = row.split()[0]


    for key in genedict1.keys():
        if genedict1[key] in genedict2.keys():
            reciprocal[key] = genedict2[genedict1[key]]

    outputfile = open(outputname,'a')
    for key, value in reciprocal.items():
        outputfile.write(key+value)
    outputfile.close()

    readme = open("readme.txt",'a')
    readme.write("number of initial hit in input1 is " + str(len(genedict1)))
    readme.write("number of initial hit in input2 is " + str(len(genedict2)))
    readme.write("number of orthologous genes is " + str(len(reciprocal)))
    readme.close()

    output1.close()
    output2.close()
            
        
    
        

    
    
        
