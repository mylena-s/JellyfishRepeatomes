#!/usr/bin/env/python3
import sys, os
import argparse
from Bio import SeqIO

def filter_cov(original_file, coverage):
    new_file= original_file.split(".")[0]+"_"+ str(coverage) + ".fasta"

    with open(original_file) as original, open(new_file, "w") as new1:
        for record in SeqIO.parse(original_file, "fasta"):
            #Removes sequecens if they are smaller than the lenght and are not annotated as one of the exceptions
            if float(record.description.split(";")[4].split("=")[1]) >= coverage: 
                SeqIO.write(record, new1, "fasta")
    return new_file 

def count_sequences(file):
    count=0
    with open(file):
        for record in SeqIO.parse(file, "fasta"):
            count+=1
    return count

def calc_diference_n(original_file, new_file):
    return count_sequences(original_file) - count_sequences(new_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script removes sequences containing incomplete domains")
    parser.add_argument("-f", "--file", type=str, required=True, help="specify the fasta file", metavar="")
    parser.add_argument("-c", "--cov", type=float, default="0.8", help="specify the minimum domain coverage to keep", metavar="")
    args = parser.parse_args()
    original_file= args.file
    coverage = args.cov
    
    file_name= filter_cov(original_file, coverage)
    deleted_sequences= calc_diference_n(original_file, file_name)
    print("Files: \n"
          +str(file_name) +": filtered by domain coverage\n"
          + "have been created and " + str(deleted_sequences) + " sequences have been deleted")

