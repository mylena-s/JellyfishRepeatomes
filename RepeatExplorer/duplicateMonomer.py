#!/usr/bin/env python3

import sys, os
from Bio import SeqIO
import argparse

def manipulate_files(fasta):
    return fasta, "duplicated_" + fasta

def duplicate_monomer(fasta, n):
    original_file, duplicated_file = manipulate_files(fasta)

    with open(original_file), open(duplicated_file, "w") as duplicated:
        for record in SeqIO.parse(original_file, "fasta"):        
            record.seq = record.seq * n
            SeqIO.write(record, duplicated, "fasta")

def achieve_lenght(fasta, lenght):
    original_file, duplicated_file = manipulate_files(fasta)
    with open(original_file), open(duplicated_file, "w") as duplicated:
        for record in SeqIO.parse(original_file, "fasta"):        
            monomer= record.seq
            if len(monomer) >= lenght:
                record.seq=record.seq+record.seq
            else:
                while len(record.seq)< lenght:
                    record.seq = record.seq + monomer
            SeqIO.write(record, duplicated, "fasta")

def rename_headers(fastq_file, word):
    original_file = fastq_file
    corrected_file = fastq_file.split(".")[0]+"renamed"

    with open(original_file), open(corrected_file, "w") as corrected:
        for record in SeqIO.parse(original_file, "fasta"):        
            record.id = record.id + word + "/" + record.id
            record.description= ""
            SeqIO.write(record, corrected, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script duplicates N times each record inside a fasta file. If argument -w is used, this script also adds a word in each id in RepeatMasker library format")
    parser.add_argument("-f", "--file", type=str, required=True, help="specify the fasta file", metavar="")
    parser.add_argument("-n", "--number", type=int, default= 2, help="number of times to duplicate monomer (default =2)", metavar="")
    parser.add_argument("-t", "--type", type=int, default= 1, choices= [1, 2], help="running mode: 1= duplicate n times (used with -n), 2=duplicate to reach lenght (used with -l)", metavar="")
    parser.add_argument("-l", "--lenght", type=int, default=1, help="specify the minimun lenght that you want to achieve by duplicating", metavar="")
    parser.add_argument("-w", "--word", type=str, required=False, help="specify a string to be added to the sequence id. Example: #Unknown, #Satellite", metavar="")

    args = parser.parse_args()
    fasta  = args.file
    n = args.number
    mode = args.type
    lenght = args.lenght
    word= args.word
    
    if mode == 1:
        duplicate_monomer(fasta, n)
    else:
        if lenght == 1:
            print("Monomer will not be duplicated since no lenght was specified")
        else:
            achieve_lenght(fasta, lenght)

    if word != None:
        rename_headers(("duplicated_" + fasta), word)
        
    print("""Done! Thanks!
contact: mylena.santander@usp.br""")
