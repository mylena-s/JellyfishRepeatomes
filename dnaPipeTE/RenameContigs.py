#!/usr/bin/env python3

import argparse
import sys, os
import csv
from Bio import SeqIO
import pandas as pd

def read_TESorter(table):
    df = pd.read_csv(table, sep="\t")
    df.replace("TIR","DNA", inplace=True)
    df.replace("TIR","DNA", inplace=True)
    df.loc[df["Order"]=="DIRS",["Order","Superfamily"]]="LTR","DIRS"
    df.loc[df["Order"]=="Maverick",["Order","Superfamily"]]="DNA","Maverick"
    df.loc[df["Order"]=="Penelope",["Order","Superfamily"]]="LINE", "Penelope"
    df ["New_name"] =  df["#TE"]+"#"+df["Order"]+"/"+df["Superfamily"]
    df["New_name"]=df["New_name"].astype(str)
    df["Program"]="TESORTER"
    return df

def read_rmasker(table2):
    df = pd.read_csv(table2, sep="\t")
    df=pd.read_csv(table2, sep="\t", header=None, names=["#TE","qlen","qprop","Element","Classification", "hlen","hcord","hprop","other"])
    df.replace("RC/Helitron","Helitron/unknown", inplace=True)
    df ["New_name"] =  df["#TE"]+"#"+df["Classification"]
    df["New_name"]=df["New_name"].astype(str)
    df["Program"]="RMASKER"
    return df

  
def rename_sequences(fasta_file, Element, New_name, open_file, acronim):
    with open(fasta_file):
        for record in SeqIO.parse(fasta_file, "fasta"):
            if Element == record.id:
                record.id = New_name
                if acronim != None:
                    record.id= acronim +"_"+ record.id
                record.description = ""
                SeqIO.write(record, new, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""This script changes the record names inside a fasta file acording to two tables. The format of the original and resulting ids is in classical RepeatMasker format: >Sequence_name#Repeat. """)
    parser.add_argument("-f", "--file", type=str, required=True, help="specify the fasta file", metavar="")
    parser.add_argument("-ts", "--Tesorter_table", type=str, required=True, help="specify the tsv table resulting from TEsorter run", metavar="")
    parser.add_argument("-rm", "--RM_table", type=str, required=True, help="specify the dnapipete RM table: one_RM_hit_per_Trinity_contigs", metavar="")
    parser.add_argument("-sp", "--species_name", type=str, default=None, help="Sp name acronim to add as prefix", metavar="")

    #defining variables names
    args = parser.parse_args()
    fasta_file= args.file
    table=args.Tesorter_table
    table2=args.RM_table
    acronim=args.species_name
    #running program
    df1=read_TESorter(table)
    df2=read_rmasker(table2)
    #combine both tables
    df1=df1.append(df2).sort_values("#TE").drop_duplicates(subset="New_name").sort_values("#TE").reset_index()
    #remove duplicates with same annotation
    df1["Duplicated"]=df1.duplicated(subset="#TE", keep=False)
    #keep TEsorter annotation over RepeatMasker's
    indexes=df1.loc[(df1["Program"]=="TESORTER") & (df1["Duplicated"]==True)].index
    df1=df1.drop(indexes)
    df1.replace("Unknown","Unclassified", inplace=True)
    new_fasta = fasta_file.split(".")[0]+"_renamed.fasta"
    #rename fasta file
    with open(new_fasta, "a+") as new:
        df1.apply(lambda x: rename_sequences(fasta_file, x["#TE"], x.New_name, new, acronim), axis=1)      
        

