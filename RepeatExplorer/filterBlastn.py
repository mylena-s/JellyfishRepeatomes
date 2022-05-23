#!/usr/bin/env python
# coding: utf-8

# Dependecy
import pandas as pd

#load blast results and rename important columns
blastdf = pd.read_csv("all_vs_all_blast.out", sep="\t", header=None, skiprows=5)
blastdf = blastdf.rename(columns={0: 'Sequence_name_1',
                                  1:'Sequence_name_2',
                                  2:'Identity',
                                  3: 'AligmentLength',
                                  10:"evalue",
                                  12:"Seq1_len"})

def comparacion(seq1, seq2):
    "detect self-hits"
    return not seq1 == seq2

blastdf = blastdf.sort_values(["Sequence_name_1", "Sequence_name_2",
                               "AligmentLength", "Identity"],
                               ascending=False).dropna(subset=["Sequence_name_2"]).drop_duplicates(subset=["Sequence_name_1",
                                                                                                           "Sequence_name_2"],
                               keep='first')

#detect and remove self-hits
blastdf["Decision"] = blastdf.apply(lambda x: comparacion(x.Sequence_name_1,
                                                          x.Sequence_name_2),
                                    axis=1)
blastdf = blastdf.loc[blastdf["Decision"]] 
#calculate proportion of the query represented in the aligment 
blastdf["AligmentProportion_ofQuery"]=blastdf["AligmentLength"]/blastdf["Seq1_len"]
#remove reciprocal hits, keep the aligment with higher aligment length
blastdf["Both"] = blastdf[["Sequence_name_1","Sequence_name_2"]].values.tolist()
blastdf["Both"] = blastdf["Both"].apply(lambda x: str(sorted(x)))
blastdf=blastdf.sort_values("AligmentProportion_ofQuery", 
                            ascending=False).drop_duplicates(subset="Both",
                                                             keep="first")
blastdf=blastdf.loc[:,["Sequence_name_1",
                       "Sequence_name_2",
                       "Identity",
                       "AligmentLength",
                       "evalue",
                       "Seq1_len",
                       "AligmentProportion_ofQuery"]]

#remove rows that do not suprass lenght and identity thereholds 
#and save results in file
evalue=0.001
#variants
blastdf=blastdf.loc[blastdf["evalue"] <= evalue]
blastdf95=blastdf.loc[blastdf["AligmentProportion_ofQuery"] >= 0.95]
blastdf95=blastdf95.loc[blastdf95["Identity"] >= 95].to_csv("Blast95.csv",
                                                            sep="\t",
                                                            index=False)
#families
blastdf80=blastdf.loc[blastdf["AligmentProportion_ofQuery"] >= 0.8]
blastdf80=blastdf80.loc[blastdf80["Identity"] >= 80].to_csv("Blast80.csv",
                                                            sep="\t",
                                                            index=False)
#superfamilies
blastdf70=blastdf.loc[blastdf["AligmentProportion_ofQuery"] >= 0.7]
blastdf70=blastdf70.loc[blastdf70["Identity"] >= 70].to_csv("Blast70.csv",
                                                            sep="\t",
                                                            index=False)
##### end ####





