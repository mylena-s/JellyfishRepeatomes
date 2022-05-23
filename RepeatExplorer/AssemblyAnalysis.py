#!/usr/bin/env python
# coding: utf-8

#Dependencies
import pandas as pd
from pybedtools import BedTool
import os
import numpy as np


## define functions to process repeatmasker files into clusterscan input files
## edit RM bedfile to bed6 which is input for clusterscan
#read repeatmasker out file with all hits for all species

def calculate_length(start, end):
    if start > end :
        start, end = end, start
    return end-start

def add_sp(s):
    species_dict={"REGM":"Aaur","JAAB":"Acoe", "Seg":"Aaco",
                  "JABA":"Cqui","RQOL":"Smal","REGS":"Resc2",
                  "SWAQ":"Resc1","ML":"Nnom","OLM":"Cxam"}
    for key, item in species_dict.items():
        if key in s:
            return item

#concatenate all repeatmasker output files 
totaldf=pd.DataFrame()
path="./"
for file in os.listdir(path):
    if file.endswith("rm.bed"):
        df=pd.read_csv(path + file, sep="\t")
        df["family"]=df.apply(lambda x: x["family"]+"."+str(x.name), axis=1)
        df["species"]=df.chrom.apply(add_sp)
        df[["family", "species"]].to_csv(file+"all.desc", sep="\t", header=False, index=False)
        totaldf=totaldf.append(df)
totaldf.to_csv("all_repeats_rmbed.tsv", sep="\t", header=True, index=False)   
totaldf[["chrom","start","stop","family","size","strand","diverge"]].to_csv("all_repeats_forclusterscan.bed6", sep="\t", index=False, header=False) #create bed file needed for clusterscan, including all species and sequecences corresponding to satellite&like
totaldf[["family", "subclass"]].to_csv("all_repeats_forclusterscan.bed.desc", sep="\t", index=False, header=False) #creade description file needed for clusterscan

#after saving, reload
totaldf=pd.read_csv("all_repeats_rmbed.tsv", sep="\t")
totaldf["ElementLength"]=totaldf["subclass"].apply(lambda x: x.split("_")[1])
totaldf=totaldf.astype({"ElementLength":int})
totaldf["HitproportionElement"]=totaldf["size"]/totaldf["ElementLength"]

#read and merge results with tables containing info about satellites
sat=pd.read_csv("Satellites.csv", sep="\t")
satlike=pd.read_csv("satellite_like_repeats.csv", sep="\t")
multigene=pd.read_csv("Multigene_clusters.csv", sep="\t")
satandlike=pd.concat([sat, multigene, satlike], axis=0)
satandlike=satandlike[["Element","Family"]]
totaldf=pd.merge(totaldf, satandlike, how="inner", right_on="Element", left_on="subclass")

#define functions finding repeats near ends and gaps
def get_closest_byfeature(feature1, feature2, output_name):
    bed=BedTool(feature1).sort()
    bed2=BedTool(feature2).sort()
    nearby= bed.closest(bed2, stream=True, d=True, io=False)
    nearby.saveas(output_name)
    return bed, nearby
    
def get_closest(x):
    bed, nearbygaps=get_closest_byfeature(x[2], x[1], x[2]+"_vs_Gap.bed")
    bed, nearbyends=get_closest_byfeature(x[2], x[0], x[2]+"_vs_Contig_end.bed")
    return bed, nearbygaps, nearbyends

def edit_dfs(gaps_file="all_repeats_forclusterscan.bed6_vs_Gap.bed", ends_file="all_repeats_forclusterscan.bed6_vs_Contig_end.bed"):
    gaps = pd.read_csv(gaps_file, sep="\t", header = None)
    ends = pd.read_csv(ends_file, sep="\t", header = None)
    gaps = gaps.rename({0:"chr", 3:"name", 10:"Gap_end", 11:"Distance"}, axis=1)
    ends = ends.rename({0:"chr", 3:"name", 10:"Gap_end", 11:"Distance"}, axis=1)
    bedfiles = gaps.append(ends, ignore_index=True)
    bedfiles = bedfiles [bedfiles[7]!= "."]
    return bedfiles[bedfiles["Distance"]<=100]

x = ["allends.bed", "all_gaps.bed", "all_repeats_forclusterscan.bed6"]
#that is, gaps ends and repeats for alls sp. concatenated in a single file.
get_closest(x)
filtered=edit_dfs()
filtered.to_csv("all_Repeats_vs_gap_end_filtered.bed", index=False, sep="\t")

#read clusterscan results after running
#define function
def read_clusterscan_all(path="results/", features_all_file="features.tsv", bystanders_all_file="bystanders.tsv"):
    features = pd.read_csv(path + features_all_file, sep="\t")
    features["Distribution"] = "Clustered"
    features = features[["chr","name","cluster_id","category","Distribution"]]
    features = features.rename({"cluster_id":"Main_cluster"}, axis=1)
    bystanders = pd.read_csv(path + bystanders_all_file, sep="\t")
    bystanders = bystanders.groupby("name").agg({"chr":"unique", "cluster_id":"unique"}).reset_index()
    bystanders["chr"] = bystanders["chr"].apply(lambda x: x[0])
    bystanders = bystanders.rename({"cluster_id":"Bystander"}, axis=1)
    return features, bystanders

#run and merge with previously generated tables
features, bystanders = read_clusterscan_all()
totaldf=pd.merge(totaldf, features, how="left", left_on=["chrom", "family"], right_on=["chr", "name"])
totaldf["Distribution"]=totaldf.Distribution.fillna("Singleton")
totaldf=pd.merge(totaldf, bystanders, how="left", left_on=["chrom", "family"], right_on=["chr", "name"])
columns = ["chrom", "start", "stop", "family", 
          "size", "strand", "Element", "diverge","Family", 
          "linkage_id", "species", "ElementLength", "HitproportionElement", 
          "Distribution", "category", "Main_cluster", "Bystander"]
totaldf = totaldf[columns]
totaldf.to_csv("ClusteredRepeats_againstassembly_all.csv", sep="\t", index=False)

def combine_duplicates(gapsdf):
    gapsdf=gapsdf.astype({"Distance":float})
    g=gapsdf.groupby(["chr","name"])
    gapsdf=pd.concat([g.apply(lambda x: list(np.unique(x.Gap_end))),
                      g.Distance.min()],
                     keys=["Gap_end","Distance"],
                     axis=1).reset_index()
    return gapsdf

def combine_clusterscan_gaps(gapsdf, df):
    df=pd.merge(df, gapsdf, how="left",
                left_on=["chrom","family"],
                right_on=["chr","name"])
    df=df.reset_index()
    df["Cluster_distance"]=df.groupby(["Main_cluster"]).Distance.transform(min)
    return df

df=pd.read_csv("ClusteredRepeats_againstassembly.csv",
               sep="\t",
               low_memory=False)
output_name="Repeats_againstassembly_all.csv"
all_gaps=pd.read_csv("all_Repeats_vs_gap_end_filtered.bed", sep="\t")
all_gaps = combine_duplicates(all_gaps)
df = combine_clusterscan_gaps(all_gaps, df)

df.to_csv(output_name, sep="\t", index=False)
#repeatmasker annotates clusters as single hit when multiple copies are 
#in tandem without insertions in between, then, clusterscan can classify these
#as singletons. Therefore, to correct this, repeatmasker hits showing more than 
# twice the size of the dimer, is classified as clustered.
df.loc[df.HitproportionElement>2, "Distribution"] = "Clustered"
df.loc[df.HitproportionElement>2, "Main_cluster"] = "RepeatMasker"
df.loc[df.HitproportionElement>2, "Cluster_distance"] = df.loc[df.HitproportionElement>2, "Distance"]
df.to_csv(output_name+"repeatmasker_corrected", sep="\t", index=False)