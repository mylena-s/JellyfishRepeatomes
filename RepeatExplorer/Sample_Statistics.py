#!/usr/bin/env python
# coding: utf-8

#dependencies
import pandas as pd
import scipy.stats as stats
import numpy as np
import scikit_posthocs

#read table with all satellite information
#CG comparison between species
satellites=pd.read_csv("Satellites.csv", sep="\t")

#add taxonomic information to df

Discomedusae=["Aurelia", "Sanderia", "Chrysaora","Nemopilema", "Rhopilema", "Stomolophus", "Cassiopea","Lychnorhiza"]
Coronamedusae=["Nausithoe", "Linuche","Thecoscyphus"]
Order={}

for element in Discomedusae:
    Order[element]="Discomedusae"
for element in Coronamedusae:
    Order[element]="Coronamedusae"
    
satellites["Order"]=satellites.Genus.apply(lambda x: Order[x])



#statistics
# "species" means sample
#CG content distribution:  
###test to compare CG distribution of the diff species
print(stats.kruskal(*[group["CG"].values for name, group in satellites.groupby("Species")]))
res=scikit_posthocs.posthoc_dunn(satellites, val_col="CG", group_col="Species", p_adjust='fdr_bh').round(3)
res.to_csv("stats/CG_dunn_sample.csv", sep ="\t")

## compare at order level
print(stats.kruskal(*[group["CG"].values for name, group in satellites.groupby("Order")]))
res=scikit_posthocs.posthoc_dunn(satellites, val_col="CG", group_col="Order", p_adjust='fdr_bh').round(3)
res.to_csv("stats/CG_dunn_order.csv", sep ="\t")


## CG content vs genome size? 
gz=pd.read_csv("gz.txt", sep="\t", header=None)
df = pd.merge(satellites, gz, how="left", right_on=0, left_on="Species")
g=df.groupby("Species_group")
print(stats.spearmanr(g[1].unique().apply(lambda x: x[0]), g["CG"].median()))

#compare length distribution between samples
print(stats.kruskal(*[group["Consensuslength"].values for name, group in satellites.groupby("Species")]))
res2=scikit_posthocs.posthoc_dunn(satellites, val_col="Consensuslength", group_col="Species", p_adjust='fdr_bh').round(3)
res2.to_csv("stats/consensuslen_sample.csv", sep="\t")

#compare at order level
print(stats.kruskal(*[group["Consensuslength"].values for name, group in df.groupby("Order")]))

#Consensus length vs genome size? 
print(stats.spearmanr(g[1].unique().apply(lambda x: x[0]), np.average(g["Consensuslength"], weights=g["Abundance%"])))
#divergence vs genome size?
print(stats.spearmanr(g[1].unique().apply(lambda x: x[0]), g.apply(lambda x: np.average (x["Kimura%"], weights=x["Abundance%"]))))


#compare kimura distributions between samples
print(stats.kruskal(*[group["Kimura%"].values for name, group in df.groupby("Species")]))
res2=scikit_posthocs.posthoc_dunn(satellites, val_col="Kimura%", group_col="Species", p_adjust='fdr_bh').round(3)
res2.to_csv("stats/Kimura_dunn_sample.csv", sep="\t")
#compare kimura distributions between orders
print(stats.kruskal(*[group["Kimura%"].values for name, group in df.groupby("Order")]))

## correlation between kimura and consensus length in each sample
g=satellites.groupby("Species")

corr=pd.DataFrame()
for group in g:
    rho, p = stats.spearmanr(group[1].Consensuslength,group[1]["Kimura%"])
    corr.loc[group[0],"Rho"] = rho
    corr.loc[group[0],"pval"] = p
    corr.loc[group[0], "sample_size"] = len(group[1])
corr.to_csv("stats/kimura_vs_length.csv", sep="\t")

#correlation between kimura and consensus length in all samples altogether
stats.spearmanr(satellites.Consensuslength, satellites["Kimura%"])

## correlation between kimura and GC in each samples
corr=pd.DataFrame()
for group in g:
    rho, p = stats.spearmanr(group[1].CG,group[1]["Kimura%"])
    corr.loc[group[0],"Rho"] = rho
    corr.loc[group[0],"pval"] = p
    corr.loc[group[0], "sample_size"] = len(group[1])
corr.to_csv("stats/CG_vs_length.csv", sep="\t")