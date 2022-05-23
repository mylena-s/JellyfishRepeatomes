#!/usr/bin/env python
# coding: utf-8

#dependencies
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv
import seaborn as sns
from palettable.cmocean.diverging import Curl_20
import matplotlib

#define variables
cmp=Curl_20.get_mpl_colormap()

# configure matplotlib so letters are exported as text
new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'}
matplotlib.rcParams.update(new_rc_params)

#define functions

def divsum_splitter(file): 
    file_modif = file.split('.')[0]
    with open(file) as file, open( file_modif + '.csv', 'w') as modif:    
        orig_table = csv.reader(file)
        for row in orig_table:
            for index in row:
                if 'Div' in index:
                    file_csv = csv.writer(modif)
                    file_csv.writerow(row)
                    for row in orig_table:
                        csv.writer(modif)
                        file_csv.writerow(row)
    return file_modif + '.csv'

def total_values(filename, totalpb, satellites=None, filter_df=False):
    df=pd.read_csv(filename, sep="\s+").set_index("Div")
    if filter_df:
        columns=df.columns.intersection(satellites).to_list()
        df=df[columns]
    return df.sum(axis=1).div(totalpb).multiply(100)

def main(divsum_list="divsum_list_comparison.csv", element_list=None, filter_df=False, acronim="", title="Repeat landscape"):
    divsumlist=pd.read_csv(divsum_list, sep="\t", header=None).set_index(4)
    divsumlist[0]=divsumlist[0].apply(lambda x: divsum_splitter(x))
    divsumlist[2]=divsumlist[2].apply(lambda x: divsum_splitter(x))
    Raw_reads=divsumlist[[0,1]].apply(lambda x: total_values(x[0], x[1]), axis=1).T
    Assembly=divsumlist[[2,3]].apply(lambda x: total_values(x[2], x[3]), axis=1).T
    for Species in Raw_reads.columns.to_list():
        landscape_comparisonplot(Raw_reads, Assembly, Species, legendtitle="Tandem repeat landscape", acronim="tandem")
    if filter_df:
        Raw_reads=divsumlist[[0,1]].apply(lambda x: total_values(x[0], x[1], filter_df=True, satellites=element_list), axis=1).T
        Assembly=divsumlist[[2,3]].apply(lambda x: total_values(x[2], x[3], filter_df=True, satellites=element_list), axis=1).T
    resultdf=pd.DataFrame()
    for Species in Raw_reads.columns.to_list():
        tempdf=pd.DataFrame()
        landscape_comparisonplot(Raw_reads, Assembly, Species, acronim=acronim, legendtitle=title)
        tempdf.loc[Species, "Species"] = Species
        tempdf.loc[Species, "Proportion of Assembly (%)"] = Assembly.loc[:, Species].sum()
        tempdf.loc[Species, "Proportion of library (%)"] = Raw_reads.loc[:, Species].sum()
        resultdf=resultdf.append(tempdf, ignore_index=True)
    resultdf["Difference between estimates"]=resultdf["Proportion of Assembly (%)"]-resultdf["Proportion of library (%)"]
    resultdf.to_csv("against assembly/Repeat_landscape_comparisons.csv",sep="\t")
    return resultdf

def landscape_comparisonplot(readsdf, assemblydf, Species, limit=40, path="plots/landscapes/", acronim="", legendtitle="Repeat landscape"):
    cmp=Curl_20.get_mpl_colormap()
    fig,ax= plt.subplots(figsize=[6.4, 4])
    sns.barplot(y=readsdf.loc[:limit,Species], x=readsdf.loc[:limit].index,  color=cmp(0.9),  ax=ax, label="From reads")
    sns.barplot(y=assemblydf.loc[:limit,Species], x=assemblydf.loc[:limit].index, color=cmp(0.3), alpha=0.75, label="From assembly", ax=ax)
    plt.xlabel("Divergence (K2P)",fontsize=12)
    plt.ylabel("Genome proportion (%)", fontsize=12)
    plt.xticks
    fig.legend(title=legendtitle,  fontsize=10, loc='upper right', frameon=False,bbox_to_anchor=(0.4, 0.38, 0.5, 0.5))
    x_ticks = np.arange(0, limit, 5)
    ax.set_xticks(x_ticks)
    ax.tick_params(labelsize=13)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.title(label=Species, style="italic")
    plt.savefig(path+Species+"_"+acronim+"_comparisonls2.svg", dpi=300)
    return fig

def get_mean(df, level="Family", level2="family", stats="count"):
    if stats=="count":
        mean=df.groupby(["species", level])[level2].count()
        mean.index=mean.index.droplevel(1)
        mean=mean.reset_index().groupby("species")[level2].mean()
    elif stats == "sum":
        mean=df.groupby(["species", level])[level2].sum() 
        mean.index=mean.index.droplevel(1)
        mean=mean.reset_index().groupby("species")[level2].mean()
    
    if stats == "max":
        mean=df.groupby(["species", level])[level2].max()
        mean.index=mean.index.droplevel(1)
        mean=mean.reset_index().groupby("species")[level2].max()
    
    return mean

def get_stat(df, glevel, column, statistic="sum"):
    g = df.groupby(glevel)
    stat = g[column].sum()
    if statistic == "count":
        stat = g[column].count()
    if statistic == "mean":
        stat = g[column].mean()
    return stat 

#read tandem repeat tables and merge
satellites=pd.read_csv("Satellites.csv", sep="\t")
multigene=pd.read_csv("Multigene_clusters.csv", sep="\t")
satellitesandlike=pd.concat([satellites, multigene], axis=0)
element_list=satellitesandlike.Element.apply(lambda x: "Sat/"+x).to_list()

#plot repeat landscapes with both RB and AB estimates
d=main(divsum_list="divsum_list_comparison.csv", element_list=element_list, filter_df=True, acronim="SatelliteLikeMultigene", title="Multigene and satellite landscape")

### Repeat fragmentation and distribution ###
#read clusterscan + distance to gaps analysis
df = pd.read_csv("against assembly/tandem_check/Repeats_againstassembly_all.csvrepeatmasker_corrected", sep="\t")
#add unique identified to Repeatmasker clusters
df.loc[df.Main_cluster == "RepeatMasker","Main_cluster"] = df.loc[df.Main_cluster == "RepeatMasker",:].apply(lambda x: x.Main_cluster+str(x["index"]), axis=1) 
# calculate largest cluster
biggest = df.groupby(["species", "Main_cluster"])["size"].sum().reset_index()
biggest = biggest.groupby("species")["size"].max().reset_index()
pd.concat([d, biggest], axis=1).to_csv("divsum_comparison.csv", sep="\t")


#modify dataframe to rename sequences into proximal and distal, 
# into singletons and clusters
df.loc[df["Gap_end"].notna(),"Gap_end"]="Proximal"
df["Cluster_Gap_end"]=df.Gap_end
df["Distribution"] = df["Distribution"].fillna("Singleton")
df.loc[df.Cluster_distance.notna(),"Cluster_Gap_end"]=df[df.Cluster_distance.notna()].Cluster_Gap_end.fillna("Proximal")
df["Gap_end"]=df.Gap_end.fillna("Distal")
df["Cluster_Gap_end"]=df.Cluster_Gap_end.fillna("Distal")
d = dict(tuple(df.groupby('species')))

distribution=df.groupby(["species","Distribution","Cluster_Gap_end"])["size"].sum().reset_index()
distribution.to_csv("clusterscan.csv", sep="\t")
df.groupby("species").Main_cluster.nunique()
df[df.Distribution=="Singleton"].groupby("species").chrom.count()

for key, item in d.items():
    fig, ax = plt.subplots(1, 1, figsize=[8, 3], sharex=True)
    data=d[key].groupby(["Distribution","Cluster_Gap_end"])["size"].sum().reset_index()
    data["size"] = data["size"] /  1000000
    p= sns.barplot(data=data, x="Distribution", y="size", hue="Cluster_Gap_end",palette=[cmp(0.1),cmp(0.3)])
    p.set_ylabel("Total length (Mbp)", fontsize=13)
    p.set_xlabel("Distribution", fontsize=13)
    p.legend(fontsize=12, frameon=False)
    ax.tick_params(labelsize=13)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig("plots/landscapes/"+key+"Distribution_hits.svg", dpi=300)
    plt.close()



factor= 1000000 #MBP


DistrubLength=df.groupby(["species", "Distribution"])["size"].sum().reset_index().set_index("species").pivot(columns = "Distribution")
DistrubLength.columns=DistrubLength.columns.droplevel(0)

g=df.groupby("species")
summary_species=pd.concat([g.family.nunique(),#number of hits
                   g.HitproportionElement.sum().round(), #number of copies
                   (g["size"].sum()/factor).round(1), #total repeat length
                   g.HitproportionElement.mean(), #dimer coverage
                   g.Main_cluster.nunique(), #number of  clusters
                   g.Distribution.apply(lambda x: np.unique(x, return_counts=True)[1][0]), #number of clustered hits
                   g.Gap_end.apply(lambda x: np.unique(x,return_counts=True)[1][1]),
                   DistrubLength],# number of proximal hits
                   keys = ["Hit count", "Copy number","Total repeat lemgth", 
                             "monomer coverage","Cluster count","Clustered hits","Proximal hits","Length"], axis=1)


Proximaldf = df[df.Gap_end == "Proximal"]
clustered_proximal = df[(df["Distribution"]=="Clustered") & (df["Cluster_Gap_end"]=="Proximal")]
clusterdf=df.groupby(["species","Main_cluster"]).agg({"start":min, "stop":max})
clusterdf["Length"]=clusterdf.stop-clusterdf.start
clusterdf.index=clusterdf.index.droplevel(1)
clusteredf=df[df.Distribution == "Clustered"]
singletons=df[df.Distribution == "Singleton"]


other_summary = pd.concat([ get_stat(Proximaldf, "species", "HitproportionElement"), # number of proximal copies
                            get_stat(clustered_proximal, "species", "family", "count"), #number of clustered hits proximal
                            get_stat(clustered_proximal, "species", "HitproportionElement"), # number of clustered ccopies proximal
                            get_stat(clusteredf, "species", "HitproportionElement", "mean"), #mean monomer coverage
                            get_stat(clusteredf, "species", "HitproportionElement"), #number of clustered copies
                            get_mean(df), #mean hits per family
                            get_mean(df, level="Family", level2="HitproportionElement", stats="sum"), #mean copies per family
                            get_mean(df, level="Family", level2="Main_cluster", stats="count"), #mean clusters per family
                            get_mean(df, level="Family", level2="size", stats="sum")/1000, #mean family len
                            get_stat(clusterdf, "species", "Length", "mean"),
                            get_stat(singletons, "species", "HitproportionElement", "mean")], #mean cluster length per sp
                            keys = ["number of proximal copies", "number of proximal clustered hits","number of proximal clustered copies", "monomer coverage of clustered", "clustered copy number",
                                   "mean hits per family", "mean copies per family per sp", "mean clusters per family per sp", "mean family len per sp", "mean cluster length per sp",
                                   "singleton hit proportion"],
                          axis=1)




output_name = "against assembly/Results_corrected.xlsx"

summary = pd.concat([summary_species,other_summary], axis=1)
with pd.ExcelWriter(output_name) as writer:
    summary.to_excel(writer, sheet_name='Summary')



