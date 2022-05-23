#!/usr/bin/env python
# coding: utf-8

# dependencies
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from palettable.cmocean.diverging import Curl_20
import scikit_posthocs
from matplotlib_venn import venn2, venn3

#global variables
cmap = Curl_20.get_mpl_colormap()

# functions

def get_family_summary(df, column_name_ab, column_name_div, level1="Family", level2="Superfamily"):
    """Produces a summary at the specified level.
    If not indicated, standard group level is Family,
    else indicate get_family_summary(level1="Superfamily", level2="Family")"""

    g=df.groupby([level1])
    tempdf=pd.concat(
        [
            g.apply(lambda x: np.unique(x["Variant"])),
            g.apply(lambda x: np.unique(x[level2])),
            g.apply(lambda x: np.sum(x[column_name_ab])),
            g.apply(lambda x: np.average(x[column_name_div],
                                         weights=x[column_name_ab]))],
        axis=0,
        keys=["Variant",level2, column_name_ab,column_name_div],
    ).unstack().T
    return tempdf

def heatmap2(x, y, size, div, name, columns=25, lines=3):
    fig, ax = plt.subplots(figsize=(columns,lines))
    # Mapping from column names to integer coordinates
    x_labels = [v for v in x.unique()]
    y_labels = [v for v in y.unique()]
    x_to_num = {p[1]:p[0] for p in enumerate(x_labels)}
    y_to_num = {p[1]:p[0] for p in enumerate(y_labels)}
    size_scale = 1550
    newcmp=Curl_20.get_mpl_colormap()
    ax.scatter(
        x=x.map(x_to_num), y=y.map(y_to_num), s=size * size_scale,marker='s',
        c=newcmp(div))
    # Show column labels on the axes
    ax.set_xticks([x_to_num[v] for v in x_labels])
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right', fontsize=12)
    ax.set_yticks([y_to_num[v] for v in y_labels])
    ax.set_yticklabels(y_labels, fontstyle="italic", fontsize=12)
    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)
    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    plt.savefig(name)

def read_divsum(path, column_name, column_name2, sat_kept_df):
    df = pd.read_csv(path).rename({"Abundance%":column_name,
                                   "Kimura%":column_name2}, axis=1)
    df = pd.merge(df, sat_kept, how="right", left_on="Variant",
                  right_on="Variant")
    return df

## Disertation plots
### Shared satellite plots
#### Rhopilema, Aurelia and Stomolophus
# 1. Load satellite + satellite_like repeats

df1=pd.read_csv("Satellites.csv", sep="\t")
df2=pd.read_csv("satellite_like_repeats.csv", sep="\t")

#merge both tables
satellites=df1.append(df2)
# 2. Start plotting

## Rhopilema
### first, read .divsum files (harboring abundance and divergence values) and
### extract only satellites and satellite like sequences.
### to do this, merge divsum df with "satellites" df
### divsums are already processed. i.e. already harbouring "Variant"
### information normalized and on % scale.
# filter satellite df to keep only "rhopilema" satellites and avoid duplicates
sat_kept=satellites.loc[satellites.Genus=="Rhopilema",
                        ["Variant", "Family", "Superfamily"]].reset_index()
sat_kept=sat_kept.sort_values("Variant").drop_duplicates(subset="Variant")
# read divsums
# resc1 = "Rhopilema esculentum HK"
Resc1 = read_divsum("Divsum/divsums/Resc1_concatenated1238305057_processed.divsum",
                    "Resc1_ab", "Resc1_div", sat_kept)
# resc2 = "Rhopilema esculengum LI"
Resc2 = read_divsum("Divsum/divsums/Resc2_concatenated1487938001_processed.divsum",
                    "Resc2_ab", "Resc2_div", sat_kept)
# resc = "Rhopilema esculentum CL"
Resc=read_divsum("Divsum/divsums/Resc_concatenated1475847782_processed.divsum",
                 "Resc_ab", "Resc_div", sat_kept)
# merge all three datasets taking account in which level you want to plot
Resc1=get_family_summary(Resc1, "Resc1_ab", "Resc1_div", level1="Superfamily",
                         level2="Family").reset_index()
Resc2=get_family_summary(Resc2, "Resc2_ab", "Resc2_div", level1="Superfamily",
                         level2="Family").reset_index()
Resc=get_family_summary(Resc, "Resc_ab", "Resc_div", level1="Superfamily",
                        level2="Family").reset_index()
#the selection of the level will afect the comands below
#(replace superfamily by family or the inverse)
Rhopilema=pd.merge(Resc1, Resc2, how="outer", left_on="Superfamily",
                   right_on="Superfamily")
Rhopilema=pd.merge(Rhopilema, Resc, how="outer", left_on="Superfamily",
                   right_on="Superfamily")
Rhopilema.to_csv("plots/Rhopilema_abundance_comparison_SF.csv", sep="\t")
Rhop = Rhopilema[Rhopilema.Resc_ab >= 0.001]
Rhop = Rhop[Rhop.Resc1_ab >= 0.001]
Rhop = Rhop[Rhop.Resc2_ab >= 0.001]
##### shared satellites statistics #####
print("ANOVA comparison for divergence estimates of shared satellites")
print(stats.f_oneway(Rhop.Resc_div.dropna(),
                     Rhop.Resc1_div.dropna(), Rhop.Resc2_div.dropna()))
scikit_posthocs.posthoc_tukey([Rhop.Resc_div.dropna(),
                               Rhop.Resc1_div.dropna(),
                               Rhop.Resc2_div.dropna()])
print("Spearman rank correlation test between abundance values")
corr_pval=Rhop[["Resc2_ab","Resc1_ab","Resc_ab"]].corr(method=lambda x, y: stats.spearmanr(x, y)[1]) - np.eye(3)
print(Rhop[["Resc2_ab","Resc1_ab","Resc_ab"]].corr(method="spearman"))
print(corr_pval)

#### scatterplot
d = Rhopilema.loc[:,["Resc_ab","Resc1_ab","Resc2_ab"]].fillna(0)
d["max"] = d.max(axis=1)
d = d[d["max"]<=1]
sns.scatterplot(x=d.Resc1_ab, y=d.Resc_ab, color=cmap(0.3), alpha=0.9)
plt.savefig("plots/Resc1_Resc_ab_SF.svg")
plt.close()
sns.scatterplot(x=d.Resc1_ab,y=d.Resc2_ab, color=cmap(0.8), alpha=0.9)
plt.savefig("plots/Rhopilema_ab_SM.svg")

### First 15 box plots
### isolate first 15 of each species
d = Rhopilema.loc[:,["Resc1_ab","Resc2_ab","Resc_ab",
                     "Superfamily"]].set_index("Superfamily")
d_Resc = d.sort_values("Resc_ab", ascending=False).iloc[:15]
d_Resc1 = d.sort_values("Resc1_ab", ascending=False).iloc[:15]
d_Resc2 = d.sort_values("Resc2_ab", ascending=False).iloc[:15]

# merge tables
d1 = d_Resc.append(d_Resc2)
d1 = d1.append(d_Resc1)
d1 = d1.reset_index().drop_duplicates().set_index("Superfamily")
d1.to_csv("Rhopilemas.csv", sep="\t")

# modify df to plot
d1 = d1.apply(lambda x: x/x.max(), axis=1)
dd = d1.unstack().reset_index()
#set colors
col={"Resc_ab":0.3, "Resc1_ab":0.9, "Resc2_ab":0.8}
dd.columns = ['y', 'x', 'value']
dd["div"]=dd.apply(lambda x: col[x.y], axis=1)

#plot
plot_name = "plots/Rhopilema_15_abSF.svg"
heatmap2(
    x=dd['x'],
    y=dd['y'], size=dd["value"], div=dd["div"], name=plot_name)

## Aurelia
# filter satellite df to keep only "Aurelia" satellites and avoid duplicates
sat_kept=satellites.loc[satellites.Genus=="Aurelia",
                        ["Variant", "Family", "Superfamily"]].reset_index()
sat_kept=sat_kept.sort_values("Variant").drop_duplicates(subset="Variant")


# "Aaco" = Aurelia coerulea CA
Aaco=read_divsum("Divsum/divsums/Aau_concatenated992767818_processed.divsum",
                 "Aaco_ab", "Aaco_div", sat_kept)

# "Aaur" = Aurelia aurita BS
Aaur=read_divsum("Divsum/divsums/AauSRR7992476_3082509672_concatenated_processed.divsum", "Aaur_ab", "Aaur_div", sat_kept)

# "Acoe" = Aurelia coerulea IR
Acoe=read_divsum("Divsum/divsums/Acoe_1484851017_concatenated_processed.divsum",
                 "Acoe_ab","Acoe_div", sat_kept)

# merge all three datasets taking account in which level you want to plot
Aurelia=pd.merge(Aaco[["Variant","Aaco_ab","Aaco_div"]],
                 Aaur[["Variant","Aaur_ab","Aaur_div"]],
                 how="outer", left_on="Variant", right_on="Variant")
Aurelia1=pd.merge(Aurelia,
                  Acoe[["Variant","Acoe_ab","Acoe_div"]],
                  how="outer", left_on="Variant", right_on="Variant")

#the selection of the level will afect the comands below (replace superfamily by family or the inverse)
Aaco=get_family_summary(Aaco, "Aaco_ab", "Aaco_div",
                        level1="Superfamily", level2="Family").reset_index()
Aaur=get_family_summary(Aaur, "Aaur_ab", "Aaur_div",
                        level1="Superfamily", level2="Family").reset_index()
Acoe=get_family_summary(Acoe, "Acoe_ab", "Acoe_div",
                        level1="Superfamily", level2="Family").reset_index()
Aurelia=pd.merge(Aaco, Aaur, how="left", left_on="Superfamily",
                 right_on="Superfamily")
Aurelia=pd.merge(Aurelia, Acoe, how="left",
                 left_on="Superfamily", right_on="Superfamily")
Aurelia.to_csv("Aurelia_abundance_comparison.csv", sep="\t")

Aurel = Aurelia[Aurelia.Aaur_ab>=0.001]
Aurel= Aurel[Aurel.Aaco_ab>=0.001]
Aurel = Aurel[Aurel.Acoe_ab>=0.001]

##### shared satellites statistics #####

print("ANOVA divergence comparisons")
print(stats.f_oneway(Aurel.Aaur_div.dropna(),
                     Aurel.Acoe_div.dropna(),
                     Aurel.Aaco_div.dropna()))
print(scikit_posthocs.posthoc_tukey([Aurel.Aaur_div.dropna(),
                                     Aurel.Acoe_div.dropna(),
                                     Aurel.Aaco_div.dropna()]))
print("Abundance comparisons by spearman rank test")
corr = Aurel[["Aaco_ab","Acoe_ab","Aaur_ab"]].corr(method="spearman")
print(corr)
corr_pval = Aurel[["Aaco_ab","Acoe_ab",
                   "Aaur_ab"]].corr(method=lambda x, y: stats.spearmanr(x, y)[1]) - np.eye(3)
print(corr_pval)
corr.to_csv("plots/Aurelia_SF_corr.txt")
corr_pval.to_csv("plots/Aurelia_SF_corr_pval.txt")
#### scatterplot
d = Aurelia.loc[:,["Aaco_ab","Acoe_ab","Aaur_ab"]]
d["max"]=d.max(axis=1)
d = d[d["max"]<=1]
sns.scatterplot(x=d.Aaco_ab,y=d.Acoe_ab, color=cmap(0.8), alpha=0.9)
plt.savefig("Aaco_Acoe_ab.svg")
sns.scatterplot(x=d.Aaco_ab,y=d.Aaur_ab, color=cmap(0.3), alpha=0.9)
plt.savefig("plots/Aurelia_ab_SF.svg")
plt.close()

### First 15 box plots
### isolate first 15 of each species

d=Aurelia.loc[:,["Aaco_ab","Acoe_ab","Aaur_ab",
                 "Superfamily"]].set_index("Superfamily")
d_Aaur=d.sort_values("Aaur_ab", ascending=False).iloc[:15]
d_Acoe=d.sort_values("Acoe_ab", ascending=False).iloc[:15]
d_Aaco=d.sort_values("Aaco_ab", ascending=False).iloc[:15]
d1 = d_Aaur.append(d_Acoe)
d1 = d1.append(d_Aaco).drop_duplicates()
d1.max(axis=1)


#modify df to plot
d1=d1.apply(lambda x: x/x.max(), axis=1)
d1.to_csv("plots/Aurelias_SF.csv", sep="\t")
dd=d1.unstack().reset_index()
#define colors
col={"Aaur_ab":0.3, "Aaco_ab":0.9, "Acoe_ab":0.8}
dd.columns = ['y', 'x', 'value']
dd["div"]=dd.apply(lambda x: col[x.y], axis=1)

#plot
heatmap2(
     x=dd['x'],
     y=dd['y'], size=dd["value"], div=dd["div"],
     name="plots/Aurelia_15_ab_SF.svg", columns=len(d1))

#### Stomomolophus #####
# filter satellite df to keep only "Aurelia" satellites and avoid duplicates
sat_kept=satellites.loc[satellites.Genus=="Stomolophus",
                        ["Variant", "Family", "Superfamily"]].reset_index()
sat_kept=sat_kept.sort_values("Variant").drop_duplicates(subset="Variant")

# "Ssp2" = Stomolophus sp. 2 MX
Ssp2 = read_divsum("Divsum/divsums/Smel_concatenated1475752148_processed.divsum",
                   "Ssp2_ab", "Ssp2_div",sat_kept)

# "Ssp2_2" = Stomolophus sp. 2 JP
Ssp2_2 = read_divsum("Divsum/divsums/Ssp2_1482714203_concatenated_processed.divsum",
                     "Ssp2_2_ab", "Ssp2_2_div", sat_kept)

# merge both df
Ssp2 = get_family_summary(Ssp2, "Ssp2_ab", "Ssp2_div",level1="Superfamily",
                          level2="Family").reset_index()
Ssp2_2 = get_family_summary(Ssp2_2, "Ssp2_2_ab", "Ssp2_2_div",
                            level1="Superfamily", level2="Family").reset_index()
Stomolophus = pd.merge(Ssp2, Ssp2_2, how="left", left_on="Superfamily",
                       right_on="Superfamily")
Stomolophus.to_csv("Stomolophus_abundance_SFcomparison.csv", sep="\t")

## Statistics for shared satellites
# divergence comparison by T-test
print(stats.ttest_ind(Stomolophus.Ssp2_div, Stomolophus.Ssp2_2_div))
# abundance comparison by spearmanr
stats.spearmanr(Stomolophus.Ssp2_ab, Stomolophus.Ssp2_2_ab)

#### scatterplot

d = Stomolophus.loc[:,["Ssp2_ab","Ssp2_2_ab"]]
d["max"] = d.max(axis=1)
d = d[d["max"]<=1]

sns.scatterplot(x=d.Ssp2_ab,y=d.Ssp2_2_ab, color=cmap(0.3), alpha=0.9)
plt.savefig("plots/Ssp2_SFab.svg")
plt.close()

####first 15 plot
### isolate first 15 of each species and merge
d=Stomolophus.loc[:,["Ssp2_ab","Ssp2_2_ab",
                     "Superfamily"]].set_index("Superfamily")
d["max"]=d.max(axis=1)
d_Ssp2 = d.sort_values("Ssp2_ab", ascending=False).iloc[:15]
d_Ssp2_2 = d.sort_values("Ssp2_2_ab", ascending=False).iloc[:15]
d1 = d_Ssp2_2.append(d_Ssp2)

#modify df to plot
d1.drop(columns="max", inplace=True)
d1 = d1.drop_duplicates()
d1=d1.apply(lambda x: x/x.max(), axis=1)
d1.to_csv("plots/Superfamily/Stomolophus_15.csv", sep="\t")
dd = d1.unstack().reset_index()
#define colors
col = {"Ssp2_ab":0.8, "Ssp2_2_ab":0.9}

dd.columns = ['y', 'x', 'value']
dd["div"]=dd.apply(lambda x: col[x.y], axis=1)

heatmap2(
    x=dd['x'],
    y=dd['y'], size=dd["value"], name="plots/Stomolophus_10firstSF.svg",
    div=dd["div"], lines=2, columns=len(d1))

####### venn graphs ########
#Aurelias

Aaur=11
Aacoe=27
Acoe=9
Acoerul=38+1 #1 disco sat +38 

aureliavenn=venn3(subsets = (Aaur, Aacoe,2, Acoe, 2, Acoerul,2), 
      set_labels=["Aurelia aurita\nBaltic Sea", "Aurelia aurita\nCalifornia", "Aurelia coerulea"],
     set_colors=("#b7526fff","#751a5cff","#6bae92ff"), alpha = 0.7)
plt.savefig("plots/Aurelia_superfamily_Venn.svg", dpi=300)
plt.close()
#Stomolophus
Stomolophus=30
SSp2=4
SSp22=7
Stomolophusvenn=venn2(subsets= (SSp2, SSp22, Stomolophus),
                     set_labels=["Stomolophus sp.2 Mexico", "Stomolophus sp.2 Japan"],
                     set_colors=("#b7526fff","#751a5cff"), alpha=0.7)
plt.savefig("plots/Stomolophus_superfamily_venn.svg", dpi=300)
plt.close()

#Rhopilema
Resc=7
Resc1=12
Resc2=6
Rhopilemas=25
Resc12=1+29+Rhopilemas#1 disco sat
Resc01=1+Rhopilemas
Resc02=1+Rhopilemas
rhopilemavenn=venn3(subsets = (Resc, Resc2,Resc02, Resc1, Resc01, Resc12,Rhopilemas), 
      set_labels=["R. cf. esculentum", "R. esculentum (Li et al. 2020)", "R. esculentum (Nong et al. 2020)"],
     set_colors=("#6bae92ff","#b7526fff","#751a5cff"), alpha = 0.7)
plt.savefig("plots/rhopilema_superfamily_venn.svg", dpi=300)
plt.close()
######## END ########
