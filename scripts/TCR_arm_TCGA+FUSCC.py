#!/usr/bin/python -u
from __future__ import division
#import sys, getopt
import re
import numpy as np
from scipy import stats
import argparse, os, sys
import math
import matplotlib
import matplotlib.pyplot as plt
matplotlib.pyplot.switch_backend('agg')
import seaborn as sns
import pandas as pd
import statsmodels.formula.api as sm
import patsy
import statsmodels
import numpy.linalg
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from scipy.stats import mannwhitneyu
sns.set(font_scale=2)
sns.set_style("white")

### add aneuploidy, arm CNA and immuno scores. ###
info=pd.read_csv("Analysis/aneuo_arm_immu_scores.txt", delimiter='\t')
clinical=pd.read_csv("nsclc_tcga_broad_2016/data_clinical.txt", delimiter='\t')
results = pd.merge(clinical, info, on="Sample")

INFO=results.loc[results['Type'] == "LUAD"]
smk = INFO["SMOKING_HISTORY"].values
IMS = INFO["Leuk"].values
AS = INFO["AneuploidyScore(AS)"].values
samplename=INFO["Sample"]
TMB=INFO["Non-silentMutationsperMb"].values
arms=["1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p",  "6q",  "7p",  "7q",  "8p",  "8q",  "9p",  "9q",  "10p",  "10q",  "11p", "11q", "12p", "12q", "13q", "14q", "15q", "16p", "16q", "17p", "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21q", "22q"]
Pu = INFO["Purity"].values

df = pd.DataFrame({"AS": AS, "IMS": IMS, "TMB": TMB})
df2=df.dropna()
f = "AS ~ IMS + TMB"
y, X = patsy.dmatrices(f, df2, return_type='dataframe')
result = sm.OLS(y, X).fit()
print "AS"
print result.summary()
rho, pval = stats.spearmanr(df2["AS"], df2["IMS"])
print "AS spearman", rho, pval

p_1={}; p_2={}; p_3={}

for i in arms:
        ARM=INFO[i]
        ARM1 = ARM.map({-1:-1, 1:1, 0:0})
        df = pd.DataFrame({"ARM": ARM1.values, "IMS": IMS, "TMB": TMB, "AS": AS, "SampleID": samplename})
        df2=df.dropna()
        c={}
        f = "ARM ~ IMS + TMB + AS"
        y, X = patsy.dmatrices(f, df2, return_type='dataframe')
        result = sm.OLS(y, X).fit()
        rho, pval = stats.spearmanr(df2["ARM"], df2["IMS"])

        print "TCGA spearman", i, rho, pval
        print result.summary()

        c["TCGA Leuk. Frac."] = -math.log(float(result.pvalues[1]), 10)
        p_1[i] = c

        df3 = df2[df2["TMB"].astype(float) > 0]
        if i in "6p":
                gain = df3[df3["ARM"] == 1]
                lost = df3[df3["ARM"] == -1]
                neu = df3[df3["ARM"] == 0]

                t, prob = mannwhitneyu(gain["IMS"].values, neu["IMS"].values)
                print "gain", i, t, prob
                t, prob = mannwhitneyu(lost["IMS"].values, neu["IMS"].values)
                print "loss", i, t, prob

                df4 = df3.sort(['ARM'])
                df4["ARM"] = df4["ARM"].map({-1:"loss", 1:"gain", 0:"none"})
                df4.to_csv("Analysis/TCR/6p_sm_IMS.source.txt", sep="\t")

                sns.boxplot(x="ARM", y="IMS", data=df4, palette="Set2")
                sns.swarmplot(x="ARM", y="IMS", data=df4, color=".25")
                x1, x2 = 1, 2
                y, h, col = 0.7, 0.025, 'k'
                plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
                plt.text((x1+x2)*.5, y+h, "p<0.001", ha='center', va='bottom', color=col, fontsize=24)
                x1, x2 = 0, 1
                y, h, col = 0.75, 0.025, 'k'
                plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
                plt.text((x1+x2)*.5, y+h, "p="+str(round(prob,3)), ha='center', va='bottom', color=col, fontsize=24)
                plt.xlabel(str(i)+" CNA")
                plt.title("TCGA LUAD")
                plt.ylim(0, 0.9)
                plt.ylabel("Leukocyte fraction")
                plt.savefig('Analysis/TCR/'+str(i)+'_sm_IMS.pdf', bbox_inches='tight')
                plt.show()
                plt.close()
                
TCR=open("Analysis/TCR/TCR_count_SE.txt")
tcr_total={}
tcr_se={}
tcr_type={}
for line in (raw.strip().split('\t') for raw in TCR):
        if ("sample" not in line[1]):
                tcr_total[line[1]]=line[3]
                tcr_se[line[1]]=line[4]
                tcr_type[line[1]]=line[5]

AIS_broad=pd.read_csv("zyCNA3/AISMIA_segfiles/broad_values_by_arm.txt", delimiter='\t')
INV_broad=pd.read_csv("zyCNA3/invasive_segfiles/broad_values_by_arm.txt", delimiter='\t')
TMB=pd.read_csv("Analysis/TMB/sample_TMB_all.txt", delimiter='\t')

arms_ais={}; arms_inv={}
for index, row in AIS_broad.iterrows():
        arms_ais[row.values[0]]=row.to_dict()
for index, row in INV_broad.iterrows():
        arms_inv[row.values[0]]=row.to_dict()

### only test the correlation in ADC. ###
revise=["CHG019139", "CHG020088", "CHG020926", "CHG021150", "D3244LC"]
for i in arms_inv:
        cna_inv=[]
        cna_ais=[]
        sam=[]
        tcr_t=[]
        tcr_h=[]
        tcr_tp=[]
        tmb=[]
        for l in arms_inv[i]:
                if l in tcr_total:
                        if float(tcr_total[l]) > 0:
                                if l not in revise:
                                        sam.append(l)
                                        if (float(arms_inv[i][l]) > 0.1):
                                                cna_inv.append(1)
                                        elif (float(arms_inv[i][l]) < -0.1):
                                                 cna_inv.append(-1)
                                        else:
                                                cna_inv.append(0)
                                        tcr_t.append(tcr_total[l])
                                        tcr_tp.append(tcr_se[l])
                                        sample=TMB.loc[TMB['DNA-LC'] == l]
                                        tmb.append(sample["TMB"].values[0])
        for l in arms_ais[i]:
                if l in tcr_total:
                        if float(tcr_total[l]) > 0:
                                sample=TMB.loc[TMB['DNA-LC'] == l]
                                if len(sample["TMB"].values) > 0:
                                        #tmb.append(sample["TMB"].values[0])
                                        if (float(arms_ais[i][l]) > 0.1):
                                                cna_ais.append(1)
                                        elif (float(arms_ais[i][l]) < -0.1):
                                                cna_ais.append(-1)
                                        else:
                                                cna_ais.append(0)
        f = "CNA ~ TCR + TMB"
        df = pd.DataFrame({"CNA": np.array(cna_inv).astype(float), "TCR": np.array(tcr_t).astype(float), "TMB": np.array(tmb).astype(float), "SampleID": np.array(sam)})
        df2 = df[df["TCR"].astype(float) > 0]
        df3 = df2[df2["TMB"].astype(float) < 1500000]
        y, X = patsy.dmatrices(f, df3, return_type='dataframe')
        model = LogisticRegression(fit_intercept = False, C = 1e9)
        mdl = model.fit(X, y)
        result = sm.OLS(y, X).fit()
        rho, pval = stats.spearmanr(df["CNA"], df["TCR"])

        print "FUSCC", "chr", i, rho, pval, result.pvalues[1]

        c={}

        fraction1=(len(cna_ais)-cna_ais.count(0))/len(cna_ais)
        fraction2=(len(cna_inv)-cna_inv.count(0))/len(cna_inv)
        ais_gain=(cna_ais.count(1))/len(cna_ais)
        ais_loss=(cna_ais.count(-1))/len(cna_ais)
        inv_gain=(cna_inv.count(1))/len(cna_inv)
        inv_loss=(cna_inv.count(-1))/len(cna_inv)

        c["Gain"] = inv_gain - ais_gain
        c["Lost"] = inv_loss - ais_loss


        t, pval = stats.fisher_exact([[len(cna_ais)-cna_ais.count(0), cna_ais.count(0)], [len(cna_inv)-cna_inv.count(0), cna_inv.count(0)]])
        print i, cna_ais.count(0), cna_inv.count(0), c, t, pval

        t, pval_gain = stats.fisher_exact([[ais_gain, cna_ais.count(0)], [inv_gain, cna_inv.count(0)]])
        print i, pval_gain

        t, pval_loss = stats.fisher_exact([[ais_loss, cna_ais.count(0)], [inv_loss, cna_inv.count(0)]])
        print i, pval_loss

        p_3[i]=c
        if i in "6p":
                gain = df3[df3["CNA"] == 1]
                lost = df3[df3["CNA"] == -1]
                neu = df3[df3["CNA"] == 0]
                t, prob_a = mannwhitneyu(gain["TCR"].values, neu["TCR"].values)
                print "gain FUSCC", i, t, prob_a
                t, prob_b = mannwhitneyu(lost["TCR"].values, neu["TCR"].values)
                print "loss FUSCC", i, t, prob_b

                df4 = df3.sort(['CNA'])
                df4["CNA"] = df4["CNA"].map({-1:"loss", 1:"gain", 0:"none"})
                df4.to_csv("Analysis/TCR/FUSCC_6p_sm_IMS.source.txt", sep="\t")

                sns.boxplot(x="CNA", y="TCR", data=df4, palette="Set2")
                sns.swarmplot(x="CNA", y="TCR", data=df4, color=".25")

                x1, x2 = 1, 2
                y, h, col = neu["TCR"].values.max()+0.5, 0.5, 'k'
                plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
                plt.text((x1+x2)*.5, y+h, "p="+str(round(prob_a,3)), ha='center', va='bottom', color=col, fontsize=24)
                x1, x2 = 0, 1
                y, h, col = neu["TCR"].values.max()+1.5, 0.5, 'k'
                plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
                plt.text((x1+x2)*.5, y+h, "p="+str(round(prob_b,3)), ha='center', va='bottom', color=col, fontsize=24)

                plt.xlabel(str(i)+" CNA")
                plt.ylim(0, 12)
                plt.title("FUSCC LUAD")
                plt.ylabel("TCR count")
                plt.savefig('Analysis/TCR/FUSCC_'+str(i)+'_sm_IMS.pdf', bbox_inches='tight')
                plt.show()
                plt.close()
                
fraction = pd.read_csv("Analysis/aneuo_arm_immu_scores_fraction.txt", delimiter='\t')
df = pd.DataFrame({"Gained": fraction["Gained"].values.astype(float), "Lost": fraction["Lost"].values.astype(float)})
df = df.T
df.columns = ['1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q', '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q']
print df
plt.subplots(figsize=(20,8))
sns.heatmap(df, linewidth=1.2, cbar_kws={"shrink": .4, 'label': 'arm CNA fraction', "orientation": "horizontal"})
plt.savefig('Analysis/TCR/heatmap_arm_fraction_TCGA.pdf', bbox_inches='tight')
plt.show()
plt.close()

data1=pd.DataFrame.from_dict(p_1, orient='columns')
data3=pd.DataFrame.from_dict(p_3, orient='columns')
df2 = pd.DataFrame(data1, columns = ['1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q', '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q'])
df3 = pd.DataFrame(data3, columns = [' ', '1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q', '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q'])

plt.subplots(figsize=(18,7))
sns.heatmap(df2, linewidth=1.2, cbar_kws={"shrink": .4, 'label': '-log10(p)', "orientation": "horizontal"})
plt.savefig('Analysis/TCR/heatmap_arm_tcr_TCGA+FUSCC.pdf', bbox_inches='tight')
plt.show()
plt.close()

plt.subplots(figsize=(20,7))
sns.heatmap(df3, linewidth=1.2, cbar_kws={"shrink": .4, 'label': 'increased frequency', "orientation": "horizontal"})
df3.to_csv("Analysis/TCR/heatmap_fraction_FUSCC.source.txt", sep="\t")
plt.savefig('Analysis/TCR/heatmap_fraction_FUSCC.pdf', bbox_inches='tight')
plt.show()
plt.close()

df2 = df2.T
df2.to_csv("Analysis/TCR/arm_TCR_output_TCGA+FUSCC.txt", sep='\t')

