#!/usr/bin/python -u
from __future__ import division
import sys, getopt
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
import scipy.spatial
from scipy.stats import mannwhitneyu
sns.set(font_scale=2)
sns.set_style("white")

sig_sum=open("Analysis/signature/all_sample_summary.txt")
samplelist=open("Analysis/partial_clinical_info_to_plot.txt")
output=open("Analysis/signature/comparison_source.txt", 'w')

sm_sig={}
apo_sig={}
AIS_sample=[];AIS_tmb=[]
MIA_sample=[];MIA_tmb=[]
INV_sample=[];INV_tmb=[]
smoking_stats={}

for line in (raw.strip().split('\t') for raw in samplelist):
        sample=line[0]
        subtype=line[2]
        if "IIIA" not in line[2]:
                if "AIS" in line[2]:
                        AIS_sample.append(sample)
                elif "MIA" in line[2]:
                        MIA_sample.append(sample)
                else:
                        INV_sample.append(sample)

                if "Never" in line[5]:
                        smoking_stats[sample]="non-smoker"
                else:
                        smoking_stats[sample]="smoker"


for line in (raw.strip().split('\t') for raw in sig_sum):
        if "name06" not in line[0]:
                sample=line[1]
                apobec=line[7]
                smoking=line[8]
                apo_sig[sample]=round(float(apobec), 5)
                sm_sig[sample]=round(float(smoking), 5)

for sample in smoking_stats:
        if sample in apo_sig:
                if sample in AIS_sample:
                        print >>output, "AIS/MIA\t", sample, "\t", apo_sig[sample], "\t", sm_sig[sample], "\t", smoking_stats[sample]
                elif sample in MIA_sample:
                        print >>output, "AIS/MIA\t", sample, "\t", apo_sig[sample], "\t", sm_sig[sample], "\t", smoking_stats[sample]
                else:
                        print >>output, "LUAD\t", sample, "\t", apo_sig[sample], "\t", sm_sig[sample], "\t", smoking_stats[sample]
        else:
                print sample
ais_apo=[]
mia_apo=[]
inv_apo=[]
sr_ais_smk=[]
sr_mia_smk=[]
sr_inv_smk=[]

for i in apo_sig:
        if i in AIS_sample:
                ais_apo.append(apo_sig[i])
                sr_ais_smk.append(smoking_stats[i])
        if i in MIA_sample:
                mia_apo.append(apo_sig[i])
                sr_mia_smk.append(smoking_stats[i])
        if i in INV_sample:
                inv_apo.append(apo_sig[i])
                sr_inv_smk.append(smoking_stats[i])

ais_apo=np.array(ais_apo).astype(float)
mia_apo=np.array(mia_apo).astype(float)
inv_apo=np.array(inv_apo).astype(float)

a_t, a_prob = mannwhitneyu(ais_apo, mia_apo)
m_t, m_prob = mannwhitneyu(ais_apo, inv_apo)
i_t, i_prob = mannwhitneyu(mia_apo, inv_apo)
print "apobec", a_t, a_prob, m_t, m_prob, i_t, i_prob
### apobec excluding smokers. ###

ais_apo=[]
mia_apo=[]
inv_apo=[]
smoking_stats[sample]="smoker"
for i in apo_sig:
        if i in smoking_stats:
                if "non-smoker" in smoking_stats[i]:
                        if float(sm_sig[i]) > 0.8:
                                print "removed two high smoking signature", sm_sig[i]
                        if i in AIS_sample:
                                ais_apo.append(apo_sig[i])
                        if i in MIA_sample:
                                mia_apo.append(apo_sig[i])
                        if i in INV_sample:
                                inv_apo.append(apo_sig[i])

ais_apo=np.array(ais_apo).astype(float)
mia_apo=np.array(mia_apo).astype(float)
inv_apo=np.array(inv_apo).astype(float)

a_t, a_prob = mannwhitneyu(ais_apo, mia_apo)
m_t, m_prob = mannwhitneyu(ais_apo, inv_apo)
i_t, i_prob = mannwhitneyu(mia_apo, inv_apo)
print "apobec excluding smokers", a_t, a_prob, m_t, m_prob, i_t, i_prob

### test smoking signature. ###

ais_smk=[]
mia_smk=[]
inv_smk=[]
sr_ais_smk=[]
sr_mia_smk=[]
sr_inv_smk=[]

for i in sm_sig:
        if i in AIS_sample:
                ais_smk.append(sm_sig[i])
                sr_ais_smk.append(smoking_stats[i])
        if i in MIA_sample:
                mia_smk.append(sm_sig[i])
                sr_mia_smk.append(smoking_stats[i])
        if i in INV_sample:
                inv_smk.append(sm_sig[i])
                sr_inv_smk.append(smoking_stats[i])

ais_smk=np.array(ais_smk).astype(float)
mia_smk=np.array(mia_smk).astype(float)
inv_smk=np.array(inv_smk).astype(float)

a_t, a_prob = mannwhitneyu(ais_smk, mia_smk)
m_t, m_prob = mannwhitneyu(ais_smk, inv_smk)
i_t, i_prob = mannwhitneyu(mia_smk, inv_smk)
print "smoking", a_t, a_prob, m_t, m_prob, i_t, i_prob

### mutation correlation with CNA. ###

genelist=open("Analysis/mutation/sig_dif_AIS_INV_h.csv")
Genes=[]
genelist=pd.read_csv("Analysis/mutation/sig_dif_AIS_INV_h.csv", delimiter='\t')
pvals=genelist['pvalue'].astype(float)
df = pd.DataFrame({"A": genelist['gene'], "B":genelist['pvalue']})
Genes = df['A'].loc[df['B'].astype(float) < 0.05]
#Genes = df['A'].loc[df['B'].astype(float) < 1]
apo_MUT={}; apo_name={}
tmb_MUT={}; tmb_name={}
TMB={}
### find TMB in mutatn samples. ###

tmb=open("Analysis/TMB/sample_TMB_all.txt")
for line in (raw.strip().split('\t') for raw in tmb):
        if "TMB" not in line[2]:
                sample=line[1]
                #TMB[sample]=round(float(line[2]), 3)
                TMB[sample]=line[2]

maf=open("Analysis/mutation/all_merged_alteration.txt")
for line in (raw.strip().split('\t') for raw in maf):
        if line[0] in Genes.values:  ### line[0] is gene name. line[9] is sample name. ###
                if ("Splice_Site" in line[4]) or ("Missense_Mutation" in line[4]) or ("Nonsense_Mutation" in line[4]) or ("Frame_Shift_Ins" in line[4]) or ("Frame_Shift_Del" in line[4]) or ("In_Frame_Del" in line[4]) or ("In_Frame_Ins" in line[4]) or ("CNA" in line[4]):
                        if line[9] in apo_sig:
                                if line[0] in apo_name:
                                        if str(line[9]) not in apo_name[line[0]]:
                                                apo_name[line[0]]=apo_name[line[0]]+";"+str(line[9])
                                else:
                                        apo_name[line[0]]=str(line[9])
                        if line[9] in TMB:
                                if line[0] in tmb_name:
                                        if str(line[9]) not in tmb_name[line[0]]:
                                                tmb_name[line[0]]=tmb_name[line[0]]+";"+str(line[9])
                                else:
                                        tmb_name[line[0]]=str(line[9])

### find apobec in mutant samples. ###
apo_WT={}
tmb_WT={}

for gene in apo_name:
        sample_array=apo_name[gene].split(";")
        for sample in apo_sig:                       ### exclude smokers and stage III ### 
                if sample in TMB:
                        if sample not in sample_array:
                                if gene in apo_WT:
                                        apo_WT[gene]=apo_WT[gene]+";"+str(apo_sig[sample])
                                else:
                                        apo_WT[gene]=str(apo_sig[sample])
                        else:
                                if gene in apo_MUT:
                                        apo_MUT[gene]=apo_MUT[gene]+";"+str(apo_sig[sample])
                                else:
                                        apo_MUT[gene]=str(apo_sig[sample])

for gene in tmb_name:
        sample_array=tmb_name[gene].split(";")
        for sample in TMB:
                if sample in apo_sig:
                        if sample not in sample_array:
                                if gene in tmb_WT:
                                        tmb_WT[gene]=tmb_WT[gene]+";"+str(TMB[sample])
                                else:
                                        tmb_WT[gene]=str(TMB[sample])
                        else:
                                if gene in tmb_MUT:
                                        tmb_MUT[gene]=tmb_MUT[gene]+";"+str(TMB[sample])
                                else:
                                        tmb_MUT[gene]=str(TMB[sample])

OUT1=open("Analysis/signature/sig_signature.csv", 'w')
OUT2=open("Analysis/TMB/sig_tmb.csv", 'w')
print >>OUT1, "gene\tpvalue\tp\tWT\tMUT\tmean_MUT\tmf"
print >>OUT2, "gene\tpvalue\tp\tWT\tMUT\tmean_MUT\tmf"

for i in apo_MUT:
        array_f_m= apo_MUT[i].split(";")
        if i in apo_WT:
                array_f_w= apo_WT[i].split(";")
        total=len(array_f_m)+len(array_f_w)
        mf=len(array_f_m)/total
        t, prob = stats.ttest_ind(np.array(array_f_m).astype(float), np.array(array_f_w).astype(float), equal_var=False)
        #if t > 0:
        output=str(i)+"\t"+str(-math.log(prob, 10))+"\t"+str(prob)+"\t"+str(apo_WT[i])+"\t"+str(apo_MUT[i])+"\t"+str(np.mean(np.array(array_f_m).astype(float)))+"\t"+str(mf)
        print >>OUT1, output
        #else:
        #       output=str(i)+"\t"+str(0)+"\t"+str(apo_WT[i])+"\t"+str(apo_MUT[i])+"\t"+str(np.mean(np.array(array_f_m).astype(float)))+"\t"+str(mf)
        #       print >>OUT1, output

for i in tmb_MUT:
        array_f_m= tmb_MUT[i].split(";")
        if i in tmb_WT:
                array_f_w= tmb_WT[i].split(";")
        total=len(array_f_m)+len(array_f_w)
        mf=len(array_f_m)/total
        t, prob = stats.ttest_ind(np.array(array_f_m).astype(float), np.array(array_f_w).astype(float), equal_var=False)
        #if t > 0:
        output=str(i)+"\t"+str(-math.log(prob, 10))+"\t"+str(prob)+"\t"+str(tmb_WT[i])+"\t"+str(tmb_MUT[i])+"\t"+str(np.mean(np.array(array_f_m).astype(float)))+"\t"+str(mf)
        print >>OUT2, output
        #else:
        #       output=str(i)+"\t"+str(0)+"\t"+str(tmb_WT[i])+"\t"+str(tmb_MUT[i])+"\t"+str(np.mean(np.array(array_f_m).astype(float)))+"\t"+str(mf)
        #       print >>OUT2, output
OUT1.close()
OUT2.close()
