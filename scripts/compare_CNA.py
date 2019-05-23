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
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
import pandas as pd
import statsmodels.formula.api as sm
import statsmodels.stats.multitest as ssm
import patsy
import scipy.spatial
from scipy.stats import mannwhitneyu
sns.set(font_scale=2)
sns.set_style("white")

samplelist=open("Analysis/partial_clinical_info_to_plot.txt")
GD=open("Analysis/GD/absolute_GD.txt") ### changed to absolute called GD and purity. ###
purity={}
gd={}

AIS=[]
MIA=[]
INV=[]
for line in (raw.strip().split('\t') for raw in samplelist):
        if ("IIIA" not in line[2]):
                if "AIS" in line[2]:
                        AIS.append(line[0])
                if "MIA" in line[2]:
                        MIA.append(line[0])
                else:
                        INV.append(line[0])

for line in (raw.strip().split('\t') for raw in GD):
                if "sample" not in line[0]:
                        sample=line[0]
                        purity[sample]=line[3]
                        if "0" in line[1]:
                                gd[sample]="noWGD"
                        else:
                                gd[sample]="WGD"

### arm events. ###
AIS_broad=open("zyCNA3/AISMIA_segfiles/broad_values_by_arm.txt")
INV_broad=open("zyCNA3/invasive_segfiles/broad_values_by_arm.txt")
arm_event_AIS={}
arm_event_MIA={}
arm_event_INV={}
sampleName=[]

p_inv=0; p_ais=0
for line in (raw.strip().split('\t') for raw in INV_broad):
        if "Chromosome Arm" in line[0]:
                for i in line[1:len(line)]:
                        sampleName.append(i)
                        arm_event_INV[i]=0
        else:
                for i in range(0, len(sampleName)):
                        if float(line[i+1]) > 0.1 or float(line[i+1]) < -0.1:
                                arm_event_INV[sampleName[i]]=arm_event_INV[sampleName[i]]+1

        if line[0] in "6p":
                for i in range(0, len(sampleName)):
                        if float(line[i+1]) > 0.1:
                                p_inv = p_inv + 1
sampleName=[]
for line in (raw.strip().split('\t') for raw in AIS_broad):
        if "Chromosome Arm" in line[0]:
                for i in line[1:len(line)]:
                        sampleName.append(i)
                        if i in AIS:
                                arm_event_AIS[i]=0
                        if i in MIA:
                                arm_event_MIA[i]=0
        else:
                for i in range(0, len(sampleName)):
                        if float(line[i+1]) > 0.1 or float(line[i+1]) < -0.1:
                                if sampleName[i] in AIS:
                                        arm_event_AIS[sampleName[i]]=arm_event_AIS[sampleName[i]]+1
                                if sampleName[i] in MIA:
                                        arm_event_MIA[sampleName[i]]=arm_event_MIA[sampleName[i]]+1
        if line[0] in "6p":
                for i in range(0, len(sampleName)):
                        if float(line[i+1]) > 0.1:
                                p_ais = p_ais + 1

num_arm_AIS=[]
num_arm_MIA=[]
num_arm_INV=[]
num_arm_INV_n=[]
num_arm_INV_r=[]
AIS_p=[]
MIA_p=[]
INV_p=[]
AIS_WGS=[]
MIA_WGS=[]
INV_WGS=[]
sample_arm_AIS=[]
sample_arm_MIA=[]
sample_arm_INV=[]

for i in arm_event_AIS:
        if (i in purity) and (i in gd):
                sample_arm_AIS.append(i)
                num_arm_AIS.append(arm_event_AIS[i])
                AIS_WGS.append(gd[i])
                AIS_p.append(purity[i])
        else:
                print i, "GD"
for i in arm_event_MIA:
        if (i in purity) and (i in gd):
                sample_arm_MIA.append(i)
                num_arm_MIA.append(arm_event_MIA[i])
                MIA_WGS.append(gd[i])
                MIA_p.append(purity[i])
        else:
                print i, "GD"
for i in arm_event_INV:
        if (i in purity) and (i in gd):
                sample_arm_INV.append(i)
                num_arm_INV.append(arm_event_INV[i])
                INV_WGS.append(gd[i])
                INV_p.append(purity[i])
                if i in Recurrent:
                        if "Recurrence" in Recurrent[i]:
                                num_arm_INV_n.append(arm_event_INV[i])
                        if "No_recurrence" in Recurrent[i]:
                                num_arm_INV_r.append(arm_event_INV[i])

### change from here to combine three groups and use a single linear regression. ###

Ais=np.zeros(len(num_arm_AIS), dtype=int)
Mia=np.ones(len(num_arm_MIA), dtype=int)
Inv=np.full(len(num_arm_INV), 2, dtype=int)
print Ais, Mia, Inv
num_arm_AIS=np.array(num_arm_AIS).astype(float)
num_arm_MIA=np.array(num_arm_MIA).astype(float)
num_arm_INV=np.array(num_arm_INV).astype(float)
num_arm_INV_n=np.array(num_arm_INV_n).astype(float)
num_arm_INV_r=np.array(num_arm_INV_r).astype(float)
sample_arm_AIS=np.array(sample_arm_AIS)
sample_arm_MIA=np.array(sample_arm_MIA)
sample_arm_INV=np.array(sample_arm_INV)

mu=np.concatenate((Ais, Mia, Inv), axis=0)
cn=np.concatenate((num_arm_AIS, num_arm_MIA, num_arm_INV), axis=0)
pu=np.concatenate((np.array(AIS_p).astype(float), np.array(MIA_p).astype(float), np.array(INV_p).astype(float)), axis=0)
df = pd.DataFrame({"A": mu, "B": cn, "C":pu})
f = "A ~ B + C"
y, X = patsy.dmatrices(f, df, return_type='dataframe')
print "arm AIS vs MIA vs INV"
result = sm.OLS(y, X).fit()
print result.params
print result.summary()

a_t, a_prob = mannwhitneyu(num_arm_AIS, num_arm_MIA)
m_t, m_prob = mannwhitneyu(num_arm_MIA, num_arm_INV)
i_t, i_prob = mannwhitneyu(num_arm_AIS, num_arm_INV)
r_t, r_prob = mannwhitneyu(num_arm_INV_n, num_arm_INV_r)
print "arm mannwhitneyu", a_prob, m_prob, i_prob, r_prob
print "arm median number", np.median(num_arm_AIS), np.median(num_arm_MIA), np.median(num_arm_INV), np.mean(num_arm_INV_n), np.mean(num_arm_INV_r)

df_arm_gd.to_csv("Analysis/CNA/comparison_broad.source.txt", sep="\t", index=False)

### focal events. ###

AIS_focal=open("zyCNA3/AISMIA_segfiles/focal_data_by_genes.txt")
INV_focal=open("zyCNA3/invasive_segfiles/focal_data_by_genes.txt")
focal_event_AIS={}
focal_event_MIA={}
focal_event_INV={}

sampleName=[]
for line in (raw.strip().split('\t') for raw in INV_focal):
        if "Gene Symbol" in line[0]:
                for i in line[3:len(line)]:
                        sampleName.append(i)
                        focal_event_INV[i]=1
        else:
                for i in range(0, len(sampleName)):
                        if float(line[i+3]) > 1 or float(line[i+3]) < -1:
                                        focal_event_INV[sampleName[i]]=focal_event_INV[sampleName[i]]+1
                                        if (line[0] in genes):
                                                print line[0], "INV", sampleName[i]
                                                
sampleName=[]
for line in (raw.strip().split('\t') for raw in AIS_focal):
        if "Gene Symbol" in line[0]:
                for i in line[3:len(line)]:
                        sampleName.append(i)
                        if i in AIS:
                                focal_event_AIS[i]=0
                        if i in MIA:
                                focal_event_MIA[i]=0
        else:
                for i in range(0, len(sampleName)):
                        if float(line[i+3]) > 1 or float(line[i+3]) < -1:
                                #if sampleName[i] not in samples_to_exclude:                    
                                        if sampleName[i] in AIS:
                                                focal_event_AIS[sampleName[i]]=focal_event_AIS[sampleName[i]]+1

                                        if sampleName[i] in MIA:
                                                focal_event_MIA[sampleName[i]]=focal_event_MIA[sampleName[i]]+1

                                        if (line[0] in genes):
                                                print line[0], "AIS", sampleName[i]

num_focal_AIS=[]
num_focal_MIA=[]
num_focal_INV=[]
num_focal_INV_n=[]
num_focal_INV_r=[]
AIS_p=[]
MIA_p=[]
INV_p=[]
sample_focal_AIS=[]
sample_focal_MIA=[]
sample_focal_INV=[]

for i in focal_event_AIS:
        if i in purity:
                sample_focal_AIS.append(i)
                num_focal_AIS.append(focal_event_AIS[i])
                AIS_p.append(purity[i])
for i in focal_event_MIA:
        if i in purity:
                sample_focal_MIA.append(i)
                num_focal_MIA.append(focal_event_MIA[i])
                MIA_p.append(purity[i])
for i in focal_event_INV:
        if i in purity:
                sample_focal_AIS.append(i)
                num_focal_INV.append(focal_event_INV[i])
                INV_p.append(purity[i])
                if i in Recurrent:
                        if "Recurrence" in Recurrent[i]:
                                num_focal_INV_n.append(focal_event_INV[i])
                        if "No_recurrence" in Recurrent[i]:
                                num_focal_INV_r.append(focal_event_INV[i])

Ais=np.zeros(len(num_focal_AIS), dtype=int)
Mia=np.ones(len(num_focal_MIA), dtype=int)
Inv=np.full(len(num_focal_INV), 2, dtype=int)
num_focal_AIS=np.array(num_focal_AIS).astype(float)
num_focal_MIA=np.array(num_focal_MIA).astype(float)
num_focal_INV=np.array(num_focal_INV).astype(float)
num_focal_INV_n=np.array(num_focal_INV_n).astype(float)
num_focal_INV_r=np.array(num_focal_INV_r).astype(float)
sample_focal_AIS=np.array(sample_focal_AIS)
sample_focal_MIA=np.array(sample_focal_MIA)
sample_focal_INV=np.array(sample_focal_INV)

a_t, a_prob = mannwhitneyu(num_focal_AIS, num_focal_MIA)
m_t, m_prob = mannwhitneyu(num_focal_MIA, num_focal_INV)
i_t, i_prob = mannwhitneyu(num_focal_AIS, num_focal_INV)
r_t, r_prob = mannwhitneyu(num_focal_INV_n, num_focal_INV_r)
print "focal mannwhitneyu", a_prob, m_prob, i_prob, r_prob
print "focal median number", np.median(num_focal_AIS), np.median(num_focal_MIA), np.median(num_focal_INV), np.median(num_focal_INV_n), np.median(num_focal_INV_r)

df2.to_csv("Analysis/CNA/comparison_focal.source.txt", sep="\t", index=False)

mu=np.concatenate((Ais, Mia, Inv), axis=0)
cn=np.concatenate((num_focal_AIS.astype(float), num_focal_MIA.astype(float), num_focal_INV.astype(float)), axis=0)
pu=np.concatenate((np.array(AIS_p).astype(float), np.array(MIA_p).astype(float), np.array(INV_p).astype(float)), axis=0)
df = pd.DataFrame({"A": mu, "B": cn, "C":pu})
df2 = df[df["B"].astype(float) < 2500]
print "samples left", df2
f = "A ~ B + C"
y, X = patsy.dmatrices(f, df2, return_type='dataframe')
print "focal AIS vs MIA vs ADC"
result = sm.OLS(y, X).fit()
print result.summary()

### plotting correlation of focal CNA with purity. ###

plt.scatter(df2["C"][df2["A"].astype(float) == 2], df2["B"][df2["A"].astype(float) == 2], color='slateblue', alpha=1, s=50, label='LUAD')
plt.scatter(df2["C"][df2["A"].astype(float) == 1], df2["B"][df2["A"].astype(float) == 1], color='orange', alpha=1, s=50, label='MIA')
plt.scatter(df2["C"][df2["A"].astype(float) == 0], df2["B"][df2["A"].astype(float) == 0], color='lightseagreen', alpha=1, s=50, label='AIS')

### add regression line. ###
slope, intercept, r_value, p_value, std_err = stats.linregress(df2["B"], df2["C"])
print "r value, p value", r_value, p_value

sns.regplot(x=df2["C"], y=df2["B"], scatter=False)
plt.legend(loc='upper center', bbox_to_anchor=(1.15, 0.8), prop={'size': 20}, frameon=True)
plt.text(0.8, 1700, "r="+str(round(r_value, 3)), ha='center', va='bottom', fontsize=20)
plt.ylabel("Number of genes affected by focal CNA", fontsize=20)
plt.xlabel("Tumor purity", fontsize=20)
plt.savefig('Analysis/CNA/comparison_scatter.pdf', bbox_inches='tight')
plt.show()
plt.close()

df2.to_csv("Analysis/CNA/comparison_scatter.source.txt", sep="\t", index=False)

### mutation correlation with CNA. ###

genelist=open("Analysis/mutation/sig_dif_AIS_INV_h.csv")

genelist=pd.read_csv("Analysis/mutation/sig_dif_AIS_INV_h.csv", delimiter='\t')
pvals=genelist['pvalue'].astype(float)
louis, pc = ssm.fdrcorrection(pvals, alpha=0.05)
df = pd.DataFrame({"A": genelist['gene'], "B":genelist['pvalue'], "C": pc})
Genes = df['A'].loc[df['B'].astype(float) < 0.05]
focal_MUT={}; focal_name={}; gd_name={}
arm_MUT={}; arm_name={}a

maf=open("Analysis/mutation/all_merged_alteration.txt")
for line in (raw.strip().split('\t') for raw in maf):
        if line[0] in Genes.values:  ### line[0] is gene name. line[9] is sample name. ###
                if ("Splice_Site" in line[4]) or ("Missense_Mutation" in line[4]) or ("Nonsense_Mutation" in line[4]) or ("Frame_Shift_Ins" in line[4]) or ("Frame_Shift_Del" in line[4]) or ("In_Frame_Del" in line[4]) or ("In_Frame_Ins" in line[4]) or ("CNA" in line[4]):

                        if (line[9] in focal_event_AIS): ### check number of focal events in mutated samples.  ###
                                if line[0] in focal_name:
                                        if str(line[9]) not in focal_name[line[0]]:
                                                focal_name[line[0]]=focal_name[line[0]]+";"+str(line[9])
                                else:
                                        focal_name[line[0]]=str(line[9])

                        if (line[9] in focal_event_MIA):
                                if line[0] in focal_name:
                                        if str(line[9]) not in focal_name[line[0]]:
                                                focal_name[line[0]]=focal_name[line[0]]+";"+str(line[9])
                                else:
                                        focal_name[line[0]]=str(line[9])

                        if (line[9] in focal_event_INV):
                                if line[0] in focal_name:
                                        if str(line[9]) not in focal_name[line[0]]:
                                                focal_name[line[0]]=focal_name[line[0]]+";"+str(line[9])
                                else:
                                        focal_name[line[0]]=str(line[9])
                        if (line[9] in arm_event_AIS): ### check number of arm events in mutated samples. ###
                                if line[0] in arm_name:
                                        arm_name[line[0]]=arm_name[line[0]]+";"+str(line[9])
                                else:
                                        arm_name[line[0]]=str(line[9])

                        if (line[9] in arm_event_MIA):
                                if line[0] in arm_name:
                                        arm_name[line[0]]=arm_name[line[0]]+";"+str(line[9])
                                else:
                                        arm_name[line[0]]=str(line[9])

                        if (line[9] in arm_event_INV):
                                if line[0] in arm_name:
                                        arm_name[line[0]]=arm_name[line[0]]+";"+str(line[9])
                                else:
                                        arm_name[line[0]]=str(line[9])

                        if line[9] in gd:
                                if line[0] in gd_name:
                                        gd_name[line[0]]=gd_name[line[0]]+";"+str(line[9])
                                else:
                                        gd_name[line[0]]=str(line[9])
                                        
### find focal and arm numbers in wildtype samples. ###
focal_WT={}
arm_WT={}
for gene in focal_name:
        sample_array=focal_name[gene].split(";")
        for sample in focal_event_INV:
                if int(focal_event_INV[sample]) < 2500:
                        if sample not in sample_array:
                                if gene in focal_WT:
                                        focal_WT[gene]=focal_WT[gene]+";"+str(focal_event_INV[sample])
                                else:
                                        focal_WT[gene]=str(focal_event_INV[sample])
                        else:
                                if gene in focal_MUT:
                                        focal_MUT[gene]=focal_MUT[gene]+";"+str(focal_event_INV[sample])
                                else:
                                        focal_MUT[gene]=str(focal_event_INV[sample])

          
        for sample in focal_event_MIA:
                if int(focal_event_MIA[sample]) < 2500:
                        if sample not in sample_array:
                                if gene in focal_WT:
                                        focal_WT[gene]=focal_WT[gene]+";"+str(focal_event_MIA[sample])
                                else:
                                        focal_WT[gene]=str(focal_event_MIA[sample])
                        else:
                                if gene in focal_MUT:
                                        focal_MUT[gene]=focal_MUT[gene]+";"+str(focal_event_MIA[sample])
                                else:
                                        focal_MUT[gene]=str(focal_event_MIA[sample])
          
        for sample in focal_event_AIS:
                if int(focal_event_AIS[sample]) < 2500:
                        if sample not in sample_array:
                                if gene in focal_WT:
                                        focal_WT[gene]=focal_WT[gene]+";"+str(focal_event_AIS[sample])
                                else:
                                        focal_WT[gene]=str(focal_event_AIS[sample])
                        else:
                                if gene in focal_MUT:
                                        focal_MUT[gene]=focal_MUT[gene]+";"+str(focal_event_AIS[sample])
                                else:
                                        focal_MUT[gene]=str(focal_event_AIS[sample])
              
        for sample in arm_event_INV:
                if sample not in sample_array:
                        if gene in arm_WT:
                                arm_WT[gene]=arm_WT[gene]+";"+str(arm_event_INV[sample])
                        else:
                                arm_WT[gene]=str(arm_event_INV[sample])
                else:
                        if gene in arm_MUT:
                                arm_MUT[gene]=arm_MUT[gene]+";"+str(arm_event_INV[sample])
                        else:
                                arm_MUT[gene]=str(arm_event_INV[sample])

        for sample in arm_event_MIA:
                if sample not in sample_array:
                        if gene in arm_WT:
                                arm_WT[gene]=arm_WT[gene]+";"+str(arm_event_MIA[sample])
                        else:
                                arm_WT[gene]=str(arm_event_MIA[sample])
                else:
                        if gene in arm_MUT:
                                arm_MUT[gene]=arm_MUT[gene]+";"+str(arm_event_MIA[sample])
                        else:
                                arm_MUT[gene]=str(arm_event_MIA[sample])
        for sample in arm_event_AIS:
                if sample not in sample_array:
                        if gene in arm_WT:
                                arm_WT[gene]=arm_WT[gene]+";"+str(arm_event_AIS[sample])
                        else:
                                arm_WT[gene]=str(arm_event_AIS[sample])
                else:
                        if gene in arm_MUT:
                                arm_MUT[gene]=arm_MUT[gene]+";"+str(arm_event_AIS[sample])
                        else:
                                arm_MUT[gene]=str(arm_event_AIS[sample])

OUT1=open("Analysis/CNA/sig_focal.csv", 'w')
OUT2=open("Analysis/CNA/sig_arm.csv", 'w')
print >>OUT1, "gene\tpvalue\tp\tWT\tMUT\tmean_MUT\tmf"
print >>OUT2, "gene\tpvalue\tp\tWT\tMUT\tmean_MUT\tmf"
for i in focal_MUT:
        array_f_m= focal_MUT[i].split(";")
        if i in focal_WT:
                array_f_w= focal_WT[i].split(";")
        total=len(array_f_m)+len(array_f_w)
        mf=len(array_f_m)/total

        t, prob = stats.ttest_ind(np.array(array_f_m).astype(int), np.array(array_f_w).astype(int), equal_var=False)

        output=str(i)+"\t"+str(-math.log(prob, 10))+"\t"+str(prob)+"\t"+str(focal_WT[i])+"\t"+str(focal_MUT[i])+"\t"+str(np.mean(np.array(array_f_m).astype(int)))+"\t"+str(mf)
        print >>OUT1, output



for i in arm_MUT:
        array_a_m= arm_MUT[i].split(";")
        if i in arm_WT:
                array_a_w= arm_WT[i].split(";")
        total=len(array_a_m)+len(array_a_w)
        mf=len(array_a_m)/total

        t, prob = stats.ttest_ind(np.array(array_a_m).astype(int), np.array(array_a_w).astype(int), equal_var=False)

        output=str(i)+"\t"+str(-math.log(prob, 10))+"\t"+str(prob)+"\t"+str(arm_WT[i])+"\t"+str(arm_MUT[i])+"\t"+str(np.mean(np.array(array_a_m).astype(int)))+"\t"+str(mf)
        print str(i), np.mean(np.array(array_a_m).astype(int)), np.mean(np.array(array_a_w).astype(int))
        print >>OUT2, output

OUT1.close()
OUT2.close()

