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
import pandas as pd
sns.set(font_scale=2)
sns.set_style("white")

### read in clincial information, TMB, purity, mutaiton, RNA-seq read count. ###
to_exclude_bach=["CHG019733", "CHG019540","CHG019681","CHG018352","CHG019764","CHG019781","CHG020285","CHG021801","CHG021772","CHG021050","CHG019679","CHG020845", "CHG019539", "CHG019682", "CHG018351", "CHG019765", "CHG019780", "CHG020284", "CHG021802", "CHG021771", "CHG021049", "CHG019680", "CHG020844"]

batch_info=open("batch_info.txt")
for line in (raw.strip().split('\t') for raw in batch_info):
        if line[4] == "3":
                to_exclude_bach.append(line[0])

mn={};sub={};label={}
C={};DNA={}
RNAlist=open("Analysis/samplelist_Recurrence_status.txt")
for line in (raw.strip().split('\t') for raw in RNAlist):
        if ("RNA-LC" not in line[0]) and (line[0] not in to_exclude_bach):
                mn[line[0]]=line[1]
                label[line[0]]="Tumor"
                label[line[1]]="Normal"
                if "NA" in line[5]:
                        sub[line[0]]="AIS/MIA"
                        sub[line[1]]="AIS/MIA"
                else:
                        sub[line[0]]="LUAD"
                        sub[line[1]]="LUAD"

                DNA[line[0]]=line[2]

rna_count=open("MIXCR/count_RNA_all.txt")
for line in (raw.strip().split('\t') for raw in rna_count):
        sample=line[0]
        count=line[1]
        C[sample]=count
 
 
SE={}
total_RPM={}
tcr_se=open("Analysis/TCR/TCR_count_SE.txt")
for line in (raw.strip().split('\t') for raw in tcr_se):
        total_RPM[line[0]]=line[3]
        if "NA" not in line[4]:
                SE[line[0]]=line[4]

expression=pd.read_csv("/cga/meyerson/Data/FUSCC_LUAD/RNA_seq/RSEM_results/TPM_RSEM_without_gene_id_duplications_deleted.csv", delimiter=',')
#genes=["CD274", "CD8A", "CD8B", "CD3E", "CD2", "CD3D", "CD4"]
RRR_P={}
#for gene in genes:
cd8b= expression.loc[expression['gene_name'] == "CD8A"]
cd4 = expression.loc[expression['gene_name'] == "CD3E"]
CD8B=cd8b.to_dict()
CD4=cd4.to_dict()
### find clone count per type. ###

tcr_seen={}
top_clone={}
n_top_clone={}
m_SE={}
for sample in mn:
        if (sample in total_RPM) and (mn[sample] not in to_exclude_bach):
                print sample
                tcr='MIXCR/'+str(sample)+'.clones.txt'
                tcr_file=open(tcr)
                for line in (raw.strip().split('\t') for raw in tcr_file):
                        if "cloneId" not in line[0]:
                                if ("TRA" in line[5]) or ("TRB" in line[5]):
                                        if sample not in top_clone:
                                                #RPM=float(line[1])
                                                RPM=float(line[1])/float(C[sample])*1000000
                                                top_clone[sample]=RPM/float(total_RPM[sample])
                                                top_seq=line[3]
                                                print RPM, float(total_RPM[sample])
                n_tcr='MIXCR/'+str(mn[sample])+'.clones.txt'
                tcr_file=open(n_tcr)
                if sample in top_clone:
                        for line in (raw.strip().split('\t') for raw in tcr_file):
                                if ("TRA" in line[5]) or ("TRB" in line[5]):
                                        if line[3] == top_seq:
                                                #RPM=float(line[1])
                                                RPM=float(line[1])/float(C[mn[sample]])*1000000
                                                n_top_clone[mn[sample]]=RPM/float(total_RPM[mn[sample]])

output=open("Analysis/TCR/normal_ais_increase.source.txt",'w')
for sample in top_clone:
        if mn[sample] in n_top_clone:
                print >>output, DNA[sample], "\t", label[sample], "\t", sub[sample], "\t", top_clone[sample], "\t", n_top_clone[mn[sample]]

print "new", len(top_clone), len(n_top_clone)
OUT=open("Analysis/TCR/normal_LUAD_SE.source.txt", 'w')

for sample in sub:
        if (sub[sample]=="LUAD") and (sample in mn):
                if (sample in SE) and (mn[sample] in SE):
                                m_SE[sample]=float(SE[sample])
                                m_SE[mn[sample]]=float(SE[mn[sample]])
                                print >>OUT, DNA[sample], "\t", m_SE[sample], "\t", mn[sample], "\t", m_SE[mn[sample]]
                                df = pd.DataFrame({" ": ["Normal", "LUAD"], "se":[m_SE[mn[sample]], m_SE[sample]]})
                                if m_SE[sample] < m_SE[mn[sample]]:
                                        sns.pointplot(x=" ", y="se", color="purple", data=df, alpha=0.7)
                                        #plt.plot(names, values, linestyle='-', color='red', marker='o')
                                if m_SE[sample] > m_SE[mn[sample]]:
                                        sns.pointplot(x=" ", y="se", color="grey", data=df, alpha=0.7)
                                        #plt.plot(names, values, linestyle='-', color='blue', marker='o')

plt.text(0.5, 2.32, "Paired t test p=0.90", ha='center', va='bottom', fontsize=24)
plt.ylabel("Entropy score of T cell clones")
plt.savefig('Analysis/TCR/normal_LUAD_SE.pdf', bbox_inches='tight')
plt.show()
plt.close()

OUT=open("Analysis/TCR/normal_AIS_SE.source.txt", 'w')
for sample in sub:
        if (sub[sample]=="AIS/MIA") and (sample in mn):
                if (sample in SE) and (mn[sample] in SE):
                        #if float(SE[mn[sample]]) > 0.9:
                                m_SE[sample]=float(SE[sample])
                                m_SE[mn[sample]]=float(SE[mn[sample]])
                                df = pd.DataFrame({" ": ["Normal", "AIS/MIA"], "se":[m_SE[mn[sample]], m_SE[sample]]})
                                print >>OUT, DNA[sample], "\t", m_SE[sample], "\t", mn[sample], "\t", m_SE[mn[sample]]
                                if m_SE[sample] < m_SE[mn[sample]]:
                                        g = sns.pointplot(x=" ", y="se", color="purple", data=df, alpha=0.7, legend_out = True)
                                        #plt.plot('x', 'se', data=df, linestyle='-', color='red', marker='o')
                                if m_SE[sample] > m_SE[mn[sample]]:
                                        g = sns.pointplot(x=" ", y="se", color="grey", data=df, alpha=0.7)
                                #plt.plot('x', 'se', data=df, linestyle='-', color='blue', marker='o')

plt.text(0.5, 2.31, "Paired t test p=0.13", ha='center', va='bottom', fontsize=24)
plt.ylabel("Entropy score of T cell clones")
plt.show()
plt.savefig('Analysis/TCR/normal_AIS_SE.pdf', bbox_inches='tight')
plt.show()
plt.close()

n=[];t=[]
for sample in sub:
        if (sub[sample]=="AIS/MIA") and (sample in mn):
                if (sample in top_clone) and (mn[sample] in n_top_clone):
                        t.append(float(top_clone[sample]))
                        n.append(float(n_top_clone[mn[sample]]))
                        print sample, float(top_clone[sample]), float(n_top_clone[mn[sample]])
nn=np.array(n)
tt=np.array(t)
a_t, a_prob = stats.ttest_rel(nn, tt)
print "top clone ais", a_t, a_prob, len(nn), len(tt), np.mean(nn), np.mean(tt)

n=[];t=[]
for sample in sub:
        if (sub[sample]=="LUAD") and (sample in mn):
                if (sample in top_clone) and (mn[sample] in n_top_clone):
                        t.append(float(top_clone[sample]))
                        n.append(float(n_top_clone[mn[sample]]))
                        print sample, float(top_clone[sample]), float(n_top_clone[mn[sample]])

nn=np.array(n)
tt=np.array(t)
a_t, a_prob = stats.ttest_rel(nn, tt)
print "top clone adc", a_t, a_prob, len(nn), len(tt), np.mean(nn), np.mean(tt)



