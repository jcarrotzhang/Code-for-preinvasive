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
AIS=[];MIA=[];RNA_AIS=[];RNA_normal=[];RNA_rec=[];RNA_norec=[]
subclonal_samples=[]; clonal_samples=[]; sub_samples=[]; c_samples=[]
hlaloh=open("Analysis/HLA/HLA_analysis_output_2.txt")
batch_info=open("batch_info.txt")

TMB={}
tmb=open("Analysis/TMB/sample_TMB_all.txt")
for line in (raw.strip().split("\t") for raw in tmb):
        TMB[line[1]]=line[2]
        
to_exclude_bach=["CHG019733", "CHG019540","CHG019681","CHG018352","CHG019764","CHG019781","CHG020285","CHG021801","CHG021772","CHG021050","CHG019679","CHG020845", "CHG019539", "CHG019682", "CHG018351", "CHG019765", "CHG019780", "CHG020284", "CHG021802", "CHG021771", "CHG021049", "CHG019680", "CHG020844"]
for line in (raw.strip().split('\t') for raw in batch_info):
        if line[4] == "3":
                to_exclude_bach.append(line[0])

rna={}
RNAlist=open("Analysis/samplelist_Recurrence_status.txt")
for line in (raw.strip().split('\t') for raw in RNAlist):
        RNA=line[0]
        DNA=line[2]
        rna[RNA]=DNA
        rna[DNA]=RNA
        rna[line[1]]=line[3]

        if RNA not in to_exclude_bach:
                if "NA" in line[5]:
                        RNA_AIS.append(RNA)
                        RNA_normal.append(line[1])
                if "1" in line[5]:
                        RNA_rec.append(RNA)
                        RNA_normal.append(line[1])
                if "0" in line[5]:
                        RNA_norec.append(RNA)
                        RNA_normal.append(line[1])

Purity=open("Analysis/GD/absolute_GD.txt")
puri={}
LOH={}
for line in (raw.strip().split('\t') for raw in Purity):
        if "sample" not in line[0]:
                sample=line[0]
                purity = line[3]
                puri[rna[sample]]=purity
                if sample in loh:
                        LOH[rna[sample]]="HLALOH"
                else:
                        LOH[rna[sample]]="HLAWT"

mutation={}
maf=open("m2/consensus_wStrelka.lite.maf")
for line in (raw.strip().split('\t') for raw in maf):
        if ("#" not in line[0]) and ("Hugo_Symbol" not in line[0]):
                sample=line[9]
                gene=line[0]
                if "EGFR" in gene:
                        if sample in rna:
                                mutation[rna[sample]]=gene

C={}
rna_count=open("MIXCR/count_RNA_all.txt")
for line in (raw.strip().split('\t') for raw in rna_count):
        sample=line[0]
        count=line[1]
        #if sample in rna:
        C[sample]=count
### find clone count per type. ###
def find_clone(clone_count):
        a=[]
        type_count={}
        for sample in clone_count:
                count=clone_count[sample].split(";")
                a.append(len(count))
                for i in range(0, len(count)):
                        key = sample+"_"+str(i)
                        type_count[key]=count[i]
        num = np.array(a).max()
        return type_count, num

def creat_array(type_count, num, samplearray):
        array={}
        for i in range(0, num):
                for s in samplearray:
                        key=s+"_"+str(i)
                        if i in array:
                                array[i]=array[i]+";"+type_count[key]
                        else:
                                array[i]=type_count[key]
        return array

def add_type(type_count, num,  samplearray):
        for i in range(0, num):
                for s in samplearray:
                        k = s+"_"+str(i)
                        if k not in type_count:
                                type_count[k]="0"
        return type_count

def total_count(clone_count):
        c_tcr={} ### dictionary for count. ###

        for i in clone_count:
                if len(clone_count[i]) > 0:
                        clone=clone_count[i].split(";")
                        if (i in puri):
                                if (float(puri[i]) < 0.8) and (float(puri[i]) > 0.2):
                                        clone=np.array(clone).astype(int)
                                        #RPM = clone.sum()/float(C[i])*1000000
                                        #j = [n for n in clone if n >= 5]
                                        #print j

                                        #if j >0:
                                        #p = 1-float(puri[i])
                                        #fr = clone.sum()/p
                                        #RPM = fr/float(C[i])*1000000
                                                #t_RPM=RPM/(1-float(puri[i]))   ### not just add RPM but also normalize with tumor purity. ###
                                        RPM = clone.sum()/float(C[i])*1000000
                                        c_tcr[i]=RPM
                        else:
                                if i in RNA_normal:
                                        clone=np.array(clone).astype(int)
                                        RPM = clone.sum()/float(C[i])*1000000
                                        c_tcr[i]=RPM
                                #else:
                                #       c_tcr[i]=np.nan
                                        #clone=np.array(clone).astype(int)
                                        #RPM = clone.sum()/float(C[i])*1000000
                                        #c_tcr[i]=RPM
                        #else:
                        #       clone=np.array(clone).astype(int)
                        #       RPM = clone.sum()/float(C[i])*1000000
                        #       c_tcr[i]=RPM

                else:
                        print "nothing", i
                #       c_tcr[i]=0
                #if c_tcr[i] == 0:
                #       print i

        return c_tcr
def entropy(clone_count):
        SE_tcr={}  ### dictionary for entropy inex. ###
        for i in clone_count:
                clone=clone_count[i].split(";")
                if (len(clone) > 1):
                        #if (len(clone) > 9):
                        #       clone=np.array(clone).astype(int)[0:10]
                        #else:
                        clone=np.array(clone).astype(int)
                        RPM=[]
                        for c in clone:
                                if c > 4:
                                        #RPM.append(c/float(C[i])*1000000)
                                        RPM.append(c)
                                        #RPM = [c/float(C[i])*1000000 for c in clone]
                        if (len(RPM) > 9):
                                #if (np.array(RPM).astype(int)).sum() > 50:
                                        #if len(RPM) > 9:
                                        entropy=stats.entropy(RPM[0:10])
                                        #entropy=stats.entropy(RPM)
                                        Hn=entropy
                                        #else:
                                        #entropy=stats.entropy(RPM)  ### change to use the RPM value for SE. ### 
                                        #Hn=entropy/np.log(len(RPM))
                                        SE_tcr[i]=Hn

                                        if i in RNA_normal:
                                                print "normal", i, RPM, entropy, Hn
                                        if i in RNA_AIS:
                                                print "AIS", i, RPM, entropy, Hn
                                        if i in RNA_rec:
                                                print "ADC", i, RPM, entropy, Hn
                                        if i in RNA_norec:
                                                print "ADC", i, RPM, entropy, Hn
                        #else:
                        #       SE_tcr[i]=1
        return SE_tcr
def read_tcr(tcr):
        tcr_file=open(tcr)
        m=re.search('MIXCR/(.+?).clones.txt', tcr)
        sample=m.group(1)
        index=1
        Count=[]
        for line in (raw.strip().split('\t') for raw in tcr_file):
                if "cloneId" not in line[0]:
                        if ("TRA" in line[5]) or ("TRB" in line[5]):
                                #if "TGCAGTGCTAGAGAGTCGACTAGCGATCCAAAAAATGAGCAGTTCTTC" not in line[3]:
                                if int(line[1]) > 0:
                                        Count.append(str(line[1]))

        Count=";".join(Count)
        return Count

### run analysis for AISMIA and ADC group. ####
TCR=open("Analysis/TCR/samplelist")
ind=0
samplearray_AISMIA=[]; samplearray_rec=[]; samplearray_norec=[]; samplearray_normal=[]
clone_count_AISMIA={}; clone_count_rec={}; clone_count_norec={}; clone_count_normal={}

for line in (raw.strip().split('\t') for raw in TCR):
        ind=ind+1
        tcr=line[0]
        m=re.search('MIXCR/(.+?).clones.txt', line[0])
        sample=m.group(1)
        if (sample in RNA_AIS):
                samplearray_AISMIA.append(sample)
                clone_count_AISMIA[sample]=read_tcr(tcr)

        if (sample in RNA_rec):
                samplearray_rec.append(sample)
                clone_count_rec[sample]=read_tcr(tcr)

        if (sample in RNA_norec):
                samplearray_norec.append(sample)
                clone_count_norec[sample]=read_tcr(tcr)

        if (sample in RNA_normal):
                samplearray_normal.append(sample)
                clone_count_normal[sample]=read_tcr(tcr)

### perform enstropy test and output results. ###
c_tcr_AISMIA = total_count(clone_count_AISMIA)
c_tcr_norec = total_count(clone_count_norec)
c_tcr_rec = total_count(clone_count_rec)
c_tcr_normal = total_count(clone_count_normal)

SE_tcr_AISMIA  = entropy(clone_count_AISMIA)
SE_tcr_norec = entropy(clone_count_norec)
SE_tcr_rec = entropy(clone_count_rec)  ###  entropy array, total count array, entropy dictionary, total count dictionary 
SE_tcr_normal = entropy(clone_count_normal)

type_count, num = find_clone(clone_count_AISMIA)
type_count = add_type(type_count, num, samplearray_AISMIA)
array = creat_array(type_count, num, samplearray_AISMIA)
output=open("Analysis/TCR/AISMIA_output_stats.csv", 'w')
out="type_id"+"\t"+"\t".join(samplearray_AISMIA)
print >>output, out
for i in array:
        out=str(i)+"\t"+"\t".join(array[i].split(";"))
        print >>output, out
output.close()
type_count, num = find_clone(clone_count_normal)
type_count = add_type(type_count, num, samplearray_normal)
array = creat_array(type_count, num, samplearray_normal)
output=open("Analysis/TCR/normal_output_stats.csv", 'w')
out="type_id"+"\t"+"\t".join(samplearray_normal)
print >>output, out
for i in array:
        out=str(i)+"\t"+"\t".join(array[i].split(";"))
        print >>output, out
output.close()

### TCR count and diversity in correlation with TMB and H index. ###

sn=[]
se_a=[]; se_b=[]; se_c=[]
for i in clone_count_AISMIA:
        se_a.append(len(clone_count_AISMIA[i].split(";"))) ### TCR clone type.

for i in clone_count_rec:
        se_b.append(len(clone_count_rec[i].split(";")))

for i in clone_count_norec:
        se_b.append(len(clone_count_norec[i].split(";")))

for i in clone_count_normal:
        se_c.append(len(clone_count_normal[i].split(";")))

output=open("Analysis/TCR/TCR_count_SE.txt", 'w')
print >>output, "RNA_ID\tsample\ttype\ttotal_TCR_count\tSE_H\ttotal_TCR_types\tpurity\tHLA"
for i in c_tcr_normal:
        if i in SE_tcr_normal:
                out=str(i)+"\t"+str(rna[i])+"\tNormal\t"+str(float(c_tcr_normal[i]))+"\t"+str(SE_tcr_normal[i])+"\t"+str(len(clone_count_normal[i].split(";")))+"\tNA\tNA"
                print >>output, out
        else:
                out=str(i)+"\t"+str(rna[i])+"\tNormal\t"+str(float(c_tcr_normal[i]))+"\t"+str("NA")+"\t"+str(len(clone_count_normal[i].split(";")))+"\tNA\tNA"
                print >>output, out


for i in c_tcr_AISMIA:
        if i in SE_tcr_AISMIA:
                        out=str(i)+"\t"+str(rna[i])+"\tAIS/MIA\t"+str(float(c_tcr_AISMIA[i]))+"\t"+str(SE_tcr_AISMIA[i])+"\t"+str(len(clone_count_AISMIA[i].split(";")))+"\t"+puri[i]+"\t"+LOH[i]
                        print >>output, out
        else:
        #       if c_tcr_norec[i] < 15:
                        out=str(i)+"\t"+str(rna[i])+"\tAIS/MIA\t"+str(float(c_tcr_AISMIA[i]))+"\t"+str("NA")+"\t"+str(len(clone_count_AISMIA[i].split(";")))+"\t"+puri[i]+"\t"+LOH[i]
                        print >>output, out


for i in c_tcr_norec:
        if i in SE_tcr_norec:
                        out=str(i)+"\t"+str(rna[i])+"\tADC noREC\t"+str(float(c_tcr_norec[i]))+"\t"+str(SE_tcr_norec[i])+"\t"+str(len(clone_count_norec[i].split(";")))+"\t"+puri[i]+"\t"+LOH[i]
                        print >>output, out
        else:
        #       if c_tcr_norec[i] < 15:
                        out=str(i)+"\t"+str(rna[i])+"\tADC noREC\t"+str(float(c_tcr_norec[i]))+"\t"+str("NA")+"\t"+str(len(clone_count_norec[i].split(";")))+"\t"+puri[i]+"\t"+LOH[i]
                        print >>output, out
for i in c_tcr_rec:
        if i in SE_tcr_rec:
                        out=str(i)+"\t"+str(rna[i])+"\tADC REC\t"+str(float(c_tcr_rec[i]))+"\t"+str(SE_tcr_rec[i])+"\t"+str(len(clone_count_rec[i].split(";")))+"\t"+puri[i]+"\t"+LOH[i]
                        print >>output, out
        else:
        #       if c_tcr_norec[i] < 15:
                        out=str(i)+"\t"+str(rna[i])+"\tADC REC\t"+str(float(c_tcr_rec[i]))+"\t"+str("NA")+"\t"+str(len(clone_count_rec[i].split(";")))+"\t"+puri[i]+"\t"+LOH[i]
                        print >>output, out


#se_a=np.array(se_a).astype(int)
#se_b=np.array(se_b).astype(int)
#se_c=np.array(se_c).astype(int)
#a_t, a_prob = mannwhitneyu(se_a, se_b)
#b_t, b_prob = mannwhitneyu(se_a, se_c)
#print "T cell clone types", a_t, a_prob, b_t, b_prob, np.mean(se_a), np.mean(se_b), np.mean(se_c)

output.close()

### test TCR richness and clonality ###
df1=pd.read_csv('Analysis/TCR/TCR_count_SE.txt', delimiter="\t")
df2=pd.read_csv('Analysis/TCR/TCGA_TCR_count_SE.txt', delimiter="\t")
normal = df1[df1["type"] == "Normal"]
AISMIA = df1[df1["type"] == "AIS/MIA"]
ADC_norec = df1[df1["type"] == "ADC noREC"]
ADC_rec = df1[df1["type"] == "ADC REC"]
ADC =  df1[(df1["type"] == "ADC noREC") | (df1["type"] == "ADC REC")]

#Htest_inv=np.concatenate((ADC_norec["SE_H"].dropna().values, ADC_rec["SE_H"].dropna().values), axis=0)
#total_inv=np.concatenate((ADC_norec["total_TCR_count"].dropna().values, ADC_rec["total_TCR_count"].dropna().values), axis=0)

h_t, h_prob_a = mannwhitneyu(AISMIA["SE_H"].dropna().values, ADC["SE_H"].dropna().values)
s_t, s_prob_a = mannwhitneyu(AISMIA["total_TCR_count"].values, ADC["total_TCR_count"].values)
print "total RPM AIS ADC", s_t, s_prob_a
print "Hn test AIS ADC", h_t, h_prob_a

h_t, h_prob_b = mannwhitneyu(normal["SE_H"].dropna().values, ADC["SE_H"].dropna().values)
s_t, s_prob_b = mannwhitneyu(normal["total_TCR_count"].dropna().values, ADC["total_TCR_count"].dropna().values)
print "Hn test normal AIS", h_t, h_prob_b
print "total RPM normal AIS", s_t, s_prob_b

h_t, h_prob_c = mannwhitneyu(ADC_norec["SE_H"].dropna().values, ADC_rec["SE_H"].dropna().values)
s_t, s_prob_c = mannwhitneyu(ADC_norec["total_TCR_count"].dropna().values, ADC_rec["total_TCR_count"].dropna().values)
print "Hn test noRec Rec", h_t, h_prob_c
print "total RPM noRec Rec", s_t, s_prob_c

loh=ADC[ADC["HLA"]=="HLALOH"]
neu=ADC[ADC["HLA"]=="HLAWT"]

h_t, h_prob_c = mannwhitneyu(loh["SE_H"].dropna().values, neu["SE_H"].dropna().values)
s_t, s_prob_c = mannwhitneyu(loh["total_TCR_count"].dropna().values, neu["total_TCR_count"].dropna().values)
print "Hn test LOH WT", h_t, h_prob_c
print "total RPM LOH WT", s_t, s_prob_c
