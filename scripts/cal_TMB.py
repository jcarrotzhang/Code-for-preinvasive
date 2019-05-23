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
import patsy
import scipy.spatial
from scipy.stats import mannwhitneyu
sns.set(font_scale=2)
sns.set_style("white")

basecov = 30
maf=open("m2/consensus_wStrelka.lite.maf")
samplelist=open("Analysis/partial_clinical_info_to_plot.txt")
purity=open("Analysis/GD/absolute_GD.txt") ### changed to absolute called GD and purity. ###
types=["Intron", "IGR", "RNA", "Silent"]

TMB={}
AIS_sample=[];AIS_tmb=[];ais_smoking={}
MIA_sample=[];MIA_tmb=[];mia_smoking={}
INV_sample=[];INV_tmb=[];INV_tmb_r=[];INV_tmb_n=[];inv_smoking={}
missense={}
AIS_tmb2=[]
MIA_tmb2=[]
INV_tmb2=[]

Recurrent={}
recurrent=open("Analysis/samplelist_Recurrence_status_new.txt")
for line in (raw.strip().split('\t') for raw in recurrent):
        if "Patient_ID" not in line:
                sample=line[2]
                recurrent=line[5]
                Recurrent[sample]=recurrent
                
for line in (raw.strip().split('\t') for raw in maf):
        for line in (raw.strip().split('\t') for raw in maf):
                if ("#" not in line[0]) and ("Hugo_Symbol" not in line[0]):
                        sample=line[9]
                        if line[4] not in types:
                                if sample in TMB:
                                        TMB[sample]=TMB[sample]+1
                                else:
                                        TMB[sample]=1
  
for line in (raw.strip().split('\t') for raw in samplelist):
       sample=line[0]
       subtype=line[2]
       #if ("Never" in line[5]) and ("CHG021156" not in sample):
       if "III" not in line[2]:
               if "AIS" in line[2]:
                       AIS_sample.append(sample)
                       ais_smoking[sample]="0"
               elif "MIA" in line[2]:
                       MIA_sample.append(sample)
                       mia_smoking[sample]="0"
               else:
                       INV_sample.append(sample)
                       inv_smoking[sample]="0"
               if sample not in TMB:
                        TMB[sample]=0

ais_purity={}
mia_purity={}
inv_purity={}
output=open("Analysis/TMB/sample_TMB.txt", 'w')
for line in (raw.strip().split('\t') for raw in purity):
        sample=line[0]
        if "sample" not in line[0]:
                if sample in AIS_sample:
                        ais_purity[sample]=line[3]
                elif sample in MIA_sample:
                        mia_purity[sample]=line[3]
                else:
                        inv_purity[sample]=line[3]
                        
print >>output, "group\tDNA-LC\tTMB\tmutation"
AIS_p=[];MIA_p=[];INV_p=[];ais_sm=[];mia_sm=[];inv_sm=[]
ais_sam=[]
mia_sam=[]
inv_sam=[]
for k, v in sorted(TMB.iteritems(), key=lambda (k,v): (v,k)):
        if int(v) > 0:
                tmb2=math.log(int(v)/basecov, 2)
                tmb=int(v)/basecov
        else:
                tmb2=float('NaN')
                tmb=0
                print tmb2

        if k in AIS_sample:
                
                if k in ais_purity:
                        ais_sam.append(k)
                        AIS_tmb.append(tmb)
                        AIS_tmb2.append(tmb2)
                        AIS_p.append(ais_purity[k])
                        ais_sm.append(ais_smoking[k])
                    
                        out="AIS"+"\t"+str(k)+"\t"+str(tmb)+"\t"+str(v)
                        print >>output, out
        if k in MIA_sample:
                
                if k in mia_purity:
                        mia_sam.append(k)
                        MIA_p.append(mia_purity[k])
                        mia_sm.append(mia_smoking[k])
                        MIA_tmb.append(tmb)
                        MIA_tmb2.append(tmb2)
                        
                        out="MIA"+"\t"+str(k)+"\t"+str(tmb)+"\t"+str(v)
                        print >>output, out
        if k in INV_sample:
                
                if k in inv_purity:
                        inv_sam.append(k)
                        INV_tmb.append(tmb)
                        INV_tmb2.append(tmb2)
                        INV_p.append(inv_purity[k])
                        inv_sm.append(inv_smoking[k])
                        #out="INV"+"\t"+str(k)+"\t"+str(tmb)+"\t"+str(v)+"\t"+str(missense[k])
                        out="INV"+"\t"+str(k)+"\t"+str(tmb)+"\t"+str(v)
                        print >>output, out

        if (k not in ais_purity) and (k not in mia_purity) and (k not in inv_purity):
                print "purity", k

### add linear regression. ###
AIS_tmb=np.array(AIS_tmb).astype(float)
MIA_tmb=np.array(MIA_tmb).astype(float)
INV_tmb=np.array(INV_tmb).astype(float)
AIS_tmb2=np.array(AIS_tmb2).astype(float)
MIA_tmb2=np.array(MIA_tmb2).astype(float)
INV_tmb2=np.array(INV_tmb2).astype(float)

AIS_p=np.array(AIS_p).astype(float)
MIA_p=np.array(MIA_p).astype(float)
INV_p=np.array(INV_p).astype(float)
ais_sam=np.array(ais_sam)
mia_sam=np.array(mia_sam)
inv_sam=np.array(inv_sam)
Ais=np.zeros(len(AIS_tmb), dtype=int)
Mia=np.ones(len(MIA_tmb), dtype=int)
Inv=np.full(len(INV_tmb), 2, dtype=int)

snm=np.concatenate((ais_sam, mia_sam, inv_sam), axis=0)
mu=np.concatenate((Ais, Mia, Inv), axis=0)
mb=np.concatenate((AIS_tmb, MIA_tmb, INV_tmb), axis=0)
mb2=np.concatenate((AIS_tmb2, MIA_tmb2, INV_tmb2), axis=0)
pu=np.concatenate((AIS_p, MIA_p, INV_p), axis=0)
A=np.array(["AIS"]*len(AIS_tmb))
B=np.array(["MIA"]*len(MIA_tmb))
C=np.array(["LUAD"]*len(INV_tmb))
t=np.concatenate((A, B, C), axis=0)
df = pd.DataFrame({"A": mu, "B": mb, "C":pu, "D":snm, "E":t, "B2": mb2})

df.to_csv("Analysis/TMB/plot.source.txt", sep="\t", index=False)
