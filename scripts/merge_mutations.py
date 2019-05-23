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

choice = str(sys.argv[1])

genes_mutation=["MET", "RET", "RAF1", "VAV1", "SOS1", "KRAS", "RIT1", "EGFR", "ARAF", "BRAF", "ERBB2", "HRAS", "KIT", "NRAS", "RASA1", "RBM10", "MGA", "RB1", "SMARCA4", "STK11", "TP53", "MAP2K1"]
genes_CNA=["EGFR", "ARAF", "KRAS", "MAP2K1", "FGFR1", "FGFR2", "FGFR3", "RIT1", "MYC", "TERT", "MDM2", "CDK4", "SFTA3"]

maf=open("m2/consensus_wStrelka_lite.wheader.maf")
#AIS_MUT=open("/cga/meyerson/Data/FUSCC_LUAD/m2/after_exclusion/consensus_AISMIA_with_header.lite.maf")
#INV_MUT=open("/cga/meyerson/Data/FUSCC_LUAD/m2/after_exclusion/consensus_invasive_with_header.lite.maf")
CNA_1=open("/cga/meyerson/Data/FUSCC_LUAD/zyCNA3/AISMIA_segfiles/focal_data_by_genes.txt")
CNA_2=open("/cga/meyerson/Data/FUSCC_LUAD/zyCNA3/invasive_segfiles/focal_data_by_genes.txt")
Fusion=open("/cga/meyerson/Data/FUSCC_LUAD/RNA_seq/fusion_tmp/Fusions_without_normal_tissue.txt")
AIS_samplelist=open("/cga/meyerson/Data/FUSCC_LUAD/samplelist_AISMIA_with_exclusion.txt")
INV_samplelist=open("/cga/meyerson/Data/FUSCC_LUAD/samplelist_invasive_with_exclusion.txt")
AIS_OUT=open("Analysis/mutation/AIS_merged_alteration.txt", 'w')
INV_OUT=open("Analysis/mutation/INV_merged_alteration.txt", 'w')
clinical_info=open("Analysis/partial_clinical_info.txt")

AIS_sample=[]
INV_sample=[]
for line in (raw.strip().split('\t') for raw in AIS_samplelist):
        sample=line[0]
        AIS_sample.append(sample)

for line in (raw.strip().split('\t') for raw in INV_samplelist):
        sample=line[0]
        INV_sample.append(sample)

mutation_count_AIS={}; mutation_count_INV={}
CNA_count_AIS={}; CNA_count_INV={}
Fusion_count_AIS={}; Fusion_count_INV={}

clin_info={}
for line in (raw.strip().split('\t') for raw in clinical_info):
        sample=line[3]
        clin_info[sample]="\t".join(line[5:len(line)-1])

print >>AIS_OUT, "gene\tchr\tstart\tend\tvariant_class\tVariant_type\tref\talt\talt_2\tsample\tNormal_Sample_Barcode\tProtein_change\tGender\tAge\tSmoking_status\tPathologic_stage"
print >>INV_OUT, "gene\tchr\tstart\tend\tvariant_class\tVariant_type\tref\talt\talt_2\tsample\tNormal_Sample_Barcode\tProtein_change\tGender\tAge\tSmoking_status\tPathologic_stage"
for line in (raw.strip().split('\t') for raw in maf):
        if "Hugo_Symbol" not in line[0]:
                gene=line[0]
                sample=line[9]
                if ("Splice_Site" in line[4]) or ("Missense_Mutation" in line[4]) or ("Nonsense_Mutation" in line[4]) or ("Frame_Shift_Ins" in line[4]) or ("Frame_Shift_Del" in line[4]) or ("In_Frame_Del" in line[4]) or ("In_Frame_Ins" in line[4]):
                        if ("EGFR" in gene) and (len(line) == 12):
                                if ("L62R" in line[11]) or ("S768I" in line[11]) or ("R776C" in line[11]) or ("T263I" in line[11]) or ("p.965_971VVDADEY>D" in line[11]):
                                        next
                                else:
                                        if sample in AIS_sample:
                                                if gene in mutation_count_AIS:
                                                        mutation_count_AIS[gene] = mutation_count_AIS[gene] + 1
                                                else:
                                                        mutation_count_AIS[gene] = 1
                                                print >>AIS_OUT, "\t".join(line[0:11])+"\t"+clin_info[sample]
                                        if sample in INV_sample:
                                                if gene in mutation_count_INV:
                                                        mutation_count_INV[gene] = mutation_count_INV[gene] + 1
                                                else:
                                                         mutation_count_INV[gene] = 1
                                                print >>INV_OUT, "\t".join(line[0:11])+"\t"+clin_info[sample]
                        elif sample in AIS_sample:
                                if gene in mutation_count_AIS:
                                        mutation_count_AIS[gene] = mutation_count_AIS[gene] + 1
                                else:
                                        mutation_count_AIS[gene] = 1
                                print >>AIS_OUT, "\t".join(line[0:11])+"\t"+clin_info[sample]
                        elif sample in INV_sample:
                                if gene in mutation_count_INV:
                                        mutation_count_INV[gene] = mutation_count_INV[gene] + 1
                                else:
                                        mutation_count_INV[gene] = 1
                                print >>INV_OUT, "\t".join(line[0:11])+"\t"+clin_info[sample]
                        else:
                                next

samples_to_exclude=["CHG019582", "CHG021148", "CHG020069", "CHG020900", "CHG019599", "CHG021188", "CHG020136", "CHG021865", "CHG019582", "CHG021148", "CHG020069", "CHG020900", "CHG019599", "CHG021188", "CHG021188", "CHG020136", "CHG021865", "CHG019582", "CHG021148", "CHG020069", "CHG020900", "CHG019599", "CHG021188", "CHG020136", "CHG021865", "CHG019582", "CHG021148", "CHG020069", "CHG020900", "CHG019599", "CHG021188", "CHG020136", "CHG021865", "D4912LC", "D4894LC"]

sampleName=[]
for line in (raw.strip().split('\t') for raw in CNA_1):
        if "Gene Symbol" in line[0]:
                for i in line[3:len(line)]:
                        sampleName.append(i)
        else:
                gene=line[0]
                if gene in genes_CNA:
                        for i in range(0, len(sampleName)):
                                if float(line[i+3]) > 1:
                                        if (sampleName[i] in AIS_sample) and (sampleName[i] not in samples_to_exclude):
                                                if gene+"_CNA" in mutation_count_AIS:
                                                        mutation_count_AIS[gene+"_CNA"] = mutation_count_AIS[gene+"_CNA"] + 1
                                                else:
                                                        mutation_count_AIS[gene+"_CNA"] = 1
                                                print >>AIS_OUT, gene+"_CNA\t1\t0\t0\tCNA\tCNA\tN\tN\tN\t"+sampleName[i]+"\t"+sampleName[i]+"\t"+line[i+3]+"\t"+clin_info[sampleName[i]]

sampleName=[]
for line in (raw.strip().split('\t') for raw in CNA_2):
        if "Gene Symbol" in line[0]:
                for i in line[3:len(line)]:
                        sampleName.append(i)
        else:
                gene=line[0]
                if gene in genes_CNA:
                        for i in range(0, len(sampleName)):
                                if float(line[i+3]) > 1:
                                        if (sampleName[i] in INV_sample) and (sampleName[i] not in samples_to_exclude):
                                                if gene+"_CNA" in mutation_count_INV:
                                                        mutation_count_INV[gene+"_CNA"] = mutation_count_INV[gene+"_CNA"] + 1
                                                else:
                                                        mutation_count_INV[gene+"_CNA"] = 1
                                                print >>INV_OUT, gene+"_CNA\t1\t0\t0\tCNA\tCNA\tN\tN\tN\t"+sampleName[i]+"\t"+sampleName[i]+"\t"+line[i+3]+"\t"+clin_info[sampleName[i]]
for line in (raw.strip().split('\t') for raw in Fusion):
        gene=line[1]
        sample=line[0]
        if sample in AIS_sample:
                if gene+"_Fusion" in mutation_count_AIS:
                        mutation_count_AIS[gene+"_Fusion"] = mutation_count_AIS[gene+"_Fusion"] + 1
                else:
                        mutation_count_AIS[gene+"_Fusion"] = 1
                print >>AIS_OUT, gene+"_Fusion\t1\t0\t0\tFusion\tFusion\tN\tN\tN\t"+sample+"\t"+sample+"\t"+line[2]+"\t"+clin_info[sample]
        if sample in INV_sample:
                if gene+"_Fusion" in mutation_count_INV:
                        mutation_count_INV[gene+"_Fusion"] = mutation_count_INV[gene+"_Fusion"] + 1
                else:
                        mutation_count_INV[gene+"_Fusion"] = 1
                print >>INV_OUT, gene+"_Fusion\t1\t0\t0\tFusion\tFusion\tN\tN\tN\t"+sample+"\t"+sample+"\t"+line[2]+"\t"+clin_info[sample]

o=[]; p=[]
OUT1=open("Analysis/mutation/sig_dif_AIS_INV_l_all_gene.csv", 'w')
OUT2=open("Analysis/mutation/sig_dif_AIS_INV_h_all_gene.csv", 'w')
print >>OUT1, "gene\tAIS_count\tINV_count\tdiff_freq_AIS\tdiff_freq_INV\toddsratio\tlog10_pvalue\tpvalue\tdiff_frequency"
print >>OUT2, "gene\tAIS_count\tINV_count\tdiff_freq_AIS\tdiff_freq_INV\toddsratio\tlog10_pvalue\tpvalue\tdiff_frequency"
sns.set(font_scale=2)
sns.set_style("white")

for i in mutation_count_INV:
        AIS_size=(len(AIS_sample)-1)/2
        INV_size=len(INV_sample)/2

        if i in mutation_count_AIS:
                wt_count_AIS=AIS_size-int(mutation_count_AIS[i])
                wt_count_INV=INV_size-int(mutation_count_INV[i])
                diff_freq_AIS=int(mutation_count_AIS[i])/AIS_size
                diff_freq_INV=int(mutation_count_INV[i])/INV_size
                diff_freq=diff_freq_INV-diff_freq_AIS
        else:
                mutation_count_AIS[i]=0
                wt_count_AIS=AIS_size
                wt_count_INV=INV_size-int(mutation_count_INV[i])
                diff_freq_AIS=0
                diff_freq_INV=int(mutation_count_INV[i])/INV_size
                diff_freq=diff_freq_INV-diff_freq_AIS

        oddsratio, pvalue  = stats.fisher_exact([[int(mutation_count_INV[i]), wt_count_INV], [int(mutation_count_AIS[i]), wt_count_AIS]])
        o.append(oddsratio)
        o.append(oddsratio)
        p.append(-math.log(pvalue, 10))

        if "yes" in choice:
                if (i in genes_mutation) or ("_CNA" in i) or ("_Fusion" in i):
                        if diff_freq < 0:
                                output=str(i)+"\t"+str(mutation_count_AIS[i])+"\t"+str(mutation_count_INV[i])+"\t"+str(diff_freq_AIS)+"\t"+str(diff_freq_INV)+"\t"+str(oddsratio)+"\t"+str(-math.log(pvalue, 10))+"\t"+str(pvalue)+"\t"+str(diff_freq)
                                print >>OUT1, output
                        else:
                                output=str(i)+"\t"+str(mutation_count_AIS[i])+"\t"+str(mutation_count_INV[i])+"\t"+str(diff_freq_AIS)+"\t"+str(diff_freq_INV)+"\t"+str(oddsratio)+"\t"+str(-math.log(pvalue, 10))+"\t"+str(pvalue)+"\t"+str(diff_freq)
                                print >>OUT2, output
        else:
                if diff_freq < 0:
                        output=str(i)+"\t"+str(mutation_count_AIS[i])+"\t"+str(mutation_count_INV[i])+"\t"+str(diff_freq_AIS)+"\t"+str(diff_freq_INV)+"\t"+str(oddsratio)+"\t"+str(-math.log(pvalue, 10))+"\t"+str(pvalue)+"\t"+str(diff_freq)
                        print >>OUT1, output
                else:
                        output=str(i)+"\t"+str(mutation_count_AIS[i])+"\t"+str(mutation_count_INV[i])+"\t"+str(diff_freq_AIS)+"\t"+str(diff_freq_INV)+"\t"+str(oddsratio)+"\t"+str(-math.log(pvalue, 10))+"\t"+str(pvalue)+"\t"+str(diff_freq)
                        print >>OUT2, output

for i in mutation_count_AIS:
        AIS_size=(len(AIS_sample)-1)/2
        INV_size=len(INV_sample)/2

        if i not in mutation_count_INV:
                mutation_count_INV[i]=0
                wt_count_INV=INV_size
                wt_count_AIS=AIS_size-int(mutation_count_AIS[i])
                diff_freq_INV=0
                diff_freq_AIS=int(mutation_count_AIS[i])/AIS_size
                diff_freq=diff_freq_INV-diff_freq_AIS
                oddsratio, pvalue  = stats.fisher_exact([[int(mutation_count_INV[i]), wt_count_INV], [int(mutation_count_AIS[i]), wt_count_AIS]])
                o.append(oddsratio)
                p.append(-math.log(pvalue, 10))
                output=str(i)+"\t"+str(mutation_count_AIS[i])+"\t"+str(mutation_count_INV[i])+"\t"+str(diff_freq_AIS)+"\t"+str(diff_freq_INV)+"\t"+str(oddsratio)+"\t"+str(-math.log(pvalue, 10))+"\t"+str(pvalue)+"\t"+str(diff_freq)
                print >>OUT1, output
OUT1.close()
OUT2.close()
