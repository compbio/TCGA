#!/usr/bin/env python

# File: sample_info_table.py
# Name: HoJoon Lee
# Desc: Script produces the table show the avaiable data (snv, cnv, WES, RNA-seq) for each sample
# Created in: 12/09/2016

import re, os, sys, pandas as pd

tumor_type = '01' ### 01 -> primary tumor, 02 -> Recurrent Solid Tumor, 06 -> Metastatic   code table from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
cancer = 'STAD'

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Connetct the downlaoded folder names from GDC to TCGA ID and file size      #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
vcf= {}
tcga_id={}
infile1 = "./VCF_mutect/file_manifest.txt"
df1 = pd.read_table(infile1, sep="\t")
for index1,row in df1.iterrows():
    sampleId = row["TCGA Entity ID"].split("-")
    ID = sampleId[0] + '-' + sampleId[1] + '-' + sampleId[2]
    
    if ID in tcga_id:
        tcga_id[ID]+=1
    else:
        tcga_id[ID]=1
    
    if ID in vcf:
        vcf[ID].append(row["TCGA Entity ID"])
    else:
        vcf.setdefault(ID,[])
        vcf[ID].append(row["TCGA Entity ID"])

            
tcnv={}
ncnv={}
infile2 = "./CN/file_manifest.txt"
df2 = pd.read_table(infile2, sep="\t")
print df2.columns
for index2,row2 in df2.iterrows():
    sampleId2 = row2["TCGA Entity ID"].split("-")
    Id2 = sampleId2[0] + '-' + sampleId2[1] + '-' + sampleId2[2]
    
    if Id2 in tcga_id:
        tcga_id[Id2]+=1
    else:
        tcga_id[Id2]=1
    
    match2 = re.match("0",sampleId2[3],flags=0)
    if match2:
        if Id2 in tcnv:
            tcnv[Id2].append(row2["TCGA Entity ID"])
        else:
            tcnv.setdefault(Id2,[])
            tcnv[Id2].append(row2["TCGA Entity ID"])
    else:
        if Id2 in ncnv:
            ncnv[Id2].append(row2["TCGA Entity ID"])
        else:
            ncnv.setdefault(Id2,[])
            ncnv[Id2].append(row2["TCGA Entity ID"])

twes={}
nwes={}
infile3 = "/mnt/tcga_data/bam/" + cancer + "/WES/file_manifest.txt"
df3 = pd.read_table(infile3, sep="\t")
print df3.columns
for index3,row3 in df3.iterrows():
    sampleId3 = row3["TCGA Entity ID"].split("-")
    Id3 = sampleId3[0] + '-' + sampleId3[1] + '-' + sampleId3[2]
    
    if Id3 in tcga_id:
        tcga_id[Id3]+=1
    else:
        tcga_id[Id3]=1
    
    match3 = re.match("0",sampleId3[3],flags=0)
    if match3:
        if Id3 in twes:
            twes[Id3].append(row3["TCGA Entity ID"])
        else:
            twes.setdefault(Id3,[])
            twes[Id3].append(row3["TCGA Entity ID"])
    else:
        if Id3 in nwes:
            nwes[Id3].append(row3["TCGA Entity ID"])
        else:
            nwes.setdefault(Id3,[])
            nwes[Id3].append(row3["TCGA Entity ID"])

trna={}
nrna={}
infile4 = "/mnt/tcga_data/bam/" + cancer + "/RNA_seq/file_manifest.txt"
df4 = pd.read_table(infile4, sep="\t")
for index4,row4 in df4.iterrows():
    sampleId4 = row4["TCGA Entity ID"].split("-")
    Id4 = sampleId4[0] + '-' + sampleId4[1] + '-' + sampleId4[2]
    
    if Id4 in tcga_id:
        tcga_id[Id4]+=1
    else:
        tcga_id[Id4]=1
    
    match4 = re.match("0",sampleId4[3],flags=0)
    if match4:
        if Id4 in trna:
            trna[Id4].append(row4["TCGA Entity ID"])
        else:
            trna.setdefault(Id4,[])
            trna[Id4].append(row4["TCGA Entity ID"])
    else:
        if Id4 in nrna:
            nrna[Id4].append(row4["TCGA Entity ID"])
        else:
            nrna.setdefault(Id4,[])
            nrna[Id4].append(row4["TCGA Entity ID"])

batch_size = 100
batch_num = 1
batch_count= 1
fout = open('gdc_sample_status.txt', 'w')
fout.write('sample' + '\t' + '#Files' + '\t' + 'snv' + '\t' + 'cnv_t' + '\t' + 'wes_t' + '\t' + 'rnaseq_T' + '\t' + '#avaibale_T' + '\t' + 'cnv_N' + '\t' + 'wes_N' + '\t' + 'rnaseq_N' + '\t' + '#avaibale_N' + '\t' + 'batch_num' + '\n')
for key, value in tcga_id.iteritems():
    fout.write(key + '\t' + str(value) + '\t')
    tcount=0
    if key in vcf:
        vcfs = set(vcf[key])
        vcfList = ",".join(vcfs)
        fout.write(vcfList)
        tcount+=1
    fout.write('\t')
    
    if key in tcnv:
        cnvts = set(tcnv[key])
        cnvtList = ",".join(cnvts)
        fout.write(cnvtList)
        tcount+=1
    fout.write('\t')
    
    if key in twes:
        wests = set(twes[key])
        westList = ",".join(wests)
        fout.write(westList)
        tcount+=1
    fout.write('\t')
    
    if key in trna:
        rnats = set(trna[key])
        rnatList = ",".join(rnats)
        fout.write(rnatList)
        tcount+=1
    fout.write('\t' + str(tcount) + '\t')
    
    ncount=0
    if key in ncnv:
        cnvns = set(ncnv[key])
        cnvnList = ",".join(cnvns)
        fout.write(cnvnList)
        ncount+=1
    fout.write('\t')
    
    if key in nwes:
        wesns = set(nwes[key])
        wesnList = ",".join(wesns)
        fout.write(wesnList)
        ncount+=1
    fout.write('\t')
    
    if key in nrna:
        rnans = set(nrna[key])
        rnanList = ",".join(rnans)
        fout.write(rnanList)
        ncount+=1
    fout.write('\t' + str(ncount) + '\t')
    
    if(tcount ==4):
        fout.write(str(batch_num))
        if(batch_count % batch_size == 0):
            batch_num+=1
        batch_count+=1
    else:
        fout.write('nan')
    fout.write('\n')

