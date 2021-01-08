# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import re
import glob
import pandas as pd
#os.chdir("/home/naji/Downloads") #specified folder with slurm-log files for laptop
os.chdir("F:\maulana\MetaData") #specified folder with slurm-log files for desktop
files = glob.glob('slurm*.out') #listing slurm files within directory 
a={'SRA': [], 'REF': [], 'clean_reads':[], 'filtered_reads':[], 'total_reads':[]}
df=pd.DataFrame(a)   
for i,k in enumerate(files):
    with open(k, 'r') as file:
        slurm_out = file.read().replace('\n', '')
    #finding individual SRA    
    individuRegex=re.compile(r'\bSRR\d{7}') #creating regex pattern of SRA(boundary SRR followed by 7 digits)
    sample=individuRegex.search(slurm_out) #serch regex pattern in the slurm-out
    print("individu in this file is " + sample.group())
    df.loc[i,"SRA"]=sample.group()
    #finding clean reads
    clean_readsRegex=re.compile(r'(Processed) (\d+) (total reads)') #brackets automatically group the pattern
    clean_reads=clean_readsRegex.search(slurm_out).group(2) #search and take the second group of matching strings
    print("clean reads are: " + clean_reads)
    df.loc[i,"clean_reads"]=int(clean_reads)
    #finding filtered reads
    filtered_readsRegex=re.compile(r'(BaseRecalibrator -) (\d+) (read\(s\) filtered by)')
    filtered_reads=filtered_readsRegex.search(slurm_out).group(2)
    print("filtered reads are: " + filtered_reads)
    df.loc[i,"filtered_reads"]=int(filtered_reads)
    #finding refseq works!
    refseqRegex=re.compile(r'.{10}\.fa ') # | ARS_UCD1.2.fa  creating regex pattern of either brahman or ars refseq
    refseq=refseqRegex.search(slurm_out)
    print("the refseq is " + refseq.group())
    df.loc[i,"REF"]=refseq.group()
    #trial
    #refseqRegex=re.compile(r'UOA') #creating regex pattern of either brahman or ars refseq
    #refseq=refseqRegex.search(slurm_out)
    #print("the refseq is " + refseq.group(1))
    #df.loc[i,"REF"]=refseq.group()    
    #adding total reads
    df.loc[i,"total_reads"]=df.loc[i,"filtered_reads"]+df.loc[i,"clean_reads"]
##!!Notes: apparently total_reads above are the clean reads only, thus filtered reads must be added

