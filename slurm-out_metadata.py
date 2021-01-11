# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import re
import glob
import pandas as pd
from datetime import datetime
#os.chdir("/home/naji/Downloads") #specified folder with slurm-log files for laptop
os.chdir("F:\maulana\MetaData") #specified folder with slurm-log files for desktop
files = glob.glob('slurm*.out') #listing slurm files within directory 
a={'SRA': [], 'REF': [], 'clean_reads':[], 'filtered_reads':[], 'total_reads':[], 'time':[]}
df=pd.DataFrame(a)   
for i,k in enumerate(files):
    with open(k, 'r') as file:
        slurm_out = file.read().replace('\n', '')
    if re.search(r'\bBaseRecalibrator\b',slurm_out): #if find string continue extract the data
        #finding individual SRA    
        individuRegex=re.compile(r'SRR\d{7}') #creating regex pattern of SRA(boundary SRR followed by 7 digits)
        sample=individuRegex.search(slurm_out) #serch regex pattern in the slurm-out
        df.loc[i,"SRA"]=str(sample.group())
        #finding clean reads
        clean_readsRegex=re.compile(r'(Processed) (\d+) (total reads)') #brackets automatically group the output pattern
        clean_reads=clean_readsRegex.search(slurm_out).group(2) #search and take the second group of matching strings
        df.loc[i,"clean_reads"]=int(clean_reads)
        #finding filtered reads
        filtered_readsRegex=re.compile(r'(BaseRecalibrator -) (\d+) (read\(s\) filtered by)')
        filtered_reads=filtered_readsRegex.search(slurm_out).group(2)
        df.loc[i,"filtered_reads"]=int(filtered_reads)
        #finding refseq works!
        refseqRegex=re.compile(r'.{10}\.fa') # | ARS_UCD1.2.fa  creating regex pattern of either brahman or ars refseq
        refseq=refseqRegex.search(slurm_out)
        df.loc[i,"REF"]=str(refseq.group())
        #total_reads 
        df.loc[i,"total_reads"]=df.loc[i,"filtered_reads"]+df.loc[i,"clean_reads"]
        #finding time stamp of finished baserecalibrator; square brackets must be escaped by backslash
        time=re.search('(engine\[)(.*) CES*T\] org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator', slurm_out)
        datetime_object=datetime.strptime(time.group(2), '%B %d, %Y %H:%M:%S %p' ) #details on https://strftime.org/       
        df.loc[i,"time"]=datetime_object
    else:
        print("nggak ketemu ee di" + str(k))

#sort df by columns 'SRA','REF','time'; drop rows containing same 'SRA','REF'; keeping the last row by its 'time'
fe=df.sort_values(by=['SRA','REF','time']).drop_duplicates(['SRA','REF'],keep='last') #this is ideal

#write to table
fe.to_csv(r'F:\maulana\MetaData\alignment_metadata.txt', sep="\t",index=False, doublequote=False)
