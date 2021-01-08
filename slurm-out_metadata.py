# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import re
import glob
os.chdir("/home/naji/Downloads") #specified folder with slurm-log files
files = glob.glob('slurm*.out') #listing slurm files within directory 

for i,k in enumerate(files):
    with open(k, 'r') as file:
        slurm_out = file.read().replace('\n', '')
    individuRegex=re.compile(r'\bSRR\d{7}') #creating regex pattern of SRA(boundary SRR followed by 7 digits)
    sample=individuRegex.search(slurm_out) #serch regex pattern in the slurm-out
    print("individu in this file is " + sample.group())
    refseqRegex=re.compile(r'UOA_Brahman_1.fa | ARS_UCD1.2.fa ') #creating regex pattern of either brahman or ars refseq
    refseq=refseqRegex.search(slurm_out)
    print("the refseq is " + refseq.group())
    total_readsRegex=re.compile(r'(Processed) (\d+) (total reads)') #brackets automatically group the pattern
    total_reads=total_readsRegex.search(slurm_out).group(2) #search and take the second group of matching strings
    print("total reads are: " + total_reads)
    filtered_readsRegex=re.compile(r'(BaseRecalibrator -) (\d+) (read\(s\) filtered by)')
    filtered_reads=filtered_readsRegex.search(slurm_out).group(2)
    print("filtered reads are: " + filtered_reads)
##!!Notes: apparently total_reads above are the clean reads only, thus filtered reads must be added

