# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 09:32:05 2021

@author: NUWI_352019
"""

import pandas as pd
import numpy as np
import re

#Reading segmented taurine 
segment=pd.read_csv(r'F:\maulana\adjusted_dataset\segmented_putative_taurus_regions.txt',
                       sep="\t", usecols=[3,4], skiprows=1, names=["Total_Delta", "Size"])

#Reading compacted_annotated_taurine_segments
BED=pd.read_csv(r'F:\maulana\adjusted_dataset\compacted_taurine_segments_annotated.bed',
                sep="\t", skiprows=3,usecols=[0,1,2,3], names=["Chr","Start","End","Effect"])

#Taking only coding protein genes from 'Effect' Column and assign it to 'kamus'
kamus={}
tick=0
for i in range(len(BED)):
    a=BED.iloc[i,3]
    mylist = re.findall( r'Gene:(.*?):protein_coding', a)
    myset = list(set(mylist))
    myset= [i for i in myset if 'U6' not in i] #removing U6 snRNA
    myset= [i for i in myset if 'ENSBIX' not in i] #removing ENSBIX transcript
    kamus[tick.real]=myset
    tick+=1

#Counting how many genes are in the taurine regions
i=0
for key, valu in kamus.items():
    i+=len(valu)
i #there are 2226 genes 

#Overlapped with positive genes reported by Koufarioutis
positive_gen = pd.read_csv(r'F:\maulana\adjusted_dataset\koufarioutis_genes_taurus.txt', header=None, sep="\t").to_numpy()
clean_kamus={}
for k,v in kamus.items():
    #print("urutan ke", k, "gennya", v)
    clean_kamus[k] = np.intersect1d(v, positive_gen)
clean_kamus

#Counting how many genes are in taurine regions that overlapped with genes reported by Koufarioutis
i=0
for key, valu in clean_kamus.items():
    i+=len(valu)
i #there are 66 genes 

#Preparing columns with positive and neutral genes
df=pd.DataFrame({'positive_genes':pd.Series(clean_kamus),'genes':pd.Series(kamus)})
df["neutral_genes"]=" "
for i in range(100):
    a=df.iloc[i,0]
    b=df.iloc[i,1]
    c=np.setdiff1d(b,a)
    df.iloc[i,2]=c
df['positive_genes'] = [', '.join(map(str, l)) for l in df['positive_genes']] #Column of lists, convert list to string as a new column
df['neutral_genes'] = [', '.join(map(str, l)) for l in df['neutral_genes']] #Column of lists, convert list to string as a new column
del df['genes']

#Adding Columns with positive and neutral genes to taurine segments
del BED['Effect']
BED = BED.merge(segment, how="inner", on=BED.index)
del BED['key_0']
BED = BED.merge(df, how="inner", on=BED.index)
del BED['key_0']
BED.to_csv(r'F:\maulana\adjusted_dataset\compacted_taurine_segments_with_positive&neutral_genes.txt',
           index=False,sep=";")
