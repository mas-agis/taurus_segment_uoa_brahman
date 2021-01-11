# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 11:20:42 2021

@author: NUWI_352019
"""


import pandas as pd

#read alignment table
align=pd.read_csv(r'F:\maulana\MetaData\alignment_metadata.txt',delimiter="\t",header=0,usecols=[0,1,2,3,4])
#read used sample table
sample=pd.read_excel(r'F:\maulana\adjusted_dataset\adjusted_dataset.xlsx', sheet_name='Sheet3',usecols=[0,1,3] )
#read SRA list taurus
list_taurus=pd.read_excel(r"D:\maulana\Cattle\Archive\SRA_List_From_Ben_Nov'19\Taurus.xlsx", usecols=[0,2,21,22,23,24] )
#read SRA list indicus
list_indicus=pd.read_excel(r"D:\maulana\Cattle\Archive\SRA_List_From_Ben_Nov'19\Indicus.xlsx", usecols=[0,16,2,17,3,25] )
#concatened sra list taurus and indicus
frames=[list_taurus,list_indicus]
result=pd.concat(frames, ignore_index=True)
#separate align by REF columns
align.REF.unique()
ars=align[align["REF"]=="Btau5.0.1Y.fa"]
uoa=align[align["REF"]=="_Brahman_1.fa"]
#combine tables of ARS and UOA by pd.merge/join
comb=pd.merge(ars, uoa, how='outer', on='SRA', suffixes=("_ars", "_uoa"))
#inner join combination table with result table
comb=pd.merge(comb, result, how='inner', left_on='SRA', right_on='Run')
#inner join sample table with combination table
final=pd.merge(sample, comb, how='inner', left_on='SRA', right_on='SRA')

#write to csv
final.to_excel(r'F:\maulana\adjusted_dataset\adjusted_dataset_full.xlsx')
