# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 10:32:07 2021

@author: NUWI_352019
"""

# Import packages
import os
import glob
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn import preprocessing

#Working directory
os.getcwd()
os.chdir(r"D:\maulana\Analysis\Analysis_freq") #folder must contain only vcftools-freq2 output files

#function for tabulation of frequency output
def fetchdata(ref="ARS",filter=0.5,scan=1000000):
    summary=pd.DataFrame(columns=("breed","chr","window","jumlah")) #empyty dataframe
    def chr_length(reference='ARS'):
        index=list(range(1,30))
        if reference=="ARS":
            length=[158534110,136231102,121005158,120000601,120089316,117806340,
                    110682743,113319770,105454467,103308737,106982474,87216183,
                    83472345,82403003,85007780,81013979,73167244,65820629,
                    63449741,71974595,69862954,60773035,52498615,62317253,
                    42350435,51992305,45612108,45940150,51098607]
            refseq=pd.Series(data=length,index=index)
        else:
            length=[156943404,135844732,120786946,119487020,120179969,117220384,
                    110433717,113109275,104466507,103156416,106190332,86744025,
                    83353338,82466886,84270218,80235084,72619036,62583388,64467025,
                    71581568,69837837,59995071,53510974,62414538,42605158,51310239,
                    45108716,45896848,51254507]
            refseq=pd.Series(data=length,index=index)
        return refseq
    refseq=chr_length(reference=ref)  
    def daftar(reference='ARS'): #list files in the current directory matching the pattern name
        list_of_files = glob.glob('*.txt')  # create the list of file
        list_of_files3 = glob.glob('freq1_sr*.txt') #create the list of file with UOA
        if reference=="ARS":
            z=list(set(list_of_files) - set(list_of_files3))
        else:
            z=list_of_files3
        return z  
    files=daftar(reference=ref)
    for nomer, file in enumerate(files):
        ras=re.search('freq1_(.+?).txt', file).group(1) #extract string between "freq_" and ".txt"
        data=pd.read_csv(file,delimiter="\t",header=None,usecols=[0,1,3,5])
        data.columns=["chr","pos","n_al","alt"]
        data=data.loc[data['alt'] > filter ] #filter alt allel more than 0.5
        chr=data.chr.unique().tolist()
        for i in chr:
            chr_size=refseq[i] #getting size of the current chr
            window=np.arange(1,chr_size,scan)#creating interval from 1 to length of chr1 by 1M space
            window=np.append(window,chr_size) #adding closing for the interval
            data1=data.loc[data['chr'] == i]
            data1["terbagi"]=pd.cut(data1.pos,bins=window,labels=False,include_lowest=True)
            dipisah=data1.groupby(by="terbagi").size() #counting how many variants within each scanning window
            ditambah=pd.DataFrame(dipisah)
            ditambah["window"]=ditambah.index
            ditambah["breed"]=ras
            ditambah["chr"]=i
            ditambah.columns=["jumlah","window","breed","chr"]
            ditambah = ditambah[['breed', 'chr', 'window', 'jumlah']]
            summary=summary.append(ditambah)
    summary.chr=summary["chr"].astype(int) #changing column from object to integer
    summary.window=summary["window"].astype(int) #changing column from object to integer
    summary.jumlah=summary["jumlah"].astype(int) #changing column from object to integer            
    return summary        

###For alignment to ARS_UCD1.2
coba=fetchdata(ref="ARS", filter=0.95)
list_all=list(coba.breed.unique())
list_taur=["jersey","simmental","holstein","angus"]
list_taur=["jersey","simmental","holstein","angus","hereford","shorthorn"]
list_zeb=list(set(list_all)-set(list_taur))
#function to group breeds into taurus/indicus
def race(x):
    if x in list_taur:
        y="taurus"
    else:
        y="indicus"
    return y
coba["rase"]=coba.breed.str[:].apply(race) #new column based on string of another column

#Lineplot with 29 subplots, hue on rase (fig.1)
g = sns.FacetGrid(data=coba, col="chr", col_wrap=6, height=2, hue="rase", 
                  hue_order=["indicus","taurus"], sharex=False)
g.map(sns.lineplot, "window", "jumlah", alpha=.7)
g.set_titles(col_template="chr{col_name}")#, row_template="{row_name}")
g.set_axis_labels(" ", " ") #suppress x&y-label in each of grid
#g.add_legend(title='',labels=['Bos indicus', 'Bos taurus'])
g.fig.text(x=0, y=0.5, 
           verticalalignment='center', #make sure it's aligned at center vertically
           s='NFAA sites', #this is the text in the ylabel
           size=12, #customize the fontsize if you will
           rotation=90) #vertical text - overall ylabel
g.fig.text(x=0.5, y=0, 
           horizontalalignment='center', #make sure it's aligned at center horizontally
           s='Genome position in Mb', #this is the text in the xlabel
           size=12) #overall xlabel

#transform number of NFAA sites to Z-score separately for indicus and taurus
indicus=coba.loc[coba['rase'] == "indicus" ]
taurus=coba.loc[coba['rase'] == "taurus" ]
indicus["z"] = stats.zscore(indicus["jumlah"], nan_policy="omit")
taurus["z"] = stats.zscore(taurus["jumlah"], nan_policy="omit")
combined = pd.concat([indicus, taurus], axis=0, join="outer")

#Lineplot with 29 subplots, hue on rase (fig.1)
g = sns.FacetGrid(data=combined, col="chr", col_wrap=6, height=2, hue="rase", 
                  hue_order=["indicus","taurus"], sharex=False)
g.map(sns.lineplot, "window", "z", alpha=.7)
g.set_titles(col_template="chr{col_name}")#, row_template="{row_name}")
g.set_axis_labels(" ", " ") #suppress x&y-label in each of grid
#g.add_legend(title='',labels=['Bos indicus', 'Bos taurus'])
g.fig.text(x=0, y=0.5, 
           verticalalignment='center', #make sure it's aligned at center vertically
           s='Z score of NFAA sites', #this is the text in the ylabel
           size=12, #customize the fontsize if you will
           rotation=90) #vertical text - overall ylabel
g.fig.text(x=0.5, y=0, 
           horizontalalignment='center', #make sure it's aligned at center horizontally
           s='Genome position in Mb', #this is the text in the xlabel
           size=12) #overall xlabel

#Boxplots with 29 subplots, hue on rase(supplementary fig.1 )
fig, axes = plt.subplots(5, 6,squeeze=True, sharey=True, figsize=(16,14))
u=list(range(1,30))
pos=axes.flatten() #transform axes format from [a,k] to only single number
for i, k in enumerate(u):
    data_temp=coba.loc[coba['chr'] == k]
    sns.boxplot(x="chr", y="jumlah",hue="rase",data=data_temp, linewidth=1, ax=pos[i],hue_order=["indicus","taurus"])
    pos[i].get_legend().remove()
    handles, labels = pos[i].get_legend_handles_labels()
    pos[i].set_title('chr' + str(k)) #set title for subplots
    pos[i].set_xticks([])
    pos[i].xaxis.label.set_visible(False)
    labels=['Bos indicus', 'Bos taurus']
    pos[i].set_ylabel("")
    #pos[i].tick_params(top='off', bottom='off', left='on', right='off', labelleft='on', labelbottom='off')
fig.legend(handles, labels, loc='center right')
#fig.text(0.5, 0.04, 'common X', ha='center')
fig.text(0.07, 0.5, 'NFAA sites', va='center', rotation='vertical')
fig.delaxes(pos[29]) #deleting subplot number 29


###For alignment to UOA_Brahman_1
test=fetchdata(ref="UOA", filter=0.95)
#test=fetchdata(ref="UOA", filter=0.95, scan=100000) #trial
listall=list(test.breed.unique())
listtaur=["sr_jersey","sr_simmental","sr_holstein","sr_angus"]
listtaur=["sr_jersey","sr_simmental","sr_holstein","sr_angus","sr_hereford",
          "sr_shorthorn"]
listzeb=list(set(listall)-set(listtaur))
#function to group breeds into taurus/indicus
def ras(x):
    if x in listtaur:
        y="taurus"
    else:
        y="indicus"
    return y
test["rase"]=test.breed.str[:].apply(ras) #new column based on string of another column

#transform number of NFAA sites to Z-score separately for indicus and taurus
indicus=test.loc[test['rase'] == "indicus" ]
taurus=test.loc[test['rase'] == "taurus" ]
indicus["z"] = stats.zscore(indicus["jumlah"], nan_policy="omit")
taurus["z"] = stats.zscore(taurus["jumlah"], nan_policy="omit")
combined1 = pd.concat([indicus, taurus], axis=0, join="outer")

#Lineplot with 29 subplots, hue on rase (fig.1)
g = sns.FacetGrid(data=combined1, col="chr", col_wrap=6, height=2, hue="rase", 
                  hue_order=["indicus","taurus"], sharex=False)
g.map(sns.lineplot, "window", "z", alpha=.7)
g.set_titles(col_template="chr{col_name}")#, row_template="{row_name}")
g.set_axis_labels(" ", " ") #suppress x&y-label in each of grid
#g.add_legend(title='',labels=['Bos indicus', 'Bos taurus'])
g.fig.text(x=0, y=0.5, 
           verticalalignment='center', #make sure it's aligned at center vertically
           s='Z score of NFAA sites', #this is the text in the ylabel
           size=12, #customize the fontsize if you will
           rotation=90) #vertical text - overall ylabel
g.fig.text(x=0.5, y=0, 
           horizontalalignment='center', #make sure it's aligned at center horizontally
           s='Genome position in Mb', #this is the text in the xlabel
           size=12) #overall xlabel

#Lineplot with 29 subplots, hue on rase (fig.2)
g = sns.FacetGrid(data=test, col="chr", height=2, hue="rase", 
                  hue_order=["indicus","taurus"],col_wrap=6, sharex=False)
g.map(sns.lineplot, "window", "jumlah", alpha=.7)
g.set_titles(col_template="chr{col_name}")#, row_template="{row_name}")
g.set_axis_labels(" ", " ")
g.add_legend(title='',labels=['Bos taurus indicus', 'Bos taurus taurus'])
g.fig.text(x=0, y=0.5, 
           verticalalignment='center', #make sure it's aligned at center vertically
           s='NFAA sites', #this is the text in the ylabel
           size=12, #customize the fontsize if you will
           rotation=90) #vertical text - overall ylabel
g.fig.text(x=0.5, y=0, 
           horizontalalignment='center', #make sure it's aligned at center horizontally
           s='Genome position in Mb', #this is the text in the xlabel
           size=12) #overall xlabel

#Boxplots with 29 subplots, hue on rase(supplementary fig.2)
fig, axes = plt.subplots(5, 6,squeeze=True, sharey=True, figsize=(16,14))
u=list(range(1,30))
pos=axes.flatten() #transform axes format from [a,k] to only single number
for i, k in enumerate(u):
    data_temp=test.loc[test['chr'] == k]
    sns.boxplot(x="chr", y="jumlah",hue="rase",data=data_temp, linewidth=1, ax=pos[i],hue_order=["indicus","taurus"])
    pos[i].get_legend().remove()
    handles, labels = pos[i].get_legend_handles_labels()
    pos[i].set_title('chr' + str(k)) #set title for subplots
    pos[i].set_xticks([])
    pos[i].xaxis.label.set_visible(False)
    labels=['Bos indicus', 'Bos taurus']
    pos[i].set_ylabel("")
    #pos[i].tick_params(top='off', bottom='off', left='on', right='off', labelleft='on', labelbottom='off')
fig.legend(handles, labels, loc='center right')
#fig.text(0.5, 0.04, 'common X', ha='center')
fig.text(0.07, 0.5, 'NFAA sites', va='center', rotation='vertical')
fig.delaxes(pos[29]) #deleting subplot number 29


###Delta values part
#pivot dataframe, using multiple index 'chr and window'
test_lagi=pd.pivot_table(test, values='jumlah', index=['chr', 'window'], columns=['rase'])
test_lagi.head()
#delta1 as different values of indicus-mean minus taurus sites
test_lagi = test_lagi.assign(delta1=lambda x: test_lagi.mean().taurus - x.taurus )
test_lagi.head()
#Plotting histogram (fig.18)
g=sns.distplot(test_lagi['delta1'], hist=True, kde=True, norm_hist=False,
             bins=int(180/5), color = 'blue', 
             hist_kws={'edgecolor':'black'}, 
             axlabel=("Number of corrected taurus SNPs (Delta)"))
g.set(yticks=[])
g.set(xlabel="Number of corrected taurus SNPs (Delta)", ylabel="Density function")

#plot using general mean and std of delta1
rata=test_lagi.mean().delta1
st=test_lagi.std().delta1
batas=st*1.5 + rata
batas1=rata-st*1.5

#plot delta1 value
fig, axes = plt.subplots(5, 6,squeeze=False, sharey=True, sharex=False, figsize=(16,14))
u=list(range(1,30))
pos=axes.flatten() #transform axes format from [a,k] to only single number
for i, k in enumerate(u):
    data_temp=test_lagi.loc[test_lagi.index.get_level_values(0).isin({k}),:].reset_index()
    sns.lineplot(x="window", y="delta1",data=data_temp, linewidth=1,ax=pos[i])
    #pos[i].set(ylim=(rata, None)) #set ylim
    pos[i].axhline(rata, ls='dotted', c="firebrick", label="mean")
    pos[i].axhspan(ymin=batas1, ymax=batas, color="lightgrey", label="1.5 sd")
    handles, labels = pos[i].get_legend_handles_labels()
    pos[i].set_title('Chr ' + str(k)) #set title for subplots
    #pos[i].xaxis.label.set_visible(False)
    #if i in [23,24,25,26,27,28]:
     #   pos[i].set_xticks([0,30,60,90,120,150])
    #else:
     #    pos[i].set_xticks([])
    pos[i].set_ylabel("")
    pos[i].set_xlabel("")
    #pos[i].tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='off')
#fig.legend(handles, labels, loc='center right')
#fig.text(0.5, 0.04, 'Genome position in Mb', ha='center')
fig.text(0.5, 0, 'Genome position in Mb', ha='center')
#fig.text(0.07, 0.5, 'Delta', va='center', rotation='vertical')
fig.text(0, 0.5, 'Delta', va='center', rotation='vertical')
fig.tight_layout()
fig.delaxes(pos[29]) #deleting subplot number 29

#filter dataset above 1.5 of standard deviation
lagilagi=test_lagi[test_lagi["delta1"]>batas].reset_index()
lagilagi=lagilagi.assign(start= lambda x: x.window*1000000)
lagilagi=lagilagi.assign(end= lambda x: x.window*1000000 + 1000000-1)
lagilagi = lagilagi[["chr","start","end"]]
lagilagi.tail()

#output regions to bed file
np.savetxt(r'F:\maulana\adjusted_dataset\region_of_interest.bed', lagilagi.values, 
           fmt='%s', delimiter="\t") #all windows
np.savetxt(r'F:\maulana\adjusted_dataset\sites_count.txt',
           lagilagi,fmt='%s',
           header='rase chr window start end indicus taurus delta', 
           delimiter="\t") 

#Compacting/segmenting dataset and counting the mean of taurus introgressed regions
df=test_lagi[test_lagi["delta1"]>batas].reset_index() #extracting regions passing treshold of batas
df=df.assign(start= lambda x: x.window*1000000)
df=df.assign(end= lambda x: x.window*1000000 + 1000000-1)
h=list(set(df.chr.tolist())) #extracting chr number in df.chr, take the unique number, and set in a list
df_segment=pd.DataFrame(columns=("chr","start","end","delta1", "size")) #empyty dataframe
for i in h:
    temp=df[df["chr"]==i] #subset df per chromosome inquired
    mask = temp['start'] != temp['end'].shift()+1 #checking whether current row of 'start' doesn't match previous row 'end'
    temp1=temp.groupby(mask.cumsum()).agg({'chr':'first','start':'first', 'end':'last', 'delta1':'sum'}) #extract four columns of 'temp' based on grouping of 'mask'
    temp1=temp1.assign(size= lambda x: (x.end+1 - x.start)/1000000) #assign 'size' column to the temp1
    df_segment=df_segment.append(temp1) #append each temporary chr dataframe to final dataframe 

#saving the whole segments with taurus introgression
np.savetxt(r'F:\maulana\adjusted_dataset\segmented_putative_taurus_regions.txt',
           df_segment,fmt='%s',
           header='Chr Start End Delta Size', 
           delimiter="\t")  

#Overlapped genes in adjusted dataset with reported by Koufarioutis
gene1 = pd.read_csv(r'F:\maulana\adjusted_dataset\coding_genes.txt', header=None, sep="\t").to_numpy()
gene3 = pd.read_csv(r'F:\maulana\adjusted_dataset\koufarioutis_genes_taurus.txt', header=None, sep="\t").to_numpy()
inter=np.intersect1d(gene1, gene3) #intersect annotated genes with taurus_genes reported by Koufarioutis
gene1_filtered=np.setdiff1d(gene1, inter) #takes only elements of 'gene1' that not found in 'inter'
np.savetxt(r'F:\maulana\adjusted_dataset\gene_intersect_taurus.txt',inter,fmt='%s')#save intersect
np.savetxt(r'F:\maulana\adjusted_dataset\gene_non_intersect_taurus.txt',gene1_filtered,fmt='%s')#save non_intersect
gene4=pd.read_csv(r'D:\maulana\Cattle\ResourceBundle\ncbi-genomes-2020-08-03\UOA_genes_full.txt',
                  header=None, sep="\t").to_numpy() #all genes in UOA list
gene4_filtered=np.setdiff1d(gene4, gene1)
np.savetxt(r'F:\maulana\adjusted_dataset\gene_putative_indicus.txt',gene4_filtered,fmt='%s')#save intersect


##Comparing neutral taurus genes to default genes GO functions
GO=pd.read_csv(r'F:\maulana\adjusted_dataset\gene_non_intersect_taurus_panther_chart.txt',
               delimiter="\t",header=None, usecols=[1,2])
GO=GO.rename(columns={1: "go", 2: "count"})
GO=GO.set_index('go')
GO_default=pd.read_csv(r'F:\maulana\adjusted_dataset\default_genes_pantherChart.txt',
               delimiter="\t",header=None, usecols=[1,2])
GO_default=GO_default.rename(columns={1: "go", 2: "count"})
GO_default=GO_default.set_index('go')
GO_comb=GO_default.join(GO, lsuffix="_def", rsuffix="_gen")
GO_comb=GO_comb.assign(def_per=lambda x: GO_comb.count_def/GO_comb.count_def.sum())
GO_comb=GO_comb.assign(gen_per=lambda x: GO_comb.count_gen/GO_comb.count_gen.sum())
GO_comb=GO_comb.assign(changes=lambda x: (GO_comb.gen_per-GO_comb.def_per)*100)
#Save to txt GO_comb sorted by changes
dat=GO_comb.sort_values(by=['changes'],ascending=False).reset_index()
np.savetxt(r'F:\maulana\adjusted_dataset\panther_comparison.txt',
           dat,fmt='%s',delimiter="\t")  
#subplotting for each GO term
x = GO_comb.index.to_numpy()
xi=range(len(x))
fig, ax = plt.subplots(figsize=(16,14))
for i,k in enumerate(xi):
    lab=x[i]
    y=GO_comb.changes[i]
    ax.bar(x=i, height=y, width=0.8, label=lab)
    #Attach a text label above each bar
    ax.annotate(x[i], xy=(i,y),xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')        
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Difference to default in percent')
ax.set_xlabel("GO_Slim Biological Process")
#ax.set.axhline(0,ls='dotted', c="firebrick")
ax.set_xticks([])
ax.set_xticklabels(labels)
ax.legend()
fig.tight_layout()
plt.show()  
#subplotting for only positive value
xi=GO_comb.index[GO_comb['changes'] > 0 ].tolist() #getting index of rows where 'changes' column is positive
fig, ax = plt.subplots(figsize=(16,14))
for i,k in enumerate(xi):
    lab=x[i]
    y=GO_comb.changes[k]
    ax.bar(x=i, height=y, width=0.8, label=lab)
    #Attach a text label above each bar
    #ax.annotate(k, xy=(i,y),xytext=(0, 3),  # 3 points vertical offset
     #               textcoords="offset points",
      #              ha='center', va='bottom')        
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Difference to default in percent')
#ax.set_xlabel("GO_Slim Biological Process")
#ax.set.axhline(0,ls='dotted', c="firebrick")
ax.set_xticks([])
ax.set_xticklabels(labels)
ax.legend(title="GO biological process")
fig.tight_layout()
plt.show()  



###individual breed plotting for showing effect of taurus introgressed in UOA_Brahman
###For alignment to UOA_Brahman_1
test=fetchdata(ref="UOA", filter=0.95)
#test=fetchdata(ref="UOA", filter=0.95, scan=100000) #trial
test.head()
#adding columns rase
listall=list(test.breed.unique())
listtaur=["sr_jersey","sr_simmental","sr_holstein","sr_angus"]
listtaur=["sr_jersey","sr_simmental","sr_holstein","sr_angus","sr_hereford",
          "sr_shorthorn"]
listzeb=list(set(listall)-set(listtaur))
#function to group breeds into taurus/indicus
def ras(x):
    if x in listtaur:
        y="taurus"
    else:
        y="indicus"
    return y
test["rase"]=test.breed.str[:].apply(ras) #new column based on string of another column
#putative taurus-introgressed regions 
intro_reg = pd.read_csv(r"D:\maulana\adjusted_dataset\draft\Old_additional_file_2.txt",
                        sep= "\t", usecols=[0,1,2])

#keeping only indicus and plot each breed separately
indicus = test.loc[test['rase'] == "indicus"]
for sp in indicus.breed.unique():
  
    subset = test.loc[test['breed'] == sp]
    #plot delta1 value
    fig, axes = plt.subplots(5, 6,squeeze=False, sharey=True, sharex=False, figsize=(16,14))
    u=list(range(1,30))
    pos=axes.flatten() #transform axes format from [a,k] to only single number
    for i, k in enumerate(u):
        data_temp = subset.loc[subset["chr"] == k]
        #data_temp=subset.loc[subset.index.get_level_values(0).isin({k}),:].reset_index()
        sns.lineplot(x="window", y="jumlah",data=data_temp, linewidth=1, ax=pos[i])
        #temporary data for intro_reg
        data_temp1 = intro_reg.loc[intro_reg["Chr"] == k]
        for l in range(len(data_temp1)):
            pos[i].axvspan(data_temp1.iloc[l,1], data_temp1.iloc[l,2], color="lightgreen")
            #pos[i].axhline(rata, ls='dotted', c="firebrick", label="mean")
            #   pos[i].axhspan(ymin=batas1, ymax=batas, color="lightgrey", label="1.5 sd")
            #handles, labels = pos[i].get_legend_handles_labels()
            pos[i].set_title('Chr ' + str(k)) #set title for subplots
            #pos[i].xaxis.label.set_visible(False)
            #if i in [23,24,25,26,27,28]:
                #   pos[i].set_xticks([0,30,60,90,120,150])
                #else:
                    #    pos[i].set_xticks([])
            pos[i].set_ylabel("")
            pos[i].set_xlabel("")
            #pos[i].tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='off')
    fig.legend(title = sp, loc='lower right', title_fontsize=20)
    #fig.text(0.5, 0.04, 'Genome position in Mb', ha='center')
    fig.text(0.5, 0, 'Genome position in Mb', ha='center')
    #fig.text(0.07, 0.5, 'Delta', va='center', rotation='vertical')
    fig.text(0, 0.5, 'NFAA sites', va='center', rotation='vertical')
    fig.tight_layout()
    fig.delaxes(pos[29]) #deleting subplot number 

#keeping only taurus  and plot each breed separately
taurus = test.loc[test['rase'] == "taurus"]
for sp in taurus.breed.unique():
    subset = test.loc[test['breed'] == sp]
    #plot delta1 value
    fig, axes = plt.subplots(5, 6,squeeze=False, sharey=True, sharex=False, figsize=(16,14))
    u=list(range(1,30))
    pos=axes.flatten() #transform axes format from [a,k] to only single number
    for i, k in enumerate(u):
        data_temp = subset.loc[subset["chr"] == k]
        #data_temp=subset.loc[subset.index.get_level_values(0).isin({k}),:].reset_index()
        sns.lineplot(x="window", y="jumlah",data=data_temp, linewidth=1, ax=pos[i])
        #temporary data for intro_reg
        data_temp1 = intro_reg.loc[intro_reg["Chr"] == k]
        for l in range(len(data_temp1)):
            pos[i].axvspan(data_temp1.iloc[l,1], data_temp1.iloc[l,2], color="lightblue")
            #pos[i].axhline(rata, ls='dotted', c="firebrick", label="mean")
            #   pos[i].axhspan(ymin=batas1, ymax=batas, color="lightgrey", label="1.5 sd")
            #handles, labels = pos[i].get_legend_handles_labels()
            pos[i].set_title('Chr ' + str(k)) #set title for subplots
            #pos[i].xaxis.label.set_visible(False)
            #if i in [23,24,25,26,27,28]:
                #   pos[i].set_xticks([0,30,60,90,120,150])
                #else:
                    #    pos[i].set_xticks([])
            pos[i].set_ylabel("")
            pos[i].set_xlabel("")
            #pos[i].tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='off')
    fig.legend(title = sp, loc='lower right', title_fontsize=20)
    #fig.text(0.5, 0.04, 'Genome position in Mb', ha='center')
    fig.text(0.5, 0, 'Genome position in Mb', ha='center')
    #fig.text(0.07, 0.5, 'Delta', va='center', rotation='vertical')
    fig.text(0, 0.5, 'NFAA sites', va='center', rotation='vertical')
    fig.tight_layout()
    fig.delaxes(pos[29]) #deleting subplot number 
    
    

#multibreeds in each plots - boran1, brahman2, gir3, indianzebu4, mangshi6, nelore7, 
breeds = [indicus.breed.unique()[1]]
breeds.append(indicus.breed.unique()[7]) #2,3,4,6,7
less_breeds = [breeds[2]]
less_breeds.append(breeds[1])
less_breeds.append(breeds[3])
    
    fig, axes = plt.subplots(5, 6,squeeze=False, sharey=True, sharex=False, figsize=(16,14))
    u=list(range(1,30))
    pos=axes.flatten() #transform axes format from [a,k] to only single number
    for i, k in enumerate(u):
        subset = test.loc[test['chr'] == k]
        for sp in less_breeds:
            data_temp = subset.loc[subset["breed"] == sp]
            mean_X = data_temp.jumlah.mean()
            data_temp = data_temp.assign(norm = lambda x: data_temp.jumlah - mean_X)
            #data_temp=subset.loc[subset.index.get_level_values(0).isin({k}),:].reset_index()
            sns.lineplot(x="window", y="norm",data=data_temp, linewidth=1, ax=pos[i])
            #temporary data for intro_reg
        data_temp1 = intro_reg.loc[intro_reg["Chr"] == k]
        for l in range(len(data_temp1)):
            pos[i].axvspan(data_temp1.iloc[l,1], data_temp1.iloc[l,2], color="lightgreen")
            #pos[i].axhline(rata, ls='dotted', c="firebrick", label="mean")
            #   pos[i].axhspan(ymin=batas1, ymax=batas, color="lightgrey", label="1.5 sd")
            handles, labels = pos[i].get_legend_handles_labels()
            pos[i].set_title('Chr ' + str(k)) #set title for subplots
            #pos[i].xaxis.label.set_visible(False)
            #if i in [23,24,25,26,27,28]:
                #   pos[i].set_xticks([0,30,60,90,120,150])
                #else:
                    #    pos[i].set_xticks([])
            pos[i].set_ylabel("")
            pos[i].set_xlabel("")
            #pos[i].tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='off')
    fig.legend(handles, labels, title = sp, loc='lower right', title_fontsize=20)
    #fig.text(0.5, 0.04, 'Genome position in Mb', ha='center')
    fig.text(0.5, 0, 'Genome position in Mb', ha='center')
    #fig.text(0.07, 0.5, 'Delta', va='center', rotation='vertical')
    fig.text(0, 0.5, 'NFAA sites', va='center', rotation='vertical')
    fig.tight_layout()
    fig.delaxes(pos[29]) #deleting subplot number 


    fig, axes = plt.subplots(5, 6,squeeze=False, sharey=True, sharex=False, figsize=(16,14))
    u=list(range(1,30))
    pos=axes.flatten() #transform axes format from [a,k] to only single number  
    test_subset = test[test["breed"].isin(less_breeds)]
    for i, k in enumerate(u):
        subset = test_subset.loc[test_subset['chr'] == k]
        #for sp in less_breeds:
         #   data_temp = subset.loc[subset["breed"] == sp]
          #  mean_X = data_temp.jumlah.mean()
           # data_temp = data_temp.assign(norm = lambda x: data_temp.jumlah - mean_X)
            #data_temp=subset.loc[subset.index.get_level_values(0).isin({k}),:].reset_index()
        sns.lineplot(x="window", y="jumlah", data=subset, hue= "breed", linewidth=1, ax=pos[i])
            #temporary data for intro_reg
        data_temp1 = intro_reg.loc[intro_reg["Chr"] == k]
        for l in range(len(data_temp1)):
            pos[i].axvspan(data_temp1.iloc[l,1], data_temp1.iloc[l,2], color="lightgreen")
            #pos[i].axhline(rata, ls='dotted', c="firebrick", label="mean")
            #   pos[i].axhspan(ymin=batas1, ymax=batas, color="lightgrey", label="1.5 sd")
            #handles, labels = pos[i].get_legend_handles_labels()
            pos[i].set_title('Chr ' + str(k)) #set title for subplots
            #pos[i].xaxis.label.set_visible(False)
            #if i in [23,24,25,26,27,28]:
                #   pos[i].set_xticks([0,30,60,90,120,150])
                #else:
                    #    pos[i].set_xticks([])
            pos[i].legend().remove()        
            pos[i].set_ylabel("")
            pos[i].set_xlabel("")
            #pos[i].tick_params(top='off', bottom='on', left='on', right='off', labelleft='on', labelbottom='off')
    fig.legend(handles, labels, title = sp, loc='lower right', title_fontsize=20)
    #fig.text(0.5, 0.04, 'Genome position in Mb', ha='center')
    fig.text(0.5, 0, 'Genome position in Mb', ha='center')
    #fig.text(0.07, 0.5, 'Delta', va='center', rotation='vertical')
    fig.text(0, 0.5, 'NFAA sites', va='center', rotation='vertical')
    fig.tight_layout()
    fig.delaxes(pos[29]) #deleting subplot number 

g = sns.FacetGrid(data=combined, col="chr", col_wrap=6, height=2, hue="rase", 
                  hue_order=["indicus","taurus"], sharex=False)
g.map(sns.lineplot, "window", "z", alpha=.7)
g.set_titles(col_template="chr{col_name}")#, row_template="{row_name}")
g.set_axis_labels(" ", " ") #suppress x&y-label in each of grid
#g.add_legend(title='',labels=['Bos indicus', 'Bos taurus'])
g.fig.text(x=0, y=0.5, 
           verticalalignment='center', #make sure it's aligned at center vertically
           s='Z score of NFAA sites', #this is the text in the ylabel
           size=12, #customize the fontsize if you will
           rotation=90) #vertical text - overall ylabel
g.fig.text(x=0.5, y=0, 
           horizontalalignment='center', #make sure it's aligned at center horizontally
           s='Genome position in Mb', #this is the text in the xlabel
           size=12) #overall xlabel

multi = pd.DataFrame({})
min_max_scaler = preprocessing.MinMaxScaler()
for breed in breeds:
    temp = test.loc[test['breed'] == breed]
    mean_X = temp.jumlah.mean()
    min_X = temp.jumlah.min()
    max_X = temp.jumlah.max()
    std_X = temp.jumlah.std()
    temp = temp.assign(norm = lambda x: (temp.jumlah - min_X) / (max_X - min_X))
    multi = pd.concat([multi, temp], axis = 0, join = "outer")
   
g = sns.FacetGrid(data=multi, col="chr", col_wrap=6, height=2, hue="breed", sharex=False)
g.map(sns.lineplot, "window", "jumlah", alpha=.7)
g.set_titles(col_template="chr{col_name}")#, row_template="{row_name}")
g.set_axis_labels(" ", " ") #suppress x&y-label in each of grid
#g.add_legend(title='',labels=['Bos indicus', 'Bos taurus'])
g.fig.text(x=0, y=0.5, 
           verticalalignment='center', #make sure it's aligned at center vertically
           s='NFAA sites', #this is the text in the ylabel
           size=12, #customize the fontsize if you will
           rotation=90) #vertical text - overall ylabel
g.fig.text(x=0.5, y=0, 
           horizontalalignment='center', #make sure it's aligned at center horizontally
           s='Genome position in Mb', #this is the text in the xlabel
           size=12) #overall xlabel        
temp_array_scaled = min_max_scaler.fit_transform(temp_array)
temp["z"] = temp_array_scaled
multi = pd.concat([multi, temp], axis = 0, join = "outer")
    
multi = test.loc[test['breed'] == "indicus"]
indicus["z"] = stats.zscore(indicus["jumlah"], nan_policy="omit")
taurus["z"] = stats.zscore(taurus["jumlah"], nan_policy="omit")
combined1 = pd.concat([indicus, taurus], axis=0, join="outer")

(df-df.min())/(df.max()-df.min())

#smoothing the graph using scipy.interpolate
from scipy.interpolate import interp1d

f_cubic = interp1d(x, y, kind= "cubic")

help(interp1d)
# Plot.
plt.plot(x, y, 'o', label='data')
plt.plot(xnew, f_linear(xnew), '-', label='linear')
plt.plot(xnew, f_cubic(xnew), '--', label='cubic')
plt.legend(loc='best')
plt.show()