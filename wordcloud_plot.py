# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 10:09:08 2021

@author: NUWI_352019
"""
# Import packages
import re
import pandas as pd 
import matplotlib.pyplot as plt
from wordcloud import WordCloud, STOPWORDS# Generate word cloud

# Define a function to plot word cloud
def plot_cloud(wordcloud):
    # Set figure size
    plt.figure(figsize=(40, 30))
    # Display image
    plt.imshow(wordcloud) 
    # No axis details
    plt.axis("off");

#text wrangling, importing from panther overrepresentation non_intersect taurus 
mtext=pd.read_csv(r'F:\maulana\adjusted_dataset\gene_non_intersect_taurus_overrepresentation_test.txt',
               delimiter="\t",header=None, skiprows=(13))
plain=mtext[0].str.split("(", n = 1, expand=True) #split GO terms by '(' and expand to new column
plain=plain[0].to_numpy() #take first column to numpy
plain=set(plain) #take elemnent of numpy array
plain=' '.join([''.join(row) for row in plain]) #joining all elements with space
plain = re.sub('process', '', plain) #removing word 'process'
plain = re.sub('regulation', '', plain) #removing word 'regulation'
#plot the wrangled text
wordcloud = WordCloud(width = 3000, height = 2000, random_state=1, background_color='salmon', colormap='Pastel1', 
                      collocations=True, stopwords = STOPWORDS).generate(plain)# Plot
plot_cloud(wordcloud)

#text wrangling, importing from panther overrepresentation putative indicus
mtext=pd.read_csv(r'F:\maulana\adjusted_dataset\gene_putative_indicus_overrepresentation_test.txt',
               delimiter="\t",header=None, skiprows=(13))
plain=mtext[0].str.split("(", n = 1, expand=True) #split GO terms by '(' and expand to new column
plain=plain[0].to_numpy() #take first column to numpy
plain=set(plain) #take elemnent of numpy array
plain=' '.join([''.join(row) for row in plain]) #joining all elements with space
plain = re.sub('process', '', plain) #removing word 'process'
plain = re.sub('regulation', '', plain) #removing word 'regulation'
#plot the wrangled text
wordcloud = WordCloud(width = 3000, height = 2000, random_state=1, background_color='salmon', colormap='Pastel1', collocations=True, stopwords = STOPWORDS).generate(plain)# Plot
plot_cloud(wordcloud)

#finding keywords for manuscript
import docx2txt
my_text = docx2txt.process("/home/naji/Documents/Second_paper/Draft_Taurus_introgression_in_UOA_Bos_indicus_assembly_GM.docx")
print(my_text)
wordcloud = WordCloud(width = 3000, height = 2000, random_state=1, background_color='salmon', colormap='Pastel1', 
                      collocations=True, stopwords = STOPWORDS).generate(my_text)# Plot
plot_cloud(wordcloud)