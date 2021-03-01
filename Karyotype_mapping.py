# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 10:05:24 2021

@author: NUWI_352019
"""

import os
import matplotlib
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

def karyoplot(karyo_filename, metadata={}, part=1):
	'''
	To create a karyo_filename go to: http://genome.ucsc.edu/cgi-bin/hgTables 
	group: Mapping and Sequencing
	track: Chromosome Band 

	An example of an output (hg19, Human) is here: http://pastebin.com/6nBX6sdE 

	The script will plot dots next to loci defined in metadata as:
	metadata = {
		'1' : [2300000, 125000000, 249250621],
	}

	'''


	karyo_dict={}
	with open(karyo_filename) as karyo_f:
		lines = [x.replace(os.linesep, '').split() for x in karyo_f.readlines()]

		for chromosome in [str(x) for x in range(1,23)] + ['X', 'Y']:
			karyo_dict[chromosome] = [[y[0], int(y[1]), int(y[2]), y[3], y[4]] for y in [x for x in lines if x[0] == 'chr' + chromosome]]

	fig, ax = plt.subplots()

	DIM = 1.0

	ax.set_xlim([0.0, DIM * (1.3)])
	ax.set_ylim([0.0, DIM])

	def get_chromosome_length(chromosome):
		chromosome_start = float(min([x[1] for x in karyo_dict[chromosome]]))
		chromosome_end = float(max(x[2] for x in karyo_dict[chromosome]))
		chromosome_length = chromosome_end - chromosome_start

		return chromosome_length

	def plot_chromosome(chromosome, order):

		chromosome_length = get_chromosome_length(chromosome)
		chromosome_length_1 = get_chromosome_length('1')

		x_start = order * DIM * 0.1 
		x_end = x_start + (DIM * 0.04)
		y_start = DIM * 0.8 * (chromosome_length/chromosome_length_1)
		y_end = DIM * 0.1


		# We use the same colors as: http://circos.ca/tutorials/lessons/2d_tracks/connectors/configuration 
		colors = {
			'gpos100' : (0/255.0,0/255.0,0/255.0),
			'gpos'    : (0/255.0,0/255.0,0/255.0),
			'gpos75'  : (130/255.0,130/255.0,130/255.0),
			'gpos66'  : (160/255.0,160/255.0,160/255.0),
			'gpos50'  : (200/255.0,200/255.0,200/255.0),
			'gpos33'  : (210/255.0,210/255.0,210/255.0),
			'gpos25'  : (200/255.0,200/255.0,200/255.0),
			'gvar'    : (220/255.0,220/255.0,220/255.0),
			'gneg'    : (255/255.0,255/255.0,255/255.0),
			'acen'    : (217/255.0,47/255.0,39/255.0),
			'stalk'   : (100/255.0,127/255.0,164/255.0),
		}

		for index, piece in enumerate(karyo_dict[chromosome]):

			current_height = piece[2] - piece[1]
			current_height_sc = ((y_end - y_start) / chromosome_length) * current_height
			if index == 0:
				y_previous = y_start

			y_next = y_previous + current_height_sc

			color = colors[piece[4]]

			#plot the caryotypes
			r = Rectangle((x_start, y_previous), x_end-x_start, current_height_sc, color = color)
			ax.add_patch(r)

			y_previous = y_next

		#Plot semicircles at the beginning and end of the chromosomes
		center_x = x_start + (x_end-x_start)/2.0
		radius = (x_end-x_start)/2.0
		theta1 = 0.0
		theta2 = 180.0
		w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
		w2 = Wedge((center_x, y_end), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
		ax.add_patch(w1)
		ax.add_patch(w2)
		ax.plot([x_start, x_start], [y_start, y_end], ls='-', color='black')
		ax.plot([x_end, x_end], [y_start, y_end], ls='-', color='black')

		#Plot metadata
		if chromosome in metadata:
			for md in metadata[chromosome]:
				ax.plot([x_end + (DIM*0.015)], [y_start + (y_end-y_start) * (md/chromosome_length)], '.', color='black')

		ax.text(center_x, y_end - (DIM * 0.07), chromosome)

	if part==1:
		plot_chromosome('1', 1)
		plot_chromosome('2', 2)
		plot_chromosome('3', 3)
		plot_chromosome('4', 4)
		plot_chromosome('5', 5)
		plot_chromosome('6', 6)
		plot_chromosome('7', 7)
		plot_chromosome('8', 8)
		plot_chromosome('9', 9)
		plot_chromosome('10', 10)
		plot_chromosome('11', 11)
		plot_chromosome('12', 12)


	elif part==2:
		plot_chromosome('13', 1)
		plot_chromosome('14', 2)
		plot_chromosome('15', 3)
		plot_chromosome('16', 4)
		plot_chromosome('17', 5)
		plot_chromosome('18', 6)
		plot_chromosome('19', 7)
		plot_chromosome('20', 8)
		plot_chromosome('21', 9)
		plot_chromosome('22', 10)
		plot_chromosome('X', 11)
		plot_chromosome('Y', 12)

	else:
		raise Exception('plot argument should be either "1" or "2"')

	plt.axis('off')
	plt.show()

if __name__ == '__main__':
	import urllib
	fn = 'karyotype_hg19.txt'
	url = 'http://pastebin.com/raw.php?i=6nBX6sdE'
	if not os.path.exists(fn):
		print 'Downloading %s to local file: %s' % (url, fn)
		with open(fn, 'w') as k_file:
			f = urllib.urlopen(url)
			k_file.write(f.read())

	print 'plotting..'
	karyoplot(fn, part=1)
	karyoplot(fn, part=2)


from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas

def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print (chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']

# Height of each ideogram
chrom_height = 1

# Spacing between consecutive ideograms
chrom_spacing = 1

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 0.4

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.1

# Width, height (in inches)
figsize = (6, 8)

# Decide which chromosomes to use
chromosome_list = ['Chr%s' % i for i in range(1, 30) ]

# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase - gene_height - gene_padding
    ybase += chrom_height + chrom_spacing

# Read in ideogram.txt, downloaded from UCSC Table Browser
ideo = pandas.read_table(
    'ideogram.txt',
    skiprows=1,
    names=['chrom', 'start', 'end', 'name', 'gieStain']
)

# Filter out chromosomes not in our list
ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

# Add a new column for width
ideo['width'] = ideo.end - ideo.start

# Colors for different chromosome stains
color_lookup = {
    'gneg': (1., 1., 1.),
    'gpos25': (.6, .6, .6),
    'gpos50': (.4, .4, .4),
    'gpos75': (.2, .2, .2),
    'gpos100': (0., 0., 0.),
    'acen': (.8, .4, .4),
    'gvar': (.8, .8, .8),
    'stalk': (.9, .9, .9),
}

# Add a new column for colors
ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])


# Same thing for genes
genes = pandas.read_table(
    'ucsc_genes.txt',
    names=['chrom', 'start', 'end', 'name'],
    usecols=range(4))
genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
genes['width'] = genes.end - genes.start
genes['colors'] = '#2243a8'


fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
    ax.add_collection(collection)

# ...and the gene data
print("adding genes...")
for collection in chromosome_collections(
    genes, gene_ybase, gene_height, alpha=0.5, linewidths=0
):
    ax.add_collection(collection)

# Axes tweaking
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list)
ax.axis('tight')
plt.show()