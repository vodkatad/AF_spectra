import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.colors as colors
import matplotlib 
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.transforms as transforms
from matplotlib.patches import Rectangle

cnvs = pd.read_csv(snakemake.input.cn, sep="\t") 
boundaries = cnvs[['chr','b','e']].copy()
cnvs = cnvs.drop(columns=['chr','b','e']).transpose()
names = list(cnvs.index.values)

# we had code to sort based on model - clone ready, we now
# have an external ordering for m dictated by input.order
orderd = pd.read_csv(snakemake.input.order, sep="\t")
orderd['n'] = range(0, len(orderd))

models = [x.split('-')[0] for x in names]
clones = [x.split('-')[1] for x in names]
times = [x.split('-')[2] for x in names]
letter = [x.split('-')[3] if len(x.split('-'))==4 else '' for x in names]
p = pd.DataFrame(data={'model':models,'c':clones, 't': times, 'l':letter})
m = pd.merge(p, orderd, on="model", how="inner")
assert m.shape[0] == p.shape[0], 'Not all samples with CN have info on how to order'

m = m.sort_values(by=['n','c','t','l'])

sorted_names= []
for i in range(0, len(m)):
    if m.iloc[i]['l'] != "":
        sorted_names.append(str(m.iloc[i]['model']) +'-'+ str(m.iloc[i]['c'])+ '-' +str(m.iloc[i]['t']) +'-'+ str(m.iloc[i]['l']))
    else: 
        sorted_names.append(str(m.iloc[i]['model']) +'-'+ str(m.iloc[i]['c'])+ '-' +str(m.iloc[i]['t']))

cnvs = cnvs.reindex(sorted_names)
with open(snakemake.log.log, 'w') as f:
    f.write(str(sorted_names))
    f.write(str(cnvs.shape))

samples = m['model'].tolist()
# Create a categorical palette to identify the samples
#sample_pal = sns.color_palette("bright", len(np.unique(samples)))
sample_pal = sns.color_palette(["#cc3300","#f607b9","#9900ff","#155d00","#77a003","#0829fc","#ff9900","#ffcc33"])
#sample_pal = sns.husl_palette(len(samples), h=.5)
sample_lut = dict(zip(map(str, np.unique(samples)), sample_pal))
# Convert the palette to vectors that will be drawn on the side of the matrix
allsamples = cnvs.index.get_level_values(None)
sample_colors = pd.Series(samples, index=allsamples).map(sample_lut)

chr_limits = boundaries.index[boundaries['e'].isin(boundaries.groupby('chr', sort=False)['e'].max().values)].tolist()
chr_boundaries = np.append(0, chr_limits)
chr_list = boundaries['chr'].unique().tolist()
chrN_list = []

for x in chr_list:
    x = x[3:] #remove 'chr' for readability
    chrN_list.append(x)

#compute the position where chromosome labels will be placed on the plots
start = 0
pos_list = []
for end in chr_limits:
    pos_list.append((start+end)/2)
    start = end+1

yticklabels = False

cbar_kws={"ticks":np.arange(0,13,1)}
#import sys
#sys.setrecursionlimit(10000)
#h = sns.clustermap(cnvs, method=method, metric=metric, col_cluster=False, yticklabels = yticklabels,  cmap='RdBu_r', vmin=0, vmax=12,norm=divnorm, cbar_kws=cbar_kws)
h = sns.clustermap(cnvs, row_colors=sample_colors, col_cluster=False, row_cluster=False, yticklabels = True, cmap='RdBu_r', vmin=0, vmax=12,center=2, cbar_kws=cbar_kws, rasterized=True)
#Z = h.dendrogram_row.linkage
ax = h.ax_heatmap
for pos in chr_limits:
    ax.axvline(x=pos, color='black')

#place chromosome ticks at the right position
ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
ax.tick_params(axis='x', rotation=0, labelsize=5)
#ax.tick_params(axis='y', rotation=125, labelsize=15)

ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
ax.tick_params(axis='x', length=20, which='minor')

ax.set_xlabel("Chromosomes", fontweight='bold', fontsize=7)
ax.set_ylabel("Clones", fontsize=7, fontweight='bold')

# legend is not a leged but a ax_cbar
cax = plt.gcf().axes[-1]
cax.tick_params(labelsize=5) # nope

# A4 is 8-1/4 x 11-3/4 in
plt.gcf().set_size_inches(7, 3.8) # w, h
# cannot find a way to get tolerable linewidth cause linewidth parameter seem to be ignored, will
# scale by hand
plt.savefig(snakemake.output.plot, dpi=300)
plt.clf()
