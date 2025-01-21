#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:49:48 2022

@author: xuanzhang
"""


# draw rmsd
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import glob

def changeFontSize(ax, size):
    ax.title.set_fontsize(size)
    for item in ([ax.xaxis.label, ax.yaxis.label] +
                  ax.get_xticklabels() + ax.get_yticklabels() + 
                  ax.get_legend().get_texts()):
         item.set_fontsize(size)

#colors = sns.color_palette('cubehelix', 14)
colors = sns.color_palette("Paired", 14)
sns.palplot(colors)
sns.set_context('paper')

rundir = sys.argv[1]
output = sys.argv[2]
figfile = sys.argv[3]

rmsd_dat_file = glob.glob('{}/*/rmsd.xvg'.format(rundir))

rmsd_array = []
df_save = pd.DataFrame()

sel_variants = []
for dat in rmsd_dat_file:
    variant = dat.split("_")[0].replace('./', '')

    sel_variants.append(variant)
    times = []
    distance = []
    with open(dat, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            line = line.lstrip()
            items = line.split()
            if int(float(items[0])*1000) % 200 == 0:
                times.append(float(items[0]))
                distance.append(float(items[1]))
    rmsd_dict = dict()
    rmsd_dict['time'] = times
    rmsd_dict['distance'] = distance
    df = pd.DataFrame.from_dict(rmsd_dict, orient='columns')
    df4save = df.loc[:, ['time', 'distance']]
    df4save.columns = ['time', variant]
    df['variant'] = variant
    rmsd_array.append(df)
    if len(df_save):
        df_save = pd.merge(left=df_save, right=df4save, how='inner', left_on='time', right_on='time')
    else:
        df_save = df4save
df_merge = pd.concat(rmsd_array, axis=0)
df_merge.index = list(range(df_merge.shape[0]))

df_save.to_excel(output)
df_merge.to_csv(figfile+'.csv')

Counter(df_merge.distance)
Counter(df_merge.time)

sel_variants_sort = ['WT', 'B.1.1.529', 'BA.2', 'XBB.1.5', 'BA.2.86', 'JN.1', 'KP.2', 'KP.3', 'KP.3.1.1']
#sel_variants_sort = ['KP.3', 'KP.3.1.1', 'KP.3.R1', 'KP.3.R2', 'KP.3.1.1.R2']
#sel_variants_sort = ['WT', 'B.1.1.529', 'B.1.1.529.R1', 'BA.2', 'KP.3.1.1']
#sel_variants_sort = sel_variants
# 窄型数据
#plt.figure(figsize=(3.4, 1.7))
plt.figure(figsize=(16, 10))
ax = sns.lineplot(x='time', y='distance', hue='variant', #style=True,
             data=df_merge, alpha=.6, linewidth=2, markers='+', \
                 markersize=2.5, palette=colors, hue_order=sel_variants_sort)
#plt.ylim(-1, 20)
#plt.xlim(-1, 51)
#plt.yticks([0, 5, 10, 15, 20], ['0', '5', '10', '15', '20'])
#plt.xticks([0, 10, 20, 30, 40, 50], ['0', '10', '20', '30', '40', '50'])
plt.ylabel("RMSD(nm)")
plt.xlabel("Time(ns)")

#sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, 1), ncol=2, title=None, frameon=False,)
ax.legend().set_title(None)

sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))
changeFontSize(plt.gca(), 20)
plt.savefig(figfile, format='pdf', \
            bbox_inches='tight')
#plt.savefig("../gromacs_ace2-rbd_50ns_md/merge_statistics/merge_traj_rmsd_lineplot.png", format='png', \
#            bbox_inches='tight', dpi=600)
#plt.show()
plt.close()

