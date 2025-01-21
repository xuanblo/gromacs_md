#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:03:19 2022

@author: xuanzhang
"""

# draw rmsd
import re
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

mmpbsa_dat_file = glob.glob('{}/*/_gmxtools~MMPBSA.dat'.format(rundir))


mmpbsa_array = []
df_save = pd.DataFrame()
sel_variants = []

for dat in mmpbsa_dat_file:
    variant = dat.split("_")[0].replace('./', '')

    sel_variants.append(variant)
    times = []
    frames = []
    mmpbsa = []

    with open(dat, 'r') as f:
        for line in f:

            matObj = re.match(".*~(\d0)ns\s+(-?\S+)\(", line)
            if matObj:
                times.append(int(matObj.group(1)))
                mmpbsa.append(round(float(matObj.group(2))/4.18, 2))
    frames = [int(i/10) for i in times]
    mmpbsa_dict = dict()
    mmpbsa_dict['frames'] = frames
    mmpbsa_dict['mmpbsa'] = mmpbsa
    df = pd.DataFrame.from_dict(mmpbsa_dict, orient='columns')
    df['variant'] = variant
    mmpbsa_array.append(df)

    df4save = df.loc[:, ['frames', 'mmpbsa']]
    df4save.columns = ['frames', variant]
    if len(df_save):
        df_save = pd.merge(left=df_save, right=df4save, how='inner', left_on='frames', right_on='frames')
    else:
        df_save = df4save

df_save.to_excel(output)

print(sel_variants)
df_merge = pd.concat(mmpbsa_array, axis=0)
df_merge.index = list(range(df_merge.shape[0]))
df_merge.to_csv(figfile+'.csv')

Counter(df_merge.frames)
Counter(df_merge.mmpbsa)
#sel_variants_sort = ['WT', 'BA.1', 'BA.2.86', 'JN.1', 'KP.2.3', 'KP.2', 'KP.3', 'EG.5.1', 'XBB.1.5', 'BA.2', 'BA.5']
sel_variants_sort = ['WT', 'B.1.1.529', 'BA.2', 'XBB.1.5', 'BA.2.86', 'JN.1', 'KP.2', 'KP.3', 'KP.3.1.1']
sel_variants_sort = ['KP.3', 'KP.3.1.1', 'KP.3.R1', 'KP.3.R2', 'KP.3.1.1.R2']
# 窄型数据
# plt.figure(figsize=(3.4, 1.7))
plt.figure(figsize=(16, 10))
#sns.lineplot(x='time', y='num', hue='variant', palette=color,\
#             data=df_merge, alpha=.6, linewidth=1)
ax = sns.barplot(x='frames', y='mmpbsa', hue='variant',\
             data=df_merge, alpha=1, linewidth=2, palette=colors, hue_order=sel_variants_sort)

plt.ylabel("MM/PBSA(kcal/mol)")
plt.xlabel("Frames")
# Put the legend out of the figure
#sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, 1), ncol=2, title=None, frameon=False,)
ax.legend().set_title(None)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))

changeFontSize(plt.gca(), 20)
plt.savefig(figfile, format='pdf', \
            bbox_inches='tight')

plt.show()
plt.close()