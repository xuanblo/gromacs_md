#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 08:59:40 2022

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

hbond_dat_file = glob.glob('{}/*/hbond_num.xvg'.format(rundir))


hbond_num_array = []
df_save = pd.DataFrame()

sel_variants = []
for dat in hbond_dat_file:
    variant = dat.split("_")[0].replace('./', '')

    sel_variants.append(variant)
    times = []
    hydrogen_bonds = []
    with open(dat, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            line = line.lstrip()
            items = line.split()
            if int(float(items[0])*1000) % 200 == 0:
                times.append(float(items[0]))
                hydrogen_bonds.append(int(items[1]))
    hbond_num_dict = dict()
    hbond_num_dict['time'] = times
    #hbond_num_dict[variant] = hydrogen_bonds
    hbond_num_dict['num'] = hydrogen_bonds
    df = pd.DataFrame.from_dict(hbond_num_dict, orient='columns')
    df['variant'] = variant
    hbond_num_array.append(df)

    df4save = df.loc[:, ['time', 'num']]
    df4save.columns = ['time', variant]
    if len(df_save):
        df_save = pd.merge(left=df_save, right=df4save, how='inner', left_on='time', right_on='time')
    else:
        df_save = df4save
df_save.to_excel(output)

df_merge = pd.concat(hbond_num_array, axis=0)
df_merge.index = list(range(df_merge.shape[0]))
df_merge.to_csv(figfile+'.csv')

Counter(df_merge.num)
Counter(df_merge.time)
#sel_variants_sort = ['WT', 'BA.1', 'BA.2.86', 'JN.1', 'KP.2.3', 'KP.2', 'KP.3', 'EG.5.1', 'XBB.1.5', 'BA.2', 'BA.5']
sel_variants_sort = ['WT', 'B.1.1.529', 'BA.2', 'XBB.1.5', 'BA.2.86', 'JN.1', 'KP.2', 'KP.3', 'KP.3.1.1']
# 窄型数据
# plt.figure(figsize=(3.4, 1.7))
plt.figure(figsize=(16, 10))
#sns.lineplot(x='time', y='num', hue='variant', palette=color,\
#             data=df_merge, alpha=.6, linewidth=1)
ax = sns.lineplot(x='time', y='num', hue='variant', #style='variant',\
             data=df_merge, alpha=.6, linewidth=2, palette=colors, hue_order=sel_variants_sort)
plt.ylim(-1, 20)
plt.xlim(-1, 51)
plt.yticks([0, 5, 10, 15, 20], ['0', '5', '10', '15', '20'])
plt.xticks([0, 10, 20, 30, 40, 50], ['0', '10', '20', '30', '40', '50'])
plt.ylabel("No. of H_bonds")
plt.xlabel("Time(ns)")

# Put the legend out of the figure
#sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, 1), ncol=2, title=None, frameon=False,)
ax.legend().set_title(None)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))

changeFontSize(plt.gca(), 20)

plt.savefig(figfile, format='pdf', \
            bbox_inches='tight')

plt.show()
plt.close()

