#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:01:52 2022

@author: xuanzhang
"""

# draw rmsd
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import re
import numpy as np
import glob

sup = re.compile(r"\(.*\)")

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
output1 = sys.argv[2]
output2 = sys.argv[3]
figfile = sys.argv[4]

mmpbsa_dat_file = glob.glob('{}/*/_gmxtools~res_MMPBSA.dat'.format(rundir))

hbond_num_array = []
df_save = pd.DataFrame()

def read_mmpbsa(sel_variants, rmsd_dat_file):

    mmpbsa_dict = dict()
    mmpbsa_df = pd.DataFrame()
    with open(rmsd_dat_file, 'r') as f:
        for line in f.readlines():
            if "------" in line or "~0ns" in line or "mean" in line:
                continue
            aa = line.split()
            
            if line.startswith("#"):
                header = []
                it = re.finditer(r"[L]~\S+\(\S+\s+\S+\)",line) 
                for i in it:
                    header.append(i.group())
                resn_sel = [sup.sub("", i).replace("L~", "") for i in header if i.startswith('L')]
                mmpbsa_dict['pos'] = [i[0:3] for i in resn_sel]
                mmpbsa_dict['res'] = [i[3:] for i in resn_sel]
            else:
                time = aa[0].split("~")[1]
                items = []
                it = re.finditer(r"\s+\S+\([^\(]+\)",line) 
                for i in it:
                    items.append(i.group())
                mmpbsa = [i.lstrip() for i in items[0: len(resn_sel)]]
                en_kg = [sup.sub("", i) for i in mmpbsa]
                en_kmol = ["{:.2f}".format(float(i)/4.18) for i in en_kg]
                mmpbsa_dict[time] = en_kmol
        mmpbsa_df = pd.DataFrame.from_dict(mmpbsa_dict, orient='index').T
        mean_mmpbsa = list(mmpbsa_df.loc[:, ['10ns', '20ns', '30ns', '40ns', '50ns']].astype('float').mean(1))

        mmpbsa_df.index = mmpbsa_df.pos
        mmpbsa_df = mmpbsa_df.loc[:, ['res']]
        mmpbsa_df['variant'] = sel_variants
        mmpbsa_df['mmpbsa'] = mean_mmpbsa
    return mmpbsa_df
            
merged_mmpbsa4csv = pd.DataFrame()
merged_res4csv = pd.DataFrame()
merged_mmpbsa4plot = []
sel_variants = []
sel_variants_sort = ['WT', 'B.1.1.529', 'BA.2', 'XBB.1.5', 'BA.2.86', 'JN.1', 'KP.2', 'KP.3', 'KP.3.1.1']
for dat in mmpbsa_dat_file:
    variant = dat.split("_")[0].replace('./', '')
    if variant not in sel_variants_sort:
        continue
    sel_variants.append(variant)
    print(variant)
    print(dat)
    if not os.path.exists(dat):
        sys.exit(-1)
    df4save_mmpbsa = read_mmpbsa(variant, dat).loc[:, ['mmpbsa']]
    df4save_mmpbsa.columns = [variant]
    df4save_res = read_mmpbsa(variant, dat).loc[:, ['res']]
    df4save_res.columns = [variant]
    
    merged_mmpbsa4plot.append(read_mmpbsa(variant, dat).loc[:, ['res', 'variant', 'mmpbsa']])
    
    
    if len(merged_mmpbsa4csv):
        merged_mmpbsa4csv = pd.merge(left=merged_mmpbsa4csv, right=df4save_mmpbsa, how='outer', left_index=True, right_index=True)
    else:
        merged_mmpbsa4csv = df4save_mmpbsa
        
    if len(merged_res4csv):
        merged_res4csv = pd.merge(left=merged_res4csv, right=df4save_res, how='outer', left_index=True, right_index=True)
    else:
        merged_res4csv = df4save_res

        
merged_mmpbsa4csv.to_csv(output1)
merged_res4csv.to_csv(output2)


df_merge = pd.concat(merged_mmpbsa4plot, axis=0)
df_merge['pos'] = df_merge.index

df_merge.index = list(range(df_merge.shape[0]))

df_merge.to_csv('merged_mmpbsa4plot.csv')

mut_pos = [339, 346, 356, 368, 371, 373, 375, 376, 403, 408, 417, 440, 445, 446, 450, 452, 455, 456, 460, 477, 478, \
           481, 483, 484, 486, 490, 493, 496, 498, 501, 505]
mut_pos = np.array(mut_pos) -333
#mut_res = ['408ARG', '417LYS', '452LEU', '478THR', '484GLU', '493GLN', '498GLN', '501ASN', '505TYR']
mut_res = mut_pos

# 窄型数据
# plt.figure(figsize=(3.4, 1.7))
plt.figure(figsize=(16, 10))
ax = sns.lineplot(x='pos', y='mmpbsa', hue='variant', #style="variant",\
             data=df_merge, alpha=.6, linewidth=2, palette=colors, hue_order=sel_variants_sort)

plt.ylabel("MM/PBSA(kcal/mol)")
plt.xlabel("Residues")
new_labels = [i.upper() if i == 'WT' else i.title() for i in sel_variants]

# Put the legend out of the figure
#sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, 1), ncol=2, title=None, frameon=False,)
ax.legend().set_title(None)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1.02, 1))
for t, l in zip(ax.get_legend().get_texts(), new_labels):
    t.set_text(l)
changeFontSize(plt.gca(), 20)
plt.xticks(mut_pos, mut_res)
plt.xticks(rotation=45)
plt.ylim(-70, 70)
for i in mut_pos:
    sns.lineplot(x=[i, i], y=[-69, 0], color='grey', linewidth=2)
plt.savefig(figfile, format='pdf', \
            bbox_inches='tight')

plt.show()
plt.close()