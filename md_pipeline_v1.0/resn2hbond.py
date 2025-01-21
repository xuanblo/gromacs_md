#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 11:32:56 2022

@author: xuanzhang
"""

import re
import sys
# giving resn

resn = ['N501Y','E484K','L452R','T478K','Q493R','Q498R']
resi = [501, 484, 452, 478, 493, 498]


# read .gro file
#variant_gro = "delta/md_0_1.gro"
variant_gro = sys.argv[1]
atom2resi_dict = dict()
with open(variant_gro, "r") as f:
    for line in f.readlines():
        items = line.split()
        if len(items) > 5:
            matObj = re.match("^\s+\d+\w+\s+", line)
            if matObj and ("SOL" not in line) and ("NA" not in line):
                #print(line[0:9], line[15:21])
                atom_index = line[15:20].lstrip()
                resn = line[0:9].lstrip()
                atom2resi_dict[atom_index] = resn
      
# read h_bond result
#hbdat = "delta/hbdat.dat"
hbdat = sys.argv[2]
variant = sys.argv[3]
# 1-9506 ace2 Protein_chain_A
# 9507-12520 rbd Protein_chain_B
print("{}\tace2_resi\tace2_resn\tace2_atom\trbd_resi\trbd_resn\trbd_atom\toccupancy%".format(variant))
with open(hbdat, "r") as f:
    for line in f.readlines():
        if line.startswith("#"):
            continue
        else:
            line = line.rstrip()
            items = line.split()
            items_1_matObj = re.match('(\d+)#"(.*)"_(\w+)@.*', items[0])
            if items_1_matObj:
                donor_atom_index = items_1_matObj.group(1)
                donor_chain = items_1_matObj.group(2)
                donor_resn = items_1_matObj.group(3)
            items_2_matObj = re.match('(\d+)#"(.*)"_(\w+)@.*', items[2])
            if items_1_matObj:
                acceptor_atom_index = items_2_matObj.group(1)
                acceptor_chain = items_2_matObj.group(2)
                acceptor_resn = items_2_matObj.group(3)
            occupancy = items[3]
            if donor_chain == "Protein_chain_A":
                rbd_resi = atom2resi_dict[acceptor_atom_index][0:3]
            else:
                rbd_resi = atom2resi_dict[donor_atom_index][0:3]
                
            if int(rbd_resi) not in resi:
                continue
            
            if donor_chain == "Protein_chain_A":
                print(variant, atom2resi_dict[donor_atom_index], donor_resn, donor_atom_index, \
                      atom2resi_dict[acceptor_atom_index], acceptor_resn, acceptor_atom_index, \
                          occupancy, sep="\t")
            else:
                print(variant, atom2resi_dict[acceptor_atom_index], acceptor_resn, acceptor_atom_index, \
                      atom2resi_dict[donor_atom_index], donor_resn, donor_atom_index, \
                          occupancy, sep="\t")
            
        

