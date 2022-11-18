from contact_map import ContactFrequency
import pdbfile as p
import argparse
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import pickle
import time
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
startTime = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("-f","--files", help="contacts files, \"filename1.p filename2.p ...\"",type=str)
parser.add_argument("-p","--pdb",help="pdb of residues to check",type=str)
parser.add_argument("-o","--outfile",help="file to save output (contactmap_ovito.xyz)",default="contactmap_ovito.xyz",type=str)
args = parser.parse_args()

###
# Does the following for a contact file created with contactmap.py:
# 1. For each amino acid in the pdb, finds the average number of contacts with residues in given ranges (IDRs, SH3, hardcoded)
#       -> To do this, it just sums the contact frequencies of this amino acid with all resiudes in the given list
# 2. Saves an .xyz file with the residue coordinates, puts the contact averages in the velocity columns to be visualized in ovito
###

files = args.files.split()
pdbfile = args.pdb
if args.outfile is not None:
    outfile = open(args.outfile,'w')
else:
    outfile = sys.stdout 

residues = []
count = 1
with open(pdbfile,'r') as f:  
    delta = 0
    for i,r in enumerate(f):
        this_row = p.split(r.strip())
        if this_row[p.ATOM].strip() != "ATOM":
            delta = delta + 1
            continue
        residues.append(int(this_row[p.serial]))
    count = count + i - delta

#These are regions of interest that are hardcoded 
IDR1 = [*range(301,870)]
IDR2 = [*range(1226,1797)]
IDRs = IDR1+IDR2
SH31 = [*range(870,928)]
SH32 = [*range(1797,1855)]
SH3s = SH31+SH32

dfs = []
threshold_check = 0.01
threshold_conv = 0.1
convergence = True

for file in files:
    this_traj = ContactFrequency.from_file(file)
    this_df = this_traj.residue_contacts.df
    this_df = this_df.sparse.to_dense()
    this_df = this_df.applymap(lambda x: 0.0 if math.isnan(x) else x)
    dfs.append(this_df)
col_dim = len(dfs[0])
row_dim = len(dfs[0][0])

out_tot = np.zeros((len(residues),10))
out_IDR = np.zeros((len(residues),10))
out_SH3 = np.zeros((len(residues),10))
for i,r in enumerate(residues):
    frame = 0
    this_i = r-1
    for df in dfs:
        val_tot = 0.0
        val_IDR = 0.0
        val_SH3 = 0.0
        for j in range(0,row_dim):
            if j not in IDRs and j not in SH3s:
                continue
            this_val = df[this_i][j]
            val_tot = val_tot + this_val
            if j in IDRs: val_IDR = val_IDR + this_val
            if j in SH3s: val_SH3 = val_SH3 + this_val
        out_tot[i][frame] = val_tot
        out_IDR[i][frame] = val_IDR
        out_SH3[i][frame] = val_SH3
        frame = frame + 1
print(count,file=outfile)
print("",file=outfile)
with open(pdbfile,'r') as f:  
    delta = 0
    for i,r in enumerate(f):
        this_row = p.split(r.strip())
        if this_row[p.ATOM].strip() != "ATOM":
            delta = delta + 1
            continue
        this_i = i-delta
        this_tot = np.mean(out_tot[this_i])
        this_IDR = np.mean(out_IDR[this_i])
        this_SH3 = np.mean(out_SH3[this_i])
        print(this_i,this_row[p.x],this_row[p.y],this_row[p.z],this_IDR,this_SH3,this_tot,file=outfile)
