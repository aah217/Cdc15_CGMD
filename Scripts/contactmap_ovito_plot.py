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
#parser.add_argument("-o","--outfile",help="file to save output (fbar_heat.pdb)",default="fbar_heat.pdb",type=str)
args = parser.parse_args()
files = args.files.split()
pdbfile = args.pdb
#residues1 = [int(i) for i in args.residues.split()]
#residues2 = [int(i)+927 for i in args.residues.split()]
#residues2.reverse()
#residues = residues1+residues2
residues = []
with open(pdbfile,'r') as f:  
    delta = 0
    for i,r in enumerate(f):
        this_row = p.split(r.strip())
        if this_row[p.ATOM].strip() != "ATOM":
            delta = delta + 1
            continue
        residues.append(int(this_row[p.serial]))

IDR1 = [*range(301,870)]
IDR2 = [*range(1226,1797)]
IDRs = IDR1+IDR2
SH31 = [*range(870,928)]
SH32 = [*range(1797,1855)]
SH3s = SH31+SH32
#outfile = args.outfile
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
#out_df = pd.DataFrame().reindex_like(dfs[0])
out_tot = np.zeros((len(residues),5))
out_IDR = np.zeros((len(residues),5))
out_SH3 = np.zeros((len(residues),5))
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
file_tot = open("fbar_heat_tot.pdb",'w')
file_IDR = open("fbar_heat_IDR.pdb",'w')
file_SH3 = open("fbar_heat_SH3.pdb",'w')
with open(pdbfile,'r') as f:  
    delta = 0
    for i,r in enumerate(f):
        this_row = p.split(r.strip())
        if this_row[p.ATOM].strip() != "ATOM":
            delta = delta + 1
            continue
        this_i = i-delta
        this_row[p.serial] = p.put(this_row[p.serial],str(int(np.mean(out_tot[this_i])*1000)))
        print(p.join(this_row),file=file_tot)
        this_row[p.serial] = p.put(this_row[p.serial],str(int(np.mean(out_IDR[this_i])*1000)))
        print(p.join(this_row),file=file_IDR)
        this_row[p.serial] = p.put(this_row[p.serial],str(int(np.mean(out_SH3[this_i])*1000)))
        print(p.join(this_row),file=file_SH3)
