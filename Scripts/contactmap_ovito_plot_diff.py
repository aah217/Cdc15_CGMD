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
parser.add_argument("-f","--file", help="Ovito file of interest",type=str)
parser.add_argument("-r","--reffile",help="Ovito file to serve as reference",type=str)
parser.add_argument("-o","--outfile",help="file to save output (contactmap_ovito.xyz)",default="contactmap_diff_ovito.xyz",type=str)
args = parser.parse_args()

###
# Takes two xyz files made with contactmap_ovito_plot.py, assumes they used the same pdb/are the same length
# and takes the difference of the "velocity columns" which store contact averages and makes a new .xyz for visualization
###


file = args.file
reffile = args.reffile
if args.outfile is not None:
    outfile = open(args.outfile,'w')
else:
    outfile = sys.stdout 

file1 = []
with open(reffile,'r') as f:  
    delta = 0
    for i,r in enumerate(f):
        this_row = r.strip().split()
        file1.append(this_row)

with open(file,'r') as f:  
    delta = 0
    for i,r in enumerate(f):
        this_row = r.strip().split()
        if len(this_row) < 2:
            print(" ".join(this_row),file=outfile)
            continue
        this_IDR = float(this_row[4]) - float(file1[i][4])
        this_SH3 = float(this_row[5]) - float(file1[i][5])
        this_tot = float(this_row[6]) - float(file1[i][6])
        print(this_row[0],this_row[1],this_row[2],this_row[3],this_IDR,this_SH3,this_tot,file=outfile)
