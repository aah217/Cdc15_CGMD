import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mdtraj as md
from contact_map import ContactFrequency, ContactDifference
from adjustText import adjust_text
import sys
import time
import argparse
import pandas as pd
startTime = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="contacts file",type=str)
parser.add_argument("-t","--ticks", help="ticks to draw",type=str)
parser.add_argument("-p","--pdbfile", help="pdb file (optional)",type=str)
parser.add_argument("-o","--outfile",help="file to save output, otherwise prints to screen",type=str)
args = parser.parse_args()

###
# This creates a simple contact map using a pickle file made from contactmap.py
###

startTime = time.time()
file = args.file
use_pdb = False
if args.pdbfile is not None:
    use_pdb = True
    pdb_file = args.pdbfile
    my_topology = md.load(pdb_file).topology
if args.outfile is not None:
    outfile = open(args.outfile,'w')
else:
    outfile = sys.stdout 
my_ticks = []
if args.ticks is not None:
    my_ticks = [int(i) for i in args.ticks.split()]
trajectory_contacts = ContactFrequency.from_file(file)
print("Contact frequency loaded")
print(trajectory_contacts.n_frames,"total frames")
print("Generating plot")
data = trajectory_contacts.residue_contacts.df
with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.expand_frame_repr', False):
    print(data,file=outfile)
fig, ax = plt.subplots(1, 1)
if use_pdb:
    new_axes = [res.resSeq for res in my_topology.residues]
    data.columns = new_axes
    data.index = new_axes
    c = ax.pcolorfast(new_axes,new_axes,data, cmap='gist_heat_r', norm=colors.SymLogNorm(linthresh=0.001, linscale=1.0, vmin=0.0, vmax=1.0))
else:
    c = ax.pcolorfast(data, cmap='gist_heat_r', norm=colors.SymLogNorm(linthresh=0.001, linscale=1.0, vmin=0.0, vmax=1.0))

ax.grid(b=True,which='both',axis='both',linestyle='-')

ax.set_aspect(1)
if len(my_ticks)>0:
    ax.set_xticks(my_ticks)
    ax.set_yticks(my_ticks)
ax.tick_params(length=6,width=2)
cbar = fig.colorbar(c, ax=ax)
tick_font_size = 16
cbar.ax.tick_params(labelsize=tick_font_size)
fig.tight_layout()
plt.show()
print("Done")
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))