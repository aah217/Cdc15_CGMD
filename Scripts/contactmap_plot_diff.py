import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mdtraj as md
from contact_map import ContactFrequency, ContactDifference, AtomMismatchedContactDifference
from adjustText import adjust_text
import sys
import time
import argparse
import pandas as pd
startTime = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("-f1","--file1", help="contacts file source",type=str)
parser.add_argument("-f2","--file2", help="contacts file target",type=str)
parser.add_argument("-t","--ticks", help="ticks to draw",type=str)
parser.add_argument("-p","--pdbfile", help="pdb file (optional)",type=str)
parser.add_argument("-o","--outfile",help="file to save output, otherwise prints to screen",type=str)
args = parser.parse_args()
###
# This creates a contact map that is the difference of two pickle files made using contactmap.py
###

startTime = time.time()
file1 = args.file1
file2 = args.file2
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
trajectory_contacts1 = ContactFrequency.from_file(file1)
trajectory_contacts2 = ContactFrequency.from_file(file2)
diff = AtomMismatchedContactDifference(trajectory_contacts2, trajectory_contacts1)
print("Contact frequencies loaded")
print("1:",trajectory_contacts1.n_frames,"total frames")
print("2:",trajectory_contacts2.n_frames,"total frames")
print("Generating plot")
data = diff.residue_contacts.df
with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.expand_frame_repr', False):
    print(data,file=outfile)
fig, ax = plt.subplots(1, 1)
if use_pdb:
    new_axes = [res.resSeq for res in my_topology.residues]
    data.columns = new_axes
    data.index = new_axes
    c = ax.pcolorfast(new_axes,new_axes,data, cmap='RdBu', norm=colors.SymLogNorm(linthresh=0.001, linscale=1.0, vmin=-1.0, vmax=1.0))
else:
    c = ax.pcolorfast(data, cmap='RdBu', norm=colors.SymLogNorm(linthresh=0.001, linscale=1.0, vmin=-1.0, vmax=1.0))
ax.grid(b=True,which='both',axis='both',linestyle='-')
ax.set_aspect(1)
if len(my_ticks)>0:
    ax.set_xticks(my_ticks)
    ax.set_yticks(my_ticks)
ax.tick_params(length=6,width=2)
cbar = fig.colorbar(c, ax=ax)
tick_font_size = 16
cbar.ax.tick_params(labelsize=tick_font_size)
#fig.colorbar(c, ax=ax)
fig.tight_layout()
plt.show()

print("Done")
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))