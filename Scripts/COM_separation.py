import MDAnalysis as mda
import argparse
import time
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f","--xtcfile", help="trajectory file, xtc format",type=str)
parser.add_argument("-s1","--sel1", help="selection 1",type=str)
parser.add_argument("-s2","--sel2", help="selection 2",type=str)
parser.add_argument("-s","--strfile", help="pdb file",type=str)
parser.add_argument("-o","--outfile",help="file to save output (COM_separation.dat)",default="COM_separation.dat",type=str)
args = parser.parse_args()

###
# Does the following for an xtc trajectory:
# 1. Calculates center of mass of two selections versus time
# 2. Outputs a time series with the xyz difference of the COM locations and distance magnitude
###

xtc_file = args.xtcfile
str_file = args.strfile
outfile = open(args.outfile,'w')
sel1 = args.sel1
sel2 = args.sel2
u = mda.Universe(str_file,xtc_file)
chain1 = u.select_atoms(sel1)
chain2 = u.select_atoms(sel2)
startTime = time.time()

print("Calculating COM separation...",flush=True)
for ts in u.trajectory:
        separation = chain1.center_of_mass()-chain2.center_of_mass()
        mag = np.sqrt(separation[0]**2+separation[1]**2+separation[2]**2)
        print(ts.time, " ".join([str(np.abs(each)) for each in separation]),mag,file=outfile)
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime),flush=True)
