import mdtraj as md
from contact_map import ContactFrequency
import sys
import argparse
import time
startTime = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="trajectory file, xtc format",type=str)
parser.add_argument("-p","--pdbfile", help="pdb file",type=str)
parser.add_argument("-o","--outfile",help="file to save output (contactmap.p)",default="contactmap.p",type=str)
args = parser.parse_args()
pdb_file = args.pdbfile
xtc_file = args.file
traj_name = args.outfile
this_cutoff = 1.0
this_chunk = 1000
firstrun = True
print("Calculating contact frequency of "+xtc_file, flush=True)
for i,chunk in enumerate(md.iterload(xtc_file,top=pdb_file,chunk=this_chunk)):
    topology = chunk.topology
    trajectory_contacts = ContactFrequency(trajectory=chunk,cutoff=this_cutoff)
    if firstrun:
        firstrun = False
    else:
        trajectory_contacts.add_contact_frequency(ContactFrequency.from_file(traj_name))
    trajectory_contacts.save_to_file(traj_name)
    print("Frame",(i+1)*this_chunk, flush=True)
print("Contact frequency saved to "+traj_name, flush=True)
print(trajectory_contacts.n_frames,"total frames", flush=True)
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime), flush=True)