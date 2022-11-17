#!/bin/bash
#SBATCH -J Job_Name
#SBATCH -o myjob.%j.out
#SBATCH --partition=PNAME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 72:00:00
#SBATCH -A ###########

module load lammps

lmp -log my.log -nocite -in run.lmp >& out.txt
