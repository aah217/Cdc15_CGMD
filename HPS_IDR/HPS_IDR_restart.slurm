#!/bin/bash
#SBATCH -J Job_Name
#SBATCH -o myjob.%j.out
#SBATCH --partition=PNAME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 24:00:00
#SBATCH -A ################

module load lammps

mpirun -np 16 lmp -log my.log -partition 16x1 -nocite -in restart.lmp >& out.txt
