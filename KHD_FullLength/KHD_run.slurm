#!/bin/bash
#SBATCH -J Job_Name
#SBATCH -o myjob.%j.out
#SBATCH --partition=PNAME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 72:00:00
#SBATCH -A #########

/folder/lmp_serial -log my."$prefix".log -nocite -in run.lmp >& out."$prefix".txt & 

wait
