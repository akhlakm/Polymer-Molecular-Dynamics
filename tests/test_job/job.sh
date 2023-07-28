#!/bin/bash
#SBATCH --job-name=PE_equilibration
#SBATCH --account=GT-rramprasad3-CODA20
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --output=out.o%j
#SBATCH --error=err.e%j
module load intel/20.0.4 mvapich2/2.3.6-z2duuy lammps/20220107-mva2-dukitd
srun -n 24 lmp -in lmp.in
