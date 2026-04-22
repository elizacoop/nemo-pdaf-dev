#!/bin/bash
#SBATCH --job-name=pdaf_test
#SBATCH --time=00:19:00
#SBATCH --account=n01-CRISP
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=95

export OMP_NUM_THREADS=1

# Load environmen
module swap PrgEnv-cray/8.4.0 PrgEnv-gnu/8.4.0
module load cray-mpich/8.1.27
module load cray-hdf5-parallel/1.12.2.7
module load cray-netcdf-hdf5parallel/4.9.0.7


# Run the job
srun PDAF_offline
