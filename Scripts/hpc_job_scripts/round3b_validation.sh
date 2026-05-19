#!/bin/bash -l
#PBS -A UCIE0003
#PBS -N round3b_validation
#PBS -q main
#PBS -j oe
#PBS -m abe
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=120:ompthreads=120:mem=235GB

### Set temp to scratch
export TMPDIR=${SCRATCH}/${USER}/temp && mkdir -p $TMPDIR

# Load modules to match compile-time environment
module --force purge
module load ncarenv/23.06 intel/2023.0.0 ncarcompilers/1.0.0 netcdf/4.9.2 proj/8.2.1 openmpi/main conda/latest geos/3.9.1 hdf5/1.12.2

# Activate conda environment
conda activate mgsim

# Run script
Rscript /glade/u/home/pilowskyj/mgsim/Scripts/round3b_validation.R