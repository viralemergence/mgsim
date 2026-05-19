#!/bin/bash
## Project code
#PBS -A UCIE0003
## Job name
#PBS -N derecho_mgsim_round3a
## Queue
#PBS -q main
## Output log format
#PBS -j oe
## Email notifications
#PBS -m abe
## How long do I want on the nodes
#PBS -l walltime=12:00:00
## Computational allocation
#PBS -l select=1:ncpus=128:mpiprocs=128:mem=235gb

### Set temp to scratch
export TMPDIR=${SCRATCH}/${USER}/temp && mkdir -p $TMPDIR


# Load modules to match compile-time environment
module --force purge
module load ncarenv/23.09 craype/2.7.23 gcc/13.2.0 hdf5/1.14.3 cray-mpich/8.1.27 conda/latest proj/9.2.1 ncarcompilers/1.0.0 geos/3.12.1 gdal/3.8.1


# Activate conda environment
conda activate mgsim

### Run script
Rscript /glade/u/home/pilowskyj/mgsim/Scripts/simulation_round3a.R
