#!/bin/bash
#
#SBATCH --job-name=surfanalysis
#
#SBATCH --ntasks=8
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=16000

module load gfortran/6.3.0 hdf5/1.10.2

# used to track the z=0 haloes back in time and determine their asymmetry parameters with respect to earlier redshifts
for sn in {70..198}
do
   srun -n 1 /home/dobreschkow/surfanalysis/surfanalysis asymmetries -snapshot 199 -initial $sn -subhalos 0 -parameterfile /home/dobreschkow/surfanalysis/parameters.txt -parameterset L210_N1024-Hydro6D-new-hyades -logfile /home/dobreschkow/surfanalysis/log_sub0_199_$sn.txt &
   srun -n 1 /home/dobreschkow/surfanalysis/surfanalysis asymmetries -snapshot 199 -initial $sn -subhalos 1 -parameterfile /home/dobreschkow/surfanalysis/parameters.txt -parameterset L210_N1024-Hydro6D-new-hyades -logfile /home/dobreschkow/surfanalysis/log_sub1_199_$sn.txt &
done
wait

# used to compute lambda of haloes identified at discret redshifts
for sn in {156,131,114,100,88}
do
   srun -n 1 /home/dobreschkow/surfanalysis/surfanalysis asymmetries -snapshot $sn -initial $((sn-1)) -subhalos 0 -parameterfile /home/dobreschkow/surfanalysis/parameters.txt -parameterset L210_N1024-Hydro6D-new-hyades -logfile /home/dobreschkow/surfanalysis/log_sub0_$sn.txt &
   srun -n 1 /home/dobreschkow/surfanalysis/surfanalysis asymmetries -snapshot $sn -initial $((sn-1)) -subhalos 1 -parameterfile /home/dobreschkow/surfanalysis/parameters.txt -parameterset L210_N1024-Hydro6D-new-hyades -logfile /home/dobreschkow/surfanalysis/log_sub1_$sn.txt &
done
wait