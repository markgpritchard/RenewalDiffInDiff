#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=1
module load Julia/1.11.3-linux-x86_64 

for n in {1..22}
do
	julia scripts/analysis_sim.jl "$n" "2" "10000" "600" "2000" & 
done

wait

