#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=8
module load Julia/1.8.5-linux-x86_64

n_rounds=12

for n in {1..4}
do
	julia scripts/analysis.jl "$n" "$n_rounds" &
done

wait

