#!/bin/bash
##SBATCH --job-name=IU_QUAL_FILT
##SBATCH --output=IU_QUAL.out
#SBATCH --error=IU_QUAL.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
module load python
conda activate anvio-7


for sample in `cat LIST4`;
do
	iu-filter-quality-minoche 02_QC/$sample.ini --ignore-deflines
done
