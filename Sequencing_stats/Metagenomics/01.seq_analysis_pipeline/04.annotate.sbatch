#!/bin/bash
##SBATCH --job-name=anvio-annotate_CO
##SBATCH --output=anvi-annotate_sbatch.out
#SBATCH --error=anvi-annotate_sbatch.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --mem-per-cpu=2100
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

source ~/.bashrc
module load python
conda activate anvio-7

        anvi-gen-contigs-database -f 04_CONTIGS/contigs.fixed.fa -o 04_CONTIGS/contigs.db -n "WD_mg_all" -T 22

        anvi-run-hmms -c 04_CONTIGS/contigs.db --num-threads 22

        anvi-run-pfams -c 04_CONTIGS/contigs.db --num-threads 22 --pfam-data-dir /project2/ebchang/REF_databases/PFAM/

        anvi-run-kegg-kofams -c 04_CONTIGS/contigs.db -T 22 --kegg-data-dir /project2/ebchang/REF_databases/KOfams --just-do-it

        anvi-run-cazymes -c 04_CONTIGS/contigs.db -T 22

