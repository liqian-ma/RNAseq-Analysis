#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 32
#SBATCH --mem=300g
#SBATCH -N 1
#SBATCH --mail-user=lma11@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -J kallisto_index
#SBATCH --output=kallisto_index-%A.out
#SBATCH -D /home/a-m/lma11/scRNAseq/src/slurm-out
# ----------------Load Modules--------------------
module load Python/3.7.2-IGB-gcc-8.2.0
# ----------------Commands------------------------
cd ~/scRNAseq/genome

kb ref -i gencode.vM26.kallisto.idx -g transcripts_to_genes.kallisto.txt \
-f1 gencode.vM26.cdna.kallisto.fa -f2 gencode.vM26.intron.kallisto.fa \
-c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno \
GRCm39.primary_assembly.genome.fa \
gencode.vM26.primary_assembly.annotation.filtered_cellranger_5.0.0.gtf
