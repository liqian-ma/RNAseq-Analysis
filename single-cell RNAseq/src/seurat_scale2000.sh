#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 64
#SBATCH --mem=400g
#SBATCH -N 1
#SBATCH --mail-user=lma11@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -J seurat_scale
#SBATCH -D /home/a-m/lma11/scRNAseq/src/slurm-out
# ----------------Load Modules--------------------
module load R/4.0.3-IGB-gcc-8.2.0
# ----------------Commands------------------------
cd ../

Rscript seurat_scale2000.R

