#!/bin/bash
#SBATCH -n 24
#SBATCH --mem=120G
#SBATCH -p normal
#SBATCH --mail-user=lma11@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J cell-ranger-index
#SBATCH -D /home/a-m/lma11/scRNAseq/src/slurm-out

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load cellranger/5.0.0


### Filter gtf per 10X's recommendation at https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_2020A

### Make cellranger reference

cd /home/a-m/lma11/scRNAseq/genome

cellranger mkref \
    --genome=cellranger_5.0.0_gencode_M26_filtered \
    --fasta=GRCm39.primary_assembly.genome.fa \
    --genes=gencode.vM26.primary_assembly.annotation.filtered_cellranger_5.0.0.gtf \
    --ref-version='gencode_M26_filtered_cellranger_5.0.0' \
    --nthreads=$SLURM_NTASKS \
    --memgb=120

