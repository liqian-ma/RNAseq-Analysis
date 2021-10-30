#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=200G
#SBATCH -p normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lma11@illinois.edu
#SBATCH -J kb_count_intron
#SBATCH --output=kb-count-intron-%A.out
#SBATCH --array=1-12%3
#SBATCH -D /home/a-m/lma11/scRNAseq/src/slurm-out


### Load Modules
module load Python/3.7.2-IGB-gcc-8.2.0

cd /home/a-m/lma11/scRNAseq/src/

### point to a cellranger reference

idx=/home/a-m/lma11/scRNAseq/genome
Sample_ID=`head files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
File_folder=`head files.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 3`
R1=`head files.address.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
R2=`head files.address.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`
R3=`head files.address.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 3`
R4=`head files.address.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 4`
R5=`head files.address.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 5`
R6=`head files.address.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 6`


kb count -i $idx/gencode.vM26.kallisto.idx -g $idx/transcripts_to_genes.kallisto.txt \
-t 32 -m 200 -x 10xv3 \
--workflow lamanno --loom -c1 $idx/cdna_t2c.txt -c2 $idx/intron_t2c.txt \
-o /home/a-m/lma11/scRNAseq/results/kb_count/$Sample_ID \
$R1 $R4 \
$R2 $R5 \
$R3 $R6
