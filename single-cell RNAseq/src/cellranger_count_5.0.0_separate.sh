#!/bin/bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=520G
#SBATCH -p normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=lma11@illinois.edu
#SBATCH -J cellranger_count
#SBATCH --output=count-%A.out
#SBATCH -D /home/a-m/lma11/scRNAseq/src/slurm-out


### Load Modules
module load cellranger/5.0.0

### specifying the output folder
cd /home/a-m/lma11/scRNAseq/results/cellranger

### point to a cellranger reference

DB=/home/a-m/lma11/scRNAseq/genome/cellranger_5.0.0_gencode_M26_filtered

### Run on file
cellranger count --id=792\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool1/792/\
    --sample=792\
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=780\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool1/780/\
    --sample=780\
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=794 \
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool1/794/\
    --sample=794\
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=779\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool2/779/\
    --sample=779\
    --expect-cells=6000 \
    --no-bam


### Run on file
cellranger count --id=778\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool2/778/\
    --sample=778\
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=793\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool2/793/\
    --sample=793\
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=791\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool3/791/\
    --sample=791\
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=789\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool3/789/\
    --sample=789 \
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=786\
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool3/786/\
    --sample=786\
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=800 \
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool4/800/\
    --sample=800 \
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=795 \
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool4/795/ \
    --sample=795 \
    --expect-cells=6000 \
    --no-bam

### Run on file
cellranger count --id=797 \
    --localcores=$SLURM_CPUS_ON_NODE  \
    --localmem 500 \
    --transcriptome=$DB \
    --fastqs=/home/a-m/lma11/scRNAseq/rawseq/Project_Nelson_SC_10X_Pool4/797/ \
    --sample=797 \
    --expect-cells=6000 \
    --no-bam


