RNAseq processing
1.	Download fastq data using Globus (Globus was the way RNA sequencing core facility used to transfer data in U of I)
2.	Unzip bz2 file using command $tar -xvf
3.	Make in one file the names of all unzipped fastq.gz files, using command line $ ls ~/RNAseq/raw_seq_1/*.fastq.gz > mouse_filenames.txt
4.	Remove .fastq.gz from each name in nano using ctrl \ (find and replace), then rename file names.txt using command line $ mv
5.	Download cDNA fasta file of mouse using command $wget ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
6.	Build index using command line $module load kallisto $ kallisto index -i Mus_musculus.GRC.index Mus_musculus.GRCm38.cdna.all.fa.gz
7.	Run pseudoalignment and abundance quantification using Kallisto on all fastq.gz files and save the output to invidual folders by sample $ xargs -I{} kallisto quant -i ../index/Mus_musculus.GRC.index -o ../kall_out_1/{}/ --single -l 140 -s 41.44 {}.fastq.gz < names.txt
