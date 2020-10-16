#! /bin/bash

#SBATCH --partition=general
#SBATCH --job-name=download_suryavanshi-data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

# covert SRA data from GEO to fastq format.
# data: Suryavanshi et al Science Advances 2018, RNA-seq data: BioProject PRJNA492902
# data are 10x single cell data from decidua
# SRA: SRR7895954

# decidua
for Sample in SRR7895954
do
	/home/arc78/applications/sratoolkit.2.9.4-centos_linux64/bin/fastq-dump --split-files -I --gzip --outdir /home/arc78/scratch60/suryavanshi/data/01_fastq/decidua $Sample

        # validate
        cd /home/arc78/scratch60/suryavanshi/data/01_fastq/decidua
	
	/home/arc78/applications/sratoolkit.2.9.4-centos_linux64/bin/vdb-validate $Sample &>/home/arc78/scratch60/suryavanshi/data/01_fastq/decidua/vdb-validate/$Sample"_validation.txt"

done

# villi
for Sample in SRR7895962
do
        /home/arc78/applications/sratoolkit.2.9.4-centos_linux64/bin/fastq-dump --split-files -I --gzip --outdir /home/arc78/scratch60/suryavanshi/data/01_fastq/villi $Sample

        # validate
        cd /home/arc78/scratch60/suryavanshi/data/01_fastq/villi

        /home/arc78/applications/sratoolkit.2.9.4-centos_linux64/bin/vdb-validate $Sample &>/home/arc78/scratch60/suryavanshi/data/01_fastq/villi/vdb-validate/$Sample"_validation.txt"

done

