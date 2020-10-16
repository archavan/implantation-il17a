#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=d0_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load TopHat
module load Bowtie2
module load Python
module load HTSeq
module load SAMtools
module load FastQC
module load cutadapt

##########################################################################################
# This script aligns raw reads to the genome with tophat2, and counts the reads mapping to gene features.
# v2017.11.15
# Arun Chavan
##########################################################################################

# trim reads
cd /home/arc78/scratch60/tcells/trimming_mapping/d0_2

/home/arc78/applications/trim_galore/trim_galore --gzip --paired --output_dir ./trimmed /home/arc78/rna_seq_data/hsap/tcells/data/d0_2/d0_2_S2_L001_R1_001.fastq.gz /home/arc78/rna_seq_data/hsap/tcells/data/d0_2/d0_2_S2_L001_R2_001.fastq.gz


##########################################################################################
# 1 alignment with tophat2

# 1.1 change to tophat2 directory
cd /home/arc78/scratch60/tcells/trimming_mapping/d0_2/d0_2_tophat

# 1.2 run the tophat2 command
tophat2 -p 20 --b2-very-sensitive --no-coverage-search -G /home/arc78/genomic_data/human_genomic_data/GRCh37/Homo_sapiens.GRCh37.87.gtf --library-type fr-unstranded /home/arc78/genomic_data/human_genomic_data/GRCh37/index/hsap_GRCh37.75 /home/arc78/scratch60/tcells/trimming_mapping/d0_2/trimmed/d0_2_S2_L001_R1_001_val_1.fq.gz /home/arc78/scratch60/tcells/trimming_mapping/d0_2/trimmed/d0_2_S2_L001_R2_001_val_2.fq.gz &> d0_2_tophat2.log

# -p 8 specifies usage of 8 threads to accelerate mapping 
# --no-coverage-search disables coverage search for junctions. It is disabled by default for reads 75bp or longer. The use of --coverage-search is recommended for short reads (less than 45bp) and small number of reads (less than 10 million). The reads for naiveT cells are 50bp and are more than 10 million. 

# sort bam file
# sorting is essential for paired end reads because htseq-counts expects matepairs to be one after the other in bam file. The default of samtools sort is by coordinates, but for matepairs to be close to each other we need sorting by read name, done by -n option. 

cd tophat_out
samtools sort -n accepted_hits.bam -o accepted_hits_sorted.bam

##########################################################################################
# 2 counting with htseq-count

# 2.1 change to htseq directory
cd /home/arc78/scratch60/tcells/trimming_mapping/d0_2/d0_2_htseq

# 2.2 run the htseq command
htseq-count -m intersection-nonempty -s reverse -i gene_id -f bam /home/arc78/scratch60/tcells/trimming_mapping/d0_2/d0_2_tophat/tophat_out/accepted_hits_sorted.bam /home/arc78/genomic_data/human_genomic_data/GRCh37/Homo_sapiens.GRCh37.87.gtf > d0_2_htseq_counts.txt 2> d0_2_htseq_stderr.txt

# -m mode, -s strandedness, -i feaure attribute to be used from gtf file, -f format of the alignment file
