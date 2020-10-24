#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=ocun_Liu2016_ImpSite1
#SBATCH --ntasks=8 
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load TopHat
module load Bowtie2
module load Python
module load HTSeq
module load SAMtools

##########################################################################################
# This script aligns raw reads to the genome with tophat2, and counts the reads mapping to gene features.
# v2017.01.16
# Arun Chavan
##########################################################################################

##########################################################################################
# 1 alignment with tophat2

# 1.1 change to tophat2 directory
cd /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/ImpSite/imp1_SRR3029220/imp1_tophat

# 1.2 run the tophat2 command
tophat2 -p 8 --b2-very-sensitive  --library-type fr-unstranded /home/arc78/genomic_data/ocun2_genomic_data/index/ocun2 /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/ImpSite/SRR3029220_1.fastq.gz &> ocun_imp1_SRR3029220_tophat2.log

# -p 8 specifies usage of 8 threads to accelerate mapping 

##########################################################################################
# 2 counting with htseq-count

# 2.1 change to htseq directory
cd /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/ImpSite/imp1_SRR3029220/imp1_htseq

# 2.2 run the htseq command
htseq-count -m intersection-nonempty -s no -i gene_id -f bam /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/ImpSite/imp1_SRR3029220/imp1_tophat/tophat_out/accepted_hits.bam /home/arc78/genomic_data/ocun2_genomic_data/Oryctolagus_cuniculus.OryCun2.0.75.gtf > ocun_imp1_htseq_counts.txt 2> ocun_imp1_htseq_stderr.txt

# -m mode, -s strandedness, -i feaure attribute to be used from gtf file, -f format of the alignment file
