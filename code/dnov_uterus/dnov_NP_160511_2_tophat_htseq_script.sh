#PBS -l nodes=1:ppn=8
#PBS -m abe -M arun.chavan@yale.edu
cd $PBS_O_WORKDIR

##########################################################################################
# This script aligns raw reads to the genome with tophat2, and counts the reads mapping to gene features.
# v2016.06.08
# Arun Chavan
##########################################################################################

##########################################################################################
# 1 alignment with tophat2

# 1.1 change to tophat2 directory
cd /home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_mapping_counting/tophat2

# 1.2 run the tophat2 command
tophat2 -p 8 --b2-very-sensitive  --library-type fr-unstranded /home/arc78/genomic_data/dnov3_genomic_data/dnov3_index/dnov3 /home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_001.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_002.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_003.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_004.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_005.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_006.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_007.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_008.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_009.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_010.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_011.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_012.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_013.fastq.gz,/home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_data/160511_2_004_TGACCA_L001_R1_014.fastq.gz &> dnov_NP_160511_2_tophat2.log

# -p 8 specifies usage of 8 threads to accelerate mapping
# library type: low input sample library prepared with SMARTer kit, which works in a non-strand specific manner. 

##########################################################################################
# 2 counting with htseq-count

# 2.1 change to htseq directory
cd /home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_mapping_counting/htseq

# 2.2 run the htseq command
htseq-count -m intersection-nonempty -s no -i gene_id -f bam /home/arc78/rna_seq_data/LMD_low_input/uterus_lmd/dnov/dnov_NP_160511_2/dnov_NP_160511_2_mapping_counting/tophat2/tophat_out/accepted_hits.bam /home/arc78/genomic_data/dnov3_genomic_data/Dasypus_novemcinctus.Dasnov3.0.77.gtf > dnov_NP_160511_2_htseq_counts.txt 2> dnov_NP_160511_2_htseq_stderr.txt

# -m mode, -s strandedness, -i feaure attribute to be used from gtf file, -f format of the alignment file
