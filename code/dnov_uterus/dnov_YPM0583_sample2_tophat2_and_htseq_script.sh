#PBS -l nodes=1:ppn=8
#PBS -m abe -M arun.chavan@yale.edu
cd $PBS_O_WORKDIR

##########################################################################################
# This script aligns raw reads to the genome with tophat2, and counts the reads mapping to gene features.
# v2017.01.16
# Arun Chavan
##########################################################################################

##########################################################################################
# 1 alignment with tophat2

# 1.1 change to tophat2 directory
cd /home/arc78/rna_seq_data/dnov/placentome_YPM0583/sample_2/YPM0583_2_tophat2

# 1.2 run the tophat2 command
tophat2 -p 8 --b2-very-sensitive  --library-type fr-firststrand /home/arc78/genomic_data/dnov3_genomic_data/dnov3_index/dnov3 /home/arc78/rna_seq_data/dnov/placentome_YPM0583/sample_2/data/run1842_lane2_read1_indexD712-D501=ARM2.fastq.gz &> dnov_YPM0583_2_tophat2.log

# -p 8 specifies usage of 8 threads to accelerate mapping 

##########################################################################################
# 2 counting with htseq-count

# 2.1 change to htseq directory
cd /home/arc78/rna_seq_data/dnov/placentome_YPM0583/sample_2/YPM0583_2_htseq

# 2.2 run the htseq command
htseq-count -m intersection-nonempty -s reverse -i gene_id -f bam /home/arc78/rna_seq_data/dnov/placentome_YPM0583/sample_2/YPM0583_2_tophat2/tophat_out/accepted_hits.bam /home/arc78/genomic_data/dnov3_genomic_data/Dasypus_novemcinctus.Dasnov3.0.77.gtf > dnov_YPM0583_2_htseq_counts.txt 2> dnov_YPM0583_2_htseq_stderr.txt

# -m mode, -s strandedness, -i feaure attribute to be used from gtf file, -f format of the alignment file
