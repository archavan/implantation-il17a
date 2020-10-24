#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=ocun_sra_to_fastq
#SBATCH --ntasks=8 
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load SRA-Toolkit

# covert SRA data from GEO to fastq format.
# data: rabbit implantation stage uterus from GSE76115
# from cursory look the data appear to be single ended but I still using --split-files option to be safe.

# implantation site samples, 3 replicates
cd /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/ImpSite

fastq-dump -I --split-files --gzip --outdir /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/ImpSite SRR3029220
fastq-dump -I --split-files --gzip --outdir /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/ImpSite SRR3029221
fastq-dump -I --split-files --gzip --outdir /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/ImpSite SRR3029222

# interimplantation site samples, 3 replicates
cd /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/InterImpSite

fastq-dump -I --split-files --gzip --outdir /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/InterImpSite SRR3029223
fastq-dump -I --split-files --gzip --outdir /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/InterImpSite SRR3029224
fastq-dump -I --split-files --gzip --outdir /home/arc78/rna_seq_data/ocun/ocun_imp_Liu2016/data/InterImpSite SRR3029225
