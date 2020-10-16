#! /bin/bash
#SBATCH --partition=general
#SBATCH --job-name=suryavanshi_cellranger-count
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

### Decidua sample
## 1. rename fastq files to match 10x format 
# run only once
cd /home/arc78/scratch60/suryavanshi/data/01_fastq/decidua
mv SRR7895954_1.fastq.gz SRR7895954_S1_L001_R1_001.fastq.gz
mv SRR7895954_2.fastq.gz SRR7895954_S1_L001_R2_001.fastq.gz

### 2. run cellranger count
cd /home/arc78/scratch60/suryavanshi/data/02_cellranger-count

/home/arc78/applications/cellranger-3.1.0/cellranger count --id=decidua\
  --sample=SRR7895954\
  --fastqs=/home/arc78/scratch60/suryavanshi/data/01_fastq/decidua\
  --transcriptome=/gpfs/ysm/datasets/genomes/10xgenomics/refdata-cellranger-hg19-3.0.0

### Villi sample
## 1. rename fastq files to match 10x format 
# run only once
cd /home/arc78/scratch60/suryavanshi/data/01_fastq/villi
mv SRR7895962_1.fastq.gz SRR7895962_S1_L001_R1_001.fastq.gz
mv SRR7895962_2.fastq.gz SRR7895962_S1_L001_R2_001.fastq.gz

### 2. run cellranger count
cd /home/arc78/scratch60/suryavanshi/data/02_cellranger-count

/home/arc78/applications/cellranger-3.1.0/cellranger count --id=villi\
  --sample=SRR7895962\
  --fastqs=/home/arc78/scratch60/suryavanshi/data/01_fastq/villi\
  --transcriptome=/gpfs/ysm/datasets/genomes/10xgenomics/refdata-cellranger-hg19-3.0.0



