### packages ==================================================================
library(tidyverse)
library(edgeR)

### data ======================================================================
counts <- read.csv(file = "data/rna-seq/read-counts_hsap_esf-dsc.csv")
rownames(counts) <- counts$hsap_ensembl_gid

### edgeR =====================================================================
## prepare --------------------------------------------------------------------
data.info <- c(rep('diff', 2), rep('undiff', 2))
data.info <- factor(data.info, levels = c('undiff', 'diff'))

y <- DGEList(counts = counts[, 3:ncol(counts)], 
             genes = counts[, 1:2], 
             group = data.info) 

## calculate CPM for plotting -------------------------------------------------
# calculate scaling factors and calculate cpm data for all genes. 
# NOTE: Don't use this for testing but only for plotting. The CPM that we calculate later on don't include filtered out lowly expressed genes, and so we can't plot them using those CPM values. 
x <- calcNormFactors(y) 
cpm.data <- cpm(x)
cpm.data <- data.frame(cpm.data)
cpm.data$hsap_ensembl_gid <- rownames(cpm.data)
cpm.data <- inner_join(x = y$genes[, 1:2], y = cpm.data, 
                       by = 'hsap_ensembl_gid')

## run DE analysis ------------------------------------------------------------
keep <- rowSums(cpm(y) > 1) >= 2 # filter low count genes
y <- y[keep, , keep.lib.sizes = F]

# calculate normalization factor -- by median of read counts; 
# while the factor of whole library size is taken internally in the package
y <- calcNormFactors(y)

# estimate dispersion
design <- model.matrix( ~ data.info)
y <- estimateDisp(y, design)

plotBCV(y)

# exact test
et <- exactTest(y, pair = c('undiff', 'diff')) # exact test

pval.thre <- 1e-6
res <- topTags(et, n = dim(y)[1]) # get data for all genes without filtering on pval

# get DE genes
de <- res$table[res$table$FDR < pval.thre,] # filter by pval
de <- de[order(abs(de$logFC), decreasing = T), ] # sort

plotMD(et)
plotMA(y)

## add DE analysis result to cpm.data -----------------------------------------
et.out <- data.frame(res$table)

cpm.full <- left_join(x = cpm.data, 
                      y = et.out, 
                      by = c('hsap_ensembl_gid', "hsap_gene_name"))

cpm.full$mean_diff <- apply(cpm.full[, 3:4], MARGIN = 1, FUN = 'mean')
cpm.full$mean_undiff <- apply(cpm.full[, 5:6], MARGIN = 1, FUN = 'mean')

write.csv(
  x = cpm.full, 
  file = "results/fig_supp_03_tcell-rna-seq/edgeR_hsap_esf-dsc/hsap_esf-dsc_cpm_with-DE-data.csv", 
  row.names = F
)

### end  ======================================================================
