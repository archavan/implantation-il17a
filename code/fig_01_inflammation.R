### packages ==================================================================
library(tidyverse)
library(cowplot)

### data ======================================================================
dn.ht <- read.csv("data/rna-seq/read-counts_dnov_uterus.csv")
dn.fl <- read.csv("data/genomic-data/feature-lengths_dnov.csv")

### calculate tpm =============================================================
dn.ht <- dplyr::inner_join(dn.fl, dn.ht, by = "dnov_ensembl_gid") 

# create empty df for tpm
tpm <- data.frame(matrix(nrow = dim(dn.ht)[1], ncol = dim(dn.ht)[2]))
colnames(tpm) <- colnames(dn.ht)
tpm[, c("dnov_ensembl_gid", "dnov_feature_length", "dnov_gene_name")] <- dn.ht[, c("dnov_ensembl_gid", "dnov_feature_length", "dnov_gene_name")]

# calculate
for(i in 4:ncol(tpm)) {
  tpm[, i] <- ((dn.ht[, i] / dn.ht$dnov_feature_length) * 10 ^ 6) / sum(dn.ht[, i] / dn.ht$dnov_feature_length)
}
tpm$dnov_feature_length <- NULL

# verify that each column adds up to a million
for (i in 3:ncol(tpm)){
  print(sum(tpm[, i]))
} 

# average tpm for imp stage
tpm$dnov_imp <- rowMeans(tpm[, c("dnov_imp_1", "dnov_imp_2")])
tpm$dnov_imp_1 <- NULL
tpm$dnov_imp_2 <- NULL

# rename
tpm <- dplyr::rename(tpm, "dnov_NP" = "dnov_NP_1")

### Fig. 1, Panel B: heatmap ==================================================
# create a subset with inflammation-related genes
imp <-  tpm[tpm$dnov_gene_name %in%  c("IL1B", "IL1A", "LIF", "TNF", "IL8", "PTGS2", "PTGES", "PTGER2", "PTGER4", 'IL10', 'IL1R1'),] # without IL17
IL6 <- tpm[tpm$dnov_ensembl_gid %in% c("ENSDNOG00000017107","ENSDNOG00000023693"), ] # IL6 has two orthologs
imp <- rbind(imp, IL6)
imp <- subset(imp, dnov_ensembl_gid != "ENSDNOG00000023693") # remove un-expressed ortholog of IL6
imp[imp$dnov_ensembl_gid == 'ENSDNOG00000017107',  "dnov_gene_name"] <- 'IL6' # keep expressed ortholog of IL6
imp <- imp[order(imp$dnov_imp), ] # sort by tpm at implantation

# melt
imp.long <- pivot_longer(data = imp, 
                         cols = c("dnov_NP", "dnov_imp"), 
                         names_to = "stage", 
                         values_to = "tpm")
imp.long$stage <- factor(imp.long$stage, levels = c("dnov_NP", "dnov_imp"))
imp.long$dnov_gene_name <- factor(imp.long$dnov_gene_name, levels = imp$dnov_gene_name)

# plot
pan.b <- ggplot(data = imp.long, 
                aes(x = stage, y = dnov_gene_name)) +
  geom_tile(aes(fill = tpm), colour = "white") + 
  scale_fill_gradient(name = 'TPM', low = "white", high = "red") +
  geom_text(aes(stage, dnov_gene_name, label = round(tpm, 1)), 
            color = "black", size = 6/.pt) +
  scale_y_discrete(limits = levels(imp.long$dnov_gene_name)) +
  scale_x_discrete(breaks = c('dnov_NP', 'dnov_imp'),
                   labels = c('Non\npregnant', 'Peri\nimplantation'), 
                   position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_blank(),
        axis.ticks.y = element_line(size = 0.25),
        panel.grid = element_blank(),
        plot.margin = margin(6, 30, 6, 6),
        legend.position = c(1, 0),
        legend.justification = c(0, 0),
        legend.key.size = unit(0.6, "lines"),
        legend.title = element_text(color = "black", size = 6),
        legend.text = element_text(color = "black", size = 6))

cowplot::ggsave2(
  filename = 'results/fig_01_inflammation/fig01_panel-B_dnov-inflam-genes.pdf',
  pan.b, width = 2.1, height = 2.5, units = "in")

### Fig. 1, Panel C: revigo ===================================================
