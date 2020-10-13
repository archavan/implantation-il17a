### packages ==================================================================
library(tidyverse)

### data ======================================================================
# htseq counts
dn.ht <- read.csv("data/rna-seq/read-counts_dnov_uterus.csv") # armadillo
oc.ht <- read.csv("data/rna-seq/read-counts_ocun_uterus.csv") # rabbit
md.ht <- read.csv("data/rna-seq/read-counts_mdom_uterus.csv") # opossum

### genomic data ==============================================================
# feature lengths 
dn.fl <- read.csv("data/genomic-data/feature-lengths_dnov.csv")
oc.fl <- read.csv("data/genomic-data/feature-lengths_ocun.csv")
md.fl <- read.csv("data/genomic-data/feature-lengths_mdom.csv")

# orthologs
orth <- read.csv("data/genomic-data/one-to-one-orthologs_mdom-ocun-dnov.csv")
orth[, c("dnov_homology_type", "ocun_homology_type")] <- NULL
orth[, c("dnov_gene_name", "mdom_gene_name", "ocun_gene_name")] <- NULL

# cytokines
ck <- read.table("data/genomic-data/mdom_GO0005125_CytokineActivity_GRCh37Biomart.txt", header = TRUE, quote = "", sep = "\t")

### Fig 3A: cytokine heatmap ==================================================
## join data for all species by orthology -------------------------------------
dn.orth <- inner_join(orth, dn.ht, by = c("dnov_ensembl_gid"))
oc.orth <- inner_join(orth, oc.ht, by = c("ocun_ensembl_gid"))
md.orth <- inner_join(orth, md.ht, by = c("mdom_ensembl_gid"))

## join feature lengths -------------------------------------------------------
dn.orth <- inner_join(dn.fl, dn.orth, by = "dnov_ensembl_gid")
oc.orth <- inner_join(oc.fl, oc.orth, by = "ocun_ensembl_gid")
md.orth <- inner_join(md.fl, md.orth, by = "mdom_ensembl_gid")

## calculate tpm --------------------------------------------------------------
dn.tpm <- dn.orth
dn.tpm[, c("dnov_NP_1", "dnov_imp_1", "dnov_imp_2")] <- NA
for(i in c("dnov_NP_1", "dnov_imp_1", "dnov_imp_2")) {
  dn.tpm[, i] <- ((dn.orth[, i] / dn.orth$dnov_feature_length) * 10 ^ 6) / 
    sum(dn.orth[, i] / dn.orth$dnov_feature_length)
}
dn.tpm$dnov_feature_length <- NULL

oc.tpm <- oc.orth
oc.tpm[, c("ocun_imp_1", "ocun_imp_2", "ocun_imp_3")] <- NA
for(i in c("ocun_imp_1", "ocun_imp_2", "ocun_imp_3")) {
  oc.tpm[, i] <- ((oc.orth[, i] / oc.orth$ocun_feature_length) * 10 ^ 6) /
    sum(oc.orth[, i] / oc.orth$ocun_feature_length)
}
oc.tpm$ocun_feature_length <- NULL

md.tpm <- md.orth
md.tpm[, c("mdom_NP_1", "mdom_NP_2", "mdom_NP_3", "mdom_13.5dpc_1", "mdom_13.5dpc_2")] <- NA
for(i in c("mdom_NP_1", "mdom_NP_2", "mdom_NP_3", "mdom_13.5dpc_1", "mdom_13.5dpc_2")) {
  md.tpm[, i] <- ((md.orth[, i] / md.orth$mdom_feature_length) * 10 ^ 6) /
    sum(md.orth[, i] / md.orth$mdom_feature_length)
}
md.tpm$mdom_feature_length <- NULL

## merge ----------------------------------------------------------------------
tpm <- md.tpm %>% 
  inner_join(., dn.tpm, c("mdom_ensembl_gid", "dnov_ensembl_gid", "ocun_ensembl_gid")) %>% 
  inner_join(., oc.tpm,  c("mdom_ensembl_gid", "dnov_ensembl_gid", "ocun_ensembl_gid"))

## calculate mean tpm ---------------------------------------------------------
tpm$mdom_NP <- apply(tpm[, c("mdom_NP_1", "mdom_NP_2", "mdom_NP_3")], MARGIN = 1, FUN = mean)
tpm$mdom_13.5dpc <- apply(tpm[, c("mdom_13.5dpc_1", "mdom_13.5dpc_2")], MARGIN = 1, FUN = mean)
tpm$dnov_NP <- tpm$dnov_NP_1
tpm$dnov_imp <- apply(tpm[, c("dnov_imp_1", "dnov_imp_2")], MARGIN = 1, FUN = mean)
tpm$ocun_imp <- apply(tpm[, c("ocun_imp_1", "ocun_imp_2", "ocun_imp_3")], MARGIN = 1, FUN = mean)

## clean up -------------------------------------------------------------------
tpm <- tpm[, c("mdom_ensembl_gid", "dnov_ensembl_gid", "ocun_ensembl_gid",
               "mdom_gene_name", "dnov_gene_name", "ocun_gene_name",
               "mdom_NP", "mdom_13.5dpc",
               "dnov_NP", "dnov_imp",
               "ocun_imp")]

## subset to cytokines --------------------------------------------------------
tpm.ck <- filter(tpm, mdom_ensembl_gid %in% ck$Gene.stable.ID)

## subset to pregnant samples -------------------------------------------------
tpm.ck <- tpm.ck[, !(names(tpm.ck) %in% c("mdom_NP", "dnov_NP"))]

## keep genes expressed in at least one species -------------------------------
tpm.ck <- filter(tpm.ck, 
                 mdom_13.5dpc > 3 | dnov_imp > 3 | ocun_imp > 3) # 3 tpm cutoff

## prepare data for plotting --------------------------------------------------
# melt
tpm.ck.long <- pivot_longer(tpm.ck, 
                            cols = c("mdom_13.5dpc", "dnov_imp", "ocun_imp"),
                            names_sep = "_",
                            names_to = c("species", "stage"), 
                            values_to = "TPM")

# reorder genes based on expression pattern: opossum specific, eutheria specific, etc
geneorder <- c("IL17A", "TNF", "CCL20", "CSF3", 
               "CXCL9", "LTA",
               "IL10", "IL1A", 
               "IL20", "IL17B", "FGF2", "WNT2", "CCL25", "CCL17",    
               "INHBA", "CCL21", "IL17D","IL16","CCL19", "TNFSF15", "KITLG", "CTF1",  
               "WNT7A", "LIF", "NDP", "CCL5", "CXCL10", 
               "CXCL12", "TNFSF10", "CSF1", "TNFSF13", "BMP4", "CX3CL1", "WNT5A", "EDN1", "IL15", "LTB", "CXCL8", "IL1B")

tpm.ck.long$mdom_gene_name <- factor(tpm.ck.long$mdom_gene_name, levels = geneorder)

tpm.ck.long$species <- factor(tpm.ck.long$species, levels = c("mdom", "dnov", "ocun"))

# fontface vector for bolding opossum and eutheria specific genes
face <- ifelse(levels(tpm.ck.long$mdom_gene_name) %in% 
                 c("IL17A", "TNF", "CCL20", "CSF3", # opossum specifc
                   "WNT7A", "LIF", "NDP", "CCL5", "CXCL10"), # eutheria specific
               yes = 2, 
               no = 1
               )

ann.line.col <- "grey50" # for annotations on right side
ann.line.size <- 2

## plot -----------------------------------------------------------------------
ck.heat <- ggplot(tpm.ck.long, 
                  aes(x = species, y = mdom_gene_name, label = round(TPM, 1))) +
  geom_tile(aes(fill = sqrt(TPM)), 
            colour='white') +
  scale_fill_gradient2(low='#FFFFFF', mid = '#FF6666', high = '#FF0000', 
                       midpoint = 20, 
                       guide = 'colourbar',
                       breaks = c(0, 10, 20, 30),
                       labels = c(0, 100, 400, 900),
                       name = "TPM") +
  scale_x_discrete(breaks = c("mdom", "dnov", "ocun"),
                   labels = c("Opossum", "Armadillo", "Rabbit")) +
  coord_cartesian(expand = FALSE, clip = "off") +
  geom_text(size = 6/.pt) +
  theme_minimal() +
  theme(
    plot.background = element_blank(),
    plot.margin = margin(6, 100, 6, 6),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks.length.y = unit(1, "pt"),
    axis.text.y = element_text(size = 7, colour = "black", face = face),
    axis.text.x = element_text(size = 7, colour = "black"),
    legend.title = element_text(size = 6, colour = "black"),
    legend.text = element_text(size = 6, colour = "black"),
    legend.key.size = unit(0.75, "lines"),
    legend.position = c(1.005, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank()
    ) +
  annotate(geom = "segment", x = 3.65, xend = 3.65, y = 0.7, yend = 4.3,
           color = ann.line.col, size = ann.line.size) +
  annotate(geom = "text", x = 3.8, y = 4.3, label = "Opossum\nspecific",
           color = "black", size = 7/.pt, fontface = 2,
           hjust = 0, vjust = 1) +
  annotate(geom = "segment", x = 3.65, xend = 3.65, y = 4.7, yend = 6.3,
           color = ann.line.col, size = ann.line.size) +
  annotate(geom = "segment", x = 3.65, xend = 3.65, y = 6.7, yend = 8.3,
           color = ann.line.col, size = ann.line.size) +
  annotate(geom = "segment", x = 3.65, xend = 3.65, y = 8.7, yend = 14.3,
           color = ann.line.col, size = ann.line.size) +
  annotate(geom = "segment", x = 3.65, xend = 3.65, y = 14.7, yend = 22.3,
           color = ann.line.col, size = ann.line.size) +
  annotate(geom = "segment", x = 3.65, xend = 3.65, y = 22.7, yend = 27.3,
           color = ann.line.col, size = ann.line.size) +
  annotate(geom = "text", x = 3.8, y = 27.3, label = "Eutheria\nspecific",
           color = "black", size = 7/.pt, fontface = 2,
           hjust = 0, vjust = 1) +
  annotate(geom = "segment", x = 3.65, xend = 3.65, y = 27.7, yend = 39.3,
           color = ann.line.col, size = ann.line.size)

cowplot::ggsave2(
  filename = "results/fig_03_cytokines/fig_03_panel-A_cytokine-heatmap.pdf",
  ck.heat, width = 3.75, height = 6, units = "in"
)




