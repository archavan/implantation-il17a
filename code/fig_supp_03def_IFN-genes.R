### packages ==================================================================
library(tidyverse)
library(edgeR)
library(ggrepel)
library(ggtext)

### data ======================================================================
cpm.tcell <- read.csv("results/fig_supp_03_tcell-rna-seq/edgeR_hsap_tcells/hsap_tcells_cpm_protein_coding_with-DE-data.csv")
cpm.esf <- read.csv("results/fig_supp_03_tcell-rna-seq/edgeR_hsap_esf-dsc/hsap_esf-dsc_cpm_with-DE-data.csv")

cpm <- left_join(x = cpm.tcell, 
                 y = cpm.esf, 
                 by = c("hsap_ensembl_gid", "hsap_gene_name"), 
                 suffix = c(".tcell", ".esf"))

### genomic data ==============================================================
GO0060337 <- read.csv("data/genomic-data/hsap_GO0060337_type1InterferonSignalingPathway_GRCh37biomart.csv") # type1 interferon signaling pathway

### Panel D: heatmap ==========================================================
# heatmap of gene expression values for type 1 IFN genes and their receptors
p.dat <- cpm %>% 
  filter(hsap_ensembl_gid %in% GO0060337$Gene.stable.ID) %>% 
  .[grep("IFN", .$hsap_gene_name), ] %>% 
  select(hsap_ensembl_gid, 
         hsap_gene_name, 
         mean_control, mean_conditioned, 
         mean_undiff, mean_diff) %>% 
  pivot_longer(cols = 3:6, 
               names_to = "sample", 
               values_to = "cpm", 
               names_prefix = "mean_") %>% 
  mutate(celltype = case_when(
    sample %in% c("control", "conditioned") ~ "T cells",
    sample %in% c("undiff", "diff") ~ "Endometrial Stromal Cells"
  )) %>% 
  mutate(sample = case_when(
    sample == "control" ~ "Control",
    sample == "conditioned" ~ "Conditioned",
    sample == "undiff" ~ "ESF",
    sample == "diff" ~ "DSC"
  )) %>% 
  mutate(sample = factor(
    sample, 
    levels = c("Control", "Conditioned", "ESF", "DSC")
  )) %>% 
  mutate(celltype = factor(
    celltype, 
    levels = c("T cells", "Endometrial Stromal Cells")
  )) %>% 
  mutate(hsap_gene_name = factor(
    hsap_gene_name,
    levels = rev(c("IFNAR2", "IFNAR1", "IFNB1", "IFNA21", "IFNA17", "IFNA16", "IFNA14", "IFNA10", "IFNA8", "IFNA7", "IFNA6", "IFNA5", "IFNA4", "IFNA2", "IFNA1"))
  ))

pan.d <- ggplot(p.dat, aes(sample, hsap_gene_name, label = round(cpm, 1))) +
  geom_tile(aes(fill = cpm), color = "white") +
  facet_grid(. ~ celltype, space = "free", scales = "free") +
  scale_fill_gradient(low = "white", high = "red", name = "CPM") +
  geom_text(size = 6/.pt) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.length.y = unit(0.0, "lines"),
    axis.text.x = element_text(size = 6, color = "black", face = 2),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title = element_blank(),
    strip.text = element_text(size = 7, color = "black", face = 2),
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    plot.margin = margin(6, 30, 6, 6),
    legend.position = c(1, 0),
    legend.justification = c(0, 0)
  )

ggsave(
  filename = "results/fig_supp_03_tcell-rna-seq/fig_supp_03_panel-D_IFN-heatmap.pdf",
  pan.d, 
  width = 3.5, height = 3.5, units = "in"
  )

### Panel E: OAS genes scatterplot ============================================
oas <- cpm.tcell %>% 
  filter(grepl("^OAS", .$hsap_gene_name)) %>% 
  mutate(label = ifelse(
    hsap_gene_name %in% c("OAS1", "OAS2", "OAS3"), hsap_gene_name, NA
  )) %>% 
  ggplot(aes(sqrt(mean_control), sqrt(mean_conditioned), label = label)) +
  geom_point(aes(color = FDR < 10^-6), shape = 19, size  = 2) +
  scale_color_manual(values = c('grey', 'red'), 
                     breaks = c(FALSE, TRUE), 
                     labels = c("p > 10<sup>-6</sup>", "p < 10<sup>-6</sup>")) +
  scale_x_continuous(limits = c(0,20)) +
  scale_y_continuous(limits = c(0,20)) +
  geom_text_repel(size = 7/.pt) +
  geom_abline(slope = 1, intercept = 0, size = 0.25, colour='grey') + 
  labs(x = "sqrt(CPM<sub>control</sub>)", 
       y = "sqrt(CPM<sub>conditioned</sub>)") +
  theme_classic() +
  theme(
    axis.line = element_line(size = 0.25, color = "black"),
    axis.ticks = element_line(size = 0.25, color = "black"),
    axis.title.x = element_markdown(size = 7),
    axis.title.y = element_markdown(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 7),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank()
  )

cowplot::ggsave2(
  filename = "results/fig_supp_03_tcell-rna-seq/fig_supp_03_panel-E_OAS.pdf",
  oas, 
  width = 3, height = 3, units = "in"
)

### Panel F: EIF genes scatterplot ============================================
eif <- cpm.tcell %>% 
  filter(grepl("^EIF", .$hsap_gene_name)) %>% 
  mutate(label = ifelse(
    hsap_gene_name %in% c("EIF1", "EIF4B", "EIF4A2", "EIF2A", "EIF2B2", "EIF4G3", "EIF2AK2", "EIFEBP2"), 
    hsap_gene_name, NA
  )) %>% 
  ggplot(aes(sqrt(mean_control), sqrt(mean_conditioned), label = label)) +
  geom_point(aes(color = FDR < 10^-6), shape = 19, size  = 2) +
  scale_color_manual(values = c('grey', 'red'), 
                     breaks = c(FALSE, TRUE), 
                     labels = c("p > 10<sup>-6</sup>", "p < 10<sup>-6</sup>")) +
  scale_x_continuous(limits = c(0,45)) +
  scale_y_continuous(limits = c(0,45)) +
  geom_text_repel(size = 7/.pt) +
  geom_abline(slope = 1, intercept = 0, size = 0.25, colour='grey') + 
  labs(x = "sqrt(CPM<sub>control</sub>)", 
       y = "sqrt(CPM<sub>conditioned</sub>)") +
  theme_classic() +
  theme(
    axis.line = element_line(size = 0.25, color = "black"),
    axis.ticks = element_line(size = 0.25, color = "black"),
    axis.title.x = element_markdown(size = 7),
    axis.title.y = element_markdown(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    legend.title = element_blank(),
    legend.text = element_markdown(size = 7),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank()
  )

cowplot::ggsave2(
  filename = "results/fig_supp_03_tcell-rna-seq/fig_supp_03_panel-F_EIF.pdf",
  eif, 
  width = 3, height = 3, units = "in"
)

### end =======================================================================
