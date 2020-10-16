### pacakges ==================================================================
library(tidyverse)

### data ======================================================================
dec <- read.csv("data/hsap_single-cell/decidua_cpm_mean.csv")
vil <- read.csv("data/hsap_single-cell/villi_cpm_mean.csv")

### prepare data for plotting =================================================
# join
dat <- inner_join(dec, vil, by = "hsap_gene_name")

# subset
geneorder <- c("IL17A", "TNF", "CCL20", "CSF3", 
               "CXCL9", "LTA",
               "IL10", "IL1A", 
               "IL20", "IL17B", "FGF2", "WNT2", "CCL25", "CCL17",    
               "INHBA", "CCL21", "IL17D","IL16","CCL19", "TNFSF15", "KITLG", "CTF1",  
               "WNT7A", "LIF", "NDP", "CCL5", "CXCL10", 
               "CXCL12", "TNFSF10", "CSF1", "TNFSF13", "BMP4", "CX3CL1", "WNT5A", "EDN1", "IL15", "LTB", "IL8", "IL1B") # Same gene order as in Fig 3A.

dat <- dat[which(dat$hsap_gene_name %in% geneorder), ]

# melt
dat.long <- pivot_longer(data = dat, 
                         cols = c("decidua", "villi"), 
                         names_to = "tissue", 
                         values_to = "cpm")

dat.long$hsap_gene_name <- factor(dat.long$hsap_gene_name, levels = geneorder)

### plot ======================================================================
p <- ggplot(dat = dat.long,
            aes(x = tissue, y = hsap_gene_name, label = round(cpm, 1))) +
  geom_tile(aes(fill = sqrt(cpm)), color = "white") +
  scale_fill_gradient2(low='#FFFFFF', mid = '#FF6666', high = '#FF0000', 
                       midpoint = 22, 
                       guide = 'colourbar',
                       name = "CPM",
                       breaks = c(0, 10, 20, 30),
                       labels = c(0, 100, 400, 900)) +
  geom_text(size = 6/.pt) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 6, colour = 'black'),
    axis.text.y = element_text(size = 6, colour = 'black'),
    axis.title.x = element_text(size = 7, colour = 'black'),
    axis.ticks = element_line(size = 0.25, colour = "black"),
    legend.key.size = unit(0.75, "lines"),
    legend.title = element_text(size = 6, colour = "black"),
    legend.text = element_text(size = 6, colour = "black")
  ) 

ggsave(
  filename = "results/fig_supp_01_hsap-single-cell/fig_supp_01_hsap-cytokines.pdf",
  p,
  width = 3, height = 6, units = "in"
)
  
### end =======================================================================
  
