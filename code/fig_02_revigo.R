### packages ==================================================================
library(tidyverse)
library(scales)
library(ggrepel)

### data ======================================================================
# htseq counts for panal A
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

### tpm calculation for orthologs =============================================
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

### prepare gene lists for GO enrichment ======================================
# background
tpm %>% 
  select(mdom_gene_name) %>% 
  write.table(
    ., 
    "results/fig_02_go/genelist_background.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )

# opossum specific & 2x up in opossum NP to 13.5dpc
tpm %>% 
  filter(mdom_13.5dpc > 3 & 
           dnov_imp < 3 &
           ocun_imp < 3 &
           mdom_13.5dpc / mdom_NP > 2) %>% 
  select(mdom_gene_name) %>% 
  write.table(
    .,
    "results/fig_02_go/genelist_Md-over3tpm_DnOc-below3tpm_Md-2xUP.txt",
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )

# eutheria specific
tpm %>% 
  filter(mdom_13.5dpc < 3 &
           dnov_imp > 3 &
           ocun_imp > 3) %>% 
  select(mdom_gene_name) %>% 
  write.table(
    .,
    "results/fig_02_go/genelist_Md_below3tpm_DnOc-over3tpm.txt",
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )

# Performed GO enrichment on GOrilla GO enrichment tool online (http://cbl-gorilla.cs.technion.ac.il/), and the output was processed on revigo (http://revigo.irb.hr/). The code for plotting the dotplot below is from the revigo output, and was modified for aesthetic adjustments. 

### revigo plot: opossum specific =============================================
# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.
revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability", "show");
revigo.data <- rbind(c("GO:0010951","negative regulation of endopeptidase activity", 0.157, 6.381,-0.996, 4.304,-8.3893,0.675,0.000, 'yes'),
                     c("GO:0015711","organic anion transport", 1.192, 3.957, 5.071, 5.184,-8.8297,0.775,0.000, 'no'),
                     c("GO:0032501","multicellular organismal process", 2.373, 1.655,-1.146, 5.483,-4.6498,0.975,0.000, 'no'),
                     c("GO:0060577","pulmonary vein morphogenesis", 0.000,-1.825, 6.790, 1.756,-3.4510,0.875,0.000, 'no'),
                     c("GO:0071332","cellular response to fructose stimulus", 0.000, 1.366, 1.883, 1.653,-3.4510,0.939,0.017, 'no'),
                     c("GO:0071825","protein-lipid complex subunit organization", 0.010,-4.498, 5.650, 3.105,-4.7721,0.904,0.021, 'no'),
                     c("GO:0042737","drug catabolic process", 0.001,-0.746,-6.389, 2.286,-3.6364,0.895,0.034, 'no'),
                     c("GO:1901615","organic hydroxy compound metabolic process", 0.831, 0.615,-7.830, 5.028,-3.6289,0.943,0.040, 'no'),
                     c("GO:0006641","triglyceride metabolic process", 0.038,-4.416,-4.268, 3.687,-5.7852,0.704,0.041, 'yes'),
                     c("GO:0017144","drug metabolic process", 0.058,-0.187,-2.434, 3.868,-3.1002,0.941,0.043, 'no'),
                     c("GO:0006629","lipid metabolic process", 3.522,-5.590,-1.041, 5.655,-5.0580,0.838,0.140, 'no'),
                     c("GO:0008284","positive regulation of cell proliferation", 0.151, 6.035,-1.501, 4.288,-3.1226,0.794,0.209, 'no'),
                     c("GO:1901605","alpha-amino acid metabolic process", 3.625,-5.987, 0.094, 5.668,-3.3883,0.762,0.236, 'no'),
                     c("GO:0065008","regulation of biological quality", 3.395, 5.810,-2.525, 5.639,-4.9830,0.823,0.256, 'no'),
                     c("GO:0010817","regulation of hormone levels", 0.161, 5.194,-3.784, 4.314,-4.7100,0.795,0.257, 'yes'),
                     c("GO:0044281","small molecule metabolic process",15.138,-6.119,-1.168, 6.288,-6.1367,0.830,0.302, 'no'),
                     c("GO:0006730","one-carbon metabolic process", 0.328,-5.492, 0.767, 4.625,-3.2441,0.827,0.355, 'no'),
                     c("GO:0071827","plasma lipoprotein particle organization", 0.007,-3.234, 5.694, 2.945,-4.8601,0.753,0.388, 'yes'),
                     c("GO:0006811","ion transport", 5.344, 3.795, 5.483, 5.836,-5.0472,0.865,0.389, 'no'),
                     c("GO:0034367","macromolecular complex remodeling", 0.007,-4.138, 5.389, 2.969,-4.3382,0.905,0.390, 'no'),
                     c("GO:0050892","intestinal absorption", 0.006, 0.282, 6.343, 2.913,-4.4437,0.802,0.408, 'no'),
                     c("GO:1904640","response to methionine", 0.000, 0.825, 1.825, 1.398,-3.4510,0.954,0.409, 'no'),
                     c("GO:0070328","triglyceride homeostasis", 0.008, 4.332,-4.786, 3.013,-3.6198,0.803,0.412, 'yes'),
                     c("GO:0001523","retinoid metabolic process", 0.013,-4.801,-4.641, 3.234,-3.5287,0.760,0.418, 'yes'),
                     c("GO:0006638","neutral lipid metabolic process", 0.042,-4.928,-3.775, 3.730,-4.9872,0.757,0.449, 'no'),
                     c("GO:1901571","fatty acid derivative transport", 0.029, 2.896, 6.004, 3.564,-3.4202,0.815,0.453, 'yes'),
                     c("GO:0008202","steroid metabolic process", 0.161,-4.827,-4.093, 4.315,-4.4510,0.760,0.454, 'no'),
                     c("GO:0071347","cellular response to interleukin-1", 0.014, 0.255, 2.058, 3.260,-3.1643,0.935,0.466, 'yes'),
                     c("GO:0043267","negative regulation of potassium ion transport", 0.004, 5.484, 2.001, 2.762,-3.1637,0.750,0.467, 'no'),
                     c("GO:0036148","phosphatidylglycerol acyl-chain remodeling", 0.000,-3.479,-5.826, 1.255,-3.0114,0.794,0.476, 'no'),
                     c("GO:0050878","regulation of body fluid levels", 0.074, 5.042,-4.182, 3.976,-3.7520,0.804,0.479, 'no'),
                     c("GO:0036149","phosphatidylinositol acyl-chain remodeling", 0.000,-3.768,-5.863, 1.398,-3.4622,0.791,0.484, 'no'),
                     c("GO:0042445","hormone metabolic process", 0.090, 4.668,-4.046, 4.064,-5.0942,0.785,0.486, 'no'),
                     c("GO:0097272","ammonia homeostasis", 0.001, 4.028,-5.194, 2.049,-3.4510,0.827,0.492, 'no'),
                     c("GO:0003008","system process", 0.660,-2.101, 6.824, 4.928,-3.3565,0.893,0.500, 'no'),
                     c("GO:0006558","L-phenylalanine metabolic process", 0.075,-6.221, 1.462, 3.984,-3.4622,0.780,0.518, 'yes'),
                     c("GO:1902221","erythrose 4-phosphate/phosphoenolpyruvate family amino acid metabolic process", 0.075,-6.571, 1.455, 3.984,-3.4622,0.809,0.518, 'no'),
                     c("GO:0016042","lipid catabolic process", 0.401,-4.184,-3.956, 4.712,-3.8447,0.718,0.540, 'no'),
                     c("GO:0006869","lipid transport", 0.270, 3.325, 5.692, 4.539,-4.1343,0.787,0.551, 'no'),
                     c("GO:0048878","chemical homeostasis", 0.543, 4.938,-3.718, 4.843,-3.6498,0.769,0.561, 'no'),
                     c("GO:0071715","icosanoid transport", 0.029, 4.365, 5.038, 3.564,-3.4202,0.769,0.613, 'yes'),
                     c("GO:0046503","glycerolipid catabolic process", 0.018,-3.801,-4.733, 3.359,-4.4437,0.711,0.617, 'no'),
                     c("GO:0055081","anion homeostasis", 0.045, 4.590,-4.511, 3.759,-3.2757,0.792,0.623, 'no'),
                     c("GO:0015849","organic acid transport", 1.024, 3.124, 5.413, 5.119,-6.4737,0.766,0.633, 'no'),
                     c("GO:0009072","aromatic amino acid family metabolic process", 0.719,-6.240, 0.530, 4.965,-3.2441,0.787,0.660, 'yes'),
                     c("GO:0006820","anion transport", 1.956, 4.479, 4.732, 5.400,-7.2596,0.823,0.663, 'no'),
                     c("GO:0030162","regulation of proteolysis", 0.304, 6.392,-1.841, 4.591,-3.4559,0.747,0.669, 'no'),
                     c("GO:0006570","tyrosine metabolic process", 0.059,-6.206, 1.731, 3.881,-3.2924,0.801,0.670, 'yes'));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below
one.data <- subset(one.data,)
one.data <- one.data[order(one.data$log10_p_value, decreasing = T),] # so more significant terms get plotted on top

# plot 
p1 <- ggplot(data = one.data) +
  geom_point(aes(plot_X, plot_Y, 
                 colour = log10_p_value, 
                 size = plot_size), alpha = I(0.8)) +
  scale_size(range = c(2, 15), guide = "none") +
  scale_colour_gradientn(
    name = "log10 (p value)",
    colours = rev(c("blue", "green", "yellow", 'red')), 
    limits = c(min(one.data$log10_p_value), 0)) +
  labs (y = "semantic space x", x = "semantic space y")
  
ex <- one.data[one.data$show %in% c('yes'), ]

p1 <- p1 + geom_text_repel(data = ex, 
                           aes(plot_X, plot_Y, label = description), 
                           colour = I(alpha("black", 1)), size = 6/.pt)

one.x_range <- max(one.data$plot_X) - min(one.data$plot_X);
one.y_range <- max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10)
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10)

p1 <- p1 + theme_classic() +
  theme(
    plot.margin = margin(6, 30, 6, 6),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.background = element_blank(),
    legend.position = c(0.9, 0), 
    legend.justification = c(0, 0),
    legend.key.size = unit(0.5, "lines")
  )

cowplot::ggsave2(
  filename = "results/fig_02_go/fig02_panel-A_reviGO_opossum-specific.pdf",
  p1, 
  width = 3.5, height = 3.5, units = "in"
)

### revigo plot: eutheria specific ============================================
# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.
revigo.names2 <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability", 'show');
revigo.data2 <- rbind(c("GO:0002376","immune system process", 0.600, 0.586,-7.055, 4.886,-3.6517,0.872,0.000, 'yes'),
                     c("GO:0042110","T cell activation", 0.081, 0.095, 5.830, 4.018,-3.9318,0.338,0.000, 'yes'),
                     c("GO:0007256","activation of JNKK activity", 0.001, 4.437,-1.080, 2.260,-3.2069,0.533,0.102, 'no'),
                     c("GO:0031023","microtubule organizing center organization", 0.061,-4.303,-2.060, 3.896,-3.6234,0.690,0.129, 'no'),
                     c("GO:0072330","monocarboxylic acid biosynthetic process", 0.940,-4.008, 2.230, 5.081,-3.8508,0.683,0.159, 'yes'),
                     c("GO:0006928","movement of cell or subcellular component", 0.973,-3.870, 4.107, 5.097,-3.2865,0.683,0.201, 'no'),
                     c("GO:0007166","cell surface receptor signaling pathway", 0.920, 6.298, 0.747, 5.072,-3.0119,0.513,0.320, 'no'),
                     c("GO:0050808","synapse organization", 0.070,-4.936,-0.962, 3.956,-3.2807,0.689,0.332, 'no'),
                     c("GO:0019221","cytokine-mediated signaling pathway", 0.093, 5.995,-0.203, 4.078,-3.1904,0.504,0.435, 'yes'),
                     c("GO:0002577","regulation of antigen processing and presentation", 0.003, 3.574, 5.224, 2.609,-3.1561,0.430,0.582, 'yes'),
                     c("GO:0090025","regulation of monocyte chemotaxis", 0.005, 3.504, 2.697, 2.763,-3.9788,0.287,0.594, 'yes'),
                     c("GO:0007165","signal transduction", 6.621, 5.421, 0.798, 5.929,-4.2190,0.474,0.679, 'no'));

one.data2 <- data.frame(revigo.data2)
names(one.data2) <- revigo.names2
one.data2 <- one.data2[which(one.data2$plot_X != "null" & one.data2$plot_Y != "null"), ]
one.data2$plot_X <- as.numeric( as.character(one.data2$plot_X) )
one.data2$plot_Y <- as.numeric( as.character(one.data2$plot_Y) )
one.data2$plot_size <- as.numeric( as.character(one.data2$plot_size) )
one.data2$log10_p_value <- as.numeric( as.character(one.data2$log10_p_value) )
one.data2$frequency <- as.numeric( as.character(one.data2$frequency) )
one.data2$uniqueness <- as.numeric( as.character(one.data2$uniqueness) )
one.data2$dispensability <- as.numeric( as.character(one.data2$dispensability) )

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below
one.data2 <- subset(one.data2,)
one.data2 <- one.data2[order(one.data2$log10_p_value, decreasing = T), ]

# plot
p2 <- ggplot(data = one.data2) +
  geom_point(aes(plot_X, plot_Y, 
                 colour = log10_p_value, 
                 size = plot_size), alpha = I(0.8)) +
  scale_size(range = c(2, 15), guide = "none") +
  scale_colour_gradientn(
    name = "log10 (p value)",
    colours = rev(c("blue", "green", "yellow", 'red')), 
    limits = c(-8.8297, 0)) + # min pval is takem from Mdover3_DnOcbelow3_Md2xUP to match the legend because both plots will be shown in the same figure
  labs (y = "semantic space x", x = "semantic space y")

ex2 <- one.data2[one.data2$show %in% c('yes'), ]

p2 <- p2 + geom_text_repel(data = ex2, 
                           aes(plot_X, plot_Y, label = description), 
                           colour = I(alpha("black", 1)), size = 6/.pt)

one.x_range2 <- max(one.data2$plot_X) - min(one.data2$plot_X);
one.y_range2 <- max(one.data2$plot_Y) - min(one.data2$plot_Y);
p2 <- p2 + xlim(min(one.data2$plot_X)-one.x_range2/10,max(one.data2$plot_X)+one.x_range2/10)
p2 <- p2 + ylim(min(one.data2$plot_Y)-one.y_range2/10,max(one.data2$plot_Y)+one.y_range2/10)

p2 <- p2 + theme_classic() +
  theme(
    plot.margin = margin(6, 30, 6, 6),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.background = element_blank(),
    legend.position = c(0.9, 0), 
    legend.justification = c(0, 0),
    legend.key.size = unit(0.5, "lines")
  )

cowplot::ggsave2(
  filename = "results/fig_02_go/fig02_panel-B_reviGO_eutheria-specific.pdf",
  p2, 
  width = 3.5, height = 3.5, units = "in"
)

### end =======================================================================



