### packages ==================================================================
library(tidyverse)
library(cowplot)
library(scales)
library(ggrepel)

### data ======================================================================
dn.ht <- read.csv("data/rna-seq/read-counts_dnov_uterus.csv") # armadillo htseq counts
dn.fl <- read.csv("data/genomic-data/feature-lengths_dnov.csv") # armadillo feature lengths

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
        axis.ticks.length.y = unit(2, "pt"),
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

### GO enrichment =============================================================
# gene lists for GO enrichment
tpm %>%  # over 3 tpm in imp and upregulated 10x from NP to imp
  filter(dnov_imp > 3 & 
           (dnov_imp / dnov_NP) > 10) %>% 
  select(dnov_gene_name) %>% 
  write.table(., file = "results/fig_01_inflammation/genelist_dnov_imp-over3tpm_imp-up-10x.txt",
              quote = FALSE, row.names = FALSE, col.names = TRUE)

tpm %>% 
  select(dnov_gene_name) %>% 
  write.table(., file = "results/fig_01_inflammation/genelist_background.txt",
              quote = FALSE, row.names = FALSE, col.names = TRUE)

# Performed GO enrichment on GOrilla GO enrichment tool online (http://cbl-gorilla.cs.technion.ac.il/), and the output was processed on revigo (http://revigo.irb.hr/). The code for plotting the dotplot below is from the revigo output, and was modified for aesthetic adjustments. 

### REVIGO plot ===============================================================
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability", "show");
revigo.data <- rbind(c("GO:0002376","immune system process", 0.600, 3.087,-2.298, 4.886,-5.7447,0.987,0.000, "yes"),
                     c("GO:0005975","carbohydrate metabolic process", 5.260, 1.496,-2.561, 5.829,-4.1198,0.966,0.000,"no"),
                     c("GO:0022603","regulation of anatomical structure morphogenesis", 0.830,-6.976, 0.421, 5.027,-6.2403,0.647,0.000,"no"),
                     c("GO:0022617","extracellular matrix disassembly", 0.007,-0.948, 2.924, 2.939,-6.4248,0.849,0.000,"yes"),
                     c("GO:0040011","locomotion", 0.997, 2.284,-4.103, 5.107,-5.1018,0.987,0.000,"no"),
                     c("GO:0042221","response to chemical", 3.071,-0.390,-6.732, 5.595,-6.5528,0.806,0.000,"no"),
                     c("GO:0051179","localization",18.495, 1.630,-1.151, 6.375,-3.2798,0.989,0.000,"no"),
                     c("GO:0030574","collagen catabolic process", 0.005,-4.046, 6.178, 2.836,-5.4134,0.795,0.042,"yes"),
                     c("GO:0030203","glycosaminoglycan metabolic process", 0.805, 5.156,-3.079, 5.014,-3.9788,0.934,0.077,"no"),
                     c("GO:0005996","monosaccharide metabolic process", 0.817, 3.308,-0.116, 5.021,-4.4012,0.880,0.107,"no"),
                     c("GO:0050900","leukocyte migration", 0.060, 3.034, 3.393, 3.890,-6.2248,0.745,0.110,"yes"),
                     c("GO:0001775","cell activation", 0.171, 1.513, 5.302, 4.341,-4.0894,0.871,0.136,"no"),
                     c("GO:0006721","terpenoid metabolic process", 0.260, 0.583, 4.776, 4.524,-4.0128,0.811,0.152,"no"),
                     c("GO:0006928","movement of cell or subcellular component", 0.973, 0.742, 4.538, 5.097,-3.3170,0.855,0.177,"no"),
                     c("GO:0032350","regulation of hormone metabolic process", 0.007,-6.473,-3.720, 2.929,-4.0195,0.793,0.184,"yes"),
                     c("GO:0016192","vesicle-mediated transport", 1.085, 4.349, 2.054, 5.144,-6.1013,0.927,0.231,"yes"),
                     c("GO:0030155","regulation of cell adhesion", 0.129,-7.426,-3.765, 4.219,-3.1203,0.777,0.231,"no"),
                     c("GO:0006898","receptor-mediated endocytosis", 0.095, 5.176, 1.470, 4.086,-3.9031,0.927,0.249,"no"),
                     c("GO:0042127","regulation of cell proliferation", 0.313,-7.386,-2.790, 4.603,-3.2503,0.724,0.251,"no"),
                     c("GO:0042981","regulation of apoptotic process", 0.313,-6.492,-0.062, 4.604,-5.6840,0.675,0.251,"yes"),
                     c("GO:0051174","regulation of phosphorus metabolic process", 0.580,-7.166,-1.553, 4.872,-3.1694,0.719,0.266,"no"),
                     c("GO:0070528","protein kinase C signaling", 0.007,-3.693,-5.386, 2.967,-4.3635,0.692,0.284,"no"),
                     c("GO:0043062","extracellular structure organization", 0.061,-0.869, 4.053, 3.894,-4.8327,0.838,0.286,"no"),
                     c("GO:0032879","regulation of localization", 0.726,-4.989,-1.240, 4.969,-4.2336,0.710,0.291,"no"),
                     c("GO:0048518","positive regulation of biological process", 1.744,-7.886,-1.660, 5.350,-5.8794,0.728,0.299,"no"),
                     c("GO:0032940","secretion by cell", 0.763, 2.676, 4.078, 4.991,-3.4763,0.818,0.304,"no"),
                     c("GO:0071526","semaphorin-plexin signaling pathway", 0.019,-3.912,-5.106, 3.381,-4.4976,0.676,0.306,"no"),
                     c("GO:0050790","regulation of catalytic activity", 1.575,-7.464,-2.304, 5.306,-3.0660,0.729,0.307,"no"),
                     c("GO:0065009","regulation of molecular function", 1.726,-7.275,-1.895, 5.345,-3.0575,0.739,0.311,"no"),
                     c("GO:0006066","alcohol metabolic process", 0.422, 5.111,-0.990, 4.734,-4.2832,0.885,0.312,"no"),
                     c("GO:0048522","positive regulation of cellular process", 1.585,-7.631,-1.098, 5.308,-5.1451,0.601,0.323,"no"),
                     c("GO:0034330","cell junction organization", 0.056,-0.868, 4.807, 3.855,-3.2815,0.839,0.327,"no"),
                     c("GO:0048523","negative regulation of cellular process", 1.830,-7.668,-2.226, 5.371,-3.5513,0.675,0.333,"no"),
                     c("GO:0048519","negative regulation of biological process", 1.984,-7.805,-1.873, 5.406,-3.1068,0.725,0.334,"no"),
                     c("GO:0070482","response to oxygen levels", 0.055, 0.608,-7.252, 3.849,-3.1864,0.855,0.336,"no"),
                     c("GO:1901700","response to oxygen-containing compound", 0.503,-0.219,-7.160, 4.810,-5.7055,0.794,0.421,"no"),
                     c("GO:0051128","regulation of cellular component organization", 1.586,-5.973,-1.322, 5.308,-3.3820,0.674,0.422,"no"),
                     c("GO:0006952","defense response", 0.568,-0.042,-6.766, 4.863,-4.0985,0.815,0.427,"no"),
                     c("GO:0009966","regulation of signal transduction", 0.857,-4.198,-4.307, 5.041,-6.4828,0.532,0.448,"no"),
                     c("GO:0006810","transport",17.616, 4.795, 2.415, 6.354,-3.9547,0.908,0.454,"no"),
                     c("GO:0048583","regulation of response to stimulus", 1.120,-4.384,-4.307, 5.158,-6.1506,0.637,0.463,"yes"),
                     c("GO:0046697","decidualization", 0.003,-4.861, 5.609, 2.539,-3.4353,0.796,0.467,"no"),
                     c("GO:0051239","regulation of multicellular organismal process", 0.628,-6.355, 1.494, 4.906,-4.4101,0.640,0.490,"no"),
                     c("GO:0006954","inflammatory response", 0.110, 0.259,-6.604, 4.151,-4.0585,0.827,0.491,"no"),
                     c("GO:0019216","regulation of lipid metabolic process", 0.095,-5.808,-0.730, 4.086,-3.9957,0.703,0.500,"no"),
                     c("GO:0034103","regulation of tissue remodeling", 0.011,-6.296, 2.628, 3.152,-3.0386,0.689,0.515,"no"),
                     c("GO:0010646","regulation of cell communication", 0.929,-6.960,-3.082, 5.076,-6.1409,0.697,0.517,"no"),
                     c("GO:0044236","multicellular organism metabolic process", 0.016,-4.045, 6.013, 3.315,-4.4045,0.810,0.529,"no"),
                     c("GO:0023051","regulation of signaling", 0.934,-7.039,-3.357, 5.079,-5.9031,0.704,0.532,"no"),
                     c("GO:0006022","aminoglycan metabolic process", 0.883, 5.292,-3.313, 5.054,-3.8894,0.943,0.542,"no"),
                     c("GO:0006665","sphingolipid metabolic process", 0.093, 1.480, 3.740, 4.077,-3.2534,0.819,0.543,"no"),
                     c("GO:0048771","tissue remodeling", 0.028,-4.393, 6.064, 3.557,-3.2233,0.815,0.551,"no"),
                     c("GO:0046903","secretion", 0.810, 4.611, 2.120, 5.017,-3.5867,0.862,0.557,"no"),
                     c("GO:0001958","endochondral ossification", 0.005,-4.804, 5.716, 2.836,-3.1013,0.783,0.562,"no"),
                     c("GO:0071363","cellular response to growth factor stimulus", 0.147, 0.062,-6.980, 4.276,-4.2388,0.769,0.572,"no"),
                     c("GO:0045603","positive regulation of endothelial cell differentiation", 0.003,-5.947, 2.653, 2.599,-3.1079,0.589,0.575,"no"),
                     c("GO:0048568","embryonic organ development", 0.113,-4.987, 5.600, 4.160,-3.5171,0.751,0.618,"no"),
                     c("GO:0007169","transmembrane receptor protein tyrosine kinase signaling pathway", 0.167,-4.065,-4.662, 4.331,-3.6990,0.627,0.620,"no"),
                     c("GO:0006720","isoprenoid metabolic process", 0.464, 0.574, 4.978, 4.775,-3.2291,0.810,0.620,"no"),
                     c("GO:0071248","cellular response to metal ion", 0.031,-0.162,-7.320, 3.601,-3.1475,0.801,0.624,"no"),
                     c("GO:0030516","regulation of axon extension", 0.017,-5.087, 2.356, 3.327,-5.1656,0.561,0.647,"no"),
                     c("GO:0001817","regulation of cytokine production", 0.108,-6.272, 1.957, 4.141,-3.4365,0.643,0.648,"yes"),
                     c("GO:0045834","positive regulation of lipid metabolic process", 0.031,-6.085,-0.162, 3.596,-4.3925,0.633,0.651,"no"),
                     c("GO:0032496","response to lipopolysaccharide", 0.043, 0.426,-6.569, 3.745,-3.5467,0.802,0.654,"no"),
                     c("GO:0031344","regulation of cell projection organization", 0.123,-4.713,-0.010, 4.199,-3.0783,0.669,0.658,"no"),
                     c("GO:0002274","myeloid leukocyte activation", 0.033, 3.021, 4.347, 3.622,-4.1146,0.851,0.659,"yes"),
                     c("GO:0002526","acute inflammatory response", 0.016, 1.098,-6.709, 3.318,-3.2534,0.846,0.662,"yes"),
                     c("GO:0051270","regulation of cellular component movement", 0.169,-3.496,-0.184, 4.337,-4.9788,0.653,0.665,"no"),
                     c("GO:0040012","regulation of locomotion", 0.183,-5.816,-2.579, 4.372,-4.1831,0.703,0.665,"no"),
                     c("GO:0021785","branchiomotor neuron axon guidance", 0.002,-1.979,-1.262, 2.318,-3.2660,0.584,0.682,"no"),
                     c("GO:0070887","cellular response to chemical stimulus", 1.007,-0.577,-7.018, 5.111,-5.5498,0.756,0.683,"no"),
                     c("GO:0006897","endocytosis", 0.235, 5.108, 1.686, 4.480,-3.1733,0.922,0.688,"no"),
                     c("GO:1903510","mucopolysaccharide metabolic process", 0.016, 5.058,-3.070, 3.323,-3.2291,0.945,0.693,"no"),
                     c("GO:0008203","cholesterol metabolic process", 0.028, 2.549, 0.812, 3.554,-3.1586,0.865,0.696,"yes"),
                     c("GO:0001525","angiogenesis", 0.096,-4.935, 5.444, 4.088,-5.0521,0.752,0.700, "yes"));

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
#head(one.data);

one.data <- subset(one.data,)
one.data <- one.data[order(one.data$log10_p_value, decreasing = T), ] # so more significant terms get plotted on top

# plot
p1 <- ggplot(data = one.data)
p1 <- p1 + 
  geom_point(aes(plot_X, plot_Y, colour = log10_p_value, size = plot_size), 
             alpha = I(0.8) ) + 
  scale_colour_gradientn(
    name = "log10 (p value)",
    colours = rev(c("blue", "green", "yellow", 'red')), 
    limits = c(min(one.data$log10_p_value), 0)) +
  scale_size(range = c(2, 15), guide = FALSE) +
  labs (y = "semantic space x", x = "semantic space y")

ex <- one.data [one.data$show %in% c('yes'), ]

p1 <- p1 + geom_text_repel(data = ex, 
                           aes(plot_X, plot_Y, label = description), 
                           colour = I(alpha("black", 1)), size = 6/.pt)


one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);

p1 <- p1 + theme_classic() +
  theme(
    plot.margin = margin(6, 30, 6, 6),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.background = element_blank(),
    legend.position = c(0.9, 0), 
    legend.justification = c(0, 0),
    legend.key.size = unit(0.5, "lines")
    )

cowplot::ggsave2(
  filename = "results/fig_01_inflammation/fig01_panel-C_dnov-reviGO.pdf",
  p1, 
  width = 3.5, height = 3.5, units = "in"
)

### end  ======================================================================

