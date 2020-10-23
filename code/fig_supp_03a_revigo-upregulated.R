### packages ==================================================================
library(tidyverse)
library(scales)
library(ggrepel)

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability", "show");
revigo.data <- rbind(c("GO:0002295","T-helper cell lineage commitment", 0.002, 5.714,-2.988, 2.470,-3.4449,0.705,0.000, "yes"),
                     c("GO:0002376","immune system process", 0.600,-1.177,-2.179, 4.886,-4.7520,0.990,0.000, "no"),
                     c("GO:0006807","nitrogen compound metabolic process",38.744,-0.461,-6.000, 6.696,-4.5214,0.962,0.000, "no"),
                     c("GO:0009987","cellular process",63.780,-1.287,-3.306, 6.913,-3.7545,0.996,0.000, "no"),
                     c("GO:0016197","endosomal transport", 0.131, 0.686, 2.165, 4.225,-4.3487,0.967,0.000, "yes"),
                     c("GO:0032502","developmental process", 2.812, 0.698,-5.383, 5.557,-3.0565,0.990,0.000, "no"),
                     c("GO:0051704","multi-organism process", 0.751,-0.894,-6.745, 4.984,-3.1568,0.990,0.000, "no"),
                     #c("GO:0065007","biological regulation",20.498,-0.591,-1.504, 6.420,-8.1675,0.992,0.000, "no"),
                     c("GO:0071840","cellular component organization or biogenesis", 8.568,-1.126,-7.206, 6.041,-4.3893,0.991,0.000, "no"),
                     c("GO:0048525","negative regulation of viral process", 0.015, 4.872, 6.549, 3.291,-6.5421,0.688,0.017, "yes"),
                     c("GO:0006470","protein dephosphorylation", 0.585,-3.653,-0.267, 4.875,-5.0851,0.863,0.023, "yes"),
                     c("GO:0008283","cell proliferation", 0.394, 0.621,-4.422, 4.704,-3.7905,0.963,0.050, "no"),
                     c("GO:0009058","biosynthetic process",31.611, 0.164,-3.450, 6.608,-3.7545,0.963,0.066, "no"),
                     c("GO:0044238","primary metabolic process",53.743,-2.135,-5.373, 6.839,-4.3737,0.959,0.089, "yes"),
                     c("GO:0000380","alternative mRNA splicing, via spliceosome", 0.021,-3.686, 1.889, 3.425,-3.3706,0.879,0.103, "no"),
                     c("GO:0071704","organic substance metabolic process",58.357,-2.564,-4.492, 6.874,-4.1656,0.958,0.120, "no"),
                     c("GO:0006928","movement of cell or subcellular component", 0.973, 1.436,-1.314, 5.097,-3.7258,0.918,0.122, "no"),
                     c("GO:0000910","cytokinesis", 0.315, 1.290,-2.268, 4.606,-3.5391,0.925,0.180, "no"),
                     c("GO:0044260","cellular macromolecule metabolic process",34.276,-5.459, 2.108, 6.643,-9.0691,0.821,0.181, "no"),
                     c("GO:0010557","positive regulation of macromolecule biosynthetic process", 0.633, 2.467, 5.538, 4.909,-5.9586,0.525,0.191, "no"),
                     c("GO:0044237","cellular metabolic process",53.061,-5.907,-1.599, 6.833,-4.3197,0.907,0.191, "no"),
                     c("GO:0006793","phosphorus metabolic process",13.507,-6.519,-1.547, 6.239,-5.6345,0.901,0.194, "no"),
                     c("GO:1901360","organic cyclic compound metabolic process",30.324,-4.220,-1.705, 6.590,-4.2197,0.927,0.198, "no"),
                     c("GO:0043170","macromolecule metabolic process",39.491,-4.786,-1.419, 6.705,-6.1785,0.923,0.224, "no"),
                     c("GO:0018130","heterocycle biosynthetic process",17.388,-5.505, 3.318, 6.348,-8.0000,0.798,0.231, "no"),
                     c("GO:0044093","positive regulation of molecular function", 0.890, 5.612, 4.856, 5.058,-4.3936,0.700,0.258, "no"),
                     c("GO:0006725","cellular aromatic compound metabolic process",29.628,-6.000, 0.445, 6.580,-4.9747,0.885,0.260, "no"),
                     c("GO:0046483","heterocycle metabolic process",29.664,-6.337, 0.147, 6.580,-4.7773,0.885,0.260, "no"),
                     c("GO:0008104","protein localization", 2.626, 0.771, 2.758, 5.527,-3.7520,0.951,0.268, "no"),
                     c("GO:0006950","response to stress", 4.575, 4.950,-4.690, 5.769,-5.5258,0.829,0.268, "no"),
                     c("GO:0034641","cellular nitrogen compound metabolic process",34.137,-5.936, 1.147, 6.641,-4.0083,0.840,0.277, "no"),
                     c("GO:0048519","negative regulation of biological process", 1.984, 6.220, 4.720, 5.406,-5.8928,0.684,0.295, "no"),
                     c("GO:0065009","regulation of molecular function", 1.726, 5.762, 4.441, 5.345,-3.6576,0.701,0.316, "no"),
                     c("GO:0008630","intrinsic apoptotic signaling pathway in response to DNA damage", 0.021, 6.062,-0.839, 3.434,-3.8508,0.673,0.320, "no"),
                     c("GO:0051128","regulation of cellular component organization", 1.586, 4.198, 4.468, 5.308,-5.7258,0.648,0.328, "no"),
                     c("GO:0043412","macromolecule modification", 9.785,-3.103, 4.940, 6.099,-3.9469,0.895,0.331, "no"),
                     c("GO:0048518","positive regulation of biological process", 1.744, 6.077, 4.923, 5.350,-3.2757,0.688,0.332, "no"),
                     c("GO:0034058","endosomal vesicle fusion", 0.006,-0.483, 3.629, 2.877,-4.5157,0.910,0.351, "no"),
                     c("GO:0016311","dephosphorylation", 1.250,-4.288,-4.445, 5.205,-4.0953,0.910,0.371, "no"),
                     c("GO:0007259","JAK-STAT cascade", 0.033, 6.573,-0.477, 3.630,-3.7570,0.683,0.375, "yes"),
                     c("GO:0033036","macromolecule localization", 3.030, 0.837, 3.863, 5.590,-3.5784,0.961,0.385, "no"),
                     c("GO:0097696","STAT cascade", 0.033, 6.446,-0.538, 3.631,-3.7570,0.683,0.386, "no"),
                     c("GO:0006978","DNA damage response, signal transduction by p53 class mediator resulting in transcription of p21 class mediator", 0.003, 6.313,-1.359, 2.621,-3.3706,0.705,0.389, "no"),
                     c("GO:0009607","response to biotic stimulus", 0.342, 4.368,-4.979, 4.643,-3.8928,0.861,0.421, "no"),
                     c("GO:0006357","regulation of transcription from RNA polymerase II promoter", 1.273, 1.412, 5.275, 5.213,-3.9626,0.568,0.422, "no"),
                     c("GO:0006952","defense response", 0.568, 5.022,-4.393, 4.863,-4.0691,0.830,0.447, "no"),
                     c("GO:0001817","regulation of cytokine production", 0.108, 5.507, 3.505, 4.141,-6.2441,0.686,0.449, "yes"),
                     c("GO:0016050","vesicle organization", 0.130,-1.236, 1.817, 4.221,-3.0491,0.935,0.451, "no"),
                     c("GO:1901362","organic cyclic compound biosynthetic process",17.871,-4.564, 4.607, 6.360,-7.1574,0.825,0.454, "no"),
                     c("GO:0034654","nucleobase-containing compound biosynthetic process",14.533,-4.911, 3.310, 6.271,-8.5331,0.756,0.468, "no"),
                     c("GO:0050794","regulation of cellular process",18.840, 5.350, 4.404, 6.383,-9.9547,0.569,0.474, "no"),
                     c("GO:0010033","response to organic substance", 0.900, 4.700,-4.954, 5.062,-4.7328,0.824,0.474, "no"),
                     c("GO:0007166","cell surface receptor signaling pathway", 0.920, 6.146, 0.249, 5.072,-4.9547,0.598,0.475, "yes"),
                     c("GO:0019438","aromatic compound biosynthetic process",16.954,-5.308, 3.368, 6.338,-7.9318,0.799,0.477, "no"),
                     c("GO:0048583","regulation of response to stimulus", 1.120, 6.359, 0.335, 5.158,-4.5346,0.627,0.487, "no"),
                     c("GO:0034645","cellular macromolecule biosynthetic process",19.291,-4.605, 3.404, 6.394,-6.2381,0.755,0.495, "no"),
                     c("GO:0009605","response to external stimulus", 1.370, 4.805,-4.912, 5.245,-3.1733,0.845,0.501, "no"),
                     c("GO:0009059","macromolecule biosynthetic process",19.548,-4.447, 4.744, 6.399,-5.4295,0.823,0.506, "no"),
                     c("GO:0016569","covalent chromatin modification", 0.424,-1.965, 3.757, 4.736,-3.1618,0.871,0.507, "no"),
                     c("GO:0010646","regulation of cell communication", 0.929, 6.577, 3.315, 5.076,-3.8356,0.667,0.521, "no"),
                     c("GO:0044271","cellular nitrogen compound biosynthetic process",22.502,-5.381, 3.348, 6.460,-5.9626,0.787,0.536, "no"),
                     c("GO:0023051","regulation of signaling", 0.934, 6.804, 3.303, 5.079,-3.8447,0.682,0.536, "no"),
                     c("GO:0061024","membrane organization", 0.759,-1.073, 1.904, 4.989,-3.5751,0.927,0.540, "no"),
                     c("GO:0080134","regulation of response to stress", 0.337, 5.997,-0.235, 4.636,-3.5768,0.606,0.542, "no"),
                     c("GO:0003093","regulation of glomerular filtration", 0.001, 4.143, 4.199, 2.196,-3.2373,0.767,0.555, "no"),
                     c("GO:0060337","type I interferon signaling pathway", 0.007, 6.048,-1.751, 2.943,-7.8729,0.615,0.560, "yes"),
                     c("GO:0090304","nucleic acid metabolic process",21.449,-5.351, 2.127, 6.440,-6.0640,0.772,0.562, "no"),
                     c("GO:0042221","response to chemical", 3.071, 4.821,-4.723, 5.595,-3.6440,0.835,0.562, "no"),
                     c("GO:0006468","protein phosphorylation", 4.137,-4.217, 0.887, 5.725,-3.5406,0.831,0.568, "no"),
                     c("GO:2001034","positive regulation of double-strand break repair via nonhomologous end joining", 0.002, 5.120,-0.408, 2.326,-3.1296,0.581,0.569, "no"),
                     c("GO:0016070","RNA metabolic process",15.951,-5.147, 2.191, 6.311,-3.9431,0.767,0.577, "yes"),
                     c("GO:0034112","positive regulation of homotypic cell-cell adhesion", 0.001, 4.644, 4.659, 1.973,-4.3233,0.713,0.578, "no"),
                     c("GO:0044249","cellular biosynthetic process",30.048,-5.098, 4.189, 6.586,-3.9547,0.806,0.585, "no"),
                     c("GO:0034097","response to cytokine", 0.136, 4.810,-5.247, 4.242,-5.0835,0.829,0.597, "no"),
                     c("GO:0006139","nucleobase-containing compound metabolic process",26.547,-5.640, 1.668, 6.532,-5.7747,0.796,0.597, "no"),
                     c("GO:0030155","regulation of cell adhesion", 0.129, 6.430, 3.427, 4.219,-5.5287,0.717,0.598, "no"),
                     c("GO:0048523","negative regulation of cellular process", 1.830, 5.992, 5.009, 5.371,-5.2426,0.611,0.601, "no"),
                     c("GO:0043900","regulation of multi-organism process", 0.079, 5.161, 6.134, 4.008,-5.3778,0.727,0.610, "no"),
                     c("GO:0051239","regulation of multicellular organismal process", 0.628, 5.929, 3.595, 4.906,-5.4597,0.668,0.614, "no"),
                     c("GO:0001816","cytokine production", 0.120, 3.183, 2.200, 4.187,-3.6180,0.900,0.619, "no"),
                     c("GO:0008286","insulin receptor signaling pathway", 0.024, 6.388,-1.211, 3.484,-3.1713,0.661,0.625, "no"),
                     c("GO:0035556","intracellular signal transduction", 4.000, 6.114, 0.398, 5.710,-4.5114,0.539,0.641, "no"),
                     c("GO:0030335","positive regulation of cell migration", 0.076, 4.023, 5.272, 3.988,-4.2668,0.612,0.650, "no"),
                     c("GO:0006996","organelle organization", 3.595,-1.125, 2.020, 5.664,-3.0191,0.915,0.652, "no"),
                     c("GO:1901576","organic substance biosynthetic process",30.365,-4.277, 5.476, 6.591,-3.9066,0.837,0.658, "no"),
                     c("GO:0006974","cellular response to DNA damage stimulus", 2.360, 5.323,-4.329, 5.481,-3.5331,0.770,0.660, "no"),
                     c("GO:0035455","response to interferon-alpha", 0.003, 4.148,-5.794, 2.606,-4.8996,0.846,0.660, "yes"),
                     c("GO:0006796","phosphate-containing compound metabolic process",13.110,-5.887,-2.168, 6.226,-5.7747,0.879,0.664, "no"),
                     c("GO:0035456","response to interferon-beta", 0.004, 4.258,-5.751, 2.692,-3.4763,0.845,0.667, "no"),
                     c("GO:0035549","positive regulation of interferon-beta secretion", 0.000, 3.349, 5.019, 1.462,-3.5129,0.682,0.676, "no"),
                     c("GO:0040012","regulation of locomotion", 0.183, 5.174, 5.429, 4.372,-3.6757,0.723,0.677, "no"),
                     c("GO:0006302","double-strand break repair", 0.211, 2.563,-4.492, 4.432,-3.8356,0.726,0.689, "no"),
                     c("GO:0006355","regulation of transcription, DNA-templated", 9.917, 1.091, 4.714, 6.105,-8.5302,0.463,0.695, "yes"),
                     c("GO:0019222","regulation of metabolic process",11.942, 6.002, 4.071, 6.185,-7.7905,0.604,0.698, "no"));

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

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below
one.data <- subset(one.data[one.data$log10_p_value < -3.9,] )
one.data <- one.data[order(one.data$log10_p_value, decreasing = T),]

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
  filename = "results/fig_supp_03_tcell-rna-seq/fig_supp_03_panel-A_reviGO_upregulated.pdf",
  p1, 
  width = 3.75, height = 3.75, units = "in"
)
