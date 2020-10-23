

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0006413","translational initiation", 0.518,-3.862, 6.557, 4.823,-91.6576,0.781,0.000),
c("GO:0008150","biological_process",100.000,-0.817, 0.923, 7.108,-4.4559,1.000,0.000),
c("GO:0008152","metabolic process",75.387, 0.206,-0.872, 6.986,-22.9872,0.999,0.000),
c("GO:0009987","cellular process",63.780, 1.036, 1.080, 6.913,-10.7055,0.998,0.000),
c("GO:0048002","antigen processing and presentation of peptide antigen", 0.013, 2.471,-0.798, 3.236,-3.4711,0.970,0.000),
c("GO:0051179","localization",18.495, 0.202, 0.208, 6.375,-11.2644,0.995,0.000),
c("GO:0051704","multi-organism process", 0.751,-2.817,-1.211, 4.984,-29.1445,0.994,0.000),
c("GO:0070972","protein localization to endoplasmic reticulum", 0.187, 0.756,-6.413, 4.381,-86.0424,0.871,0.000),
c("GO:0071840","cellular component organization or biogenesis", 8.568,-3.587,-1.467, 6.041,-6.6861,0.995,0.000),
c("GO:0009056","catabolic process", 4.820, 1.429,-0.050, 5.791,-28.6498,0.977,0.017),
c("GO:0009058","biosynthetic process",31.611,-1.003,-1.215, 6.608,-15.4559,0.971,0.033),
c("GO:0033554","cellular response to stress", 2.967, 6.402, 1.732, 5.581,-8.4498,0.860,0.037),
c("GO:0006457","protein folding", 0.903, 0.451, 1.792, 5.064,-6.5045,0.951,0.040),
c("GO:0022900","electron transport chain", 0.564, 4.828,-3.044, 4.860,-15.3372,0.822,0.057),
c("GO:0006091","generation of precursor metabolites and energy", 1.940,-2.736, 1.372, 5.396,-11.8794,0.921,0.065),
c("GO:1901360","organic cyclic compound metabolic process",30.324,-3.807, 2.090, 6.590,-16.7282,0.930,0.069),
c("GO:0044237","cellular metabolic process",53.061,-6.603,-1.661, 6.833,-27.1772,0.913,0.078),
c("GO:0006807","nitrogen compound metabolic process",38.744, 3.106,-1.203, 6.696,-15.8356,0.970,0.088),
c("GO:0000184","nuclear-transcribed mRNA catabolic process, nonsense-mediated decay", 0.030,-3.696, 4.636, 3.588,-77.6799,0.775,0.106),
c("GO:0010605","negative regulation of macromolecule metabolic process", 1.169, 2.427, 7.156, 5.176,-18.5607,0.709,0.116),
c("GO:0071704","organic substance metabolic process",58.357,-0.941,-3.179, 6.874,-12.6819,0.968,0.119),
c("GO:0044238","primary metabolic process",53.743,-0.478,-3.767, 6.839,-12.0804,0.968,0.120),
c("GO:1903332","regulation of protein folding", 0.001, 3.887, 5.512, 1.929,-3.8125,0.884,0.162),
c("GO:0061469","regulation of type B pancreatic cell proliferation", 0.002, 3.103, 5.185, 2.354,-4.3197,0.867,0.172),
c("GO:0046483","heterocycle metabolic process",29.664,-7.095,-0.275, 6.580,-19.5072,0.888,0.176),
c("GO:0007059","chromosome segregation", 0.476, 5.141,-2.147, 4.786,-7.0980,0.863,0.178),
c("GO:0043603","cellular amide metabolic process", 6.879,-7.363, 2.862, 5.946,-46.2684,0.860,0.178),
c("GO:1903047","mitotic cell cycle process", 0.514, 4.835,-2.027, 4.819,-12.3002,0.722,0.179),
c("GO:0019538","protein metabolic process",18.489,-2.578, 7.155, 6.375,-16.8665,0.885,0.183),
c("GO:0019083","viral transcription", 0.012,-4.332, 4.940, 3.203,-74.5607,0.846,0.184),
c("GO:0051301","cell division", 1.230, 5.571,-3.049, 5.198,-6.7520,0.852,0.195),
c("GO:0043170","macromolecule metabolic process",39.491,-4.873, 1.570, 6.705,-14.5560,0.927,0.211),
c("GO:0051348","negative regulation of transferase activity", 0.088, 4.669, 6.353, 4.053,-10.4559,0.837,0.219),
c("GO:0007049","cell cycle", 1.885, 5.457,-2.413, 5.384,-3.7144,0.846,0.223),
c("GO:0015833","peptide transport", 0.298, 0.166,-6.398, 4.582,-23.7620,0.922,0.227),
c("GO:0006839","mitochondrial transport", 0.182, 1.207,-6.787, 4.369,-7.2055,0.938,0.236),
c("GO:0006725","cellular aromatic compound metabolic process",29.628,-6.826,-0.042, 6.580,-18.3410,0.888,0.245),
c("GO:0042981","regulation of apoptotic process", 0.313, 4.376, 4.095, 4.604,-6.3344,0.737,0.259),
c("GO:0034660","ncRNA metabolic process", 3.407,-5.349, 5.114, 5.641,-51.1605,0.805,0.269),
c("GO:0051641","cellular localization", 2.041, 0.093,-6.277, 5.418,-22.3206,0.925,0.283),
c("GO:0045454","cell redox homeostasis", 0.861, 5.079, 3.860, 5.043,-4.6840,0.759,0.288),
c("GO:0055086","nucleobase-containing small molecule metabolic process", 4.917,-5.958, 1.351, 5.800,-4.4660,0.739,0.305),
c("GO:0034641","cellular nitrogen compound metabolic process",34.137,-6.915, 2.488, 6.641,-23.1457,0.839,0.310),
c("GO:0070498","interleukin-1-mediated signaling pathway", 0.005, 6.086, 3.948, 2.792,-3.7959,0.818,0.315),
c("GO:0048519","negative regulation of biological process", 1.984, 4.925, 6.406, 5.406,-12.0650,0.838,0.316),
c("GO:0055114","oxidation-reduction process",15.060, 1.589,-4.518, 6.286,-7.0580,0.881,0.320),
c("GO:1902600","hydrogen ion transmembrane transport", 1.015, 0.073,-6.922, 5.115,-6.7905,0.859,0.324),
c("GO:1903320","regulation of protein modification by small protein conjugation or removal", 0.078, 1.712, 7.429, 3.999,-8.8633,0.720,0.326),
c("GO:0006984","ER-nucleus signaling pathway", 0.013, 5.496, 4.181, 3.215,-3.8386,0.819,0.340),
c("GO:0071705","nitrogen compound transport", 1.767, 0.608,-6.896, 5.355,-18.1169,0.924,0.347),
c("GO:0038061","NIK/NF-kappaB signaling", 0.018, 5.422, 4.397, 3.360,-4.8996,0.810,0.350),
c("GO:1901796","regulation of signal transduction by p53 class mediator", 0.012, 6.103, 4.583, 3.189,-3.5670,0.798,0.350),
c("GO:0033036","macromolecule localization", 3.030, 0.605,-6.658, 5.590,-20.8601,0.922,0.372),
c("GO:0016071","mRNA metabolic process", 0.798,-6.206, 5.348, 5.010,-50.7258,0.832,0.376),
c("GO:0006364","rRNA processing", 0.952,-5.375, 4.538, 5.087,-55.7496,0.689,0.384),
c("GO:1901566","organonitrogen compound biosynthetic process",14.064,-7.033, 3.503, 6.256,-38.6737,0.842,0.391),
c("GO:0009894","regulation of catabolic process", 0.146, 3.072, 6.748, 4.272,-5.0857,0.768,0.394),
c("GO:0010608","posttranscriptional regulation of gene expression", 0.719, 1.106, 7.209, 4.965,-5.2757,0.769,0.403),
c("GO:0044260","cellular macromolecule metabolic process",34.276,-5.243, 5.487, 6.643,-15.9830,0.834,0.408),
c("GO:1901564","organonitrogen compound metabolic process",17.886,-6.875, 3.341, 6.361,-19.5751,0.884,0.415),
c("GO:0071702","organic substance transport", 4.980, 0.443,-6.916, 5.805,-13.0186,0.915,0.423),
c("GO:0002181","cytoplasmic translation", 0.064,-3.323, 6.878, 3.915,-46.7471,0.815,0.429),
c("GO:0010035","response to inorganic substance", 0.317, 6.331, 1.646, 4.609,-3.4134,0.928,0.430),
c("GO:0071826","ribonucleoprotein complex subunit organization", 0.377,-4.083,-4.838, 4.685,-16.4750,0.834,0.439),
c("GO:0051246","regulation of protein metabolic process", 1.551, 1.592, 7.284, 5.299,-3.7595,0.751,0.439),
c("GO:0097031","mitochondrial respiratory chain complex I biogenesis", 0.013,-3.826,-3.984, 3.208,-4.8665,0.898,0.440),
c("GO:0022411","cellular component disassembly", 0.423, 1.982,-3.166, 4.734,-6.5670,0.779,0.444),
c("GO:0007005","mitochondrion organization", 0.418,-3.333,-3.827, 4.729,-4.7033,0.839,0.445),
c("GO:0044267","cellular protein metabolic process",14.293,-3.948, 6.663, 6.263,-20.4191,0.800,0.450),
c("GO:0006396","RNA processing", 3.210,-4.817, 5.665, 5.615,-44.2020,0.796,0.453),
c("GO:0042149","cellular response to glucose starvation", 0.016, 6.651, 1.479, 3.301,-3.5654,0.896,0.479),
c("GO:0051881","regulation of mitochondrial membrane potential", 0.014, 4.662, 5.984, 3.264,-3.7375,0.884,0.480),
c("GO:0006139","nucleobase-containing compound metabolic process",26.547,-6.637, 3.086, 6.532,-19.1013,0.804,0.484),
c("GO:0044249","cellular biosynthetic process",30.048,-7.466, 1.219, 6.586,-15.7033,0.852,0.498));

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

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
