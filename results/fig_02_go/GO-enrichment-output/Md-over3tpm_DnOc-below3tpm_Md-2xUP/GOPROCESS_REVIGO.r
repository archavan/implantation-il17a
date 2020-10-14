

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
revigo.data <- rbind(c("GO:0010951","negative regulation of endopeptidase activity", 0.157, 6.381,-0.996, 4.304,-8.3893,0.675,0.000),
c("GO:0015711","organic anion transport", 1.192, 3.957, 5.071, 5.184,-8.8297,0.775,0.000),
c("GO:0032501","multicellular organismal process", 2.373, 1.655,-1.146, 5.483,-4.6498,0.975,0.000),
c("GO:0060577","pulmonary vein morphogenesis", 0.000,-1.825, 6.790, 1.756,-3.4510,0.875,0.000),
c("GO:0071332","cellular response to fructose stimulus", 0.000, 1.366, 1.883, 1.653,-3.4510,0.939,0.017),
c("GO:0071825","protein-lipid complex subunit organization", 0.010,-4.498, 5.650, 3.105,-4.7721,0.904,0.021),
c("GO:0042737","drug catabolic process", 0.001,-0.746,-6.389, 2.286,-3.6364,0.895,0.034),
c("GO:1901615","organic hydroxy compound metabolic process", 0.831, 0.615,-7.830, 5.028,-3.6289,0.943,0.040),
c("GO:0006641","triglyceride metabolic process", 0.038,-4.416,-4.268, 3.687,-5.7852,0.704,0.041),
c("GO:0017144","drug metabolic process", 0.058,-0.187,-2.434, 3.868,-3.1002,0.941,0.043),
c("GO:0006629","lipid metabolic process", 3.522,-5.590,-1.041, 5.655,-5.0580,0.838,0.140),
c("GO:0008284","positive regulation of cell proliferation", 0.151, 6.035,-1.501, 4.288,-3.1226,0.794,0.209),
c("GO:1901605","alpha-amino acid metabolic process", 3.625,-5.987, 0.094, 5.668,-3.3883,0.762,0.236),
c("GO:0065008","regulation of biological quality", 3.395, 5.810,-2.525, 5.639,-4.9830,0.823,0.256),
c("GO:0010817","regulation of hormone levels", 0.161, 5.194,-3.784, 4.314,-4.7100,0.795,0.257),
c("GO:0044281","small molecule metabolic process",15.138,-6.119,-1.168, 6.288,-6.1367,0.830,0.302),
c("GO:0006730","one-carbon metabolic process", 0.328,-5.492, 0.767, 4.625,-3.2441,0.827,0.355),
c("GO:0071827","plasma lipoprotein particle organization", 0.007,-3.234, 5.694, 2.945,-4.8601,0.753,0.388),
c("GO:0006811","ion transport", 5.344, 3.795, 5.483, 5.836,-5.0472,0.865,0.389),
c("GO:0034367","macromolecular complex remodeling", 0.007,-4.138, 5.389, 2.969,-4.3382,0.905,0.390),
c("GO:0050892","intestinal absorption", 0.006, 0.282, 6.343, 2.913,-4.4437,0.802,0.408),
c("GO:1904640","response to methionine", 0.000, 0.825, 1.825, 1.398,-3.4510,0.954,0.409),
c("GO:0070328","triglyceride homeostasis", 0.008, 4.332,-4.786, 3.013,-3.6198,0.803,0.412),
c("GO:0001523","retinoid metabolic process", 0.013,-4.801,-4.641, 3.234,-3.5287,0.760,0.418),
c("GO:0006638","neutral lipid metabolic process", 0.042,-4.928,-3.775, 3.730,-4.9872,0.757,0.449),
c("GO:1901571","fatty acid derivative transport", 0.029, 2.896, 6.004, 3.564,-3.4202,0.815,0.453),
c("GO:0008202","steroid metabolic process", 0.161,-4.827,-4.093, 4.315,-4.4510,0.760,0.454),
c("GO:0071347","cellular response to interleukin-1", 0.014, 0.255, 2.058, 3.260,-3.1643,0.935,0.466),
c("GO:0043267","negative regulation of potassium ion transport", 0.004, 5.484, 2.001, 2.762,-3.1637,0.750,0.467),
c("GO:0036148","phosphatidylglycerol acyl-chain remodeling", 0.000,-3.479,-5.826, 1.255,-3.0114,0.794,0.476),
c("GO:0050878","regulation of body fluid levels", 0.074, 5.042,-4.182, 3.976,-3.7520,0.804,0.479),
c("GO:0036149","phosphatidylinositol acyl-chain remodeling", 0.000,-3.768,-5.863, 1.398,-3.4622,0.791,0.484),
c("GO:0042445","hormone metabolic process", 0.090, 4.668,-4.046, 4.064,-5.0942,0.785,0.486),
c("GO:0097272","ammonia homeostasis", 0.001, 4.028,-5.194, 2.049,-3.4510,0.827,0.492),
c("GO:0003008","system process", 0.660,-2.101, 6.824, 4.928,-3.3565,0.893,0.500),
c("GO:0006558","L-phenylalanine metabolic process", 0.075,-6.221, 1.462, 3.984,-3.4622,0.780,0.518),
c("GO:1902221","erythrose 4-phosphate/phosphoenolpyruvate family amino acid metabolic process", 0.075,-6.571, 1.455, 3.984,-3.4622,0.809,0.518),
c("GO:0016042","lipid catabolic process", 0.401,-4.184,-3.956, 4.712,-3.8447,0.718,0.540),
c("GO:0006869","lipid transport", 0.270, 3.325, 5.692, 4.539,-4.1343,0.787,0.551),
c("GO:0048878","chemical homeostasis", 0.543, 4.938,-3.718, 4.843,-3.6498,0.769,0.561),
c("GO:0071715","icosanoid transport", 0.029, 4.365, 5.038, 3.564,-3.4202,0.769,0.613),
c("GO:0046503","glycerolipid catabolic process", 0.018,-3.801,-4.733, 3.359,-4.4437,0.711,0.617),
c("GO:0055081","anion homeostasis", 0.045, 4.590,-4.511, 3.759,-3.2757,0.792,0.623),
c("GO:0015849","organic acid transport", 1.024, 3.124, 5.413, 5.119,-6.4737,0.766,0.633),
c("GO:0009072","aromatic amino acid family metabolic process", 0.719,-6.240, 0.530, 4.965,-3.2441,0.787,0.660),
c("GO:0006820","anion transport", 1.956, 4.479, 4.732, 5.400,-7.2596,0.823,0.663),
c("GO:0030162","regulation of proteolysis", 0.304, 6.392,-1.841, 4.591,-3.4559,0.747,0.669),
c("GO:0006570","tyrosine metabolic process", 0.059,-6.206, 1.731, 3.881,-3.2924,0.801,0.670));

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
