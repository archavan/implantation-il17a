

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
revigo.data <- rbind(c("GO:0002376","immune system process", 0.600, 0.586,-7.055, 4.886,-3.6517,0.872,0.000),
c("GO:0042110","T cell activation", 0.081, 0.095, 5.830, 4.018,-3.9318,0.338,0.000),
c("GO:0007256","activation of JNKK activity", 0.001, 4.437,-1.080, 2.260,-3.2069,0.533,0.102),
c("GO:0031023","microtubule organizing center organization", 0.061,-4.303,-2.060, 3.896,-3.6234,0.690,0.129),
c("GO:0072330","monocarboxylic acid biosynthetic process", 0.940,-4.008, 2.230, 5.081,-3.8508,0.683,0.159),
c("GO:0006928","movement of cell or subcellular component", 0.973,-3.870, 4.107, 5.097,-3.2865,0.683,0.201),
c("GO:0007166","cell surface receptor signaling pathway", 0.920, 6.298, 0.747, 5.072,-3.0119,0.513,0.320),
c("GO:0050808","synapse organization", 0.070,-4.936,-0.962, 3.956,-3.2807,0.689,0.332),
c("GO:0019221","cytokine-mediated signaling pathway", 0.093, 5.995,-0.203, 4.078,-3.1904,0.504,0.435),
c("GO:0002577","regulation of antigen processing and presentation", 0.003, 3.574, 5.224, 2.609,-3.1561,0.430,0.582),
c("GO:0090025","regulation of monocyte chemotaxis", 0.005, 3.504, 2.697, 2.763,-3.9788,0.287,0.594),
c("GO:0007165","signal transduction", 6.621, 5.421, 0.798, 5.929,-4.2190,0.474,0.679));

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
