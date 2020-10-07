### packages ==================================================================
library(tidyverse)
library(cowplot)
library(ggtext)

### read absorbance data ======================================================
## IL17A
il17a <- list()
il17a[["stds"]] <- read.csv("data/elisa/20170621_ELISA_IL17A_standards.csv")
il17a[["samp"]] <- read.csv("data/elisa/20170621_ELISA_IL17A_samples.csv")
il17a$samp$sample_name <- factor(il17a$samp$sample_name, 
                                 levels = c("T_unstim", "T_stim", "d0_imm", "d2_imm", "d10_imm"))
il17a[["dilution"]] <- 1 # undiluted samples were used for ELISA

## IL6
il6 <- list()
il6[["stds"]] <- read.csv("data/elisa/20170224_ELISA_IL6_standards.csv") # same batch as samples
il6[["samp"]] <- read.csv("data/elisa/20170326_ELISA_IL6_samples.csv") 
il6$samp$sample_name <- factor(il6$samp$sample_name, 
                               levels = c("d0", "d2"))
il6[["dilution"]] <- 5 # 1:5 diluted samples were used for ELISA

### calculate concentrations ==================================================
## fit linear models
il17a$lm <- lm(data = il17a$stds, 
               formula = log10(abs_corrected[2:8]) ~ log10(conc_pgml[2:8]))
il6$lm <- lm(data = il6$stds,
             formula = log10(abs_corrected[1:7]) ~ log10(conc_pgml[1:7]))

## calculate concentrations
il17a$samp$conc_pgml <- il17a$dilution * (
  10 ^ (
    (log10(il17a$samp$abs_corrected) - il17a$lm$coefficients[[1]]) /
      il17a$lm$coefficients[[2]]
  )
)

il6$samp$conc_pgml <- il6$dilution * (
  10 ^ (
    (log10(il6$samp$abs_corrected) - il6$lm$coefficients[[1]]) /
      il6$lm$coefficients[[2]]
  )
)

### stats =====================================================================
il17a$anova <- aov(lm(data = il17a$samp, formula = conc_pgml ~ sample_name))
il17a$tukey <- TukeyHSD(il17a$anova)

il6$anova <- aov(lm(data = il6$samp, formula = conc_pgml ~ sample_name))
il6$tukey <- TukeyHSD(il6$anova)

### basic plots ===============================================================
## theme and colors
rectcol1 <- "white" # colors for rectagle annotations. see below.
rectcol2 <- "#FFF2CC"

plot_theme <-  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) 

## IL17A
il17a$p <- ggplot(data = il17a$samp, aes(x = sample_name, y = conc_pgml)) +
  geom_hline(yintercept = mean(il17a$samp$conc_pgml[il17a$samp$sample_name == "T_unstim"]),
             linetype = 2, color = "grey", size = 0.25) +
  geom_jitter(size = 1.5, shape = 21, alpha = 1,
              stroke = 0.25, width = 0.075, color = "black", fill = rectcol2) +
  coord_cartesian(ylim = c(0, 70), xlim = c(0, 6), clip = "off", expand = FALSE) +
  stat_summary(data = il17a$samp, fun.data = 'mean_se',
               size = 0.25, colour = 'black', geom = 'errorbar', width = 0.2, 
               position = position_nudge(x = 0.30)) +
  stat_summary(data = il17a$samp, fun.data = 'mean_se',
               geom = 'point', size = 0.75, shape = 21, fill="black",
               position = position_nudge(x = 0.30)) +
  ylab('IL17A (pg/ml)') +
  plot_theme

## IL6
il6$p <- ggplot(data = il6$samp, aes(x = sample_name, y = conc_pgml)) +
  geom_hline(yintercept = mean(il6$samp$conc_pgml[il6$samp$sample_name == "d0"]),
             linetype = 2, color = "grey", size = 0.25) +
  geom_jitter(size = 1.5, shape = 21, alpha = 1,
              stroke = 0.25, width = 0.075, color = "black", fill = rectcol2) +
  coord_cartesian(ylim = c(0, 500), xlim = c(0, 3), clip = "off", expand = FALSE) +
  stat_summary(data = il6$samp, fun.data = 'mean_se',
               size = 0.25, colour = 'black', geom = 'errorbar', width = 0.2, 
               position = position_nudge(x = 0.30)) +
  stat_summary(data = il6$samp, fun.data = 'mean_se',
               geom = 'point', size = 0.75, shape = 21, fill="black",
               position = position_nudge(x = 0.30)) +
  ylab('IL6 (pg/ml)') +
  plot_theme

### add p values ==============================================================
## IL17A ----------------------------------------------------------------------
# segment data
il17a.seg.adj <- 0.1
il17a.tick.height <- 1
il17a.seg.size <- 0.25

il17a.segs <- data.frame(
  seg = c("1", "2", "3"),
  sample1 = c(3, 3, 4),
  sample2 = c(5, 4, 5),
  pval = c("0.02", "2.8e-5", "0.09"), 
  y = c(8, 15, 15)
)

il17a.segs$x <- il17a.segs$sample1 + il17a.seg.adj
il17a.segs$xend <- il17a.segs$sample2 - il17a.seg.adj
il17a.segs$tick.y <- il17a.segs$y - il17a.tick.height
il17a.segs$tick.yend  <- il17a.segs$y + il17a.tick.height

# add segments and pvals to plot
for(i in 1:nrow(il17a.segs)){
  il17a$p <- il17a$p + # segments
    annotate(geom = "segment", size = il17a.seg.size,
             x = il17a.segs$x[i], xend = il17a.segs$xend[i], 
             y = il17a.segs$y[i], yend = il17a.segs$y[i])
  il17a$p <- il17a$p + # left tick
    annotate(geom = "segment", size = il17a.seg.size,
             x = il17a.segs$x[i], xend = il17a.segs$x[i],
             y = il17a.segs$tick.y[i], yend = il17a.segs$tick.yend[i])
  il17a$p <- il17a$p + # right tick
    annotate(geom = "segment", size = il17a.seg.size,
             x = il17a.segs$xend[i], xend = il17a.segs$xend[i],
             y = il17a.segs$tick.y[i], yend = il17a.segs$tick.yend[i])
  il17a$p <- il17a$p + # pvalue
    annotate(geom = "text", size = 6/.pt, vjust = -0.5,
             x = mean(c(il17a.segs$x[i], il17a.segs$xend[i])), y = il17a.segs$y[i],
             label = il17a.segs$pval[i])
}

il17a$p <- il17a$p +
  annotate(geom = "richtext", label = "Tukey test *p* values",
           x = 4, y = il17a.segs$y[1], size = 6/.pt, vjust = 1.75, 
           fill = NA, 
           label.color = NA, # remove background and outline
           label.padding = grid::unit(0, "pt"))

## IL6 ------------------------------------------------------------------------
# segment data
il6.seg.adj <- 0.1
il6.tick.height <- il17a.tick.height * (500 / 70) # have to convert (ymax of il6 over ymax of il17a)
il6.seg.size <- 0.25

il6.segs <- data.frame(
  seg = c("1"),
  sample1 = c(1),
  sample2 = c(2),
  pval = c("0.77"), 
  y = c(70)
)

il6.segs$x <- il6.segs$sample1 + il6.seg.adj
il6.segs$xend <- il6.segs$sample2 - il6.seg.adj
il6.segs$tick.y <- il6.segs$y - il6.tick.height
il6.segs$tick.yend  <- il6.segs$y + il6.tick.height

# add to plot
for(i in 1:nrow(il6.segs)){
  il6$p <- il6$p + # segments
    annotate(geom = "segment", size = il6.seg.size,
             x = il6.segs$x[i], xend = il6.segs$xend[i], 
             y = il6.segs$y[i], yend = il6.segs$y[i])
  il6$p <- il6$p + # left tick
    annotate(geom = "segment", size = il6.seg.size,
             x = il6.segs$x[i], xend = il6.segs$x[i],
             y = il6.segs$tick.y[i], yend = il6.segs$tick.yend[i])
  il6$p <- il6$p + # right tick
    annotate(geom = "segment", size = il6.seg.size,
             x = il6.segs$xend[i], xend = il6.segs$xend[i],
             y = il6.segs$tick.y[i], yend = il6.segs$tick.yend[i])
  il6$p <- il6$p + # pvalue
    annotate(geom = "text", size = 6/.pt, vjust = -0.5,
             x = mean(c(il6.segs$x[i], il6.segs$xend[i])), y = il6.segs$y[i],
             label = il6.segs$pval[i])
}

il6$p <- il6$p +
  annotate(geom = "richtext", label = "*t* test *p* value",
           x = 1.5, y = il6.segs$y[1], size = 6/.pt, vjust = 1.75, 
           fill = NA, 
           label.color = NA, # remove background and outline
           label.padding = grid::unit(0, "pt"))

### add margins ===============================================================
# add margins to accommodate annotation
il17a$p.mar <- il17a$p +
  theme(plot.margin = margin(6, 6, 125, 45, unit = "pt"))
il6$p.mar <- il6$p +
  theme(plot.margin = margin(6, 6, 125, 6, unit = "pt"))

### parameters for annotation data ============================================
# common to both IL17A and IL6
annsize <- 6/.pt
rectcol1 <- "#FFF2CC"
rectcol2 <- "white" #"#FFFAEA"

### add annotations: IL17A ====================================================
il17a.annstart.y <- -8
il17a.thickness.adj <- 0.9 # for adjusting heigh of box relative to text lines

il17a$ann <- data.frame(
  titles = c("Cells", "Treatment", "Medium ratio\n(Native : Decidual)", "Decidual medium", "Decidualization\nduration (days)"),
  s1 = c("Naive T", "none", "100 : 0", "NA", "NA"),
  s2 = c("Naive T", "IL6 +\nIL23 +\nTGFB1", "100 : 0", "NA", "NA"),
  s3 = c("Naive T", "IL6 +\nIL23 +\nTGFB1", "50 : 50", "Control", "NA"),
  s4 = c("Naive T", "IL6 +\nIL23 +\nTGFB1", "50 : 50", "Cond.", "2"),
  s5 = c("Naive T", "IL6 +\nIL23 +\nTGFB1", "50 : 50", "Cond.", "10"),
  rowpos = c(1:5),
  thickness = c(1.5, 3, 1.5, 1.5, 1.5) # ~ number of text lines
)

il17a$ann$mid <- NA
il17a$ann$mid[1] <- il17a.annstart.y

for(i in 2:5) {
  il17a$ann$mid[i] <- il17a$ann$mid[(i-1)] - 
    (il17a$ann$thickness[(i-1)] + il17a$ann$thickness[(i)]) * 
    (annsize * il17a.thickness.adj)
} # y midpoint for each box based on required height to accommodate text
il17a$ann$min <- il17a$ann$mid - -(il17a$ann$thickness * 
                                     (annsize * il17a.thickness.adj))
il17a$ann$max <- il17a$ann$mid -  (il17a$ann$thickness * 
                                     (annsize * il17a.thickness.adj))


is.even <- function(x) x %% 2 == 0

# text box rectangles
il17a$p.ann <- il17a$p.mar
for(i in il17a$ann$rowpos) {
  if (is.even(i))
    il17a$p.ann <- il17a$p.ann +
      annotate(geom = "rect", 
               xmin = 0, xmax = 6, 
               ymin = il17a$ann$min[i], 
               ymax = il17a$ann$max[i],
               fill = rectcol2)
  else
    il17a$p.ann <- il17a$p.ann +
      annotate(geom = "rect", 
               xmin = 0, xmax = 6, 
               ymin = il17a$ann$min[i], 
               ymax = il17a$ann$max[i],
               fill = rectcol1)
}

# text
for(i in il17a$ann$rowpos) {
  il17a$p.ann <- il17a$p.ann + 
    annotate(geom = "text",
             x = -0.2, y = il17a$ann$mid[i], 
             label = il17a$ann$titles[i],
             hjust = 1, 
             size = annsize, 
             fontface = 2) # titles
  
  for(j in 1:5) {
    il17a$p.ann <- il17a$p.ann +  
      annotate(geom = "text", 
               x = j, y = il17a$ann$mid[i], 
               label = il17a$ann[i, paste0("s", j)], 
               size = annsize) 
  }
}

### add annotations: IL6 ======================================================
# If we use the same numbers as above, the annotations will not be aligned with the IL17A plot's annotations. This is because we are using native coordinate units, and the coordinates are different between the plots, i.e. IL17A y goes from 0 to 70 and IL6 y goes from 0 to 500. Using the same y numbers as IL17A plot for IL6 will squish the annotations. We have to multiply all y numbers with (500/70)
il17a.ymax <- 70
il6.ymax <- 500
conv <- il6.ymax / il17a.ymax

il6.annstart.y <- il17a.annstart.y * conv
il6.thickness.adj <- il17a.thickness.adj * conv # for adjusting heigh of box relative to text lines

il6$ann <- data.frame(
  titles = c("Cells", "Treatment", "Medium ratio\n(Native : Decidual)", "Decidual medium", "Decidualization\nduration (days)"),
  s1 = c("M.phage", "LPS", "50 : 50", "Control", "NA"),
  s2 = c("M.phage", "LPS", "50 : 50", "Cond.", "2"),
  rowpos = c(1:5),
  thickness = c(1.5, 3, 1.5, 1.5, 1.5) # ~ number of text lines
)

il6$ann$mid <- NA
il6$ann$mid[1] <- il6.annstart.y

for(i in 2:5) {
  il6$ann$mid[i] <- il6$ann$mid[(i-1)] - 
    (il6$ann$thickness[(i-1)] + il6$ann$thickness[(i)]) * 
    (annsize * il6.thickness.adj)
} # y midpoint for each box based on required height to accommodate text
il6$ann$min <- il6$ann$mid - -(il6$ann$thickness * (annsize * il6.thickness.adj))
il6$ann$max <- il6$ann$mid -  (il6$ann$thickness * (annsize * il6.thickness.adj))


is.even <- function(x) x %% 2 == 0

# text box rectangles
il6$p.ann <- il6$p.mar
for(i in il6$ann$rowpos) {
  if (is.even(i))
    il6$p.ann <- il6$p.ann +
      annotate(geom = "rect", 
               xmin = 0, xmax = 3, 
               ymin = il6$ann$min[i], 
               ymax = il6$ann$max[i],
               fill = rectcol2)
  else
    il6$p.ann <- il6$p.ann +
      annotate(geom = "rect", 
               xmin = 0, xmax = 3, 
               ymin = il6$ann$min[i], 
               ymax = il6$ann$max[i],
               fill = rectcol1)
}

# text
for(i in il6$ann$rowpos) {
  for(j in 1:2) {
    il6$p.ann <- il6$p.ann +  
      annotate(geom = "text", 
               x = j, y = il6$ann$mid[i], 
               label = il6$ann[i, paste0("s", j)], 
               size = annsize) 
  }
}

### align plots ===============================================================
aligned <- align_plots(il17a$p.ann, il6$p.ann, align = "h", axis = "tb") 
(gridtest <- plot_grid(aligned[[1]], aligned[[2]], rel_widths = c(10, 5)))
cowplot::ggsave2(gridtest, filename = "~/Downloads/gridtest.pdf",
                 width = 6, height = 4.4, units = "in")
   







