### packages ==================================================================
library(tidyverse)
library(ggtext)

### data ======================================================================
ct <- read.csv("data/qpcr/20170913_qPCR_Th17_ct-data.csv")

### process data ==============================================================
# calculate dCt
ct <- ct %>% 
  group_by(sample_name) %>% 
  mutate(dCt = Ct - Ct[target_name == "TBP"]) %>% 
  ungroup()

# calculate ddCt
ct <- ct %>% 
  group_by(target_name) %>% 
  mutate(ddCt = dCt - mean(dCt[sample == "unstim"])) %>% 
  ungroup()

# calculate log 2 fold change
ct$log2fc <- -ct$ddCt

# stats
il17a.anova <- aov(lm(data = ct[ct$target_name == "IL17A", ],
                      formula = log2fc ~ sample))
il17a.tukey <- TukeyHSD(il17a.anova)

### prepare for plotting ======================================================
dat <- ct[which(ct$target_name != "TBP"), ]
dat$sample <- factor(dat$sample, levels = c("unstim", "stim", "d0", "d2"))

### plot ======================================================================
## basic plot -----------------------------------------------------------------
p <- ggplot(dat, aes(sample, log2fc)) +
  geom_hline(yintercept = mean(dat$log2fc[dat$sample == "unstim"]),
             linetype = 2, color = "black", size = 0.25) +
  geom_jitter(width = 0.1, 
              shape = 21, fill = "#ebe9dc", colour = "black", size = 1.75) +
  labs(y = "log<sub>2</sub> fold change<br>(over unstimulated)") +
  stat_summary(data = dat, fun.data = 'mean_se',
               size = 0.25, colour = 'black', geom = 'errorbar', width = 0.1, 
               position = position_nudge(x = 0.25)) +
  stat_summary(data = dat, fun.data = 'mean_se',
               geom = 'point', size = 0.75, shape = 21, fill="black",
               position = position_nudge(x = 0.25)) +
  coord_cartesian(ylim = c(-3, 4), xlim = c(0.5, 4.5), clip = "off", expand = FALSE) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7, colour = "black"),
        axis.title.y = element_markdown(size = 7),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25)) 

## add margin -----------------------------------------------------------------
p.mar <- p + theme(plot.margin = margin(25, 6, 125, 40))

## add p values ---------------------------------------------------------------
# segment data
seg.adj <- 0.1
tick.height <- 0.1
seg.size <- 0.25

segs <- data.frame(
  seg = c("1", "2"),
  sample1 = c(2, 3),
  sample2 = c(3, 4),
  pval = c("3.1e-4", "3.3e-2"), 
  y = c(4.2, 4.2)
)

segs$x <- segs$sample1 + seg.adj
segs$xend <- segs$sample2 - seg.adj
segs$tick.y <- segs$y - tick.height
segs$tick.yend  <- segs$y + tick.height

# add segments and pvals to plot
p.pval <- p.mar +
  annotate(geom = "segment", size = seg.size, # segment
           x = segs$x, xend = segs$xend,
           y = segs$y, yend = segs$y) +
  annotate(geom = "segment", size = seg.size, # left ticks
           x = segs$x, xend  = segs$x,
           y = segs$tick.y, yend = segs$tick.yend) +
  annotate(geom = "segment", size = seg.size, # right ticks
           x = segs$xend, xend  = segs$xend,
           y = segs$tick.y, yend = segs$tick.yend)

for(i in 1:nrow(segs)) { # p values
  p.pval <- p.pval +
    annotate(geom = "text", size = 5.5/.pt,
             x = mean(c(segs$sample1[i], segs$sample2[i])), y = 4.3,
             label = paste0("tukey p\n", segs$pval[i]),
             hjust = 0.5, vjust = 0)
}

### add annotation text boxes -------------------------------------------------
# parameters
annsize <- 6/.pt
rectcol1 <- "#ebe9dc"
rectcol2 <- "white" #"#FFFAEA"
annstart.y <- -3.5
thickness.adj <- 0.1 # for adjusting heigh of box relative to text lines

# dataframe with annotation data
ann <- data.frame(
  titles = c("Cells", "Treatment", "Medium ratio\n(Native : Decidual)", "Decidual medium", "Decidualization\nduration (days)"),
  s1 = c("Naive T", "none", "100 : 0", "NA", "NA"),
  s2 = c("Naive T", "IL6 +\nIL23 +\nTGFB1", "100 : 0", "NA", "NA"),
  s3 = c("Naive T", "IL6 +\nIL23 +\nTGFB1", "50 : 50", "Control", "NA"),
  s4 = c("Naive T", "IL6 +\nIL23 +\nTGFB1", "50 : 50", "Cond.", "2"),
  rowpos = c(1:5),
  thickness = c(1.5, 3, 1.5, 1.5, 1.5) # ~ number of text lines
)

ann$mid <- NA
ann$mid[1] <- annstart.y
for(i in 2:5) {
  ann$mid[i] <- ann$mid[(i-1)] - 
    (ann$thickness[(i-1)] + ann$thickness[(i)]) * 
    (annsize * thickness.adj)
} # y midpoint for each box based on required height to accommodate text

ann$min <- ann$mid - -(ann$thickness * 
                         (annsize * thickness.adj))
ann$max <- ann$mid -  (ann$thickness * 
                         (annsize * thickness.adj))


is.even <- function(x) x %% 2 == 0

# text box rectangles
p.ann <- p.pval
for(i in ann$rowpos) {
  if (is.even(i))
    p.ann <- p.ann +
      annotate(geom = "rect", 
               xmin = 0.5, xmax = 4.5, 
               ymin = ann$min[i], 
               ymax = ann$max[i],
               fill = rectcol2)
  else
    p.ann <- p.ann +
      annotate(geom = "rect", 
               xmin = 0.5, xmax = 4.5, 
               ymin = ann$min[i], 
               ymax = ann$max[i],
               fill = rectcol1)
}

# text
for(i in ann$rowpos) {
  p.ann <- p.ann + 
    annotate(geom = "text",
             x = 0.45, y = ann$mid[i], 
             label = ann$titles[i],
             hjust = 1, 
             size = annsize, 
             fontface = 2) # titles
  
  for(j in 1:5) {
    p.ann <- p.ann +  
      annotate(geom = "text", 
               x = j, y = ann$mid[i], 
               label = ann[i, paste0("s", j)], 
               size = annsize) 
  }
}

cowplot::ggsave2(
  filename = "results/fig_supp_02_qpcr/fig_supp_02_IL17A-qPCR.pdf",
  p.ann,
  width = 3.6, height = 4.25, units = "in"
)

 ## end ========================================================================
