require (ggplot2)
require (gridExtra)

#LOADS DATA

neu.data <- read.table (file = "ESTIMATE_evolution_raw_data.txt", sep = "\t", header = TRUE)

#DEFINES THE NEUTRAL EVOLUTION VARIABLE

neu.data$Group <- ifelse (neu.data$r2.all > 0.98, "Neutral", "Not neutral")

#THESE ARE THE GLOBAL ANCOVA ANALYSES. THIS IS COMMENTED SO YOU CAN RE-CREATE THE FIGURE
#BY RUNNING THE CODE THAT'S BELOW. CD274 IS PD-L1 AND PDCD1LG2 IS PD-L2

##GLOBAL NEUTRAL EVOLUTION ANCOVA
#anova (lm (ESTIMATE ~ type + Neo + CD274 + PDCD1LG2 + Group, data = sub.analyze))

#GLOBAL NEUTRAL EVOLUTION ACCOUNTING ALSO FOR MUTATION RATE (mu.m.all)
#anova (lm (ESTIMATE ~ type + Neo + CD274 + PDCD1LG2 + mu.m.all + Group, data = sub.analyze))

#MUTATION RATE AND IMMUNITY ACCOUNTING ALSO FOR NEUTRAL EVOLUTION
#anova (lm (ESTIMATE ~ type + Neo + CD274 + PDCD1LG2 + Group + mu.m.all, data = sub.analyze))

#THE CODE FROM THIS POINT UNTIL THE END RE-CREATE THE FIGURE
#PREPARES A DATA FRAME WITH THE LABELS FOR THE TISSUE-SPECIFIC PLOT

p.estimate <- data.frame (type = unique (neu.data$type), p = 1, Mean = 0, y = 3000)

for (i in 1:nrow (p.estimate)) {
	sub.t <- subset (neu.data, type == as.character (p.estimate[i,1]))
	res <- wilcox.test (sub.t$ESTIMATE ~ sub.t$Group)
	p.estimate[i,2] <- res$p.value
	p.estimate[i,3] <- mean (sub.t$ESTIMATE)
}

p.estimate$Label <- ifelse (p.estimate$p > 0.05, "", ifelse (p.estimate$p > 0.01, "*", ifelse (p.estimate$p > 0.001, "**", "***")))

#THIS SORTS THE CANCER TYPES ACCORDING TO WHETHER THERE IS A CORRELATION BETWEEN NEUTRAL EVOLUTION AND IMMUNITY OR NOT AND THEIR ESTIMATE MEAN
neu.data$type <- factor (neu.data$type, levels = p.estimate[order(p.estimate$Label, -p.estimate$Mean),]$type)

#THIS CREATES THE BOXPLOT WITH THE ESTIMATE SCORES
estimate.plot <- ggplot (neu.data, aes (x = type, y = ESTIMATE)) + geom_boxplot (outlier.shape = NA, aes(dodge = Group)) + geom_point (position = position_jitterdodge (jitter.width = 0.3), aes (color = Group)) + theme (axis.text.y = element_text (color = "black"), axis.ticks.y = element_line (color = "black"), axis.line.y = element_line (color = "black"), axis.ticks.x = element_line (color = "black"), axis.text.x = element_text (color = "black"), axis.line.x = element_line (color = "black"), legend.position = "none", panel.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text (face = "italic")) + scale_color_manual (values = c ("#F62A00", "#4CB5F5")) + geom_text (data = subset(p.estimate, p > 0.05), aes (x = type, y = y, label = Label), size = 3) + geom_text (data = subset(p.estimate, p < 0.05), aes (x = type, y = y, label = Label), size = 5)

#THIS CREATES THE CORRELATION PLOT BETWEEN MUTATION RATE AND ESTIMATE
mu.plot <- ggplot (neu.data, aes (x = mu.m.all, y = ESTIMATE)) + geom_point (aes (color = Group)) + scale_x_log10() + facet_wrap (~Group) + geom_smooth (method = "lm") + theme (axis.text = element_text (color = "black"), axis.ticks = element_line (color = "black"), axis.line.x = element_line (color = "black"), axis.line.y = element_line (color = "black"), panel.background = element_blank(), axis.title.y = element_text (face = "italic")) + scale_color_manual (values = c ("#F62A00", "#4CB5F5")) + xlab ("Mutation rate") + theme (legend.position = "bottom", strip.background = element_blank(), strip.text = element_text (face = "bold.italic"))

#PUT IT TOGETHER IN A SINGLE PANEL
gA=ggplot_gtable(ggplot_build(estimate.plot))
gB=ggplot_gtable(ggplot_build(mu.plot))
maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)
figure.1 <- grid.arrange (gA, gB, ncol = 1)

#SAVE IT
ggsave (figure.1, file = "Figure_1.pdf")