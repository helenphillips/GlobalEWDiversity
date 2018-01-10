figure <- "C:\\Users\\hp39wasi\\Google Drive\\sWorm\\sWorm_WS1\\TheorySubgroup\\Figure"

library(ggplot2)
library(reshape)

all <- c(638,1442,996,1383,26541,1293)
soil <- c(37,46,28,48,50,61)

## Boring line chart
pdf(file = file.path(figure, "lineChart.pdf"))
plot(log10(soil) ~ log10(all), pch = 19, cex =2, xlim = c(2.5,5), ylim = c(1.2, 1.9))
abline(a=log10(0.1), b=1, col='black', lty=3)  #  10% line
abline(a=log10(0.01), b=1, col='black', lty=3)  #  1% line
labs <- c("SER", "IBT", "MT", "SAR", "Niche/Neutral", "Meta")
text(log10(all), log10(soil), labels = labs, adj = c(1.2,1))
dev.off()

#####  Alternate Pie Chart

x <- data.frame(theory = labs, soil = soil, all = all)

d <- melt(x, id.vars = c("theory"))
names(d) <- c("theory", "biome", "n")
d$log10n <- log10(d$n)

soil <- d[d$biome == "soil",]
all <- d[d$biome == "all",]

pdf(file = file.path(figure, "pieChart.pdf"))
ggplot(all, aes(x=theory, y = log10n, col = theory)) +
  geom_bar(stat = "identity", fill = 'blue', alpha = 0.25) +
  geom_bar(data = soil, stat = "identity", fill= "green", alpha = 0.25,
           aes(x = theory, y = log10n, col = theory)) + 
  coord_polar(theta = "x") +
  theme_bw() +
  theme(legend.position="none") +
  labs(x = "Biodiversity Theory", y = "log10- Number of Publications", 
       title  = " ")
dev.off()
