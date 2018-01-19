# setwd("C:/Users/hp39wasi/sWorm/Theory")
library(ggplot2)


all <- read.csv("AllSoilSearchResults.csv")

byyear <- as.data.frame(table(all$PY, all$Theory))
names(byyear) <- c("year", "theory", "freq")


pdf(file = file.path(figs, "PublicationYear_theory.pdf"), width = 11)
ggplot(data = byyear, aes(x = factor(year), y = freq, color = theory)) +       
  geom_line(aes(group = theory)) + geom_point() +
  labs(x = "Year", y = "Frequency", colour = "")
dev.off()
