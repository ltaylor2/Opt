setwd("~/Homework/Optimization and Uncertainty/Opt5")

require(ggplot2)
require(dplyr)
full <- read.csv("ex_1_2.csv", header=TRUE, stringsAsFactors=TRUE)

ggplot(full, aes(x=Type, y=Time)) +
  geom_boxplot() +
  ylab("Time") +
  xlab("Algorithm Type") +
  theme_bw()
