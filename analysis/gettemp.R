setwd("~/research/tippt/META-2021/analysis")

df1 <- read.csv("orig.csv")
df2 <- read.csv("useex.csv")

library(ggplot2)

ggplot(rbind(cbind(df1, group='FAIR-only'), cbind(df2, group='RMA Forcing')), aes(time, T, colour=group)) +
    geom_line() + xlab(NULL) + ylab(NULL) + theme_bw()
