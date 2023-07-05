setwd("~/research/tippt/META-2021/analysis")

library(ggplot2)
library(reshape2)

df1 <- read.csv("orig.csv")
df2 <- read.csv("useex.csv")

ggplot(rbind(cbind(df1, group='FAIR-only'), cbind(df2, group='RMA Forcing')), aes(time, T, colour=group)) +
    geom_line() + xlab(NULL) + ylab(NULL) + theme_bw()

df <- read.csv("compare.csv")
df2 <- melt(df, 'Year')
df2$scen <- ifelse(grepl("CP.Base", df2$variable), 'CP-Base', 'RCP4.5')
df2$gas <- ifelse(grepl("CO2", df2$variable), 'CO2',
           ifelse(grepl("CH4", df2$variable), 'CH4', 'FEx'))

ggplot(df2, aes(Year, value, colour=scen)) +
    facet_wrap(~ gas, scales='free_y', ncol=1) + geom_line() +
    theme_bw() + xlab(NULL) + ylab(NULL)
