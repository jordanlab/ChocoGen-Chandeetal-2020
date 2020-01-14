setwd("D:/Dropbox/ColombianGenomeAnalysis/ChocoGen/Analysis/ChocoGen3/GlobalAncestry")
rm(list=ls())

admixture = read.table("auto.admixture", header = T, sep="\t")
mapping = read.table("sampleInfo.txt", header=T, sep="\t")
data = merge(admixture, mapping, by.x="IND", by.y="Ind")

clm = data[data$Pop == "CLM", 2:4]
chg = data[data$Pop == "CHG", 2:4]

library(ade4)

pdf("clm.triPlot.pdf", useDingbats=F)
plot = triangle.plot(clm, show=F, draw.line=F, min3 = c(0,0,0), max3 = c(1,1,1))
points(plot, col="green", cex=1, pch=19)
dev.off()

pdf("chg.triPlot.pdf", useDingbats=F)
plot = triangle.plot(chg, show=F, draw.line=F, min3 = c(0,0,0), max3 = c(1,1,1))
points(plot, col="purple", cex=1, pch=19)
dev.off()

