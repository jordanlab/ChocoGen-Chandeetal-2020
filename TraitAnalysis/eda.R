d = read.table("clm_vs_cho.txt", header=T, sep="\t", quote="")
d$polarizedFst = d$Fst
d$polarizedFst[d$Ref.Freq.CLM < d$Ref.Freq.CHO] = d$polarizedFst[d$Ref.Freq.CLM < d$Ref.Freq.CHO]*-1

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

library(Hmisc)

set = unique(d$HigherLevel)
cols = as.character(jet.colors(length(set)))

pdf("FstDist.pdf", useDingbats = F)
hist(d$polarizedFst, freq = T, las=1, breaks = 20, xlim=c(-0.6, 0.6), xlab="Polarized Fst CHO -> 0 -> CLM", main="")
dev.off()



pdf("all.pdf", useDingbats = F)
plot(d$polarizedFst[abs(d$polarizedFst) >= 0.25], log(d$OR[abs(d$polarizedFst) >= 0.25]), las=1, pch=19, cex=1, col=cols, xlab = "Polarized Fst (CHO -> 0 -> CLM)", ylab = "Log Odds Ratio")
points(d$polarizedFst[abs(d$polarizedFst) < 0.25], log(d$OR[abs(d$polarizedFst) < 0.25]), las=1, pch=19, cex=0.5, col=alpha(cols,0.5))
legend("topright", legend = unique(d$HigherLevel), fill=cols)
dev.off()


pdf("all.pdf", useDingbats = F)
plot(d$polarizedFst, log(d$OR), las=1, pch=19, cex=0.5, col=cols, xlab = "Polarized Fst (CHO -> 0 -> CLM)", ylab = "Log Odds Ratio")
legend("topright", legend = unique(d$HigherLevel), fill=cols)
dev.off()



for(i in 1:length(set)){
  s = d[d$HigherLevel == set[i],]
  #print(s[1,])
  pdf(paste0(set[i], ".pdf"), useDingbats = F)
  plot(s$polarizedFst, log(s$OR), las=1, pch=19, cex=1, col=cols[i], main=set[i], xlim = c(-0.5, 0.5), xlab = "Polarized Fst (CHO -> 0 -> CLM)", ylab = "Log Odds Ratio")
  abline(h = 0)
  # abline(v = 0.25)
  # abline(v = -0.25)
  dev.off()
}



d = read.table("highConf-v3.txt", sep="\t", header=T, quote="", na.strings = "")

pdf("bar.pdf", useDingbats = F)
barplot(d$Polarized.Fst, names.arg = d$Trait, las=2, ylim=c(-0.6, 0.4), col=c(rep("blue",10), rep("orange", 10)), ylab="Polarized FST")
dev.off()

pdf("bar-freq.pdf", useDingbats = F)
barplot(t(cbind(d$EffectAlleleCLM, d$EffectAlleleCHO)), beside=T, names.arg = d$Trait, las=2, col=c("green", "purple"), ylim=c(0,1))
dev.off()

library(Hmisc)
library(gplots)
cols=colorRampPalette(c("blue", "black","red"))

pdf("heatmap.pdf", useDingbats = F)
heatmap.2(cbind(d$EffectAlleleCLM, d$EffectAlleleCHO), col=cols(250), trace="none", density.info=c("none"), Rowv=FALSE, Colv = FALSE, dendrogram = "none")
dev.off()
