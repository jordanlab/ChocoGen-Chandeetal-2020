#!/usr/bin/env Rscript
library("optparse")
 
option_list = list(
  make_option(c("-d", "--dir"), type="character", default="./", 
              help="Directory containing PRS text files [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="ImputedMatchedSNP-PRS", 
              help="Output file prefix [default= %default]", metavar="character")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(data.table)
library(dplyr)
readGrs = function(fileName, sampleInfo){
  grs = my.read.table(fileName, comment.char = "##")
  colnames(grs) = c("IND", "GRS", "SnpCount")
  grs = merge(grs, sampleInfo, by="IND")
}

my.read.table = function(fileName, comment.char){
  clean.lines = sub(paste0(comment.char, ".*"), "", readLines(fileName))
  clean.lines = sub("SAMPLE.*", "", clean.lines)
  read.table(text = paste(clean.lines, collapse = "\n"))
}
# metadata files
sampleInfo = read.csv("Cho-CLM_sampleinfo.csv", fileEncoding="latin1")
supInfo = read.csv("super_population_info.csv",fileEncoding="latin1")
supInfo = supInfo[, c(1,3)]
colnames(supInfo) = c("POP", "SUP")
sampleInfo = merge(sampleInfo, supInfo, by="POP")
supNames = c("AFR", "EUR", "AMR", "SAS", "EAS")
c5 = c("#484496", "#F4A500", "#328A4F", "#944116", "#D92414")
names(c5) = supNames

fl = list.files(opt$dir, "*.txt")
temp = list()

for (file in fl){
  t = readGrs(file, sampleInfo)
  temp[[i]] = t %>%
    group_by(POP) %>%
    dplyr::summarize(Mean = mean(GRS, na.rm=TRUE)) %>%
    mutate(Trait = gsub(".txt", "", file)) %>%
    mutate(Diff = Mean - lag(Mean)) %>%
    filter(!is.na(Diff)) %>%
    select(-POP, -Mean)
    temp[[i]]$p.val = t %>%
      group_by(POP) %>%
      nest() %>%
      spread(key = POP, value = data) %>% 
      mutate(
        t_test = map2(CHO, CLM, ~{t.test(.x$GRS, .y$GRS) %>% tidy()}),
        CHO = map(CHO, nrow),
        CLM = map(CLM, nrow)
      ) %>% 
      unnest() %>% select(p.value) %>% unlist() %>% unname()
    i = i+1
  }

df = do.call(rbind.data.frame, temp) %>% arrange(-Diff) %>% mutate(p.adj = p.adjust(p.val))
df.filt <- df %>% filter(p.adj < 0.05)
write.table(df, paste0(opt$out, "-ImputedMatchedSNP-PRS.tsv"), sep="\t", row.names=F, quote=F)
write.table(df.filt, paste0(opt$out,"ImputedMatchedSNP-PRS.filt.tsv"), sep="\t", row.names=F, quote=F)


pdf(paste0(opt$out, "-ImputedMatchedSNP-PRS.pdf"), useDingbats=F)
df %>% 
  mutate(x = 1:length(Diff), Pop = ifelse(Diff > 0, "CLM", "CHO" )) %>% 
  ggplot(., aes(x=x, y=Diff, color=Pop, fill=Pop)) + 
  geom_bar(stat = "identity") + 
  scale_color_manual(values=c("mediumpurple", "green"), guide=FALSE) +
  scale_fill_manual(values=c("mediumpurple", "green"), guide=FALSE) +
  ggthemes::theme_few() + labs(y = "delta", x = "") + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
dev.off()

pdf(paste0(opt$out, "-ImputedMatchedSNP-PRS.filt.pdf"), useDingbats=F)
df.filt %>% 
  mutate(x = 1:length(Diff), Pop = ifelse(Diff > 0, "CLM", "CHO" )) %>% 
  ggplot(., aes(x=x, y=Diff, color=Pop, fill=Pop)) + 
  geom_bar(stat = "identity") + 
  scale_color_manual(values=c("mediumpurple", "green"), guide=FALSE) +
  scale_fill_manual(values=c("mediumpurple", "green"), guide=FALSE) +
  ggthemes::theme_few() + labs(y = "delta", x = "") + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
dev.off()



return
## Scratch space
readGrs(file, sampleInfo) %>%
  group_by(POP) %>%
  nest() %>%
  spread(key = POP, value = data) %>% 
  mutate(
    t_test = map2(CHO, CLM, ~{t.test(.x$GRS, .y$GRS) %>% tidy()}),
    CHO = map(CHO, nrow),
    CLM = map(CLM, nrow)
  ) %>% 
  unnest() %>% select(p.value)