library(data.table)
library(summarytools)
library(psych)
library("ggpubr")

path <- "~/Fig1/data/"
file <- "mocha.calls.filtered.txt"
mocha <- read.table(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
head(mocha)
mocha$p_arm[mocha$p_arm=="Y"]<-"1"
mocha$p_arm[mocha$p_arm=="T"]<-"1"
mocha$p_arm[mocha$p_arm=="N"]<-"0"
mocha$q_arm[mocha$q_arm=="Y"]<-"1"
mocha$q_arm[mocha$q_arm=="T"]<-"1"
mocha$q_arm[mocha$q_arm=="N"]<-"0"
mocha$arm<-paste(mocha$p_arm,mocha$q_arm,sep="")
mocha$arm[mocha$arm=="10"]<-"p"
mocha$arm[mocha$arm=="01"]<-"q"
mocha$arm[mocha$arm=="11"]<-"pq"
mocha$event<-paste(mocha$chrom,mocha$arm, mocha$type,sep="_")
mocha$chrom_arm<-paste(mocha$chrom,mocha$arm,sep="_")
mocha$nr<-"1"
mocha$nr<-as.numeric(mocha$nr)
###
mocha %>% nrow %>% paste('events')
mocha %>% pull(Sample_ID) %>% unique %>% length %>% paste('unique individuals with CNV')
mocha %>% 
  dplyr::count(type) %>% mutate(prop = signif(n/sum(n), 2))
summary(mocha$cf)

###now non-same-arm events count these separately
library(tidyverse)
mocha_distinct <- mocha  %>% group_by(Sample_ID) %>% distinct(event, .keep_all = TRUE)
mocha_distinct<- mocha_distinct %>% ungroup
mocha_distinct %>% nrow %>% paste('distinct events')
mocha_distinct %>% pull(Sample_ID) %>% unique %>% length %>% paste('unique individuals with CNV')
mocha_distinct %>% 
  dplyr::count(type) %>% mutate(prop = signif(n/sum(n), 2))
summary(mocha$cf)

mocha_distinct %>%
  filter(chrom != 23) %>%
  dplyr::count(Sample_ID) %>% dplyr::count(n) %>% mutate(prop = signif(nn/sum(nn), 2))

mocha$cf<-as.numeric(mocha$cf)

cat('Cell fraction ranges:', mocha$cf %>% range %>% signif(2) %>% paste(collapse = '-') %>% paste('\n'))
cat('median cell fraction:', mocha$cf %>% median %>% signif(2) %>% paste('\n'))

chr17p<- mocha %>%
  filter(chrom_arm %in% c("17_p", "17_pq")) %>%
  dplyr::count(Sample_ID) %>% dplyr::count(n) %>% mutate(prop = signif(nn/sum(nn), 2))

CN_LOH <- mocha %>%
  filter(type == "CN-LOH") %>%
  dplyr::count(Sample_ID) %>% dplyr::count(n) %>% mutate(prop = signif(nn/sum(nn), 2))

CN_LOH_17 <- CN_LOH %>%
  filter(chrom_arm == "17_p")%>%
  dplyr::count(Sample_ID) %>% dplyr::count(n) %>% mutate(prop = signif(nn/sum(nn), 2))


CN_LOH_9<- CN_LOH %>%
  filter(chrom_arm == "9_p")%>%
  dplyr::count(Sample_ID) %>% dplyr::count(n) %>% mutate(prop = signif(nn/sum(nn), 2))


chr21_cases<-mocha %>%
  filter(Sample_ID %in% c("3","4","16","17","20","55","31","37","46","54","14","70","71","72","73","63"))

####REMOVING chr21from the calculaiton of nr of CNAs
chr_21cases<-chr21_cases %>% filter(chrom!= 21)
count21<- chr21_cases %>% group_by(Sample_ID) %>% distinct(event, .keep_all = TRUE) %>% 
  dplyr::count(Sample_ID)

no_chr21_cases<-mocha %>%
  filter(!Sample_ID %in% c("3","4","16","17","20","55","31","37","46","54","14","70","71","72","73","63"))

no21_count<-no_chr21_cases %>% group_by(Sample_ID) %>% distinct(event, .keep_all = TRUE) %>% 
  dplyr::count(Sample_ID)
  
# add 10 no_event cases

Sample_ID<-c("11","13","15","33","36","47","50","64","69","86")

no_event<-as.data.table(Sample_ID)
no_event$n<-"0"
rownames(no_event)<-NULL
no21_count<-rbind(no_event,no21_count)

no21_count$chr21<-"0"
count21$chr21<-"1"
all<-rbind(no21_count,count21)
all$nrevents<-as.numeric(all$n)
group_by(all, chr21) %>%
  dplyr::summarise(
    count = n(),
    median = median(nrevents, na.rm = TRUE),
    IQR = IQR(nrevents, na.rm = TRUE),
    range=range(nrevents, na.rm = TRUE)
  )


all$chr21status<-all$chr21
all$chr21status[all$chr21status=="0"]<-"no_chr21amp"
all$chr21status[all$chr21status=="1"]<-"chr21amp"

levels(all$chr21status)<-rev(levels(all$chr21status))

plot<-ggboxplot(all, x = "chr21status", y = "nrevents", color = "black", 
          color.palette= c("black"),
          fill = "chr21status", 
          palette = c("blue", "darkred"),
          ylab = "Nr_CNAs", xlab = "Chr21 status", size=0.6)

plot+ theme(axis.line = element_line(size = 1))

res <- wilcox.test(nrevents ~ chr21, data = all,
                   exact = FALSE)
res
res$p.value

no21_count$n
summary(count21$n)
no21_count$n<-as.numeric(no21_count$n)
summary(no21_count$n)
all
table(all$chr21status)

write.table(all , file='/Users/charlotteb/Documents/DPhil/iAMP21_manuscript/data_check/github/data/Fig1_EDFig1/1F_all.tsv', sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
