
library(AllelicImbalance)
library(devtools)
install_github("parklab/ShatterSeek")

searchArea <- GRanges(seqnames = c("21"), ranges = IRanges(37508643,40196388))
searchArea <- GRanges(seqnames = c("17"), ranges = IRanges(79478301,79478361))
pathToFiles <- system.file("extdata/ERP000101_subset",
                           package = "AllelicImbalance")
reads <- impBamGAL(pathToFiles, searchArea, verbose = FALSE)
heterozygotePositions <- scanForHeterozygotes(reads, verbose = FALSE)
countList <- getAlleleCounts(reads, heterozygotePositions, verbose = FALSE)
a.simple <- ASEsetFromCountList(heterozygotePositions, countList)
a.simple

#Getting searchArea from genesymbol
library(org.Hs.eg.db)
searchArea<-getAreaFromGeneNames("ACTG1",org.Hs.eg.db)
#Getting rs-IDs
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
## Warning: replacing previous import 'utils::findMatches' by
## 'S4Vectors::findMatches' when loading 'SNPlocs.Hsapiens.dbSNP144.GRCh37'
gr <- rowRanges(a.simple)
updatedGRanges<-getSnpIdFromLocation(gr, SNPlocs.Hsapiens.dbSNP144.GRCh37)
rowRanges(a.simple)<-updatedGRanges

#simulate phenotype data
pdata <- DataFrame(
  Treatment=sample(c("ChIP", "Input"),length(reads),replace=TRUE),
  Gender=sample(c("male", "female"),length(reads),replace=TRUE),
  row.names=paste("individual",1:length(reads),sep=""))
#make new ASEset with pdata
a.new <- ASEsetFromCountList(
  heterozygotePositions,
  countList,
  colData=pdata)
#add to existing object
colData(a.simple) <- pdata

#infer and add genotype require declaration of the reference allele
ref(a.simple) <- randomRef(a.simple)
genotype(a.simple) <- inferGenotypes(a.simple)
#access to genotype information requires not only the information about the
#reference allele but also the alternative allele
alt(a.simple) <- inferAltAllele(a.simple)
genotype(a.simple)[,1:4]
## [,1] [,2] [,3] [,4]
## [1,] "G/A" "G/A" NA "G/A"
## [2,] "T/G" "T/G" "T/G" "T/G"
## [3,] "C/G" "C/G" "C/G" "C/G"


barplot(a.stranded[2], strand="+", xlab.text="", legend.interspace=2)
barplot(a.simple[2], type="fraction", cex.main = 0.7)

library(org.Hs.eg.db)
barplot(a.simple[1],OrgDb=org.Hs.eg.db,
        cex.main = 0.7,
        xlab.text="",
        ypos.annotation = 0.08,
        annotation.interspace=1.3,
        legend.interspace=1.5
)

searchArea <- GRanges(seqnames = c("chr21"), ranges = IRanges(37508643,40196388))

###upload BAM files from RNAseq dataset
reads1 <- impBamGAL(pathToFiles, searchArea, verbose = FALSE)
heterozygotePositions1 <- scanForHeterozygotes(reads1, verbose = FALSE)
name_het<-names(heterozygotePositions1)
name_het<-substring(name_het, 4)
names(heterozygotePositions1)<-name_het

countList1 <- getAlleleCounts(reads1, heterozygotePositions1, verbose = FALSE)
a.simple1 <- ASEsetFromCountList(heterozygotePositions1, countList1)
a.simple1
