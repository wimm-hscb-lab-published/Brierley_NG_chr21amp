library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(RColorBrewer)
library(viridis)
library(ggstar)
library(tidyr)
path<-"~/Fig1/data/"
file<-"mydata.coeff.tsv"

mydata.coeff<- as.matrix(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))
file <- "mydata.p.tsv"
mydata.p <-as.matrix(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))
rownames(mydata.coeff)<-colnames(mydata.coeff)
rownames(mydata.p)<-colnames(mydata.p)


corrplot(mydata.coeff, is.corr = FALSE, type=c("upper"), method = "circle", p.mat=mydata.p, sig.level = .05, insig="label_sig", col = brewer.pal(n = 11, name = "RdBu"))

# R version 4.0.4 (2021-02-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggstar_1.0.4       viridis_0.6.5      viridisLite_0.4.2  corrplot_0.92     
# [5] tibble_3.2.1       tidyr_1.3.1        cmprsk_2.2-11      haven_2.5.4       
# [9] ggsurvfit_1.1.0    reshape2_1.4.4     RColorBrewer_1.1-3 GenVisR_1.20.0    
# [13] readxl_1.4.3       readr_2.1.5        kableExtra_1.4.0   ggrepel_0.9.5     
# [17] cowplot_1.1.3      IRdisplay_1.1      survminer_0.4.9    ggpubr_0.4.0      
# [21] survival_3.6-4     patchwork_1.2.0    pec_2020.11.17     prodlim_2023.08.28
# [25] ggforce_0.4.2      glue_1.7.0         ggbeeswarm_0.7.2   ggplot2_3.5.1     
# [29] data.table_1.15.4  stringr_1.5.1      dplyr_1.1.4   