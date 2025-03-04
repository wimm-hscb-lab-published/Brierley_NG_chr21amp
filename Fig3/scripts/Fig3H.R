

library(data.table)
library(EnhancedVolcano)
library(dplyr)


###load data from github
df<- read.table("~/Fig3/Peaks_in_DA.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

df$DYRK1Apeak<-df$SYMBOL
df$DYRK1Apeak<-ifelse((df$SYMBOL=="DYRK1A" & df$log2FoldChange>1 &df$padj<0.0501), "DYRK1A sig peak", "NA" )


df$DYRK1Apeak[df$SYMBOL=="DYRK1A" & df$log2FoldChange>1&df$padj>0.05]<-"DYRK1A log2FC >1"

genetype1 <- c('DYRK1A sig peak')
genetype2 <- c('DYRK1A log2FC >1')

keyvals.shape <- ifelse(
  df$DYRK1Apeak %in% genetype1, 17,
  ifelse(df$DYRK1Apeak %in% genetype2, 15,
         1))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 1] <- 'Non DYRK1A peak'
names(keyvals.shape)[keyvals.shape == 17] <- 'DYRK1A DA peak'
names(keyvals.shape)[keyvals.shape == 15] <-'DYRK1A log2FC >1 peak'

keyvals.colour <- ifelse(
  df$DYRK1Apeak %in% genetype1, 'red4',
  ifelse(df$DYRK1Apeak %in% genetype2, 'royalblue',
         'black'))
keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <-'DYRK1A log2FC >1 peak'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'Non DYRK1A peak'
names(keyvals.colour)[keyvals.colour == 'red4']<- 'DYRK1A DA peak'


p<-  EnhancedVolcano(df,
                     lab = df$SYMBOL,
                     x = 'log2FoldChange',
                     y = 'padj',
                     #selectLab = df$SYMBOL[which(names(keyvals.shape) %in% c('DYRK1A DA peak', 'DYRK1A log2FC >1 peak'))],
                     selectLab = "",
                     title = 'Differentially accessible peaks in DE genes',
                     subtitle = bquote(italic('chr21amp vs non-chr21amp MPNAML')),
                     xlab = bquote(~Log[2]~ 'fold change'),
                     shapeCustom= keyvals.shape,
                     colCustom = keyvals.colour,
                     pCutoff = 0.05,
                     FCcutoff = 1.0,
                     pointSize = 4.0,
                     labSize = 6.0,
                     xlim = c(-4,4),
                     ylim = c(0, 3),
                     labCol = 'black',
                     labFace = 'bold',
                                     colAlpha = 4/5,
                     legendPosition = 'right',
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                          max.overlaps=30)
p

# 
# > sessionInfo()
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
#   [1] dplyr_1.1.4           EnhancedVolcano_1.6.0 ggrepel_0.9.5        
# [4] ggplot2_3.5.1         data.table_1.15.4    
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.12         rstudioapi_0.16.0   magrittr_2.0.3     
# [4] tidyselect_1.2.1    munsell_0.5.1       colorspace_2.1-0   
# [7] R6_2.5.1            rlang_1.1.3         fansi_1.0.6        
# [10] tools_4.0.4         grid_4.0.4          gtable_0.3.5       
# [13] xfun_0.43           tinytex_0.51        utf8_1.2.4         
# [16] cli_3.6.2           withr_3.0.0         tibble_3.2.1       
# [19] lifecycle_1.0.4     BiocManager_1.30.23 vctrs_0.6.5        
# [22] glue_1.7.0          compiler_4.0.4      pillar_1.9.0       
# [25] generics_0.1.3      scales_1.3.0        pkgconfig_2.0.3  
# 
# 
