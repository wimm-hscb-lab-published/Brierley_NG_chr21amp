####
library(plyr)
library(pheatmap)


############################################################
####################### GSEA HEATMAP #######################
############################################################

###Supplementary table 8 and 9, also provided in github as split data per condition
path<-"~/Fig4/data/"
# Define GSEA result files
files <- c(
  "gsea_report_for_chr21amp_1688635164519.tsv",
           "gsea_report_for_nonchr21amp_1688635164519.tsv",
           "gsea_report_for_chr21amp_MPNAML_1688635146171.tsv",
           "gsea_report_for_HC_1688635146171.tsv"
          )
labels <- c("chr21amp vs non-chr21amp",
            "non-chr21amp vs chr21amp",
            "chr21amp vs HC",
            "HC vs chr21amp "
            )

df.list <- list()

for(i in 1:length(files)) {
  
  file <- files[i]
  df <- read.table(paste(path,file,sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
  
  . <- grep("|", df$pathway, fixed=TRUE, value=TRUE)
  
  if(length(.) != 0) {
    
    df <- df[-grep("|", df$pathway, fixed=TRUE), ]
    
  }
  
  
  df.list[[i]] <- df
  
}

length(df.list)

# Tabulate padj
. <- lapply(df.list, function(x) {(x[,c("NAME", "NOM.p.val")])})
names(.[[1]])[2] <- labels[1]
names(.[[2]])[2] <- labels[2]
names(.[[3]])[2] <- labels[3]
names(.[[4]])[2] <- labels[4]

. <- Reduce(function(x,y) join(x=x, y=y, by="NAME", type="full"), .)
row.names(.) <- .$NAME
.$NAME <- NULL
padj <- .

# Tabulate nes
. <- lapply(df.list, function(x) {(x[,c("NAME", "NES")])})
names(.[[1]])[2] <- labels[1]
names(.[[2]])[2] <- labels[2]
names(.[[3]])[2] <- labels[3]
names(.[[4]])[2] <- labels[4]

. <- Reduce(function(x,y) join(x=x, y=y, by="NAME", type="full"), .)
row.names(.) <- .$pathway
.$pathway <- NULL
nes <- .


row.names(nes) <- nes$NAME
nes$NAME <- NULL
# Censor NES < 1
threshold <- 1
nes[abs(nes) < threshold] <- NA



# Remove non-sig. gene sets
. <- apply(nes, 1, function(x) {sum(!is.na(x)) != 0})
nes <- nes[.,]

# Reorder for aesthetic purpose
nes$sum <- rowSums(nes, na.rm=TRUE)
nes <- nes[order(nes$sum, decreasing=TRUE), ]
nes$sum <- NULL

# Tidy gene set names

row.names(nes) <- gsub("HALLMARK_", "", row.names(nes))
row.names(nes) <- gsub("_", " ", row.names(nes))

# Subset selected pathways

data_vhc<-nes[,c(-1,-2)]
data_vnon<-nes[,c(-3,-4)]

# Remove non-sig. gene sets
. <- apply(data_vhc, 1, function(x) {sum(!is.na(x)) != 0})
data_vhc <- data_vhc[.,]

. <- apply(data_vnon, 1, function(x) {sum(!is.na(x)) != 0})
data_vnon <- data_vnon[.,]


hc_no_pathways<-c("KRAS SIGNALING DN","COAGULATION","MITOTIC SPINDLE","DAUER STAT3 TARGETS DN")
data_vhc<-data_vhc %>% filter(!rownames(data_vhc) %in% hc_no_pathways)
xx<-pheatmap(data_vhc, cluster_rows=FALSE, cluster_cols=FALSE, scale="none", fontsize_row=9, border_color="white", legend=TRUE, angle_col="45",cellwidth = 15, cellheight = 15, # <--- changed here
             width = 15, height=11)



chr21_no_pathways<-c("KRAS SIGNALING DN","HYPOXIA","TGF BETA SIGNALING","PI3K AKT MTOR SIGNALING","KRAS SIGNALING UP","HEDGEHOG SIGNALING","MITOTIC SPINDLE","KEGG WNT SIGNALING PATHWAY")
data_vnon<-data_vnon %>% filter(!rownames(data_vnon) %in% chr21_no_pathways)

xx<-pheatmap(data_vnon, cluster_rows=FALSE, cluster_cols=FALSE, scale="none", fontsize_row=9, border_color="white", legend=TRUE, angle_col="45",cellwidth = 20, cellheight = 20, # <--- changed here
             width = 7, height=9.1)

xx


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
#   [1] pheatmap_1.0.12       plyr_1.8.9            dplyr_1.1.4          
# [4] EnhancedVolcano_1.6.0 ggrepel_0.9.5         ggplot2_3.5.1        
# [7] data.table_1.15.4    
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.12         rstudioapi_0.16.0   magrittr_2.0.3     
# [4] tidyselect_1.2.1    munsell_0.5.1       colorspace_2.1-0   
# [7] R6_2.5.1            rlang_1.1.3         fansi_1.0.6        
# [10] tools_4.0.4         grid_4.0.4          gtable_0.3.5       
# [13] xfun_0.43           tinytex_0.51        utf8_1.2.4         
# [16] cli_3.6.2           withr_3.0.0         tibble_3.2.1       
# [19] lifecycle_1.0.4     RColorBrewer_1.1-3  BiocManager_1.30.23
# [22] vctrs_0.6.5         glue_1.7.0          compiler_4.0.4     
# [25] pillar_1.9.0        generics_0.1.3      scales_1.3.0       
# [28] pkgconfig_2.0.3    