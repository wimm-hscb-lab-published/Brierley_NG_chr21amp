library(ggplot2)
library(Seurat)
library(ggpubr)
library(tidyr)
##Download seurat object from zenodo/chr21amp_seurat.rds

samples_cd34 <- readRDS(file = "~/chr21amp_seurat.rds")
##DYRK1A violins


Idents(samples_cd34)<-factor(samples_cd34$prelim_celltype_broad_gb, levels=c("HSC","MPP",   "MkEP",     "GMP_early", "GMP_late" , "GMP_cycling","EryP_early",   "EryP_late"))

c = VlnPlot(samples_cd34, features = "DYRK1A", idents = c("HSC","MPP",   "MkEP",     "GMP_early", "GMP_late" , "GMP_cycling","EryP_early",   "EryP_late"), split.by = "chr21status_id", 
            cols=c("#73a871","#0D0E47","#d8443c"), pt.size=0.2) +
   geom_boxplot(position=position_dodge(1), lwd=0.2, color="white", width =0.3)+ #position=position_dodge(1), lwd=0.05, color="white")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


c
###
table(samples_cd34@meta.data$identity,samples_cd34$HCvsAML)

DYRK1A_exp<- GetAssayData(samples_cd34, assay = "RNA")["DYRK1A",]
DYRK1A_exp<-as.data.frame(DYRK1A_exp)
Samples.BL<-AddMetaData(samples_cd34, metadata=DYRK1A_exp)

###calc nr of cells expressing 
library(magrittr)
library(dplyr)


DYRK1Adata<-data.frame("x" = Samples.BL@meta.data[,"prelim_celltype_broad_gb"],
                       "y" = Samples.BL$DYRK1A_exp,
                       "Pos" = as.logical(Samples.BL$DYRK1A_exp > 0),
                       "grp" = Samples.BL@meta.data[,"chr21status_id"])


plot.sum <- DYRK1Adata %>% 
  dplyr::group_by(x, grp) %>% 
  dplyr::summarise(M = median(y),
                   x.lab = paste0(unique(grp), 
                                  ": ", sum(Pos), "/", length(Pos), " (", 
                                  round(mean(Pos)*100, 1), "%)"),
                   Pos = mean(Pos))


x.labs <- aggregate(data.frame(plot.sum$x.lab),
                    by = list(plot.sum$x),
                    function(x) paste(x, collapse = "\n")) %>%
  mutate(x.lab = paste(Group.1, plot.sum.x.lab, sep = "\n")) %>%
  dplyr::select(x.lab)



####cell type % change


table_cd34_2<-prop.table(table(md$identity,md$prelim_celltype_broad_gb),margin=1)

table_cd34_2<-as.data.frame(table_cd34_2)

library(tidyr)

colnames(table_cd34_2)<-c("identity", "Cell_type", "Freq")
table_cd34_2$percent <- (table_cd34_2$Freq*100)
table_cd34_2<-table_cd34_2[table_cd34_2$Freq != 0, ]

table_cd34_2$chr21status<-table_cd34_2$identity
table_cd34_2$chr21status<-recode(table_cd34_2$chr21status, "HC_BAS33"="HC","HC_BAS78"="HC","HC_BAS84"="HC","HC_BAS86"="HC","HC_BAS93"="HC","TNO13_BL"="BP_MPN_no_chr21",
                                 "TNO10_BL"="BP_MPN_no_chr21","TNO17_BL"="BP_MPN_chr21","TNO19_BL"="BP_MPN_no_chr21",
                                 "TNO21_BL"="BP_MPN_no_chr21","TNO28_BL"="BP_MPN_no_chr21","TNO37_BL"="BP_MPN_chr21",
                                 "TNO42_BL"="BP_MPN_no_chr21","TNO48_BL"="BP_MPN_no_chr21","TNO49_BL"="BP_MPN_no_chr21",
                                 "TNO5_BL"="BP_MPN_no_chr21")

table_cd34_2$Cell_type<-factor(table_cd34_2$Cell_type, levels=(c("HSC","MPP",   "MkEP",     "GMP_early", "GMP_late" , "GMP_cycling","EryP_early",   "EryP_late")))
library(ggpubr)
cols=c("#73a871","#0D0E47","#d8443c")
p1<-ggplot(data=table_cd34_2, aes(x=chr21status, y=percent,  color=chr21status))+ #fill=chr21status,
  geom_boxplot() +
  scale_color_manual(values=cols) +
  facet_wrap(~Cell_type,scale="free",nrow = 2)+
  geom_jitter(color="black", size=0.4, alpha=0.2) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = TRUE)+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + ylab("percent of CD34+ cells") +
  ggtitle("Cell type distribution by chr21amp disease status") +
  xlab("") + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

p1


####cell type % change


celltype_MPP <-table_cd34_2 %>% filter(Cell_type == "MPP")
kruskal.test(celltype_MPP$Freq ~ celltype_MPP$chr21status, data = celltype_MPP)


celltype_MPP <-table_cd34_2 %>% filter(Cell_type == "EryP_late")
kruskal.test(celltype_MPP$Freq ~ celltype_MPP$chr21status, data = celltype_MPP)

# > sessionInfo()
# R version 4.3.0 (2023-04-21)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /ceph/package/u22/R-base/4.3.0/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/London
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.6.0           tidyr_1.3.1            plyr_1.8.9             
# [5] ggrepel_0.9.6          data.table_1.16.4         Matrix_1.7-1          
# [9] dplyr_1.1.4            Seurat_4.4.0           ggplot2_3.5.1          
# [13] sp_2.1-4              

# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.20              splines_4.3.0                
# [3] bitops_1.0-7                  filelock_1.0.2               
# [5] tibble_3.2.1                  polyclip_1.10-4              
# [7] lifecycle_1.0.3               globals_0.16.2               
# [9] lattice_0.21-8                MASS_7.3-60                  
# [11] magrittr_2.0.3                plotly_4.10.2                
# [13] yaml_2.3.7                    httpuv_1.6.11                
# [15] sctransform_0.4.0             sp_2.1-1                     
# [17] spatstat.sparse_3.0-2         reticulate_1.34.0            
# [19] cowplot_1.1.1                 pbapply_1.7-2                
# [21] DBI_1.1.3                     RColorBrewer_1.1-3           
# [23] abind_1.4-5                   zlibbioc_1.46.0              
# [25] Rtsne_0.16                    GenomicRanges_1.52.0         
# [27] purrr_1.0.2                   BiocGenerics_0.46.0          
# [29] RCurl_1.98-1.12               rappdirs_0.3.3               
# [31] circlize_0.4.15               GenomeInfoDbData_1.2.10      
# [33] IRanges_2.34.1                S4Vectors_0.38.1             
# [35] ggrepel_0.9.3                 irlba_2.3.5.1                
# [37] listenv_0.9.0                 spatstat.utils_3.0-3         
# [39] goftest_1.2-3                 spatstat.random_3.1-5        
# [41] fitdistrplus_1.1-11           parallelly_1.36.0            
# [43] DelayedMatrixStats_1.22.1     leiden_0.4.3                 
# [45] codetools_0.2-19              DelayedArray_0.26.3          
# [47] tidyselect_1.2.0              shape_1.4.6                  
# [49] farver_2.1.1                  viridis_0.6.3                
# [51] matrixStats_1.0.0             stats4_4.3.0                 
# [53] BiocFileCache_2.8.0           spatstat.explore_3.2-1       
# [55] jsonlite_1.8.5                ellipsis_0.3.2               
# [57] progressr_0.13.0              ggridges_0.5.4               
# [59] survival_3.5-5                tools_4.3.0                  
# [61] ica_1.0-3                     Rcpp_1.0.10                  
# [63] glue_1.6.2                    gridExtra_2.3                
# [65] MatrixGenerics_1.12.2         GenomeInfoDb_1.36.1          
# [67] withr_2.5.0                   BiocManager_1.30.21          
# [69] fastmap_1.1.1                 fansi_1.0.4                  
# [71] digest_0.6.32                 R6_2.5.1                     
# [73] mime_0.12                     colorspace_2.1-0             
# [75] scattermore_1.2               tensor_1.5                   
# [77] spatstat.data_3.0-1           RSQLite_2.3.1                
# [79] celldex_1.10.1                utf8_1.2.3                   
# [81] tidyr_1.3.1                   generics_0.1.3               
# [83] data.table_1.14.8             httr_1.4.6                   
# [85] htmlwidgets_1.6.2             S4Arrays_1.2.0               
# [87] uwot_0.1.15                   pkgconfig_2.0.3              
# [89] gtable_0.3.3                  blob_1.2.4                   
# [91] lmtest_0.9-40                 XVector_0.40.0               
# [93] htmltools_0.5.5               scales_1.3.0                 
# [95] Biobase_2.60.0                png_0.1-8                    
# [97] rstudioapi_0.16.0             reshape2_1.4.4               
# [99] nlme_3.1-162                  curl_5.0.1                   
# [101] cachem_1.0.8                  zoo_1.8-12                   
# [103] GlobalOptions_0.1.2           stringr_1.5.0                
# [105] BiocVersion_3.17.1            KernSmooth_2.23-21           
# [107] parallel_4.3.0                miniUI_0.1.1.1               
# [109] vipor_0.4.5                   AnnotationDbi_1.62.1         
# [111] ggrastr_1.0.1                 pillar_1.9.0                 
# [113] grid_4.3.0                    vctrs_0.6.5                  
# [115] RANN_2.6.1                    promises_1.2.0.1             
# [117] dbplyr_2.3.2                  xtable_1.8-4                 
# [119] cluster_2.1.4                 beeswarm_0.4.0               
# [121] cli_3.6.2                     compiler_4.3.0               
# [123] rlang_1.1.3                   crayon_1.5.2                 
# [125] future.apply_1.11.0           labeling_0.4.2               
# [127] plyr_1.8.8                    ggbeeswarm_0.7.2             
# [129] stringi_1.7.12                viridisLite_0.4.2            
# [131] deldir_1.0-9                  munsell_0.5.0                
# [133] Biostrings_2.68.1             lazyeval_0.2.2               
# [135] spatstat.geom_3.2-1           Matrix_1.6-1                 
# [137] ExperimentHub_2.8.0           patchwork_1.2.0              
# [139] sparseMatrixStats_1.12.1      bit64_4.0.5                  
# [141] future_1.33.2                 KEGGREST_1.40.0              
# [143] shiny_1.7.4                   SummarizedExperiment_1.30.2  
# [145] interactiveDisplayBase_1.38.0 AnnotationHub_3.8.0          
# [147] ROCR_1.0-11                   igraph_1.5.0                 
# [149] memoise_2.0.1                 bit_4.0.5  
