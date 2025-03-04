

#####Fig4K 

library(AUCell)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)

##read in datasets from github
targeting_DYRK1A <-read.table('~/Fig4/data/GRN_DYRK1A.tsv', sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
topRegulators<-read.table('~/Fig4/data/topRegulators_scaled_HSPC_BL_gbclass.tsv', sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
regulonActivity_byCellType_Scaled<-read.table('~/Fig4/data/regulonActivity_byCellType_Scaled.tsv', sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)




topRegulators$TF<-topRegulators$Regulon

chr21celltype<-c("HSC_MPNAML_chr21amp","MPP_MPNAML_chr21amp",  "GMP_early_MPNAML_chr21amp" ,"GMP_late_MPNAML_chr21amp", "MkEP_MPNAML_chr21amp")

topRegulators$TF<-substr(topRegulators$TF,1,nchar(topRegulators$TF)-3)

dyrk_regs<-topRegulators %>% filter(TF %in% targeting_DYRK1A$TF)

stats<-topRegulators %>% filter(Regulon %in% c("STAT5B(+)","TP53(+)") ) %>%slice(1:2)
dyrk_regs<-dyrk_regs %>% filter(CellType %in% chr21celltype)


dyrk_regs<-dyrk_regs %>% 
  arrange(desc(RelativeActivity)) %>% 
  group_by(CellType) %>% dplyr::slice(1:10)
dyrk_regs<-rbind(dyrk_regs,stats)
dyrk_regs<-unique(dyrk_regs$Regulon)

celltypesumm<-c(   
  
  "GMP_early_HC_no_chr21amp"  ,      "GMP_early_MPNAML_chr21amp"   ,   
  "GMP_early_MPNAML_no_chr21amp",    "GMP_late_HC_no_chr21amp" ,       
  "GMP_late_MPNAML_chr21amp"   ,     "GMP_late_MPNAML_no_chr21amp"  ,  
  "HSC_HC_no_chr21amp"     ,         "HSC_MPNAML_chr21amp"    ,        
  "HSC_MPNAML_no_chr21amp"    ,           
  "MkEP_HC_no_chr21amp"     ,        "MkEP_MPNAML_chr21amp" ,          
  "MkEP_MPNAML_no_chr21amp"    ,     "MPP_HC_no_chr21amp"    ,         
  "MPP_MPNAML_chr21amp"     ,        "MPP_MPNAML_no_chr21amp") 
celltype<-c(  
  "GMP_early"  ,  "GMP_early"  , "GMP_early"  , 
  "GMP_late"  , "GMP_late"  ,"GMP_late"  ,
  "HSC"     ,      "HSC"     ,    "HSC"     ,    
  "MkEP"     ,   "MkEP"     , "MkEP"     , 
  "MPP"     ,        "MPP",        "MPP" )  
chr21status<-   c(
  
  "no_chr21amp","chr21amp",  "no_chr21amp",
  "no_chr21amp","chr21amp",  "no_chr21amp",
  "no_chr21amp","chr21amp",  "no_chr21amp",
  "no_chr21amp","chr21amp",  "no_chr21amp",
  "no_chr21amp","chr21amp",  "no_chr21amp")
HCvsAML<-   c(
  
  "HC", "MPNAML",  "MPNAML", 
  "HC", "MPNAML",  "MPNAML", 
  "HC", "MPNAML",  "MPNAML", 
  "HC", "MPNAML",  "MPNAML", 
  "HC", "MPNAML",  "MPNAML")

df.pheno<-as.data.frame(celltypesumm)
df.pheno$celltype<-celltype
df.pheno$chr21status<-chr21status
df.pheno$HCvsAML<-HCvsAML
rownames(df.pheno)<-df.pheno$celltypesumm
df.pheno$celltypesumm<-NULL

#####



regulonActivity_byCellType_Scaled2<-regulonActivity_byCellType_Scaled[dyrk_regs,]
regulonActivity_byCellType_Scaled2<-as.matrix(regulonActivity_byCellType_Scaled2)



colnames<-colnames(regulonActivity_byCellType_Scaled2)

df<-df.pheno %>% arrange(factor(rownames(df.pheno), levels = colnames))


df$HCvsAML[df$HCvsAML=="MPNAML"]<-"BP_MPN"

column_ha = HeatmapAnnotation(df=df,col=list(HCvsAML=c("HC"="blue","BP_MPN"="red"), chr21status=c("no_chr21amp"="#1fa187","chr21amp"="#AE017E"),celltype=c("GMP_early"="#0d0887","GMP_late"= "#fe9b00",  "HSC"="#eee82b","MkEP"="#2beeba","MPP"="#ee2bc1")), #c("GMP"="#0d0887","HSC"="#eee82b","MkEP"="#2beeba","MPP"="#ee2bc1","EryP"="#ee2b2e"),
                              annotation_name_gp= gpar(fontsize = 16, fontface="bold" ),
                              annotation_legend_param = list(
                                HCvsAML = list(
                                  title = "BP_MPN vs control",
                                  at = c("HC","BP_MPN"),
                                  labels = c("HC","BP_MPN")
                                ),
                                chr21status = list(
                                  title = "chr21 status",
                                  at = c("no_chr21amp","chr21amp"),
                                  labels = c("no_chr21amp","chr21amp")
                                ),
                                celltype = list(
                                  title = "cell type",
                                  #   at = c("EryP","GMP","HSC","MkEP","MPP" ),
                                  #  labels = c("EryP","GMP","HSC","MkEP","MPP" )
                                  at=c(          
                                    
                                    "GMP_early"  , 
                                    "GMP_late"  , 
                                    "HSC"     ,     
                                    "MkEP"     ,    
                                    "MPP"  ) ,
                                  labels=c(        
                                    
                                    "GMP_early"  , 
                                    "GMP_late"  , 
                                    "HSC"     ,     
                                    "MkEP"     ,    
                                    "MPP"    )
                                )))



colnames(regulonActivity_byCellType_Scaled2)<-c("GMP_early_HC_no_chr21amp"   ,  "GMP_early_BP_MPN_chr21amp" ,  
                                                "GMP_early_BP_MPN_no_chr21amp" ,"GMP_late_HC_no_chr21amp"   ,  
                                                "GMP_late_BP_MPN_chr21amp"   ,  "GMP_late_BP_MPN_no_chr21amp", 
                                                "HSC_HC_no_chr21amp"          , "HSC_BP_MPN_chr21amp"  ,       
                                                "HSC_BP_MPN_no_chr21amp"   ,    "MkEP_HC_no_chr21amp",         
                                                "MkEP_BP_MPN_chr21amp"    ,     "MkEP_BP_MPN_no_chr21amp"  ,   
                                                "MPP_HC_no_chr21amp"  ,         "MPP_BP_MPN_chr21amp",         
                                                "MPP_BP_MPN_no_chr21amp" )

ht=ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled2, name="Regulon activity",#col = viridis(n=12, option = "plasma"), 
                           cluster_rows = TRUE,cluster_columns = TRUE,
                           column_names_rot = 45, top_annotation = column_ha, heatmap_legend_param = list(labels_gp = gpar(fontsize = 16), 
                                                                                                          showcolumnnames =FALSE,#column_names_gp = grid::gpar(fontsize = 12), 
                                                                                                          row_names_gp=grid::gpar(fontsize=16,fontface="bold")))
draw(ht,padding = unit(c(2, 30, 2, 2), "mm")) # row font size


#dev.off() 
# 
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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] data.table_1.14.8     ComplexHeatmap_2.16.0 BiocParallel_1.34.2   plotly_4.10.2        
# [5] RColorBrewer_1.1-3    KernSmooth_2.23-21    AUCell_1.22.0         dplyr_1.1.4          
# [9] ggplot2_3.5.1         SeuratObject_4.1.4    Seurat_4.4.0          later_1.3.2          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.20              splines_4.3.0                 bitops_1.0-7                 
# [4] filelock_1.0.2                R.oo_1.25.0                   tibble_3.2.1                 
# [7] polyclip_1.10-4               graph_1.78.0                  XML_3.99-0.14                
# [10] lifecycle_1.0.3               doParallel_1.0.17             globals_0.16.2               
# [13] lattice_0.21-8                MASS_7.3-60                   magrittr_2.0.3               
# [16] yaml_2.3.7                    httpuv_1.6.11                 sctransform_0.4.0            
# [19] sp_2.1-1                      spatstat.sparse_3.0-2         reticulate_1.34.0            
# [22] cowplot_1.1.1                 pbapply_1.7-2                 DBI_1.1.3                    
# [25] abind_1.4-5                   zlibbioc_1.46.0               Rtsne_0.16                   
# [28] GenomicRanges_1.52.0          R.utils_2.12.2                purrr_1.0.2                  
# [31] BiocGenerics_0.46.0           RCurl_1.98-1.12               rappdirs_0.3.3               
# [34] circlize_0.4.15               GenomeInfoDbData_1.2.10       IRanges_2.34.1               
# [37] S4Vectors_0.38.1              ggrepel_0.9.3                 irlba_2.3.5.1                
# [40] listenv_0.9.0                 spatstat.utils_3.0-3          goftest_1.2-3                
# [43] annotate_1.78.0               spatstat.random_3.1-5         fitdistrplus_1.1-11          
# [46] parallelly_1.36.0             DelayedMatrixStats_1.22.1     leiden_0.4.3                 
# [49] codetools_0.2-19              DelayedArray_0.26.3           tidyselect_1.2.0             
# [52] shape_1.4.6                   farver_2.1.1                  viridis_0.6.3                
# [55] matrixStats_1.0.0             stats4_4.3.0                  BiocFileCache_2.8.0          
# [58] spatstat.explore_3.2-1        jsonlite_1.8.5                GetoptLong_1.0.5             
# [61] ellipsis_0.3.2                progressr_0.13.0              iterators_1.0.14             
# [64] ggridges_0.5.4                survival_3.5-5                foreach_1.5.2                
# [67] tools_4.3.0                   ica_1.0-3                     Rcpp_1.0.10                  
# [70] glue_1.6.2                    gridExtra_2.3                 MatrixGenerics_1.12.2        
# [73] GenomeInfoDb_1.36.1           withr_2.5.0                   BiocManager_1.30.21          
# [76] fastmap_1.1.1                 fansi_1.0.4                   digest_0.6.32                
# [79] R6_2.5.1                      mime_0.12                     colorspace_2.1-0             
# [82] Cairo_1.6-0                   scattermore_1.2               tensor_1.5                   
# [85] spatstat.data_3.0-1           RSQLite_2.3.1                 R.methodsS3_1.8.2            
# [88] celldex_1.10.1                utf8_1.2.3                    tidyr_1.3.1                  
# [91] generics_0.1.3                httr_1.4.6                    htmlwidgets_1.6.2            
# [94] S4Arrays_1.2.0                uwot_0.1.15                   pkgconfig_2.0.3              
# [97] gtable_0.3.3                  blob_1.2.4                    lmtest_0.9-40                
# [100] XVector_0.40.0                htmltools_0.5.5               clue_0.3-64                  
# [103] GSEABase_1.62.0               scales_1.3.0                  Biobase_2.60.0               
# [106] png_0.1-8                     rstudioapi_0.16.0             rjson_0.2.21                 
# [109] reshape2_1.4.4                nlme_3.1-162                  curl_5.0.1                   
# [112] cachem_1.0.8                  zoo_1.8-12                    GlobalOptions_0.1.2          
# [115] stringr_1.5.0                 BiocVersion_3.17.1            parallel_4.3.0               
# [118] miniUI_0.1.1.1                vipor_0.4.5                   AnnotationDbi_1.62.1         
# [121] ggrastr_1.0.1                 pillar_1.9.0                  vctrs_0.6.5                  
# [124] RANN_2.6.1                    promises_1.2.0.1              dbplyr_2.3.2                 
# [127] xtable_1.8-4                  cluster_2.1.4                 beeswarm_0.4.0               
# [130] magick_2.7.4                  cli_3.6.2                     compiler_4.3.0               
# [133] rlang_1.1.3                   crayon_1.5.2                  future.apply_1.11.0          
# [136] labeling_0.4.2                plyr_1.8.8                    ggbeeswarm_0.7.2             
# [139] stringi_1.7.12                viridisLite_0.4.2             deldir_1.0-9                 
# [142] munsell_0.5.0                 Biostrings_2.68.1             lazyeval_0.2.2               
# [145] spatstat.geom_3.2-1           Matrix_1.6-1                  ExperimentHub_2.8.0          
# [148] patchwork_1.2.0               sparseMatrixStats_1.12.1      bit64_4.0.5                  
# [151] future_1.33.2                 KEGGREST_1.40.0               shiny_1.7.4                  
# [154] SummarizedExperiment_1.30.2   interactiveDisplayBase_1.38.0 AnnotationHub_3.8.0          
# [157] ROCR_1.0-11                   igraph_1.5.0                  memoise_2.0.1                
# [160] bit_4.0.5            
