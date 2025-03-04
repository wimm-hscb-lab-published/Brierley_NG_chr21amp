
library(ggplot2)
library(cowplot)
library(remotes)
library(reticulate)
library(SingCellaR)
library(dplyr)
library(grid)
library(gridExtra)

###download TARGETseq data from https://zenodo.org/record/8060602 : MFp53MPNAML_integration.revised.rdata
##for pre-processing: incorporate metadata incl numbat annotation for processed postQC object from github : MPN_QC_metadata.tsv
#load both files and subset
md <- read.table("/github/Fig3/data/MPN_QC_metadata.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
load(file="~/MFp53MPNAML_integration.revised.rdata")


MPN_QC<-MPN.integration.all[,as.character(md$Cell)]
MPN_QC@meta.data<-md


#FIGURE 3A CALCULATIONS

#need to save the updated R data object
cells_chr21amp<-MPN_QC@meta.data %>% filter(chr21.donor.type =="chr21amp_MPNAML")
dim(cells_chr21amp)
#2205 cells total, remove those without TP53 genotyping information
cells_chr21amp<-cells_chr21amp %>% filter(chr21.genotype2!="chr21amp_MPNAML_TP53NA")

#1903 cells data available of total 2205

#examine the pre-leukaemic cells
chr21amp_preLSC<-cells_chr21amp %>% filter(chr21.genotype2=="chr21amp_preLSC")
table(chr21amp_preLSC$Genotype_curated)

# JAK2    JAK2_TET2 JAK2HOM_TET2            TET2_JAK2         WT 
# 51            2            4                  122          107 
# so 107 WT and 122+4+2+ 51 JAK2 single
 #179 jak2 single mutant, 107 wt/non jak2/TP53

#chr21.genotype4incorporates info on chr21amp and genotyping status
table(cells_chr21amp$chr21.genotype4) #162 tp53mut, nochr21amp and 1455 ch21amp,tp53mut


########violin plots
###run function


my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
#' Plot violin-plot for gene expression per custom group of cells
#' @param  object The SingCellaR object.
#' @param  custom_group_of_cells The list of custome group of cells.
#' @param  gene_list The vector of gene names.
#' @param  take_log2 is logical. If TRUE, log2 expression will be applied.
#' @param  xlab.text.size The font size of the label on the x-axis. Default 5
#' @param  point.size The size of data point. Default 0.2
#' @param  point.alpha The alpha parameter of the data point. Default 0.1
#' @param  grid.ncol the column number of the grid of the containing plots to be displayed. Default 3
#' @param  grid.nrow the row number of the grid of the containing plots to be displayed. Default 3
#' @export


plot_violin_for_genes_per_custom_group_of_cells<-function(object,custom_group_of_cells=list(),gene_list=c(),
                                                          take_log2=T,xlab.text.size=5,ylab.text.size=12,point.size=0.2,point.alpha=0.1,
                                                          grid.ncol = 3,grid.nrow=3){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  
  ####################################
  umi.dat<-get_normalized_umi(object)
  my.p<-list()
  for(i in 1:length(gene_list)){
    my.gene<-as.character(gene_list)[i]
    my.sub.dat<-umi.dat[rownames(umi.dat)==my.gene,]
    my.new.dat<-data.frame()
    for(k in 1:length(custom_group_of_cells)){
      x.cell.names<-custom_group_of_cells[[k]]
      x.expr<-my.sub.dat[names(my.sub.dat) %in% as.character(x.cell.names)]
      ##################
      exp.val<-as.numeric(x.expr)
      n.exp<-length(exp.val[exp.val>0])
      n.total<-length(exp.val)
      my.label<-paste(names(custom_group_of_cells)[k],"\n",n.exp,"/",n.total,sep="")
      ##################
      y.legend<-""
      ##################
      if(take_log2==T){
        my.f<-data.frame(cluster=my.label,normUMI=log2(as.numeric(x.expr)+1))
        y.legend<-"log2(normalized UMI)"
      }else{
        my.f<-data.frame(cluster=my.label,normUMI=as.numeric(x.expr))
        y.legend<-"normalized UMI"
      }
      my.new.dat<-rbind(my.new.dat,my.f)
      #    compare_means(normUMI ~ cluster,  data = my.new.dat)
    }
    my.p[[i]]<-ggplot(
      data = my.new.dat,aes(x=cluster, y=normUMI,fill=cluster))+labs(y=y.legend)+
      theme(text = element_text(size = 16))+
      geom_violin(trim=TRUE,scale="width")+
      #stat_summary(fun.data=data_summary)+
      stat_summary(fun=median,  geom="point", shape=23, size=3.5, color="white", fill="white")+ #can amend to median
      geom_boxplot(outlier.shape = NA, coef=0, width=0.1, color="white", alpha=0.1) +
      geom_jitter(height = 0,size=point.size,alpha = point.alpha)+
      #    stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", hide.ns = TRUE)+
      ggtitle(my.gene)+
      theme(plot.title = element_text(face = "bold", hjust = 0.5))+
      guides(fill="none")+
      theme(axis.text=element_text(size=xlab.text.size))+
      theme(axis.text=element_text(size=ylab.text.size))+
      theme(axis.title.y=element_text(size=12), axis.title.x=element_text(size=12))+
      # theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=1))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } 
  grid.arrange(grobs = my.p,nrow=grid.nrow,ncol=grid.ncol)
}

####now select genes and metadata


WT_normal.cells<-subset(MPN_QC@meta.data,MPN_QC@meta.data$chr21.genotype4=="WT")
WT_normal.cell_ids<-WT_normal.cells$Cell

MF.cells<-subset(MPN_QC@meta.data,MPN_QC@meta.data$chr21.genotype4=="MF")
MF.cells<-MF.cells$Cell

TP53_MT_no_chr21amp.cells<-subset(MPN_QC@meta.data,MPN_QC@meta.data$chr21.genotype4=="TP53_MT_no_chr21amp")
TP53_MT_no_chr21amp.cells<-TP53_MT_no_chr21amp.cells$Cell

chr21amp_MPNAML.cells<-subset(MPN_QC@meta.data,MPN_QC@meta.data$chr21.genotype4=="chr21amp_TP53_MT")
chr21amp_MPNAML.cells<-chr21amp_MPNAML.cells$Cell

preLSC.cells<-subset(MPN_QC@meta.data,MPN_QC@meta.data$chr21.genotype4=="preLSC")
preLSC.cells<-preLSC.cells$Cell

cell_groups<-list(WT_normal.cell_ids,MF.cells,TP53_MT_no_chr21amp.cells,chr21amp_MPNAML.cells,preLSC.cells)
names(cell_groups) <- c("HC","MF","TP53_MT_no_chr21amp.cells","chr21amp_BPMPN_cell","preLSC")

####


plot_violin_for_genes_per_custom_group_of_cells(
  MPN_QC,
  custom_group_of_cells = cell_groups,
  gene_list = c("DYRK1A"),
  take_log2 = T,
  xlab.text.size = 5,
  point.size = 0.6,
  point.alpha = 0.1,
  grid.ncol = 1,
  grid.nrow = 1
)

plot_violin_for_genes_per_custom_group_of_cells(
  MPN_QC,
  custom_group_of_cells = cell_groups,
  gene_list = c("DYRK1A","MORC3","TTC3","DSCR3","PIGP"),
  take_log2 = T,
  xlab.text.size = 5,
  point.size = 0.6,
  point.alpha = 0.1,
  grid.ncol = 3,
  grid.nrow = 2
)


#####

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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] gridExtra_2.3     dplyr_1.1.4       SingCellaR_1.2.0  reticulate_1.36.1 remotes_2.5.0    
# [6] cowplot_1.1.3     ggplot2_3.5.1    
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.4                  R.utils_2.12.3              tidyselect_1.2.1           
# [4] RSQLite_2.3.6               AnnotationDbi_1.50.3        htmlwidgets_1.6.4          
# [7] BiocParallel_1.22.0         Rtsne_0.17                  munsell_0.5.1              
# [10] codetools_0.2-20            ica_1.0-3                   statmod_1.5.0              
# [13] future_1.33.2               miniUI_0.1.1.1              withr_3.0.0                
# [16] colorspace_2.1-0            Biobase_2.48.0              uuid_1.2-0                 
# [19] rstudioapi_0.16.0           Seurat_4.0.1                stats4_4.0.4               
# [22] SingleCellExperiment_1.10.1 ROCR_1.0-11                 ggsignif_0.6.4             
# [25] tensor_1.5                  listenv_0.9.1               labeling_0.4.3             
# [28] GenomeInfoDbData_1.2.3      polyclip_1.10-6             farver_2.1.1               
# [31] bit64_4.0.5                 pheatmap_1.0.12             parallelly_1.37.1          
# [34] vctrs_0.6.5                 generics_0.1.3              xfun_0.43                  
# [37] R6_2.5.1                    GenomeInfoDb_1.24.2         clue_0.3-65                
# [40] LinkedMatrix_1.4.0          cachem_1.0.8                bitops_1.0-7               
# [43] spatstat.utils_3.0-4        fgsea_1.14.0                DelayedArray_0.14.1        
# [46] promises_1.3.0              scales_1.3.0                gtable_0.3.5               
# [49] globals_0.16.3              goftest_1.2-3               rlang_1.1.3                
# [52] GlobalOptions_0.1.2         splines_4.0.4               rstatix_0.7.0              
# [55] lazyeval_0.2.2              spatstat.geom_3.2-9         broom_1.0.5                
# [58] BiocManager_1.30.23         reshape2_1.4.4              abind_1.4-5                
# [61] backports_1.4.1             httpuv_1.6.15               tools_4.0.4                
# [64] spatstat.core_2.0-0         RColorBrewer_1.1-3          proxy_0.4-27               
# [67] BiocGenerics_0.34.0         sessioninfo_1.2.2           ggridges_0.5.6             
# [70] Rcpp_1.0.12                 plyr_1.8.9                  zlibbioc_1.34.0            
# [73] purrr_1.0.2                 RCurl_1.98-1.14             ggpubr_0.4.0               
# [76] rpart_4.1.23                deldir_2.0-4                pbapply_1.7-2              
# [79] GetoptLong_1.0.5            S4Vectors_0.26.1            zoo_1.8-12                 
# [82] SeuratObject_4.0.0          SummarizedExperiment_1.18.2 haven_2.5.4                
# [85] ggrepel_0.9.5               cluster_2.1.6               tinytex_0.51               
# [88] magrittr_2.0.3              data.table_1.15.4           scattermore_1.2            
# [91] openxlsx_4.2.5.2            circlize_0.4.16             lmtest_0.9-40              
# [94] RANN_2.6.1                  fitdistrplus_1.1-11         matrixStats_1.3.0          
# [97] hms_1.1.3                   patchwork_1.2.0             mime_0.12                  
# [100] xtable_1.8-4                XML_3.99-0.16.1             rio_0.5.26                 
# [103] AUCell_1.10.0               readxl_1.4.3                IRanges_2.22.2             
# [106] shape_1.4.6.1               cccd_1.6                    compiler_4.0.4             
# [109] tibble_3.2.1                KernSmooth_2.23-22          crayon_1.5.2               
# [112] R.oo_1.26.0                 htmltools_0.5.8.1           mgcv_1.9-1                 
# [115] later_1.3.2                 tidyr_1.3.1                 RcppParallel_5.1.7         
# [118] DBI_1.2.2                   ComplexHeatmap_2.4.3        MASS_7.3-53.1              
# [121] Matrix_1.5-3                car_3.0-10                  cli_3.6.2                  
# [124] R.methodsS3_1.8.2           parallel_4.0.4              igraph_2.0.3               
# [127] GenomicRanges_1.40.0        forcats_1.0.0               pkgconfig_2.0.3            
# [130] bigmemory.sri_0.1.8         foreign_0.8-86              plotly_4.10.4              
# [133] spatstat.sparse_3.0-3       annotate_1.66.0             crochet_2.3.0              
# [136] XVector_0.28.0              stringr_1.5.1               digest_0.6.35              
# [139] sctransform_0.4.1           RcppAnnoy_0.0.22            graph_1.66.0               
# [142] spatstat.data_3.0-4         cellranger_1.1.0            leiden_0.4.3.1             
# [145] fastmatch_1.1-4             uwot_0.1.10                 GSEABase_1.50.1            
# [148] curl_5.2.1                  shiny_1.8.1.1               rjson_0.2.21               
# [151] lifecycle_1.0.4             nlme_3.1-164                jsonlite_1.8.8             
# [154] carData_3.0-5               limma_3.44.3                bigmemory_4.6.4            
# [157] viridisLite_0.4.2           fansi_1.0.6                 pillar_1.9.0               
# [160] lattice_0.22-6              fastmap_1.1.1               httr_1.4.7                 
# [163] survival_3.6-4              glue_1.7.0                  zip_2.3.1                  
# [166] FNN_1.1.4                   png_0.1-8                   bit_4.0.5                  
# [169] stringi_1.8.4               blob_1.2.4                  memoise_2.0.1              
# [172] irlba_2.3.5.1               future.apply_1.11.2    

