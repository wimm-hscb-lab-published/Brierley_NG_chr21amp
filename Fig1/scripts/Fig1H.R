###read in data from extended data table 1

library(readxl)
library(survival)
library(cmprsk)

PH_clin_mol_v2 <- read_excel("/~/Table S1.xlsx", skip =1)

View(PH_clin_mol_v2)
PH_clin_mol_v2$ttlfupdm<-PH_clin_mol_v2$`Time to last follow-up/death in months`
PH_clin_mol_v2$ttlfupdm<-as.numeric(PH_clin_mol_v2$ttlfupdm)


fit<-survfit(Surv(PH_clin_mol_v2$ttlfupdm,PH_clin_mol_v2$DeathIndicator)~ PH_clin_mol_v2$chrom_21_amp)

survdiff(Surv(PH_clin_mol_v2$ttlfupdm,PH_clin_mol_v2$DeathIndicator)~ PH_clin_mol_v2$chrom_21_amp)
print(fit)
plot(fit)
summary(fit,times=c(3,6,12))
plot(fit,col=c("blue","red"),xlim=c(0,36),xaxp=c(0,36,6),lwd=2,main="Overall Survival by chr21 amp (n=64)",xlab="Time since diagnosis (months)")
legend("topright", c("No chr21amp","chr21amp"),col=c("blue","red"),lwd=2,lty=1,bty="n")

summary(PH_clin_mol_v2$ttlfupdm)
summary(PH_clin_mol_v2$Age)
table(PH_clin_mol_v2$chromothripsis)

PH_clin_mol_v2 <- PH_clin_mol_v2[!is.na(PH_clin_mol_v2$chromothripsis), ] 
table(PH_clin_mol_v2$TP53)
test<-fisher.test(PH_clin_mol_v2$chromothripsis,PH_clin_mol_v2$TP53)
test$p.value
test

chisq.test(PH_clin_mol_v2$chromothripsis,PH_clin_mol_v2$TP53)$expected
####
table(PH_clin_mol_v2$chromothripsis,PH_clin_mol_v2$TP53)

table(PH_clin_mol_v2$chromothripsis,PH_clin_mol_v2$chrom_21)
# 
# 
# sessionInfo()
# version 4.0.4 (2021-02-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cmprsk_2.2-10      readxl_1.3.1       readr_1.4.0        kableExtra_1.3.4  
# [5] ggrepel_0.9.1      cowplot_1.1.1      IRdisplay_1.0      survminer_0.4.9   
# [9] ggpubr_0.4.0       survival_3.2-10    patchwork_1.1.1    pec_2020.11.17    
# [13] prodlim_2019.11.13 ggforce_0.3.3      glue_1.4.2         ggbeeswarm_0.6.0  
# [17] ggplot2_3.3.3      data.table_1.14.0  dplyr_1.0.5        stringr_1.4.0     
# 
# loaded via a namespace (and not attached):
#   [1] webshot_0.5.2       httr_1.4.2          repr_1.1.3          numDeriv_2016.8-1.1
# [5] tools_4.0.4         backports_1.2.1     utf8_1.2.1          R6_2.5.0           
# [9] vipor_0.4.5         DBI_1.1.1           colorspace_2.0-0    withr_2.4.1        
# [13] tidyselect_1.1.0    gridExtra_2.3       curl_4.3            compiler_4.0.4     
# [17] rvest_1.0.0         cli_2.4.0           xml2_1.3.2          labeling_0.4.2     
# [21] scales_1.1.1        survMisc_0.5.5      systemfonts_1.0.1   digest_0.6.27      
# [25] foreign_0.8-81      svglite_2.0.0       rmarkdown_2.7       rio_0.5.26         
# [29] base64enc_0.1-3     pkgconfig_2.0.3     htmltools_0.5.1.1   rlang_1.0.2        
# [33] rstudioapi_0.13     farver_2.1.0        generics_0.1.0      zoo_1.8-9          
# [37] jsonlite_1.7.2      zip_2.1.1           car_3.0-10          magrittr_2.0.1     
# [41] Matrix_1.5-3        Rcpp_1.0.6          munsell_0.5.0       fansi_0.4.2        
# [45] abind_1.4-5         lifecycle_1.0.0     stringi_1.5.3       carData_3.0-4      
# [49] MASS_7.3-53.1       plyr_1.8.6          grid_4.0.4          forcats_0.5.1      
# [53] crayon_1.4.1        lattice_0.20-41     haven_2.3.1         splines_4.0.4      
# [57] hms_1.0.0           knitr_1.31          pillar_1.5.1        timereg_1.9.8      
# [61] ggsignif_0.6.1      reshape2_1.4.4      codetools_0.2-18    evaluate_0.14      
# [65] remotes_2.3.0       BiocManager_1.30.12 vctrs_0.3.7         tweenr_1.0.2       
# [69] foreach_1.5.1       cellranger_1.1.0    gtable_0.3.0        purrr_0.3.4        
# [73] polyclip_1.10-0     tidyr_1.1.3         km.ci_0.5-2         assertthat_0.2.1   
# [77] xfun_0.22           openxlsx_4.2.3      xtable_1.8-4        broom_0.7.6        
# [81] rstatix_0.7.0       viridisLite_0.3.0   tibble_3.1.0        iterators_1.0.13   
# [85] tinytex_0.31        beeswarm_0.3.1      KMsurv_0.1-5        lava_1.6.9         
# [89] ellipsis_0.3.1
