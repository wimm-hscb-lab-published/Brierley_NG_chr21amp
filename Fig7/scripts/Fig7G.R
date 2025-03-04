
library(data.table)
library(ggplot2)
library(dplyr)
data <- read.table("~/Fig7/BCL2_expr.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
data
# Boxplot
# Definition
levels <- c("HC", "chr21amp_MPNAML", "non_chr21amp_MPNAML")
x <- factor(data$Type, levels=levels)
y <- data$exp
z <- data$Type
maintitle <- "BCL2"
ytitle <- "log2(CPM + 1)"
xtitle <- "Type"

my_pal<-c("#4C9A92","#d8443c","#0D0E47")
# Plot
plot <- ggplot() +
  geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), color="black", outlier.size=1) +
  stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
  labs(title=maintitle, x=xtitle, y=ytitle) +
  scale_fill_manual(values=my_pal)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border=element_blank(),
        plot.title=element_text(hjust = 0.5, size=12),
        plot.subtitle=element_text(hjust = 0.5, size=12),
        axis.line.y.left = element_line(color="black"),
        axis.line.x = element_line(color="black"),
        axis.title=element_text(size=12),
        axis.text.x=element_text(size=12, colour="black",angle=45,hjust=1),
        axis.text.y=element_text(size=12, colour="black"),
        legend.position="none"
  )

plot


# Statistical test
# Pair-wise comparison
stats <- pairwise.wilcox.test(data$exp, data$Type, p.adjust.method="fdr")$p.value
. <- data.frame("V1"=row.names(stats))
stats <- cbind.data.frame(., stats)
names(stats) <- ""
# stats
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
