#### Plotting pathways analysis results using ClusterProfiler/ReactomePA
# Bulk RNA-seq datasets:
# 1. CRISPR/Cas9 mTORC1 or mTORC2-depleted pHLFs stimulated with TGFbeta
# 2. mTOR inhibitors: AZD8055 (dual mTOR inhibitor) or rapamycin (partial mTORC1 inhibitor) stimulated with TGFbeta
# Input is significant gene lists (q<0.05, FC>1.5) following DESeq2 interaction analysis

#Load packages
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)

#Load the gene list of interest:
setwd("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022/Visualisation")
setwd("~/Documents/Chambers/AZD analysis YaleSummer2022/Visualisation")
genes <- read.delim("TGF-B1-RAPTOR-dep-genes-entrez.txt", header=FALSE)
genes <- read.delim("TGF-B1-RICTOR-dep-genes-entrez.txt", header=FALSE)
genes <- read.delim("rapa-interaction-entrez.txt")
genes <- read.delim("azd-interaction-entrez.txt")
genes <- read.delim("tgf-rapains-mtorc1-dep.txt", header=FALSE)

de <- as.character(genes[,1])
head(de)

#Pathway enrichment using ReactomePA
#qvalue cut off was set to 0.01, so that there are no significant mTORC2 pathways. Now set to 0.05
x <- enrichPathway(gene=de, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE, minGSSize = 10, maxGSSize = 500)
head(x)

x <- enrichKEGG(gene=de,
  organism = 'hsa',
  qvalueCutoff = 0.05,
  minGSSize = 10
)

#Save analysis file
write.csv(x, file="diffexpr_results_TGF-B1-RAPTOR-dep-pathways.csv")
write.csv(x, file="diffexpr_results_TGF-B1-RICTOR-dep-pathways.csv")
write.csv(x, file="diffexpr_results_TGF-B1-AZD-dep-pathways-Reactome.csv")
write.csv(x, file="diffexpr_results_TGF-B1-Rapa-dep-pathways-Reactome.csv")
write.csv(x, file="diffexpr_results_TGF-B1-Rapa-ins-RAPTOR-dep-pathways-KEGG.csv")

#Visualise using bar plot
barplot(x, showCategory=15,
        label_format = function(x) stringr::str_wrap(x, width=40)
        ) 

#Custom theme for dot plot for publication
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))

#Visualise using dot plot - AZD or Rapamycin-treated pHLFs:
dotplot(x, showCategory=15) +
  ylab("") +
  #scale_color_gradient("Adjusted p-value", low = "maroon1", high = "maroon4") +
  scale_color_gradient("Adjusted p-value", low = "#0C2D48", high = "#2E8BC0") +
  xlim(0,0.1) +
  scale_y_discrete(label = function(y) {
    y %>% sub("Homo sapiens\r: ", "", .)
  }) +
  mytheme


#Visualise using dot plot - mTORC1 or mTORC2 CRISPR/cas9 pHLFs:
dotplot(x, showCategory=15) + 
  ylab("") + 
  scale_color_gradient("Adjusted p-value", low = "#AC145A", high = "#560a2d") +
  #scale_color_gradient("Adjusted p-value", low = "#0097A9", high = "#004b54") +
  xlim(0,0.085) +
  scale_y_discrete(label = function(y) {
  y %>% sub("Homo sapiens\r: ", "", .)
  }) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  mytheme


#Visualise using cnet plot:
#Remove the Homo sapiens label from term labels
x@result$Description <- gsub("Homo sapiens\r: ","", as.character(x@result$Description))

#Look at the pathways and choose one to plot
selected_pathways <- sample(x@result$Description, 15)
selected_pathways

#For mTORC1-depleted pHLFs:
cnetplot(x, 
         showCategory = x$Description[c(12, 9)],
         #circular = TRUE, 
         colorEdge = TRUE,
         foldChange = FALSE,
         cex_label_category = 0.6,
         #cex_label_gene = 0.4,
         cex_category = 1.5,
         shadowtext = "none",
         #color_category = "#AC145A"
         #color_gene = "#FBC740",
         node_label = "category"
) 

#for AZD8055-treated pHLFs:
cnetplot(x, 
         showCategory = x$Description[c(1,10,9,17,7)], 
         #circular = TRUE, 
         colorEdge = TRUE,
         foldChange = FALSE,
         cex_label_category = 0.6,
         #cex_label_gene = 0.4,
         cex_category = 1.5,
         shadowtext = "none",
         #color_category = "#AC145A"
         #color_gene = "#FBC740",
         node_label = "category"
) 



