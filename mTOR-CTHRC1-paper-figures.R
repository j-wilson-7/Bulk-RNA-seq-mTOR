# Paper figures:
# Critical role for the TGFbeta1/mTORC1 signalling axis in defining the transcriptional identity of CTHRC1+ pathologic fibroblasts in IPF

#Must load TGF fold changes for each of the genes. Genes contains the list of genes (Entrez ID's)
setwd("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022")
fold.changes <- read.csv("DEseq2 interaction analysis results/diffexpr_results_crCTRL_TGF_response.csv")
mtorc1.dep <- read.csv("Visualisation/mTORC1-dep genes.csv", header=F)

fold.changes.filtered <- fold.changes %>%
  filter(SYMBOL %in% mtorc1.dep$V1)

IDfoldchange <- fold.changes.filtered[, c("ENTREZID", "log2FoldChange")]
geneList = IDfoldchange[, 2]
names(geneList) = as.character(IDfoldchange[,1])
geneList = sort(geneList, decreasing = TRUE)

#Pathway enrichment using ReactomePA
#qvalue cut off is set to 0.01
de <- as.character(IDfoldchange[,1])
x <- enrichPathway(gene=de, pvalueCutoff=0.05, qvalueCutoff=0.01, readable=TRUE,
                   minGSSize = 10, maxGSSize = 500)

#Make sure the list of fold changes has same number of genes as were used for pathway enrichment
xGENES <- x@gene


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


cnetplot(x, 
         showCategory = 5,
         #circular = TRUE, 
         colorEdge = TRUE,
         foldChange = FALSE,
         cex_label_category = 1,
         cex_label_gene = 0.8,
         cex_category = 1.5,
         shadowtext = "none"
         #color_category = "#AC145A"
         #color_gene = "#FBC740",
         #node_label = "category"
)

mytheme2 = list(
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        axis.title = element_text(face="bold"), legend.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
        axis.text=element_text(face="bold"), legend.position = "bottom", legend.box = "vertical"
  ))

heatplot(x,
         showCategory = x$Description[c(1, 3, 5, 6, 8, 9, 11, 13, 14)],
         label_format = function(x) stringr::str_wrap(x, width=31),
         foldChange = geneList
) + mytheme2 

heatplot(x,
         showCategory=15,
         foldChange = geneList)




#Barplot to show the two rapamycin-dependent pathways: Light green + dark green
barplot(x, 
        showCategory = 2,
        cex_label_category = 1,
        cex_label_gene = 0.8,
        cex_category = 1.5,
        label_format = function(x) stringr::str_wrap(x, width=40),
) + mytheme

cnetplot(x, 
         #label_format = function(x) stringr::str_wrap(x, width=10),
         #colorEdge = TRUE,
         cex_label_category = 1.2,
         cex_label_gene = 1,
         cex_category = 1.5,
         shadowtext = "none",
         color_category = "palegreen",
         color_gene = "#FBC740",
         node_label = "all"
)  + theme(legend.position = "none")



#Dotplot to show the rapa-insensitive, tgf-dep, mtorc1-dep genes:
dotplot(x, showCategory=15) + 
  ylab("") + 
  scale_color_gradient("Adjusted p-value", low = "seagreen1", high = "seagreen4") +
  #scale_color_gradient("Adjusted p-value", low = "#0097A9", high = "#004b54") +
  xlim(0,0.085) +
  scale_y_discrete(label = function(y) {
    y %>% sub("Homo sapiens\r: ", "", .)
  }) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  mytheme
