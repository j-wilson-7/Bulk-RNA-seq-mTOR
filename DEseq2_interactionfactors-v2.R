####### DEseq2 with interaction factors #######
#CRISPR mTORC1 and 2 DATASET

#### Prepare data

# Set working directory
library(DESeq2)
setwd("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022")

# Import the data
countdata <- read.delim("./Combined_counts.tabular", header=TRUE, comment.char="#", row.names=1)

#all samples (previously did each set of samples separately)
countdata1 <- countdata

# Convert to matrix
countdata1 <- as.matrix(countdata1)
head(countdata1)

# Assign condition (first three are controls, second three and third three contain two different experimental conditions)
# all samples:
(genotype1 <- factor(c(rep("crCTRL", 6), rep("crRAPTOR", 6), rep("crRICTOR", 6))))
(condition1 <- factor(rep(rep(c("MC","TGF"),each=3),3)))




#### More complex DESeq2 design with interaction factor:

# Create a coldata frame and instantiate the DESeqDataSet and set complex design
(coldata1 <- data.frame(row.names=colnames(countdata1), condition1, genotype1))
dds1 <- DESeqDataSetFromMatrix(countData=countdata1, colData=coldata1, design=~ genotype1 + condition1 + genotype1:condition1)
dds1

# Run the DESeq pipeline
dds1 <- DESeq(dds1)
resultsNames(dds1)


## Choose one of these possibilities to analyse:

## All samples
## Effect of TGF treatment in crCTRL
res1 <- results(dds1, contrast=c("condition1","TGF","MC"))
# Effect of TGF treatment in crRAPTOR
res1 <- results(dds1, list( c("condition1_TGF_vs_MC","genotype1crRAPTOR.condition1TGF") ))
# Difference between crCTRL and crRAPTOR without treatment
res1 <- results(dds1, contrast=c("genotype1","crCTRL","crRAPTOR"))
## Difference between crCTRL and crRAPTOR with TGF treatment
res1 <- results(dds1, list( c("genotype1_crRAPTOR_vs_crCTRL","genotype1crRAPTOR.condition1TGF")))
# Difference in TGF response between crCTRL samples vs crRAPTOR samples (interaction term)
res1 <- results(dds1, name="genotype1crRAPTOR.condition1TGF")
# Effect of TGF treatment in crRICTOR
res1 <- results(dds1, list( c("condition1_TGF_vs_MC","genotype1crRICTOR.condition1TGF") ))
# Difference between crCTRL and crRICTOR without treatment
res1 <- results(dds1, contrast=c("genotype1","crCTRL","crRICTOR"))
## Difference between crCTRL and crRICTOR with TGF treatment
res1 <- results(dds1, list( c("genotype1_crRICTOR_vs_crCTRL","genotype1crRICTOR.condition1TGF")))
## Difference in TGF response between crCTRL samples vs crRICTOR samples (interaction term)
res1 <- results(dds1, name="genotype1crRICTOR.condition1TGF")



## This section remains the same for all of the analyses:

# Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata1)[1] <- "Gene"
head(resdata1)

## Annotation with gene symbols and filtering

#install.packages("dplyr")
library(dplyr)
library(BiocManager)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

#Create new data frame with annotation information
anno <- AnnotationDbi::select(org.Hs.eg.db,keys=resdata1$Gene,
                              columns=c("SYMBOL"),
                              keytype="ENTREZID")

#Rename first column in resdata file from Gene to Ensembl
resdata1.labelled <- resdata1
colnames(resdata1.labelled)
names(resdata1.labelled)[names(resdata1.labelled) == "Gene"] <- "ENTREZID"

#Bind the annotation information to the results data frame
results1.annotated <- left_join(resdata1.labelled, anno,by="ENTREZID")

#Change the order of the columns (move last column to start)
results1.annotated <- results1.annotated[,c(26,1:25)]

# Remove NAs
results1.annotated <- na.omit(results1.annotated)

#Filter for abs(Log2FoldChange)>=0. and padj<0.05
results1.annotated.p0.05 <- filter(results1.annotated, padj < 0.05)
results1.annotated.p0.05.fc0.58 <- filter(results1.annotated.p0.05, abs(log2FoldChange) >= 0.58)


#If you want to make volcano plot below save the unfiltered TGF response after removing NAs
TGF.unfiltered <- results1.annotated


## Write results to one of these files:
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_crCTRL_TGF_response.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_crRAPTOR_TGF_response.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_crCTRLvscrRAPTOR_baseline.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_crCTRL_TvscrRAPTOR_T.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_interaction_crCTRLvscrRAPTOR_TGFresponse.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_crRICTOR_TGF_response.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_crCTRLvscrRICTOR_baseline.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_crCTRL_TvscrRICTOR_T.csv")
write.csv(results1.annotated.p0.05.fc0.58, file="diffexpr_results_interaction_crCTRLvscrRICTOR_TGFresponse.csv")


#To save normalised counts with symbols
norm.counts <- as.data.frame(counts(dds1, normalized=TRUE))
entrez_ids <- rownames(norm.counts)
gene_symbols <- select(org.Hs.eg.db, keys = entrez_ids, columns = "SYMBOL", keytype = "ENTREZID")
norm.counts$SYMBOL <- gene_symbols$SYMBOL
norm.counts <- norm.counts[,c(19,1:18)]

#write.csv(norm.counts, file="normalised.counts.24h")
write.csv(norm.counts, file="normalised.counts.6h")



#### Visualisation:

## PCA plot
library(ggplot2)

#Custom theme for PCA plot for publication
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))

#using rlog transformed data
rld <- rlog(dds1)
pcaData <- plotPCA(rld, intgroup = c("genotype1","condition1"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

plotPCA(rld, intgroup = c("genotype1","condition1"))

#choose colours for the plot:
col = c("darkgoldenrod1", "darkorchid", "seagreen3")
col = c("#F6BE00", "#AC145A", "#0097A9")

#Plot PCA
ggplot(pcaData, aes(PC1, PC2, shape=condition1, color=genotype1)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +

  #customisable shapes and color labels
  scale_color_manual(labels = c("Control", "crRPTOR", "crRICTOR"),
                     values = col) +
  labs(color="Treatment") +
  scale_shape_manual(labels = c("Media Control", "TGF-β1"),
                     values = c(16, 17)) +
  labs(shape="Condition") +
  geom_point(size=2) +
  mytheme +
  guides(colour = guide_legend(order = 1), 
                shape = guide_legend(order = 2))
  

## Heatmap

#Extract the top 150 most significantly TGF-B1 regulated genes
setwd("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022")

TGF.results <- read.csv("Deseq2 interaction analysis results/diffexpr_results_crCTRL_TGF_response.csv", header=TRUE, comment.char="#")
TGF.results.top150 <- head(TGF.results[order(TGF.results$padj),], 7122) #Change this to the number of genes you would like on the heatmap
TGF.results.top150 <- TGF.results.top150[, c(2,10:27)]
rownames(TGF.results.top150) <- TGF.results.top150[,1]
TGF.results.top150 <- TGF.results.top150[,c(2:19)]

#Scale the data to Z scores
heat <- t(scale(t(TGF.results.top150)))


#Prepare sample labels
col.labels.treatment <- data.frame(
  sample = c(rep("Control", 6), rep("crRAPTOR", 6), rep("crRICTOR", 6)))
col.labels.TGF <- data.frame(
  sample = c(rep(rep(c("Media Control","TGF-β1"),each=3),3)))
row.names(col.labels.TGF) <- colnames(heat)
col.labels <- cbind(col.labels.TGF, col.labels.treatment)
colnames(col.labels) <- c("  ", " ")


#Plot the pretty heatmap
pheatmap(heat,
         show_rownames=F,
         cutree_rows=2,
         treeheight_row=0,
         #legend_breaks= c(-2, 0, 2, 4),
         annotation_col = col.labels,
         show_colnames=F
         #color=hcl.colors(50,"RdBu"))
)


## Matrisome pie chart

slices <- c(74,26)
lbls <- c("TGF-β1-modulated matrisomal genes (393)", "mTORC1-regulated (136)")
pie(slices,labels = lbls, col=c("white","#AC145A"))


## Volcano plot highlighting mTORC1 and mTORC2 genes

#Load packages into R session
library(EnhancedVolcano)
pacman::p_load_gh("trinker/textshape")

setwd("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022")

#Take the unfiltered TGF results from above (TGF.unfiltered) and also take the filtered results as TGF.dep
TGF.dep <- read.csv("Deseq2 interaction analysis results/diffexpr_results_crCTRL_TGF_response.csv", header=TRUE, comment.char="#")

#Make sure all of the NA's are removed:
TGF.dep <- na.omit(TGF.dep)
TGF.unfiltered <- na.omit(TGF.unfiltered)

#Make the gene symbols into row columns
TGF.unfiltered <- column_to_rownames(TGF.unfiltered, loc=1)

#Set publication theme
mytheme2 = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
          axis.text=element_text(face="bold")
    ))

#Get the mTORC1 and 2-dep genes from the interaction term with those not in TGF control removed:
mTORC1.dep <- read.csv("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022/Visualisation/mTORC1-dep genes.csv", header=FALSE, comment.char="#")
mTORC2.dep <- read.csv("~/Documents/Chambers/mTORC1 RAPTOR RNAseq Sep-Dec2019/AnalysisYaleSummer2022/Visualisation/mTORC2-dep genes.csv", header=FALSE, comment.char="#")


#Make two separate volcano plots, one for mTORC1 and one for mTORC2:

#Set the colour scheme to label mTORC1-dep genes 
keyvals <- ifelse(
  rownames(TGF.unfiltered) %in% mTORC1.dep$V1, '#A7204C',
  ifelse(rownames(TGF.unfiltered) %in% TGF.dep$SYMBOL, '#edd2db',
         'darkgrey'))

names(keyvals)[keyvals == 'darkgrey'] <- 'Not significant'
names(keyvals)[keyvals == '#edd2db'] <- 'mTORC1-independent'
names(keyvals)[keyvals == '#A7204C'] <- 'mTORC1-dependent'

#Set the colour scheme to label mTORC2-dep genes 
keyvals <- ifelse(
  rownames(TGF.unfiltered) %in% mTORC2.dep$V1, 'yellow',
  ifelse(rownames(TGF.unfiltered) %in% TGF.dep$SYMBOL, 'orange',
         'darkgrey'))

names(keyvals)[keyvals == 'darkgrey'] <- 'Not significant'
names(keyvals)[keyvals == 'orange'] <- 'mTORC2-independent'
names(keyvals)[keyvals == 'yellow'] <- 'mTORC2-dependent'

#Plot volcano plot logFC2 and p-adjvalue0.05
p1 <- EnhancedVolcano(
  #Load data
  TGF.unfiltered,
  x = 'log2FoldChange',
  y = 'padj',
  
  #Labelling
  lab = rownames(TGF.unfiltered),
  selectLab = '',
  cutoffLineType = 'dashed',
  
  title = '', 
  subtitle = '',
  #legendlabels = c('test', 'mTOR-dependent'),
  #legendPosition = 'right',
  caption = '',
  #caption = bquote(~Log[2]~ "Fold change cutoff, 1.5; adjusted p-value cutoff, 0.05"),
  legendPosition = 'right',
  
  #Cutoffs
  pCutoff = 5*10e-03,
  FCcutoff = 0.58,
  
  #Size colour shape of points
  #col = color,
  colAlpha=1,
  colCustom = keyvals,
  pointSize = 1.5
)


#Custom axis tick marks
p1 +
  ggplot2::coord_cartesian(xlim=c(-10, 10)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-10, 10, 2)) +
  mytheme2 +
  theme(legend.title = element_blank())
#+ theme(plot.subtitle = element_text(color='black', hjust=0.5, face='italic'), plot.title = element_text(color='black', hjust=0.5))








