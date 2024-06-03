###### DEseq2 with interaction factors: Imapct of mTOR-inhibitors AZD8055 and rapamycin on TGFbeta1 effect #######

## Prepare data

# Set working directory and load packages
library(DESeq2)
#install.packages("dplyr")
library(dplyr)
library(BiocManager)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

setwd("~/Documents/Chambers/AZD analysis YaleSummer2022")

# Import the data
countdata <- read.delim("Combined_counts.tabular.txt", header=TRUE, comment.char="#", row.names=1)

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first four are controls, second four and third four contain two different experimental conditions)
(genotype <- factor(c(rep("CTRL", 8), rep("RAPA", 8), rep("AZD", 8))))
(condition <- factor(rep(rep(c("MC","TGF"),each=4),3)))

# Create a coldata frame and instantiate the DESeqDataSet and set complex design
(coldata <- data.frame(row.names=colnames(countdata), condition, genotype))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~ genotype + condition + genotype:condition)
dds

# Relevel so that CTRL is the reference, rather than AZD (DEseq2 wants to go in alphabetical order)
dds$genotype <- relevel(dds$genotype, "CTRL")

# Run the DESeq pipeline
dds <- DESeq(dds)
resultsNames(dds)


## Choose one of these possibilities to analyse:

## Effect of TGF treatment in CTRL
res <- results(dds, contrast=c("condition","TGF","MC"))

#RAPAMCIN
# Effect of TGF treatment in RAPA
res <- results(dds, list( c("condition_TGF_vs_MC","genotypeRAPA.conditionTGF") ))
# Difference between CTRL and RAPA without treatment
res <- results(dds, contrast=c("genotype","CTRL","RAPA"))
## Difference between CTRL and RAPA with TGF treatment
res <- results(dds, list( c("genotype_RAPA_vs_CTRL","genotypeRAPA.conditionTGF")))
# Difference in TGF response between CTRL samples vs RAPA samples (interaction term)
res <- results(dds, name="genotypeRAPA.conditionTGF")

#AZD8055
# Effect of TGF treatment in AZD
res <- results(dds, list( c("condition_TGF_vs_MC","genotypeAZD.conditionTGF") ))
# Difference between CTRL and AZD without treatment
res <- results(dds, contrast=c("genotype","CTRL","AZD"))
## Difference between CTRL and AZD with TGF treatment
res <- results(dds, list( c("genotype_AZD_vs_CTRL","genotypeAZD.conditionTGF")))
## Difference in TGF response between CTRL samples vs AZD samples (interaction term)
res <- results(dds, name="genotypeAZD.conditionTGF")




#Same analysis here for all scenarios:

# Merge with normalized count data
resdata<- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

## Annotation with gene symbols and filtering
#Create new data frame with annotation information
anno <- AnnotationDbi::select(org.Hs.eg.db,keys=resdata$Gene,
                              columns=c("SYMBOL"),
                              keytype="ENTREZID")

#Rename first column in resdata file from Gene to Ensembl
resdata.labelled <- resdata
colnames(resdata.labelled)
names(resdata.labelled)[names(resdata.labelled) == "Gene"] <- "ENTREZID"

#Bind the annotation information to the results data frame
results.annotated <- left_join(resdata.labelled, anno,by="ENTREZID")

#Change the order of the columns (move last column to start)
results.annotated <- results.annotated[,c(32,1:31)]

# Remove NAs
results.annotated <- na.omit(results.annotated)

#Filter for abs(Log2FoldChange)>=0. and padj<0.05
results.annotated.p0.05 <- filter(results.annotated, padj < 0.05)
results.annotated.p0.05.fc0.58 <- filter(results.annotated.p0.05, abs(log2FoldChange) >= 0.58)





## Write results to one of these files:

#CTRL
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_CTRL_TGF_response.csv")

#RAPA
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_RAPA_TGF_response.csv")
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_CTRLvsRAPA_baseline.csv")
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_CTRL_TvsRAPA_T.csv")
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_interaction_CTRLvsRAPA_TGFresponse.csv")

#AZD
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_AZD_TGF_response.csv")
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_CTRLvsAZD_baseline.csv")
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_CTRL_TvsAZD_T.csv")
write.csv(results.annotated.p0.05.fc0.58, file="diffexpr_results_interaction_CTRLvsAZD_TGFresponse.csv")


#If making heatmap with TGF response in control samples do this additional step: (Make sure NAs are removed)
TGFctrl <- results.annotated
TGF.dep <- results.annotated.p0.05.fc0.58$SYMBOL



#### Visualisation:

## Load libraries required
library(ggplot2)
library("gplots")
library(pheatmap)

#Custom theme for PCA plot for publication
mytheme = list(
        theme_classic()+
                theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                      legend.title = element_blank(),legend.position="bottom", axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
                      axis.text=element_text(face="bold")
))

#using rlog transformed data
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup = c("genotype","condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#Plot PCA
ggplot(pcaData, aes(PC1, PC2, shape=condition, color=genotype)) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        
        #customisable shapes and color labels
        scale_color_manual(labels = c("Control", "AZD8055", "Rapamycin"),
                           values = c("darkgoldenrod1", "steelblue1", "maroon3")) +
        labs(color="Treatment") +
        scale_shape_manual(labels = c("Media Control", "TGF-β1"),
                           values = c(16, 17)) +
        labs(shape="Condition") +
        geom_point(size=2) +
        mytheme +
        guides(colour = guide_legend(order = 1), 
               shape = guide_legend(order = 2))


#Heatmap
library(pheatmap)

#Extract the top 150 most significantly TGF-B1 regulated genes
setwd("~/Documents/Chambers/AZD analysis YaleSummer2022")
TGF.results <- read.csv("Diffexpr results/diffexpr_results_CTRL_TGF_response.csv", header=TRUE, comment.char="#")
TGF.results.top150 <- head(TGF.results[order(TGF.results$padj),], 150) #Change this to the number of genes you would like on the heatmap
TGF.results.top150 <- TGF.results.top150[, c(2,10:33)]
rownames(TGF.results.top150) <- TGF.results.top150[,1]
TGF.results.top150 <- TGF.results.top150[,c(2:25)]

#Scale the data to Z scores
heat <- t(scale(t(TGF.results.top150)))

#Prepare sample labels
col.labels.treatment <- data.frame(
                        sample = c(rep("Control", 8), rep("Rapamycin", 8), rep("AZD8055", 8)))
col.labels.TGF <- data.frame(
                        sample = c(rep(rep(c("Media Control","TGF-β1"),each=4),3)))
row.names(col.labels.TGF) <- colnames(heat)
col.labels <- cbind(col.labels.TGF, col.labels.treatment)
colnames(col.labels) <- c("  ", " ")

#Plot the pretty heatmap
pheatmap(heat,
         show_rownames=F,
         cutree_rows=2,
         treeheight_row=0,
         legend_breaks= c(-2, 0, 2, 4),
         annotation_col = col.labels,
         show_colnames=F
         #color=hcl.colors(50,"RdBu"))
         )
         

## Now lets plot AZD and Control only
TGF.results.top150 <- TGF.results.top150[,-c(9:20)]

#Scale the data to Z scores
heat <- t(scale(t(TGF.results.top150)))

#Prepare sample labels
col.labels <- col.labels[-c(9:20),]

#Plot the pretty heatmap
pheatmap(heat,
         show_rownames=T,
         cutree_rows=3,
         treeheight_row=0,
         annotation_col = col.labels,
         show_colnames=F,
         cluster_cols=F,
         fontsize_row = 3
         #color=hcl.colors(50,"RdBu"))
)


## Now lets plot AZD and Control only and provide it with matrisome genes
setwd("~/Documents/Chambers/AZD analysis YaleSummer2022")
matrisome <- read.csv("Visualisation/Collagen-for-heatmap.csv", header=FALSE, comment.char="#")
matrisome <- read.csv("Visualisation/Proteoglycans-for-heatmap.csv", header=FALSE, comment.char="#")
matrisome <- read.csv("Visualisation/Glycoproteins-for-heatmap.csv", header=FALSE, comment.char="#")

matrisome <- c(matrisome$V1)

matrisome.counts <- subset(TGF.results, SYMBOL %in% matrisome)
matrisome.counts <- matrisome.counts[, c(2, 10:17, 30:33)]
rownames(matrisome.counts) <- matrisome.counts[,1]
matrisome.counts <- matrisome.counts[,c(2:13)]

#Scale the data to Z scores
heat <- t(scale(t(matrisome.counts)))

#Prepare sample labels
#col.labels <- col.labels[-c(9:20),] #Don't need this if ran the previous section of code

#Set the order of the pheatmap, don't order by clustering
col_order = c("CTRL.MC.1", "CTRL.MC.2", "CTRL.MC.3", "CTRL.MC.4", "CTRL..T.1", "CTRL..T.2", "CTRL..T.3", "CTRL..T.4", "AZD..T.1", "AZD..T.2", "AZD..T.3", "AZD..T.4")

pheatmap(heat,
         show_rownames=T,
         cutree_rows=5,   #6 for collagen, 3 for proteo, 5 for glyco
         treeheight_row=0,
         #legend_breaks= c(-2, 0, 2, 4),
         annotation_col = col.labels,
         show_colnames=F,
         fontsize_row=5,
         col_order = col_order,
         cluster_cols = FALSE
         #color=hcl.colors(50,"RdBu"))
)


##Volcano plot 

#Download the packages required
#BiocManager::install('EnhancedVolcano')
#if (!require("pacman")) install.packages("pacman")

#Load packages into R session
library(EnhancedVolcano)
pacman::p_load_gh("trinker/textshape")

#Make sure you do the step above following deseq2 where you take all TGF regulated genes without p value or FC cutoffs

#Make the gene symbols into row columns
TGFctrl <- column_to_rownames(TGFctrl, loc=1)

#Set publication theme
mytheme2 = list(
        theme_classic()+
                theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                      axis.title = element_text(face="bold"), strip.text = element_text(face="bold", size=9),
                      axis.text=element_text(face="bold")
                ))

#Get the AZD-dep genes from the interaction term with those not in TGF control removed:
AZD.dep <- read.table("~/Documents/Chambers/AZD analysis YaleSummer2022/Visualisation/AZD-dep genes", header=FALSE, comment.char="#")
Rapa.dep <- read.table("~/Documents/Chambers/AZD analysis YaleSummer2022/Visualisation/Rapa-dep genes.txt", header=FALSE, comment.char="#")

#Set the colour scheme to label AZD-dep genes ONLY:
keyvals <- ifelse(
        rownames(TGFctrl) %in% AZD.dep$V1, 'royalblue',
        'lightgrey')

names(keyvals)[keyvals == 'lightgrey'] <- 'AZD8055-insensitive'
names(keyvals)[keyvals == 'royalblue'] <- 'AZD8055-sensitive'


#Set the colour scheme to label AZD-dep genes AND label AZD-indep genes too:
keyvals <- ifelse(
        rownames(TGFctrl) %in% AZD.dep$V1, '#145DA0',
        ifelse(rownames(TGFctrl) %in% TGF.dep, '#B1D4E0',
                      'darkgrey'))

names(keyvals)[keyvals == 'darkgrey'] <- 'Not significant'
names(keyvals)[keyvals == '#B1D4E0'] <- 'AZD8055-insensitive'
names(keyvals)[keyvals == '#145DA0'] <- 'AZD8055-sensitive'


#OR for rapamycin:
#Set the colour scheme to label Rapa-dep genes ONLY:
keyvals <- ifelse(
        rownames(TGFctrl) %in% Rapa.dep$V1, 'royalblue',
        'lightgrey')

names(keyvals)[keyvals == 'lightgrey'] <- 'Rapamycin-insensitive'
names(keyvals)[keyvals == 'royalblue'] <- 'Rapamycin-sensitive'


#Set the colour scheme to label Rapa-dep genes AND label Rapa-indep genes too:
keyvals <- ifelse(
        rownames(TGFctrl) %in% Rapa.dep$V1, '#FBC740',
        ifelse(rownames(TGFctrl) %in% TGF.dep, 'powderblue',
               'darkgrey'))

names(keyvals)[keyvals == 'darkgrey'] <- 'Not significant'
names(keyvals)[keyvals == 'powderblue'] <- 'Rapamycin-insensitive'
names(keyvals)[keyvals == '#FBC740'] <- 'Rapamycin-sensitive'

### TGF response with the AZD-dependent interaction term genes highlighted
#Set colour of plot
#color = c('darkgrey','darkgrey','darkgrey','#66D2D6') #change the fourth colour

#Plot volcano plot logFC2 and p-adjvalue0.05
p1 <- EnhancedVolcano(
        #Load data
        TGFctrl,
        x = 'log2FoldChange',
        y = 'padj',
        
        #Labelling
        lab = rownames(TGFctrl),
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







