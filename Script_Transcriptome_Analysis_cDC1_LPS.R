#------------------------------------------------------------------------------#
# Pre-processing Microarray Data - QC, Normalization, DEGs and Visualization
#------------------------------------------------------------------------------#
# Analysis of datasets from: GSE203450 (Microarray Affymetric Data from DCs)


# Set Working Directory.--------------------------------------------------------
setwd("") # Set the Working Directory on your computer.

# Part 1. Load required libraries.==============================================
print("loading libraries...")

# Install Bioconductor Manager (needed for some packages).----------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the required libraries.-----------------------------------------------
# BiocManager::install("GEOquery")
library(GEOquery)

# BiocManager::install("limma")
library(limma)

# Part 2. Import and explore data from GEO.=====================================
print("Importing data...")

id <- "GSE203450"
gse <- getGEO(id, GSEMatrix =TRUE, AnnotGPL=TRUE)

# Get some information about the file.----------------------------------------
list(gse)
names(pData(gse[[1]]))      # print the sample info.
length(gse)                 # check how many platforms was used.

# Select the dataset.---------------------------------------------------------
gse <- gse[[1]]

# if more than one dataset is present, you can analyse the other dataset by 
# changing the number inside the [[...]]
# e.g. gse <- gse[[2]]

# We can also use that kind of line.
# if (length(gse) > 1) idx <- grep("GPL6246", attr(gse, "names")) else idx <- 1


# Explore the expression dataframe.
pData(gse)                    # print the sample info.
fData(gse)                    # print the gene annotation.
exprs(gse)                    # print the expression data.


# Part 3. Group membership for all samples.=====================================
print("Grouping the samples...")

# Select one comparison between groups.
cdc1_lps <- "XXXXXX111000XXXXXX"          # (cDC1 LPS vs cDC1 Ctrl)

# Split Samples.----------------------------------------------------------------
sml <- strsplit(cdc1_lps, split ="")[[1]]


# Part 4. Filter out excluded samples and Create the expression matrix.=========
print("Filtering out excluded samples...")
sel <- which(sml != "X")   # excluded samples marked as "X".
sml <- sml[sel]
gse <- gse[ ,sel]

head(gse)

# Create the expression matrix.-----------------------------------------------
ex <- exprs(gse)
ex[which(ex <= 0)] <- NaN

head(ex)
hist(ex)

# Part 5. Log2 transformation and Normalization of the data.====================
print("Transforming and Normalizing data...")

# Box-and-whisker plot (all samples before normalization).--------------------
par(mar=c(7,4,2,1))
title <- paste ("GSE203450", "/", annotation(gse), sep ="")

# pdf("Microarray_data_before_transformation.pdf")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=F, las=2)
# dev.off()

# Log2 transformation.----------------------------------------------------------
exprs(gse) <- log2(ex)       # Log2 transform.

# Box-and-whisker plot (after log2 transformation).-----------------------------
par(mar=c(7,4,2,1))
title <- paste ("GSE203450", "/", annotation(gse), sep ="")

# pdf("Microarray_data_after_transformation.pdf")
boxplot(exprs(gse), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
# dev.off()

# Expression value distribution plot.-------------------------------------------
par(mar=c(4,4,2,1))
title <- paste ("GSE203450", "/", annotation(gse), " value distribution", 
                sep ="")

# pdf("Densities_after_normalization.pdf")
plotDensities(ex, main=title, legend=F)
# dev.off()

# Part 6. Normalization log-ratios with Limma.==================================
print("Normalizing data with Limma...")
exprs(gse) <- normalizeBetweenArrays(exprs(gse))        # normalize data
gse

plotDensities(gse, main=title, legend=F)


# Part 7. Assign samples to groups and set up design matrix.====================

# Define groups to compare.-----------------------------------------------------
groups <- make.names(c("cDC1_LPS","cDC1_ctrl"))

gs <- factor(sml)
levels(gs) <- groups
gse$group <- gs
design <- model.matrix(~group + 0, gse)
colnames(design) <- levels(gs)

gse <- gse[complete.cases(exprs(gse)), ] # skip missing values

# PCA and UMAP plot (dimensionality reduction).---------------------------------
# library("maptools")  # point labels without overlaps
library("umap")
library("car")

ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates

# a) PCA plot.
pca <- prcomp(t(ex), scale. = TRUE) # Transpose e scaling

par(mar=c(3,3,2,6), xpd=TRUE)
plot(pca$x[, 1:2], main="PCA plot", xlab="PC1", ylab="PC2", 
     col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.25,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
pointLabel(pca$x[, 1:2], labels = rownames(pca$x), method="SANN", cex=0.6)

# b) Plot UMAP.-----------------------------------------------------------------
ump <- umap(t(ex), n_neighbors = 3, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=3", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.25,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# Stop Point.
################################################################################


# Part 8. Differential expression with Limma-Voom.==============================
print("Calculating the differential gene expression...")

# Calculate precision weights and show plot of mean-variance trend.-------------
v <- vooma(gse, design, plot=T)

v$genes <- fData(gse) # attach gene annotations

# Fit linear model.-------------------------------------------------------------
fit  <- lmFit(v)


# Set up contrasts of interest and recalculate model coefficients.--------------
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# Compute statistics and table of top significant genes.------------------------
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)


# Part 9. Visualize adj p-values, Venn diagram and QQ plot.=====================
print("Visualizing adj p-values and venn diagram...")

# Build histogram of P-values for all genes. Normal test.-----------------------
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# Summarize test results as "up", "down" or "not expressed".--------------------
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0.5)


# Venn diagram of differential genes between groups.----------------------------
# pdf("Venn_diagram.pdf")
vennDiagram(dT, circle.col=palette())
# dev.off()


# Part 10. DEGs list.=========================================================== 
print("Creating a DEGs list...")
library(tidyverse)
DEG <- tT2[c(3,1,22,26)] %>% 
  setNames(c("gene","id","logFC","padj")) 

# DEGs selection (Define groups and cut-off points of logFC and padj).----------
DEG$Enriched <- "NS"
DEG$Enriched[DEG$logFC > 0.5 & DEG$padj < 0.05] <- "cDC1_LPS"
DEG$Enriched[DEG$logFC < -0.5 & DEG$padj < 0.05] <- "cDC1_ctrl"
# 

# Filter DEGs
DEGs <- filter(DEG, Enriched != "NS")

# Save as a CSV file in the current working directory
write.csv(DEGs, file = "degs.csv", row.names = FALSE)


# Part 11. Volcano plot Representation (logFC =|0.5|).====================
print("Creating a VP representation...")

# Obtain counts of differexpressed genes.---------------------------------------
DEG %>% 
  count(Enriched)


# Add colour, size and alpha (transparency) to volcano plot.--------------------
cols <- c("cDC1_ctrl" = "darkblue", "cDC1_LPS" = "firebrick", "NS" = "lightgrey")


# a) Select genes to highlight.=================================================

genes <- DEG %>%
  filter(gene %in% c("Il6","Il1b","Tnf","Cd86","Cd80","Clec9a","Xcr1","Cd40"))


# b) Volcano Plot with selected genes using ggplot.=============================

# install.packages("ggrepel") # OR
# devtools::install_github("slowkow/ggrepel")
library(ggrepel)

# Visualization Volcano Plot.---------------------------------------------------
ggplot(data = DEG,
       aes(x = logFC,
           y = -log10(padj))) + 
  geom_point(aes(colour = Enriched), 
             alpha = 0.5, 
             shape = 16,
             size = 1) + 
  geom_point(data = genes,
             shape = 21,
             size = 2, 
             fill = "black", 
             colour = "black") + 
  theme_classic() +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.7071), log2(1.4142)),
             linetype = "dashed") +
  geom_label_repel(data = genes, # Add labels last to appear as the top layer 
                   max.overlaps = Inf, 
                   aes(label = gene),fontface = 'italic', 
                   force = 1, nudge_y = 0.5) +
  scale_colour_manual(values = cols) +
  # scale_x_continuous(breaks = c(seq(-4, 5, 2)),     
  #                    limits = c(-4, 5),)  + 
  ggtitle("BMDCs/cDC1-ctrl versus cDC1-LPS enriched genes")



# Part 12. GO (Gene Ontology) enrichment Analysis.==============================

# BiocManager::install("DOSE")
library(DOSE)

# BiocManager::install("clusterProfiler")
library(clusterProfiler)

# BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

#------------------------------------------------------------------------------#
# # Select only group 1 for analysis.---------------------------------------------
# DEG_1 <- filter(DEG, Enriched == "cDC1_ctrl")
# DEG_1$entrez = mapIds(org.Mm.eg.db, 
#                       keys=as.character(DEG_1$gene), 
#                       column = "ENTREZID",
#                       keytype = "SYMBOL", 
#                       multiVals = "first")
# genelist_1 <- DEG_1$entrez               # Create a list with genes.
# go_cDC1 <- enrichGO(gene = genelist_1, 
#                     OrgDb = org.Mm.eg.db,
#                     pvalueCutoff = 0.05,
#                     qvalueCutoff = 0.05,
#                     ont="all",                #BP, CC, MF or all
#                     readable = T)

#------------------------------------------------------------------------------#
# Select only group 2 for analysis.---------------------------------------------
DEG_LPS <- filter(DEG, Enriched == "cDC1_LPS")
DEG_LPS$entrez = mapIds(org.Mm.eg.db, 
                        keys=as.character(DEG_LPS$gene), 
                        column = "ENTREZID",
                        keytype = "SYMBOL", 
                        multiVals = "first")
DEG_LPS <- DEG_LPS$entrez             # Create a list with genes.
go_cDC1_LPS <- enrichGO(gene = DEG_LPS, 
                        OrgDb = org.Mm.eg.db,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        ont="all",                #BP, CC, MF or all
                        readable = T)
# ## Visualization.-------------------------------------------------------------
# Dot plot.---------------------------------------------------------------------

# pdf("Dot Plot displaying GO enrichment analysis for BMDCs cDC1 stimulated with LPS.pdf")
dotplot(go_cDC1_LPS,showCategory = 7,font.size = 10) + 
  ggtitle("GO enrichment for BMDCs cDC1 stimulated with LPS")
# dev.off()

# Gene-Concept Network.---------------------------------------------------------
# pdf(("Cnetplot Top3 displaying GO enrichment analysis for BMDCs cDC1 stimulated with LPS.pdf"),20,10)
cnetplot(go_cDC1_LPS, showCategory = 3, circular = T, colorEdge = T) + 
  ggtitle("GO enrichment for BMDCs cDC1 stimulated with LPS")
# dev.off()


# Part 13. GSEA.================================================================
print("Gene Set Enrichment Analysis...")

DEG$entrez = mapIds(org.Mm.eg.db, keys=as.character(DEG$gene), column = "ENTREZID",
                    keytype = "SYMBOL", multiVals = "first")
genelist <- DEG$entrez                        # Create a list with genes.

rnk_gsea <- dplyr::bind_cols(DEG$gene,
                             as.numeric(-log10(DEG$padj) *
                                          sign(DEG$logFC)))
# rnk_gsea <- rnk_gsea[order(rnk_gsea[,2], decreasing = TRUE),]
rnk_gsea <- rnk_gsea[!is.na(rnk_gsea$...2),]
rnk_gsea_clean <- rnk_gsea[rnk_gsea$...2 !=-Inf,]
rnk_gsea_clean <- rnk_gsea_clean[rnk_gsea_clean$...2 !=Inf,]


write.table(rnk_gsea_clean, file = "BMDC_DC1_ctrl_vs_LPS.rnk",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")            # load this file into the GSEA program.

# End of script.
################################################################################