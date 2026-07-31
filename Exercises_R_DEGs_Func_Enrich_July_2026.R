############################################
# R BASICS FOR TRANSCRIPTOMICS
############################################

### SECTION 1 — Calculator Basics

# Exercise 1:
# Perform the following operations:
# 1. Add 5 + 3
# 2. Multiply 7 * 6
# 3. Divide 20 / 4
# Write your code below:


# Exercise 2:
# Compute the following expression:
# (10 + 5) * 2


############################################
### SECTION 2 — Variables

# Exercise 3:
# Create two variables:
# a = 15
# b = 4
# Then compute:
# - sum
# - difference
# - product
# Write your code below:


############################################
### SECTION 3 — Vectors (Gene Expression)

# Exercise 4:
# Create a vector representing gene expression:
gene_expr <- c(10, 50, 30, 80, 25)

# Tasks:
# - Calculate mean expression
# - Find max value
# - Find min value


# Exercise 5:
# Extract values greater than 40 (overexpressed genes)
gene_expr[gene_expr > 40]
genes[gene_expr > 40]

############################################
### SECTION 4 — Data Frames (RNA-seq style)
############################################

# Exercise 6:
# Expanded fictitious dataset (12 genes) so downstream exercises
# (volcano plot, GO, GSEA) have enough data points to be meaningful
df <- data.frame(
  gene    = c("TP53", "BRCA1", "MYC", "EGFR", "PTEN",
              "CDK1", "CCNB1", "BAX", "BCL2", "IL6",
              "TNF", "GAPDH"),
  control = c(10, 50, 30, 80, 25, 15, 20, 40, 60, 12, 18, 100),
  treated = c(20, 55, 90, 70, 30, 45, 65, 10, 58, 40, 55, 98)
)

# Tasks:
# - View dataset
# - Calculate fold change (treated / control)
# - Add fold_change column
# - Calculate log2 fold change (logFC) and add it as a column
# - Add fictitious p-values (as if from a differential expression test)
# - Adjust p-values for multiple testing (adjusted p-value / FDR)

View(df)

df$treated / df$control
df$fold_change <- df$treated / df$control

# logFC = log2(fold change) — standard in transcriptomics (e.g. DESeq2, limma)
df$logFC <- log2(df$fold_change)

# Fictitious raw p-values, one per gene
df$p_value <- c(0.001, 0.20, 0.0005, 0.03, 0.40,
                0.002, 0.001, 0.01, 0.35, 0.02,
                0.03, 0.90)

# Adjusted p-value (FDR correction, Benjamini-Hochberg method)
df$adj_p_value <- p.adjust(df$p_value, method = "BH")

df

############################################
### SECTION 5 — Volcano Plot (RNA-seq Visualization)
############################################

# Exercise 7:
# A volcano plot shows statistical significance (-log10 adj p-value)
# vs magnitude of change (logFC) for every gene.

# Define significance thresholds
logFC_cutoff <- 1        # e.g. |logFC| > 1 means at least 2-fold change
padj_cutoff  <- 0.05

# Classify each gene as Up / Down / Not significant
df$status <- "Not significant"
df$status[df$logFC > logFC_cutoff & df$adj_p_value < padj_cutoff]  <- "Upregulated"
df$status[df$logFC < -logFC_cutoff & df$adj_p_value < padj_cutoff] <- "Downregulated"

# Assign colors based on status
colors <- ifelse(df$status == "Upregulated", "red",
                 ifelse(df$status == "Downregulated", "blue", "grey"))

plot(df$logFC, -log10(df$adj_p_value),
     col = colors, pch = 19,
     xlab = "log2 Fold Change",
     ylab = "-log10(adjusted p-value)",
     main = "Volcano Plot: Treated vs Control")

abline(v = c(-logFC_cutoff, logFC_cutoff), lty = 2)
abline(h = -log10(padj_cutoff), lty = 2)

text(df$logFC, -log10(df$adj_p_value), labels = df$gene, pos = 3, cex = 0.7)

# Questions:
# - Which genes are significantly upregulated?
# - Which genes are significantly downregulated?
# - Why do we use -log10(p-value) instead of the raw p-value on the y-axis?

############################################
### SECTION 6 — Gene Ontology (GO) Basics
############################################

# Gene Ontology groups genes into shared biological categories
# (Biological Process, Molecular Function, Cellular Component).
# Here we simulate a tiny, fictitious GO annotation to understand
# the LOGIC of enrichment analysis (not real biological annotations).

# Exercise 8:
# Fictitious GO annotation table (gene -> GO term)
go_annotation <- data.frame(
  gene = c("TP53", "BRCA1", "MYC", "EGFR", "PTEN",
           "CDK1", "CCNB1", "BAX", "BCL2", "IL6",
           "TNF", "GAPDH"),
  go_term = c("Apoptosis", "DNA repair", "Cell cycle", "Signaling", "Apoptosis",
              "Cell cycle", "Cell cycle", "Apoptosis", "Apoptosis", "Immune response",
              "Immune response", "Metabolism")
)

# Task 1: Get the list of significant genes from the volcano plot classification
sig_genes <- df$gene[df$status != "Not significant"]
sig_genes

# Task 2: Merge with GO annotation to see which categories these genes belong to
merge(data.frame(gene = sig_genes), go_annotation, by = "gene")

# Task 3: Count how many significant genes fall into each GO term
table(merge(data.frame(gene = sig_genes), go_annotation, by = "gene")$go_term)

# Task 4 (conceptual): Is "Apoptosis" over-represented among significant genes
# compared to all genes? Build a 2x2 contingency table and test with fisher.test()
apoptosis_sig   <- sum(merge(data.frame(gene = sig_genes), go_annotation, by = "gene")$go_term == "Apoptosis")
apoptosis_total <- sum(go_annotation$go_term == "Apoptosis")
n_sig   <- length(sig_genes)
n_total <- nrow(df)

contingency_table <- matrix(
  c(apoptosis_sig, n_sig - apoptosis_sig,
    apoptosis_total - apoptosis_sig, n_total - n_sig - (apoptosis_total - apoptosis_sig)),
  nrow = 2,
  dimnames = list(c("In gene set", "Not in gene set"), c("Apoptosis", "Other"))
)
contingency_table
fisher.test(contingency_table)

# Questions:
# - What does the p-value from fisher.test() tell you here?
# - This is the core logic behind real GO enrichment tools
#   (e.g. clusterProfiler::enrichGO) — what do you think those tools
#   do differently/better than this toy example?

############################################
### SECTION 7 — GSEA (Gene Set Enrichment Analysis) Basics
############################################

# Unlike GO overrepresentation (which needs a fixed "significant gene" cutoff),
# GSEA uses the FULL ranked gene list (e.g. ranked by logFC or test statistic)
# and asks: does a gene set (pathway) cluster near the top or bottom of the ranking?

# Exercise 9:
# Task 1: Rank all genes by logFC (descending)
ranked_df <- df[order(-df$logFC), c("gene", "logFC")]
ranked_df

# Task 2: Define a fictitious pathway gene set (e.g. "Cell cycle")
pathway_genes <- go_annotation$gene[go_annotation$go_term == "Cell cycle"]
pathway_genes

# Task 3: Find where the pathway genes fall in the ranked list
ranked_df$in_pathway <- ranked_df$gene %in% pathway_genes
ranked_df

# Task 4 (conceptual): A simple "running score" idea —
# +1 step when we hit a pathway gene, -small step otherwise,
# walking down the ranked list (this is a simplified illustration
# of the real GSEA enrichment score, not the actual algorithm)
step_hit  <- 1
step_miss <- -1 / (nrow(ranked_df) - length(pathway_genes))

running_score <- cumsum(ifelse(ranked_df$in_pathway, step_hit, step_miss))
running_score

# Plot the running score with gene names on the x-axis, so you can see
# exactly which gene drives each jump in the curve
plot(running_score, type = "l",
     xaxt = "n",
     xlab = "",
     ylab = "Running enrichment score (toy version)",
     main = "Simplified GSEA-style Running Score: Cell Cycle pathway")

# Add gene names as x-axis labels (rotated for readability)
axis(1, at = 1:nrow(ranked_df), labels = ranked_df$gene, las = 2, cex.axis = 0.7)
mtext("Rank (by logFC, high to low)", side = 1, line = 4)

abline(h = 0, lty = 2)

# Mark the pathway genes directly on the curve (like the tick marks
# under a real GSEA enrichment plot)
points(which(ranked_df$in_pathway), running_score[ranked_df$in_pathway],
       col = "red", pch = 19, cex = 1.2)

text(which(ranked_df$in_pathway), running_score[ranked_df$in_pathway],
     labels = ranked_df$gene[ranked_df$in_pathway],
     pos = 3, col = "red", cex = 0.8, font = 2)

# Questions:
# - Does the "Cell cycle" gene set accumulate its score early (top of ranking,
#   upregulated) or late (bottom, downregulated)?
# - What would it mean biologically if a pathway peaked near the top vs the bottom?
# - In real GSEA (e.g. fgsea, clusterProfiler::GSEA), what additional step
#   is used to test whether the enrichment score is statistically significant?

############################################
# END OF EXERCISES
############################################