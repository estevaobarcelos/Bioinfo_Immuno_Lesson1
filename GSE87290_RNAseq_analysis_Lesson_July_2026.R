################################################################################
# GSE87290 — RNA-seq of human PBMCs before/after LPS (endotoxin) challenge
# Teaching script: from raw counts to differential expression
#
# Study design recap (from GEO):
#   15 healthy subjects, PBMCs sampled at BASELINE and after IV LPS (endotoxin).
#   Subjects were pre-classified as HIGH or LOW inflammatory responders based
#   on their clinical cytokine response. This gives a 2x2 factorial design:
#       group:     high / low
#       timepoint: baseline / LPS
#   with subject as a paired/blocking factor (each subject contributes both
#   a baseline and an LPS sample -> paired design, NOT independent samples).
#
# Goal of the exercise: teach the full workflow —
#   1. Import a GEO-provided raw count matrix + gene annotation
#   2. Build sample metadata correctly (matching by ID, not row order!)
#   3. QC: library size, filtering, PCA / sample clustering
#   4. Differential expression with DESeq2, accounting for pairing
#      (LPS vs baseline, and the group x timepoint interaction)
#   5. Visualization: PCA, volcano, heatmap
#   6. GSEA (ranked-list) on the LPS-vs-baseline response
#   7. High vs Low responders AT EACH timepoint (unpaired, between-subject),
#      each followed by GO over-representation analysis (ORA) to summarize
#      what distinguishes the two groups at baseline vs. after LPS
#
# We use the NCBI-generated "raw counts" matrix rather than the .diff files,
# because it lets students go through the standard DESeq2 workflow rather
# than working with pre-computed Cuffdiff output.
################################################################################


## ---- 0. Setup ---------------------------------------------------------- ##

# Install (only run once per machine) ----
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "EnhancedVolcano", "clusterProfiler",
#                         "org.Hs.eg.db", "enrichplot", "pheatmap"))
# install.packages(c("tidyverse", "ggrepel"))

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)   # nicer volcano plots
library(ggrepel)
library(clusterProfiler)   # GSEA
library(org.Hs.eg.db)      # gene ID / GO annotation for GSEA
library(enrichplot)        # GSEA plotting (gseaplot2, dotplot)

set.seed(42)

# Working directory for downloads / outputs
dir.create("GSE87290_analysis", showWarnings = FALSE)
setwd("GSE87290_analysis")


## ---- 1. Download the data ----------------------------------------------##
# GEO now hosts NCBI-harmonized, gene-level raw count matrices for many
# RNA-seq series, generated with a consistent pipeline (STAR + featureCounts
# vs Annotation Release 109.20190905, genome GRCh38.p13). These are far more
# convenient for teaching than re-aligning FASTQs, and the counts are
# comparable across GEO series processed the same way.

counts_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE87290&format=file&file=GSE87290_raw_counts_GRCh38.p13_NCBI.tsv.gz"
annot_url  <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"

download.file(counts_url, destfile = "GSE87290_raw_counts.tsv.gz", mode = "wb")
download.file(annot_url,  destfile = "Human.GRCh38.p13.annot.tsv.gz", mode = "wb")

counts_raw <- read.delim(gzfile("GSE87290_raw_counts.tsv.gz"),
                         row.names = 1, check.names = FALSE)
annot <- read.delim(gzfile("Human.GRCh38.p13.annot.tsv.gz"))

dim(counts_raw)     # genes x samples
head(counts_raw[, 1:5])
colnames(counts_raw) # should be GSM accessions


## ---- 2. Build sample metadata ------------------------------------------##
# NEVER assume the sample order in the counts matrix matches the order on
# the GEO web page — always build metadata keyed by GSM accession and then
# match() it against colnames(counts_raw). This is the single most common
# source of silent, catastrophic errors in RNA-seq analysis (mislabeled
# groups -> completely wrong DE results that still "look" plausible).

# Metadata transcribed directly from the GEO series page (Samples table).
# subject = replicate number, which is the SAME individual at both timepoints
# (paired design).
sample_info <- tribble(
  ~geo_accession,  ~subject, ~timepoint,  ~group,
  "GSM2326887",    1,        "baseline",  "high",
  "GSM2326888",    2,        "baseline",  "high",
  "GSM2326889",    3,        "baseline",  "high",
  "GSM2326890",    4,        "baseline",  "high",
  "GSM2326891",    5,        "baseline",  "high",
  "GSM2326892",    6,        "baseline",  "high",
  "GSM2326893",    7,        "baseline",  "low",
  "GSM2326894",    8,        "baseline",  "low",
  "GSM2326895",    9,        "baseline",  "low",
  "GSM2326896",    10,       "baseline",  "high",
  "GSM2326897",    11,       "baseline",  "high",
  "GSM2326898",    12,       "baseline",  "low",
  "GSM2326899",    13,       "baseline",  "low",
  "GSM2326900",    14,       "baseline",  "low",
  "GSM2326901",    15,       "baseline",  "low",
  "GSM2326902",    1,        "LPS",       "high",
  "GSM2326903",    2,        "LPS",       "high",
  "GSM2326904",    3,        "LPS",       "high",
  "GSM2326905",    4,        "LPS",       "high",
  "GSM2326906",    5,        "LPS",       "high",
  "GSM2326907",    6,        "LPS",       "high",
  "GSM2326908",    7,        "LPS",       "low",
  "GSM2326909",    8,        "LPS",       "low",
  "GSM2326910",    9,        "LPS",       "low",
  "GSM2326911",    10,       "LPS",       "high",
  "GSM2326912",    11,       "LPS",       "high",
  "GSM2326913",    12,       "LPS",       "low",
  "GSM2326914",    13,       "LPS",       "low",
  "GSM2326915",    14,       "LPS",       "low",
  "GSM2326916",    15,       "LPS",       "low"
) %>%
  mutate(
    subject   = factor(paste0("S", subject)),
    timepoint = factor(timepoint, levels = c("baseline", "LPS")),
    group     = factor(group, levels = c("low", "high"))
  )

# Match metadata to the columns actually present in the counts matrix
sample_info <- sample_info[match(colnames(counts_raw), sample_info$geo_accession), ]
rownames(sample_info) <- sample_info$geo_accession


## ---- 3. Gene annotation -------------------------------------------------##
# counts_raw rows are Entrez GeneIDs; annot maps GeneID -> Symbol/description
head(annot)
gene_map <- annot %>%
  dplyr::select(GeneID, Symbol, Description) %>%
  dplyr::mutate(GeneID = as.character(GeneID))


## ---- 4. Build the DESeq2 dataset ----------------------------------------##
# Design: ~ subject + timepoint  models the paired baseline-vs-LPS effect
# while controlling for inter-individual variability (a paired t-test
# analogue). We'll build a second dataset with an interaction term to ask
# whether the LPS response itself differs between high vs low responders.

counts_mat <- as.matrix(counts_raw)
storage.mode(counts_mat) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sample_info,
  design    = ~ subject + timepoint
)

# Annotate rowData with gene symbols for later convenience
rowData(dds)$GeneID <- rownames(dds)
rowData(dds)$Symbol <- gene_map$Symbol[match(rownames(dds), gene_map$GeneID)]


## ---- 5. QC & filtering ---------------------------------------------------##

# Library sizes
colSums(counts(dds)) %>%
  enframe(name = "geo_accession", value = "total_counts") %>%
  left_join(as_tibble(sample_info), by = "geo_accession") %>%
  ggplot(aes(x = reorder(geo_accession, total_counts), y = total_counts,
             fill = interaction(group, timepoint))) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Total mapped reads", fill = "group.timepoint",
       title = "Library size per sample") +
  theme_minimal()

# Standard low-count filter: keep genes with a reasonable count in at least
# as many samples as the smallest group (here, 7 low-responder subjects)
keep <- rowSums(counts(dds) >= 10) >= 7
dds <- dds[keep, ]
nrow(dds)   # genes retained for testing

# Variance-stabilizing transform for QC visualization (NOT for DE testing)
vsd <- vst(dds, blind = TRUE)


## ---- 6. Exploratory analysis: PCA & clustering ---------------------------##

pca_data <- plotPCA(vsd, intgroup = c("group", "timepoint"), returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = group, shape = timepoint,
                     label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.5, show.legend = FALSE) +
  labs(x = paste0("PC1: ", pct_var[1], "% variance"),
       y = paste0("PC2: ", pct_var[2], "% variance"),
       title = "PCA of PBMC transcriptomes: baseline vs LPS, high vs low responders") +
  theme_minimal()
# Expect: samples separate strongly by timepoint (LPS effect dominates),
# and within LPS, high responders should separate further from low responders.

# Sample-to-sample distance heatmap
sample_dist <- dist(t(assay(vsd)))
pheatmap(as.matrix(sample_dist),
         annotation_col = as.data.frame(colData(vsd)[, c("timepoint", "group")]),
         main = "Sample-to-sample distances")


## ---- 7. Differential expression: LPS vs baseline (paired) ----------------##

dds <- DESeq(dds)
resultsNames(dds)

res_lps_vs_baseline <- results(dds, name = "timepoint_LPS_vs_baseline",
                               alpha = 0.05)

# Shrink noisy log2FoldChange estimates for ranking/plotting. We use the
# classic normal-prior shrinkage (type = "normal") here, which is built into
# DESeq2 and needs no extra package (apeglm/ashr are alternatives but add a
# dependency). It shrinks LFCs of low-count/high-variance genes toward zero
# without touching the p-values from results() above.
res_lps_vs_baseline <- lfcShrink(dds, coef = "timepoint_LPS_vs_baseline",
                                 res = res_lps_vs_baseline, type = "normal")

res_df <- as.data.frame(res_lps_vs_baseline) %>%
  rownames_to_column("GeneID") %>%
  left_join(gene_map, by = "GeneID") %>%
  arrange(padj)

summary(res_lps_vs_baseline)
write_csv(res_df, "DE_LPS_vs_baseline_allSubjects.csv")

# Volcano plot

# Create a function to replace the function of the package EnhancedVolcano
make_volcano <- function(res_df, plot_title,
                         p_cutoff = 0.05, fc_cutoff = 1,
                         max_labels = 15) {
  
  plot_df <- res_df %>%
    mutate(
      negLog10Padj = -log10(pmax(padj, .Machine$double.xmin)),
      significance = case_when(
        !is.na(padj) & padj < p_cutoff & log2FoldChange >= fc_cutoff ~ "Up",
        !is.na(padj) & padj < p_cutoff & log2FoldChange <= -fc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  # Select the most significant genes for labeling
  label_df <- plot_df %>%
    filter(significance != "Not significant",
           !is.na(Symbol),
           Symbol != "") %>%
    arrange(padj) %>%
    slice_head(n = max_labels)
  
  ggplot(plot_df,
         aes(x = log2FoldChange,
             y = negLog10Padj)) +
    
    geom_point(aes(color = significance),
               alpha = 0.7,
               size = 1.8) +
    
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed") +
    
    geom_hline(yintercept = -log10(p_cutoff),
               linetype = "dashed") +
    
    geom_text_repel(
      data = label_df,
      aes(label = Symbol),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      show.legend = FALSE
    ) +
    
    scale_color_manual(
      values = c(
        "Up" = "firebrick",
        "Down" = "steelblue",
        "Not significant" = "grey70"
      )
    ) +
    
    labs(
      title = plot_title,
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Significance"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

# Plot the Volcano Plot
make_volcano(
  res_df,
  "LPS vs baseline (all subjects, paired design)"
)

# Heatmap of top 30 DE genes
top_genes <- res_df %>% filter(!is.na(padj)) %>% slice_min(padj, n = 30) %>% pull(GeneID)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)  # center per gene for visualization
rownames(mat) <- gene_map$Symbol[match(rownames(mat), gene_map$GeneID)]

pheatmap(mat,
         annotation_col = as.data.frame(colData(vsd)[, c("timepoint", "group")]),
         show_rownames = TRUE,
         main = "Top 30 DE genes: LPS vs baseline")


## ---- 8. Does the LPS response differ between high vs low responders? -----##
# This is the more interesting biological question (matches the study's own
# framing): test the GROUP x TIMEPOINT interaction.

dds_int <- DESeqDataSetFromMatrix(
  countData = counts_mat[rownames(dds), ],   # reuse filtered gene set
  colData   = sample_info,
  design    = ~ group + timepoint + group:timepoint
)
dds_int <- DESeq(dds_int)
resultsNames(dds_int)

res_interaction <- results(dds_int, name = "grouphigh.timepointLPS", alpha = 0.05)
res_int_df <- as.data.frame(res_interaction) %>%
  rownames_to_column("GeneID") %>%
  left_join(gene_map, by = "GeneID") %>%
  arrange(padj)

summary(res_interaction)
write_csv(res_int_df, "DE_interaction_highVsLow_LPSresponse.csv")

# Genes with a significant interaction term respond differently to LPS
# in high- vs low-responders -- these are the strongest candidates for
# driving the clinical divergence in inflammatory response.
sig_interaction <- res_int_df %>% filter(padj < 0.05)
nrow(sig_interaction)
head(sig_interaction, 20)


## ---- 9. GSEA (Gene Set Enrichment Analysis) of the LPS response ----------##
# Unlike enrichGO/enrichKEGG (over-representation on a significance cutoff),
# GSEA uses the FULL ranked gene list, so it can detect coordinated shifts
# even when no single gene passes a strict padj threshold -- well suited to
# an LPS response where many genes move modestly but consistently.
#
# Two things to teach here:
#   (a) how to build a properly ranked gene list in R and run GSEA directly
#       with clusterProfiler::GSEA / gseGO
#   (b) how to export a standalone .rnk file, in case students want to run
#       the same ranked list through the standalone GSEA desktop tool
#       (Broad Institute) or another package (e.g. fgsea) for comparison.

# 9a. Build the ranked list -------------------------------------------------
# Rank by the Wald statistic ("stat" column from results()) rather than
# log2FoldChange: it combines effect size AND its precision, which is the
# recommended ranking metric for RNA-seq GSEA (log2FC alone over-weights
# noisy, low-count genes). Use the UNSHRUNKEN results for this, since
# lfcShrink() is meant for effect-size reporting, not ranking statistics.

res_for_rank <- as.data.frame(results(dds, name = "timepoint_LPS_vs_baseline")) %>%
  rownames_to_column("GeneID") %>%
  filter(!is.na(stat)) %>%
  left_join(gene_map, by = "GeneID") %>%
  # collapse duplicate symbols by keeping the most significant entry
  group_by(Symbol) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(Symbol), Symbol != "")

ranked_stat <- res_for_rank$stat
names(ranked_stat) <- res_for_rank$GeneID           # Entrez IDs for clusterProfiler
ranked_stat <- sort(ranked_stat, decreasing = TRUE)

# Optional: export a .rnk file (SYMBOL <tab> score) for the standalone
# Broad GSEA desktop app, or for fgsea::fgsea() outside this script.
rnk_export <- res_for_rank %>% dplyr::select(Symbol, stat) %>% dplyr::arrange(dplyr::desc(stat)) 
readr::write_tsv( rnk_export, "LPS_vs_baseline.rnk", col_names = FALSE )

# 9b. Run GSEA against GO Biological Process terms ---------------------------
gsea_go <- gseGO(
  geneList     = ranked_stat,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  eps          = 0,        # get exact small p-values instead of a floor
  seed         = TRUE
)

gsea_go_df <- as.data.frame(gsea_go) %>% arrange(p.adjust)
write_csv(gsea_go_df, "GSEA_GO_BP_LPS_vs_baseline.csv")
head(gsea_go_df[, c("Description", "NES", "pvalue", "p.adjust")], 15)

# Dot plot summarizing top enriched pathways (split up/down = activated vs suppressed)
dotplot(gsea_go, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA (GO:BP) — PBMC response to LPS")

# Classic running-enrichment-score plot for the single most significant term
gseaplot2(gsea_go, geneSetID = 1,
          title = gsea_go_df$Description[1])

# Ridge plot: distribution of ranking metric within each leading gene set
# (useful for teaching what "enrichment score" is actually summarizing)
ridgeplot(gsea_go, showCategory = 15) +
  ggtitle("Distribution of rank statistic within top gene sets")

# 9c. (Optional) same idea against KEGG pathways instead of GO
# gsea_kegg <- gseKEGG(geneList = ranked_stat, organism = "hsa",
#                       pvalueCutoff = 0.05, eps = 0)
# dotplot(gsea_kegg, showCategory = 15)


## ---- 10. High vs Low responders: at baseline AND after LPS ---------------##
# So far every comparison has been WITHIN subject, across timepoint.
# Now the complementary question: at a given timepoint, how do high- and
# low-responders differ from each other? This is a BETWEEN-subject
# comparison (no pairing possible within a single timepoint), so the model
# for each subset is simply ~ group.
#
# We reuse dds's already-filtered gene set for a fair comparison across all
# three contrasts (main LPS effect, interaction, and these two).

run_group_contrast <- function(dds_full, keep_timepoint, gene_map, label) {
  # Subset samples to a single timepoint
  dds_sub <- dds_full[, dds_full$timepoint == keep_timepoint]
  dds_sub$group <- droplevels(dds_sub$group)
  design(dds_sub) <- ~ group
  
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, contrast = c("group", "high", "low"), alpha = 0.05)
  res <- lfcShrink(dds_sub, contrast = c("group", "high", "low"),
                   res = res, type = "normal")
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("GeneID") %>%
    left_join(gene_map, by = "GeneID") %>%
    arrange(padj)
  
  cat("\n===", label, "===\n")
  print(summary(res))
  write_csv(res_df, paste0("DE_highVsLow_", label, ".csv"))
  
  list(dds = dds_sub, results = res, res_df = res_df)
}

baseline_hl <- run_group_contrast(dds, "baseline", gene_map, "baseline")
lps_hl      <- run_group_contrast(dds, "LPS",      gene_map, "postLPS")

# Volcano plots for each timepoint, side by side conceptually
make_volcano(
  baseline_hl$res_df,
  "High vs Low responders — BASELINE"
)

make_volcano(
  lps_hl$res_df,
  "High vs Low responders — AFTER LPS"
)

# How many genes distinguish the groups at each timepoint? A useful number
# to discuss: the study's premise is that high/low responders look alike
# at baseline and diverge only once challenged.
n_sig_baseline <- sum(baseline_hl$res_df$padj < 0.05, na.rm = TRUE)
n_sig_lps      <- sum(lps_hl$res_df$padj < 0.05, na.rm = TRUE)
cat("Significant genes, high vs low, baseline:", n_sig_baseline, "\n")
cat("Significant genes, high vs low, post-LPS:", n_sig_lps, "\n")


## ---- 10b. Gene Ontology (over-representation) on each contrast -----------##
# Here we deliberately use ORA (enrichGO), not GSEA: the question is "what
# biological processes are over-represented among the genes that are
# SIGNIFICANTLY different between high and low responders", using a clean
# significance cutoff -- a more direct, easier-to-explain analysis for a
# first pass than the ranked-list GSEA done in section 9.
#
# universe = every gene that was actually tested (not the whole genome!) --
# this is essential for a well-calibrated ORA p-value.

run_go_ora <- function(res_df, universe_ids, label) {
  sig_genes <- res_df %>% filter(padj < 0.05) %>% pull(GeneID)
  
  if (length(sig_genes) == 0) {
    message("No significant genes for ", label, " — skipping GO ORA.")
    return(NULL)
  }
  
  ego <- enrichGO(
    gene         = sig_genes,
    universe     = universe_ids,
    OrgDb        = org.Hs.eg.db,
    keyType      = "ENTREZID",
    ont          = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable     = TRUE
  )
  
  ego_df <- as.data.frame(ego)
  write_csv(ego_df, paste0("GO_ORA_highVsLow_", label, ".csv"))
  
  if (nrow(ego_df) > 0) {
    print(dotplot(ego, showCategory = 15) +
            ggtitle(paste0("GO:BP enrichment — High vs Low, ", label)))
  } else {
    message("No significantly enriched GO terms for ", label, ".")
  }
  
  ego
}

tested_universe <- rownames(dds)   # same filtered gene set used for all DE tests

go_baseline_hl <- run_go_ora(baseline_hl$res_df, tested_universe, "baseline")
go_lps_hl      <- run_go_ora(lps_hl$res_df,      tested_universe, "postLPS")


# End of the script.
################################################################################