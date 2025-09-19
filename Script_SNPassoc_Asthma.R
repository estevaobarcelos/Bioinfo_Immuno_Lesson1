#------------------------------#
# SNPassoc Analysis Script
# Association study of asthma case-control data with 51 SNPs
#------------------------------#
# From: https://cran.r-project.org/web/packages/SNPassoc/vignettes/SNPassoc.html

# Load required library
library(SNPassoc)

# 2. Data loading
# Load the asthma dataset included in the SNPassoc package
data(asthma, package = "SNPassoc")

# Examine the structure of the dataset
str(asthma, list.len=9)
# Results: Dataset contains 1578 observations with 57 variables including
# demographic data (country, gender, age, BMI, smoking status), case-control status,
# and 51 SNP genotypes encoded as factors with 3 levels (AA, AG, GG format)

# View first few rows
asthma[1:5, 1:8]
# Results: Shows the structure of demographic variables and first few SNPs

# 3. Descriptive analysis
# Identify SNP columns (those starting with "rs")
idx <- grep("^rs", colnames(asthma))

# Set up SNP data with proper class assignment
asthma.s <- setupSNP(data=asthma, colSNPs=idx, sep="")
# Comment: sep="" indicates no separator between alleles (e.g., "GG" not "G/G")

# Examine a single SNP
head(asthma.s$rs1422993)
class(asthma.s$rs1422993)
# Results: SNP is now of class "snp" and "factor"

# Summary of a single SNP
summary(asthma.s$rs1422993)
# Results: Shows genotype frequencies (G/G: 57.2%, G/T: 36.1%, T/T: 6.7%),
# allele frequencies (G: 75.3%, T: 24.7%), and HWE p-value (0.25)

# Visualize SNP distribution
plot(asthma.s$rs1422993)
plot(asthma.s$rs1422993, type=pie)
# Results: Bar and pie charts showing genotype distribution

# Summary of all SNPs
summary(asthma.s, print=FALSE)
# Results: Table showing major allele frequencies, HWE p-values, and missing rates
# for all 51 SNPs. Most SNPs have acceptable HWE p-values (>0.05) and low missing rates

# Check missing data patterns
plotMissing(asthma.s, print.labels.SNPs = FALSE)
# Results: Missing values appear mostly random with some clustering in consecutive SNPs
# suggesting potential genotyping issues in those regions

# 4. Hardy-Weinberg equilibrium
# Test HWE for all SNPs
hwe <- tableHWE(asthma.s)
head(hwe)
# Results: First SNPs show HWE p-values > 0.05, indicating no deviation from HWE

# Test HWE separately in cases and controls
hwe2 <- tableHWE(asthma.s, casecontrol)

# Identify SNPs that are in HWE in full sample but not in controls
snpNHWE <- hwe2[,1]>0.05 & hwe2[,2]<0.05
rownames(hwe2)[snpNHWE]
hwe2[snpNHWE,]
# Results: rs1345267 shows HWE violation in controls (p=0.0496) but not in full sample

# Filter SNPs to keep only those passing HWE in controls (p >= 0.001)
snps.ok <- rownames(hwe2)[hwe2[,2]>=0.001]
pos <- which(colnames(asthma)%in%snps.ok, useNames = FALSE)
asthma.s <- setupSNP(asthma, pos, sep="")
# Comment: Using 0.001 threshold as quality control measure

# 5. SNP association analysis
# Single SNP association test
association(casecontrol ~ rs1422993, data = asthma.s)
# Results: Multiple genetic models tested. Codominant (p=0.018), dominant (p=0.008),
# and overdominant (p=0.005) models show significant association.
# OR for dominant model: 1.39 (95% CI: 1.09-1.77)

# Max-statistic test
maxstat(asthma.s$casecontrol, asthma.s$rs1422993)
# Results: MAX-statistic = 7.073, p=0.0179, confirming association

# # Test specific genetic model
# association(casecontrol ~ rs1422993, asthma.s, model="dominant")
# # Results: Only dominant model tested, p=0.0078, OR=1.39

# Adjust for covariates
association(casecontrol ~ rs1422993 + country + smoke, asthma.s)
# Results: After adjusting for country and smoking, association weakens
# Dominant model p=0.034, OR=1.33 (95% CI: 1.02-1.73)

# Stratified analysis by gender
association(casecontrol ~ rs1422993 + survival::strata(gender), asthma.s)
# Results: Association remains significant after stratification

# Subset analysis for Spain only
association(casecontrol ~ rs1422993, asthma.s, subset=country=="Spain")
# Results: No significant association in Spanish subgroup (p=0.21 for dominant model)

# Quantitative trait analysis (BMI)
association(bmi ~ rs1422993, asthma.s)
# Results: No significant association between rs1422993 and BMI (all p>0.66)

# Genome-wide association analysis
ans <- WGassociation(casecontrol, data=asthma.s)
head(ans)
# Results: P-values for all SNPs across different genetic models
# No strong associations detected in initial screening

# Adjusted genome-wide analysis
ans.adj <- WGassociation(casecontrol ~ country + smoke, asthma.s)
head(ans.adj)
# Results: Covariate-adjusted p-values for all SNPs

# Fast scanning (additive model only)
# Note: Requires development version from GitHub
# ans.fast <- scanWGassociation(casecontrol, asthma.s)

# Manhattan plot of results
plot(ans)
# Results: Visual representation of -log10(p-values) for all SNPs
# No SNPs reach Bonferroni significance threshold (0.05/51 ≈ 0.001)
# No strong evidence of association between any SNP and asthma in this dataset

# 7. Haplotype analysis

# 7.1 Haplotype estimation
library(haplo.stats)

# Select SNPs for haplotype analysis (example: three adjacent SNPs)
snpsH <- c("rs714588", "rs1023555", "rs898070")

# Prepare genotype data for haplotype analysis
genoH <- make.geno(asthma.s, snpsH)

# Estimate haplotype frequencies using expectation-maximization algorithm
em <- haplo.em(genoH, locus.label = snpsH, miss.val = c(0, NA))

# View estimated haplotype frequencies
print(em)
# Results: Shows estimated frequencies for all possible haplotype combinations
# Common haplotypes typically have frequencies >5%, while rare haplotypes are pooled

# 7.2 Haplotype association
trait <- asthma.s$casecontrol

# Fit haplotype association model using logistic regression
mod <- haplo.glm(trait ~ genoH,           
                 family = "binomial", 
                 locus.label = snpsH,
                 allele.lev = attributes(genoH)$unique.alleles,
                 control = haplo.glm.control(haplo.freq.min = 0.05))

# Display results with confidence intervals
intervals(mod)
# Results: Shows odds ratios and 95% confidence intervals for each haplotype
# compared to the reference haplotype (usually the most common one)

# Interpretation example:
# - Haplotypes with OR significantly >1 indicate increased asthma risk
# - Haplotypes with OR significantly <1 indicate protective effect
# - The global test evaluates whether haplotypes are collectively associated with the trait

# Session info
sessionInfo()
# Comment: Record package versions for reproducibility