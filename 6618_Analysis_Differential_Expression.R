## 01.0-Package Import ---------------------------------------------

library(magrittr)

# Load Bioconductor Packages
library(edgeR)
library(limma)
library(org.Mm.eg.db)

# Load Functions File
source("6618_Func.R")

## 02.0-Data Import ---------------------------------------------

load(file = "6618_Analysis_Pre_Process.RData")

## 03.0-Model Setup ---------------------------------------------

# Split sample name into factors which we put into a representative model matrix
factor.matrix = strsplit(names(read.counts), '_') %>%
  unlist() %>%
  matrix(ncol=5, byrow=TRUE)

# Pull out the factor lists from the matrix and transform into factors
# Individual factors are here
sub <- factor(factor.matrix[,1])
mut <- factor(factor.matrix[,2], levels=c("WT", "GAL", "IHH"))
pert <- factor(factor.matrix[,3], levels=c("noCNV", "CNV"))
tiss <- factor(factor.matrix[,4], levels=c("RC", "NR"))

# Joined factors are here (basically identical to what we started with)
treatfactor <- factor(paste(factor.matrix[,2], factor.matrix[,3], factor.matrix[,4], sep="."))

# Create model matrix
design <- model.matrix(~0+treatfactor)
colnames(design) <- gsub("treatfactor", "", colnames(design))

# Create contrast matrix
contr.matrix <- makeContrasts(RCvsNR_WTnoCNV = WT.noCNV.RC-WT.noCNV.NR,
                              noCNVvsCNV_WTRC = WT.noCNV.RC-WT.CNV.RC,
                              WTvsGAL_noCNVRC = WT.noCNV.RC-GAL.noCNV.RC,
                              WTvsGAL_CNVRC = WT.CNV.RC-GAL.CNV.RC,
                              noCNVvsCNV_GalRC = GAL.noCNV.RC-GAL.CNV.RC,
                              WTvsIHH_noCNVRC = WT.noCNV.RC-IHH.noCNV.RC,
                              WTvsIHH_CNVRC = WT.CNV.RC-IHH.CNV.RC,
                              noCNVvsCNV_IHHRC = IHH.noCNV.RC-IHH.CNV.RC,
                              levels=design)

## Data pre-processing ---------------------------------------------

# Create edgeR DGEList object
DE.counts <- DGEList(read.counts, samples=colnames(read.counts), 
                     group=treatfactor, remove.zeros=TRUE)

#Filter Genes with low count-per-million reads within and across exp groups
keep.exprs <- filterByExpr(DE.counts, design)
DE.counts <- DE.counts[keep.exprs, ]

# Normaliztion by the method of trimmed mean of M-value for non-biological factors
DE.counts <- calcNormFactors(DE.counts, method="TMM")

# Estimate the correlation between measurements made on the same subject
# (tissue samples from the same mouse)
corfit <- limma::duplicateCorrelation(DE.counts$counts, design, block=sub)
corfit$consensus

## 04.0-Analysis ---------------------------------------------

# Normalize and reduce heteroscedascity
# Plot mean-variance trend
par(mfrow=c(1,2))
voom_transformed <- voom(DE.counts, design, plot=TRUE)

# Fit a linear model for each gene
voomed.fitted <- lmFit(voom_transformed, design, block=sub, correlation=corfit$consensus)
voomed.fitted <- contrasts.fit(voomed.fitted, contrasts=contr.matrix)

# Compute differentian expression
voomed.fitted <- eBayes(voomed.fitted)

# Plot mean-variance trend with reduced hetroscedascity
plotSA(voomed.fitted, main="Final model: Mean-variance trend")

## 05.0-Summary ---------------------------------------------


# Extract all results to summary table
noCNVvsCNV_WTRC_results <- topTable(voomed.fitted, coef="noCNVvsCNV_GalRC", 
                          number=Inf, adjust.method="BH")

# Subset genes by DE adjusted p-values
DEgenes <- rownames(subset(noCNVvsCNV_WTRC_results, adj.P.Val < 1))

# Convert identifying scheme
anno.DE <- select(org.Mm.eg.db,
                   keys = DEgenes, # rownames
                   keytype= "ENSEMBL", # rownames are ENSEMBL identifiers
                   columns= c("SYMBOL", "GENENAME")) # what to return

anno.trim.DE <- aggregate(anno.DE['SYMBOL'], anno.DE['ENSEMBL'], 
                  FUN = function(X) paste(unique(X), collapse=", "))

# Merge dataframes so that all indentifying schemes are included
out.df <- merge(noCNVvsCNV_WTRC_results, anno.trim.DE, 
                by.x = "row.names", by.y = "ENSEMBL", all.y=TRUE)


library(xlsx)
write.xlsx(out.df, "GAL.noCNV.RC-GAL.CNV.RC.xlsx")






anno.DE[0:10,]

# Extract results as matrix representing upregulation, 
# downregulation, or not significant
dt <- decideTests(voomed.fitted)

# Find comman DE genes between different comparisons
de.common <- names(which(dt[, 'noCNVvsCNV_WTRC']!=0 & 
                           dt[, 'WTvsGAL_noCNVRC']!=0))

# Convert identifying scheme
anno.common.DE <- select(org.Mm.eg.db,
                  keys = de.common, # rownames
                  keytype= "ENSEMBL", # rownames are ENSEMBL identifiers
                  columns= c("SYMBOL", "GENENAME")) # what to return

## 06.0-Results Plots ---------------------------------------------

test_results <- topTable(voomed.fitted, coef="RCvsNR_WTnoCNV", 
                                    number=Inf, adjust.method="BH")

# Plot Adjusted p-values
ggplot(data = test_results, aes(test_results$P.Value)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05)) +
  labs(tile = "C vs. u1A P-value Distribution") +
  labs(x = "P-values", y = "Count")


# MD Plots
plotMD(voomed.fitted, column=1, status=dt[, 'RCvsNR_WTnoCNV'], 
       main=colnames(voomed.fitted)['RCvsNR_WTnoCNV'], xlim=c(-8,13))

# Uncorrected Volcano Plots
uncorr_vol(noCNVvsCNV_WTRC_results)

# Corrected Volcano Plots
corr_vol(noCNVvsCNV_WTRC_results)


