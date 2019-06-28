## 01.0-Package Import ---------------------------------------------

library(tidyverse)
library(rmarkdown)
library(magrittr)
library(fdrtool)
library(xml2)
library(readr)
library(xtable)

# Load Bioconductor Packages
library(limma)
library(DESeq2)
library(vsn)

## 02.0-Import ---------------------------------------------

# Load counts matrix output from feature counts
read.counts <- read.table("FeatureCounts Output/featureCounts_results.txt", header = TRUE)

# Store ensembl gene IDs as row names
row.names(read.counts) <- read.counts$Geneid

# Exclude the columns without read counts (columns 1 to 6 contain additional
# info such as genomic coordinates)
read.counts <- read.counts[ , -c(1:6)]

# Order and store original col names
read.counts <- read.counts[ , order(names(read.counts))]
orig_names <- names(read.counts)

# Renaming col names
names(read.counts) <- c("S01_WT_CNV_RC_1", "S02_WT_CNV_NR_2", "S04_GAL_CNV_RC_1", "S04_GAL_CNV_NR_1", 
                       "S06_GAL_CNV_RC_3", "S05_GAL_CNV_NR_2", "S05_GAL_CNV_RC_2", "S03_WT_CNV_NR_3",
                       "S06_GAL_CNV_NR_3", "S07_IHH_CNV_RC_1", "S07_IHH_CNV_NR_1", "S08_IHH_CNV_RC_2",
                       "S08_IHH_CNV_NR_2", "S09_IHH_CNV_RC_3", "S09_IHH_CNV_NR_3", "S10_IHH_noCNV_NR_1", 
                       "S10_IHH_noCNV_RC_1", "S11_IHH_noCNV_NR_2", "S11_IHH_noCNV_RC_2", "S12_IHH_noCNV_NR_3",
                       "S12_IHH_noCNV_RC_3", "S13_GAL_noCNV_NR_1", "S13_GAL_noCNV_RC_1", "S14_GAL_noCNV_NR_2",
                       "S14_GAL_noCNV_RC_2", "S15_GAL_noCNV_NR_3", "S15_GAL_noCNV_RC_3", "S16_WT_noCNV_NR_1",
                       "S16_WT_noCNV_RC_1", "S17_WT_noCNV_NR_2", "S02_WT_CNV_RC_2", "S17_WT_noCNV_RC_2",
                       "S18_WT_noCNV_NR_3", "S18_WT_noCNV_RC_3", "S01_WT_CNV_NR_1", "S03_WT_CNV_RC_3")

# Order new col names
read.counts = read.counts[ , order(names(read.counts))]

# Create colData data.frame
sample_info <- DataFrame(condition = gsub("_[0-9]+", "", names(read.counts)),
                         row.names = names(read.counts))

# Generate DESeqDataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = read.counts,
                                   colData = sample_info,
                                   design = ~ condition)
# Plot counts per sample
colSums(counts(DESeq.ds)) %>% barplot

## 03.0-Filtering and Normalizing ---------------------------------------------

# Filter genes with no reads
keep_genes <- rowSums(counts(DESeq.ds)) > 0
DESeq.ds <- DESeq.ds[keep_genes, ]

nonlog.counts <- counts(DESeq.ds, normalized = FALSE)
head(non.lo)

# Calculate size factors and plot against counts per sample
DESeq.ds <- estimateSizeFactors(DESeq.ds)
plot(sizeFactors(DESeq.ds), colSums(counts(DESeq.ds)))

# 03.1-Normalizing ---------------------------------------------

# Plot read counts and read counts normalized by size factors
par(mfrow=c(1,2)) # to plot the two box plots next to each other
## bp of non-normalized
boxplot(log2(counts(DESeq.ds)+1), notch=TRUE,
        main = "Non-normalized read counts",
        ylab="log2(read counts)")

# bp of size-factor normalized values
boxplot(log2(counts(DESeq.ds, normalize=TRUE)+1), notch=TRUE,
        main = "Size-factor-normalized read counts",
        ylab="log2(read counts)")

# non-normalized read counts plus pseudocount
log.counts <- log2(counts(DESeq.ds, normalized = FALSE) + 1)
## Assign log2 non-normalized counts to a distinct matrix within the DESeq.ds object
assay(DESeq.ds, "log.counts") <- log2(counts(DESeq.ds, normalized = FALSE) + 1)
## Assign log2 normalized counts to a distinct matrix within the DESeq.ds object
log.norm.counts <- log2(counts(DESeq.ds, normalized=TRUE) + 1)
assay(DESeq.ds, "log.norm.counts") <- log.norm.counts


# Plot log2 normalized gene counts of samples against other replicates
DESeq.ds[, c("S16_WT_noCNV_NR_1","S18_WT_noCNV_NR_3")] %>%
  assay(., "log.norm.counts") %>%
  plot(., cex=.1, main = "WT_noCNV_NR_1 vs. WT_noCNV_NR_3")

# 03.2-Reducing Variance Dependence on the mean (heteroskedasticity) ---------------------------------------------

log.norm.counts <- log2(counts(DESeq.ds, normalized=TRUE) + 1)
msd_plot <- vsn::meanSdPlot(log.norm.counts,
                            ranks=FALSE, # show the data on the original scale
                            plot = FALSE)
msd_plot$gg +
  ggtitle("Sequencing depth normalized log2(read counts)") +
  ylab("standard deviation")

# Stabilize Variance
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE)

# Plot size factor and log2-transformed vs rlog transformed data
par(mfrow=c(1,2))
plot(log.norm.counts[,1:2], cex=.1,
     main = "size factor and log2-transformed")
## the rlog-transformed counts are stored in the accessor "assay"
plot(assay(DESeq.rlog)[,7],
     assay(DESeq.rlog)[,8],
     cex=.1, main = "rlog transformed",
     xlab = colnames(assay(DESeq.rlog[,7])),
     ylab = colnames(assay(DESeq.rlog[,8])) )


rlog.norm.counts <- assay(DESeq.rlog)

msd_plot <- vsn::meanSdPlot( rlog.norm.counts, ranks=FALSE, plot = FALSE)
msd_plot$gg + ggtitle("rlog transformation")

save.image(file = "6618_Analysis_Pre_Process.RData")





