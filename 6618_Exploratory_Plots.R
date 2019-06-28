## 01.0-Package Import ---------------------------------------------

library(magrittr)
library(ggplot2)
library(pheatmap)

# Load Bioconductor Packages
library(DESeq2)
library(org.Mm.eg.db)

# Load Functions File
source("6618_Func.R")

## 02.0-Data Import ---------------------------------------------

load(file = "6618_Analysis_Pre_Process.RData")

## 03.0-Exploratory Plots ---------------------------------------------

corr_coeff <- cor(rlog.norm.counts, method = "pearson")


# Dendrograms before and after rlog
par(mfrow=c(1,2))
# Pearson corr. for rlog.norm values
as.dist(1 - corr_coeff) %>% hclust %>%
  plot( ., labels = colnames(rlog.norm.counts),
        main = "rlog transformed read counts")
# Pearson corr. for log.norm.values
as.dist( 1 - cor(log.norm.counts, method = "pearson")) %>%
  hclust %>% plot( ., labels = colnames(log.norm.counts),
                   main = "no rlog")


# Plot Heatmap
heat_plot(rlog.norm.counts)

# Plot RPE Choroid Heatmap
rlog.norm.counts.RC <- rlog.norm.counts[, grep("*RC_[0-9]", 
      colnames(rlog.norm.counts), value=TRUE)]
heat_plot(rlog.norm.counts.RC)

# Plot Heatmap
rlog.norm.counts.GAL <- rlog.norm.counts[, grep("GAL*", 
      colnames(rlog.norm.counts), value=TRUE)]
heat_plot(rlog.norm.counts.GAL)


# Plot RPE Choroid with IHH KO Heatmap
rlog.norm.counts.GAL.RC <- rlog.norm.counts[, grep("GAL*", 
      grep("*RC_[0-9]", colnames(rlog.norm.counts), value=TRUE), value=TRUE)]
heat_plot(rlog.norm.counts.GAL.RC)

rlog.norm.counts.WT <- rlog.norm.counts[, grep("WT*", 
      colnames(rlog.norm.counts), value=TRUE)]

tmp12 <- merge(rlog.norm.counts.WT, rlog.norm.counts.GAL, by=0, all=T)
rownames(tmp12) <- tmp12$Row.names; tmp12$Row.names <- NULL

rlog.norm.counts.WT.GAL <- rbind(rlog.norm.counts.WT, rlog.norm.counts.GAL)


# Calculate Principal Components and Plot

# Using the same rlog normalized count dataframes 
# created for the heatmaps above
paired_pca_plot(test)


test <- data.matrix(tmp12)

# List top 25 genes that contribute to PC1
top_genes(rlog.norm.counts)


