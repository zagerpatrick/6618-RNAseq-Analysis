heat_plot <- function(norm_counts) {
  corr_coeff.RC <- cor(norm_counts, method = "pearson")
  p <- as.dist(1-corr_coeff.RC, upper = TRUE) %>% as.matrix %>%
    pheatmap::pheatmap(., main = "Pearson correlation")
  return(p)
}

paired_pca_plot <- function(norm_counts) {
  colnames(norm_counts) <- gsub("_[0-9]", "", colnames(norm_counts))
  colnames(norm_counts) <- gsub("S[0-9][0-9]_", "", colnames(norm_counts))
  
  rv <- rowVars(norm_counts)
  top_variable <- order(rv, decreasing = TRUE)[seq_len(500)]
  pca <- prcomp(t(norm_counts[top_variable, ]))
  
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  pca.data <- data.frame(Sample=rownames(pca$x),
                         X=pca$x[,1], 
                         Y=pca$x[,2])
  
  p <- ggplot(data=pca.data, aes(x=X, y=Y, colour=Sample)) +
    geom_point(size=3) +
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
    scale_colour_brewer(palette="Paired") +
    theme_bw()
  return(p)
}


top_genes <- function(norm_counts) {
  rv <- rowVars(norm_counts)
  top_variable <- order(rv, decreasing = TRUE)[seq_len(500)]
  pca <- prcomp(t(norm_counts[top_variable, ]))

  loading_scores <- pca$rotation[,1]
  gene_scores <- abs(loading_scores) ## get the magnitudes
  gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
  top_genes <- names(gene_score_ranked[1:25])
  top_genes_hr <- AnnotationDbi::select(org.Mm.eg.db,
                       keys = top_genes, # rownames
                       keytype= "ENSEMBL", # rownames are ENSEMBL identifiers
                       columns= c("SYMBOL", "GENENAME")) # what to return
  return(top_genes_hr)
}


corr_vol <- function(df) {
  p <- ggplot(df, aes(df$logFC, -log10(df$adj.P.Val))) +
    geom_point(size = 2, alpha = 1, na.rm = T, 
               shape = 16, colour = "black") +
    geom_hline(yintercept = -log10(0.05), 
               colour = "orange", linetype = "dashed", size = 1) +
    scale_x_continuous(expand=c(0, 0), limits = c(-10, 10)) +
    ggtitle(label = "", subtitle = "Moderated t-test") +
    ylab("-log10 Adjusted p-value") + 
    xlab("Log Fold Change")
  return(p)
}

uncorr_vol <- function(df) {
  p <- ggplot(df, aes(df$logFC, -log10(df$P.Value))) +
    geom_point(size = 2, alpha = 1, na.rm = T, 
               shape = 16, colour = "black") +
    scale_x_continuous(expand=c(0, 0), limits = c(-10, 10)) +
    ggtitle(label = "", subtitle = "Moderated t-test") +
    ylab("-log10 Adjusted p-value") + 
    xlab("Log Fold Change")
  return(p)
}