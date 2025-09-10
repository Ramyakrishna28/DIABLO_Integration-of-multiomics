# ------------------------------
# General boxplot function
# ------------------------------
plot_selected_omics <- function(df, gene_col = "Gene.names", selected_genes, group_df, 
                                title_prefix = "Expression", outpath = NULL) {
  # df: omics data (first col = gene, rest = samples)
  # gene_col: column containing gene names
  # selected_genes: vector of genes to keep
  # group_df: data.frame(Sample, Group)
  # title_prefix: title for plot
  # outpath: optional path to save PDF/PNG
  
  library(reshape2)
  library(ggplot2)
  
  # clean gene names
  df[[gene_col]] <- sub("\\..*", "", df[[gene_col]])
  
  # subset by selected genes
  df_sub <- df[df[[gene_col]] %in% selected_genes, ]
  
  # melt into long format
  melted <- melt(df_sub, id.vars = gene_col)
  colnames(melted) <- c("Gene", "Sample", "Expression")
  
  # add group info
  melted$Group <- group_df$Group[match(melted$Sample, group_df$Sample)]
  
  # make plot
  p <- ggplot(melted, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Gene, scales = "free") +
    ggtitle(paste(title_prefix, "of Selected Genes")) +
    theme_minimal()
  
  # save if path given
  if (!is.null(outpath)) {
    ggsave(outpath, plot = p, width = 8, height = 6)
  }
  
  return(p)
}

# ------------------------------
# Example usage
# ------------------------------

# Metadata: adjust group labels as per your experiment
group_df_methylation <- data.frame(
  Sample = colnames(methylation)[-1],
  Group = c(rep("Flight", 7), rep("Ground", 7))
)

group_df_proteomics <- data.frame(
  Sample = colnames(proteomics)[-1],
  Group = c(rep("Flight", 7), rep("Ground", 7))
)

# Plot methylation
plot_selected_omics(
  df = methylation,
  gene_col = "Gene.names",
  selected_genes = selected_genes$methylation,
  group_df = group_df_methylation,
  title_prefix = "Methylation",
  outpath = "results/plots/methylation_boxplots.pdf"
)

# Plot proteomics
plot_selected_omics(
  df = proteomics,
  gene_col = "Gene.names",
  selected_genes = selected_genes$proteomics,
  group_df = group_df_proteomics,
  title_prefix = "Proteomics",
  outpath = "results/plots/proteomics_boxplots.pdf"
)
