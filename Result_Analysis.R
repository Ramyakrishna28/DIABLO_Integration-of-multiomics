# ---- 9. Plots: plotIndiv, circosPlot, networks ----
# 9a. Individual plot (component 1)
safe_write_plot(function() {
  plotIndiv(final_model, comp = 1, group = Y, ind.names = FALSE, legend = TRUE, title = "DIABLO - Component 1")
}, "results/plots/plotIndiv_comp1.pdf")

# 9b. Circos plot (component 1)
safe_write_plot(function() {
  circosPlot(final_model,
             comp = 1,
             cutoff = 0.8,
             Y.name = "Condition",
             var.names = list(
               transcriptomics = gene_names_trans,
               methylation     = gene_names_meth,
               proteomics      = gene_names_prot
             ))
}, "results/plots/circos_comp1.pdf", width = 10, height = 10)

# 9c. Network edges and igraph plotting per component
for (comp in 1:ncomp) {
  edges_info <- try(network(final_model, comp = comp, blocks = names(X), cutoff = 0.7)$edges, silent = TRUE)
  if (inherits(edges_info, "try-error") || is.null(edges_info) || nrow(edges_info) == 0) {
    warning(paste("No edges returned for comp", comp, "- try lowering cutoff."))
    next
 }
  edges_df <- as.data.frame(edges_info)
  edges_file <- sprintf("results/edges_comp%d.csv", comp)
  write.csv(edges_df, edges_file, row.names = FALSE)
  # first two columns are treated as node names (works with typical mixOmics output)
  el <- edges_df[, 1:2]
  colnames(el) <- c("from", "to")
  g <- igraph::graph_from_data_frame(el, directed = FALSE)
  safe_write_plot(function() {
    plot(g, vertex.label.cex = 0.7, vertex.size = 6, main = paste("DIABLO network comp", comp))
  }, sprintf("results/plots/network_comp%d.pdf", comp), width = 10, height = 8)
}

# ---- 10. Selected variables per block & component ----
for (comp in 1:ncomp) {
  sel <- selectVar(final_model, comp = comp)
  # sel is a list with entries per block
  for (block in names(sel)) {
    sel_block <- sel[[block]]
    df_sel <- NULL
    if (is.data.frame(sel_block)) {
      df_sel <- sel_block
    } else if (!is.null(sel_block$name)) {
      # try to make a data.frame
      df_sel <- data.frame(Variable = sel_block$name, sel_block$value)
    } else {
      df_sel <- as.data.frame(sel_block)
    }
    outname <- sprintf("results/selected_vars_%s_comp%d.csv", block, comp)
    write.csv(df_sel, outname, row.names = FALSE)
  }
}

# ---- 11. Save per-block selected variable lists (simple lists) ----
selected_vars_all <- lapply(1:ncomp, function(comp) selectVar(final_model, comp = comp))
saveRDS(selected_vars_all, file = "results/selected_vars_all_comps.rds")

# ---- 12. Helper: expression plots for selected genes ----
plot_selected_gene <- function(block_df, gene, group_vector, outpath = NULL) {
  # block_df: original (not transposed) df with Feature column as first col
  # gene: gene name (must match block_df$Feature)
  expr_row <- block_df[block_df$Feature == gene, -1, drop = FALSE]
  if (nrow(expr_row) == 0) {
    warning("Gene not found in block: ", gene)
    return(NULL)
  }

expr_vals <- as.numeric(expr_row[1, ])
  sample_names_local <- colnames(block_df)[-1]
  dfp <- data.frame(Sample = sample_names_local, Expression = expr_vals, Group = group_vector[sample_names_local])
  p <- ggplot(dfp, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot() + geom_jitter(width = 0.15, size = 1) +
    ggtitle(paste0("Expression - ", gene)) + theme_minimal()
  if (!is.null(outpath)) ggsave(outpath, plot = p, width = 6, height = 4)
  return(p)
}
