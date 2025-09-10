# DIABLO pipeline (generalized) with plotting and saving
# Author: [Your Name]
# Place input files in data/ and outputs will be saved to results/

# ---- 0. Packages ----
pkgs <- c("mixOmics", "dplyr", "ggplot2", "reshape2", "igraph")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

set.seed(123)

# ---- 1. File paths (edit if needed) ----
transcriptomics_file <- "data/Transcriptomics.csv"
methylation_file     <- "data/Methylation.csv"
proteomics_file      <- "data/Proteomics.csv"
meta_file            <- "data/sample_metadata.csv"   # optional: columns Sample,Group

# ---- 2. Create results folders ----
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/plots")) dir.create("results/plots", recursive = TRUE)

# ---- 3. Helper functions ----
read_and_clean <- function(path, feature_col = 1) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  # remove rows with NA in feature column and duplicates
  df <- df[!is.na(df[[feature_col]]), ]
  df <- df[!duplicated(df[[feature_col]]), ]
  colnames(df)[feature_col] <- "Feature"
  return(df)
}

safe_write_plot <- function(expr_fun, outpath, width = 8, height = 6) {
  pdf(outpath, width = width, height = height)
  try(expr_fun())
  dev.off()
}

# ---- 4. Load data ----
transcriptomics <- read_and_clean(transcriptomics_file)
methylation     <- read_and_clean(methylation_file)
proteomics      <- read_and_clean(proteomics_file)

# store feature lists
gene_names_trans <- transcriptomics$Feature
gene_names_meth  <- methylation$Feature
gene_names_prot  <- proteomics$Feature

# ---- 5. Build X list: samples x features ----
# Expect CSVs where first column = feature, remaining columns = samples
make_block_matrix <- function(df) {
  mat <- as.matrix(df[ , -1, drop = FALSE])
  mat_t <- t(mat)            # samples as rows
  colnames(mat_t) <- df$Feature
  rownames(mat_t) <- colnames(df)[-1]
  return(mat_t)
}

X <- list(
  transcriptomics = make_block_matrix(transcriptomics),
  methylation     = make_block_matrix(methylation),
  proteomics      = make_block_matrix(proteomics)
)

# check sample name consistency
sample_names <- rownames(X$transcriptomics)
if (is.null(sample_names) || length(sample_names) == 0) stop("No sample names found in transcriptomics file (check columns).")

# ensure other blocks have same samples (order)
for (b in names(X)) {
  if (!all(rownames(X[[b]]) == sample_names)) {
    # try to reorder if possible
    common <- intersect(sample_names, rownames(X[[b]]))
    if (length(common) == nrow(X[[b]])) {
      X[[b]] <- X[[b]][sample_names, , drop = FALSE]
    } else {
      warning(sprintf("Samples in block '%s' do not match transcriptomics samples. You may need to re-order or provide metadata.", b))
    }
  }
}

# ---- 6. Define sample groups (Y) ----
if (file.exists(meta_file)) {
  meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Sample", "Group") %in% colnames(meta))) {
    warning("metadata file found but should contain columns 'Sample' and 'Group'. Falling back to default grouping.")
    use_default_grouping <- TRUE
  } else {
    # create named vector of groups ordered to sample_names
    meta <- meta[match(sample_names, meta$Sample), ]
    if (any(is.na(meta$Group))) stop("Some sample names in metadata do not match sample names in the data.")
    Y <- factor(meta$Group)
    names(Y) <- meta$Sample
  }
} else {
  # default grouping: split samples 50/50 into Group1 and Group2 (user should replace)
  n <- length(sample_names)
  half <- floor(n / 2)
  groups <- c(rep("Group1", half), rep("Group2", n - half))
  Y <- factor(groups)
  names(Y) <- sample_names
  message("No metadata found at 'data/sample_metadata.csv'. Default grouping created: first half Group1, second half Group2.")
}

# reorder Y to match rows in X
Y <- factor(as.character(Y[ sample_names ]))

# ---- 7. DIABLO parameters (adjust) ----
ncomp <- 2

# recommended keepX should be <= number of variables in each block
cap_keepX <- function(x, val) { min(ncol(x), val) }
keepX <- list(
  transcriptomics = rep(cap_keepX(X$transcriptomics, 50), ncomp),
  methylation     = rep(cap_keepX(X$methylation, 50), ncomp),
  proteomics      = rep(cap_keepX(X$proteomics, 50), ncomp)
)

# design matrix: strength of links between blocks (0-1)
design <- matrix(0.5, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# ---- 8. Run DIABLO (block.splsda) ----
final_model <- block.splsda(X, Y = Y, ncomp = ncomp, keepX = keepX, design = design)
message("DIABLO model trained.")
