# Load libraries
library(VennDiagram)
library(gridExtra)
library(grid)
library(dplyr)
library(ggplotify)
library(ggplot2)

args <- commandArgs(TRUE)

top_num <- as.integer(args[3])

###########################################
##
## Venn for each top set
##
###########################################

# Put all sets into a list for easy looping
set_names <- c("pval", "accuracy", "balanced_accuracy", "f1_score", "precision", "recall", "roc_auc")

sets_list <- list(vector("list", 100),
                  vector("list", 100),
                  vector("list", 100),
                  vector("list", 100),
                  vector("list", 100),
                  vector("list", 100),
                  vector("list", 100))

for (i in 1:length(set_names)){
  for (j in 1:100){
    df_true <- read.csv(paste0(args[1], j, "_betas_causal.csv"), header=T)
    df <- read.csv(paste0(args[2], j, "_", set_names[i], ".csv"), header=T)
    sets_list[[i]][[j]] <- intersect(df$SNP[1:top_num], df_true$SNP)
  }
  
}

# ---- averaged region proportions (unchanged idea) ----
pair_avg_props_7x100 <- function(sets_list, idx1, idx2) {
  R <- length(sets_list[[1]])
  props <- matrix(0, nrow = 3, ncol = R,
                  dimnames = list(c("Aonly","Bonly","Both"), NULL))
  for (r in seq_len(R)) {
    A <- sets_list[[idx1]][[r]]
    B <- sets_list[[idx2]][[r]]
    U <- union(A, B)
    if (length(U) == 0L) {
      props[, r] <- c(0,0,0)
    } else {
      props[, r] <- c(
        length(setdiff(A, B)) / 100,
        length(setdiff(B, A)) / 100,
        length(intersect(A, B)) / 100
      )
    }
  }
  rowMeans(props, na.rm = TRUE)  # named vector: Aonly, Bonly, Both
}

# ---- draw one averaged pair with xx.xx% labels using VennDiagram only ----
avg_pair_venn_percent_VD <- function(sets_list, set_names = NULL, idx1 = 1, idx2 = 2,
                                     fill = c("skyblue", "lightpink"),
                                     cat.cex = 1.1, cex = 1.3, scale = 1000) {
  if (is.null(set_names)) {
    set_names <- if (!is.null(names(sets_list))) names(sets_list) else paste0("Set", seq_along(sets_list))
  }

  avgp <- pair_avg_props_7x100(sets_list, idx1, idx2)  # Aonly, Bonly, Both

  # Convert to areas for drawing (any common scale works)
  area1      <- (avgp["Aonly"] + avgp["Both"]) * scale
  area2      <- (avgp["Bonly"] + avgp["Both"]) * scale
  cross_area <-  avgp["Both"] * scale

  # Format % labels
  lab_Aonly <- sprintf("%.2f%%", 100 * avgp["Aonly"])
  lab_Bonly <- sprintf("%.2f%%", 100 * avgp["Bonly"])
  lab_Both  <- sprintf("%.2f%%", 100 * avgp["Both"])

  # Build the venn grobs but don't draw to device yet
  venn_out <- VennDiagram::draw.pairwise.venn(
    area1      = area1,
    area2      = area2,
    cross.area = cross_area,
    category   = c(set_names[idx1], set_names[idx2]),
    fill       = fill,
    alpha      = rep(0.6, 2),
    cat.cex    = cat.cex,
    cex        = rep(cex, 3),
    margin     = 0.05,
    ind        = FALSE  # return grobs
  )

  # Replace default numeric labels with our percentage labels.
  # In typical VennDiagram output, the first three text grobs are region labels
  # (A-only, B-only, Both), followed by two category labels.
  text_idx <- which(vapply(venn_out, function(g) inherits(g, "text"), logical(1)))

  if (length(text_idx) >= 3) {
    venn_out[[ text_idx[1] ]]$label <- lab_Aonly
    venn_out[[ text_idx[2] ]]$label <- lab_Bonly
    venn_out[[ text_idx[3] ]]$label <- lab_Both
  } else {
    # Fallback to historical 7:9 positions if text detection behaves differently
    if (length(venn_out) >= 9) {
      venn_out[[7]]$label <- lab_Aonly
      venn_out[[8]]$label <- lab_Bonly
      venn_out[[9]]$label <- lab_Both
    }
  }

  # Add a simple title as a top text grob
  ttl <- grid::textGrob(
    sprintf("%s vs %s", set_names[idx1], set_names[idx2]),
    y = grid::unit(1, "npc") - grid::unit(4, "mm"),
    gp = grid::gpar(fontface = 2)
  )

  grid::grobTree(
    grid::rectGrob(gp = grid::gpar(col = NA, fill = NA)), # noop background
    ttl,
    grid::grobTree(venn_out)
  )
}

# ---- grid of "first vs. others" (percent labels, VennDiagram-only) ----
avg_venn_first_vs_others_percent <- function(sets_list, set_names = NULL, ncol = 3) {
  if (is.null(set_names)) {
    set_names <- if (!is.null(names(sets_list))) names(sets_list) else paste0("Set", seq_along(sets_list))
  }
  K <- length(sets_list)
  grobs <- lapply(2:K, function(j) {
    avg_pair_venn_percent_VD(sets_list, set_names = set_names, idx1 = 1, idx2 = j)
  })
  nrow <- ceiling((K - 1) / ncol)
  gridExtra::arrangeGrob(grobs = grobs, ncol = ncol, nrow = nrow)
}



set_names_print <- c("P-value", "Accuracy", "Balanced Accuracy", "F1 Score", "Precision", "Recall", "AUC")

p_grid <- as.ggplot(avg_venn_first_vs_others_percent(sets_list, set_names_print))

ggsave(paste0(args[4]),
       p_grid, width = 12, height = 8, dpi = 300)