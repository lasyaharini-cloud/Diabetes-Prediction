# Dimensionality reduction (PCA / t-SNE / UMAP, 2D and 3D) for the three
# diabetes datasets, plus class-overlap metrics used to explain why class
# imbalance makes minority classes hard to predict.
#
# Outputs (in dr_results/):
#   dr_embeddings.rds          cached 2D/3D embeddings + labels
#   figures/dr_2d_<ds>.png      PCA | t-SNE | UMAP in 2D (one row, shared legend)
#   figures/dr_3d_<ds>.png      PCA | t-SNE | UMAP in 3D
#   tables/overlap_metrics.csv  kNN class-purity per class per dataset
#   tables/class_counts.csv     class sizes and imbalance ratios

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(Rtsne)
  library(umap)
  library(scatterplot3d)
  library(FNN)
})

set.seed(42)
outdir <- "dr_results"
dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

d012 <- read.csv("data/diabetes_012_health_indicators_BRFSS2021.csv")
dbal <- read.csv("data/diabetes_binary_5050split_health_indicators_BRFSS2021.csv")
preds <- setdiff(names(d012), "Diabetes_012")

datasets <- list(
  balanced = list(
    X = dbal[, preds],
    y = factor(dbal$Diabetes_binary, levels = c(0, 1),
               labels = c("No diabetes", "Diabetes")),
    palette = c("No diabetes" = "#2166AC", "Diabetes" = "#B2182B"),
    title = "Balanced binary"
  ),
  imbalanced = list(
    X = d012[, preds],
    y = factor(as.numeric(d012$Diabetes_012 > 0), levels = c(0, 1),
               labels = c("No diabetes", "Prediabetes or diabetes")),
    palette = c("No diabetes" = "#2166AC", "Prediabetes or diabetes" = "#B2182B"),
    title = "Imbalanced binary"
  ),
  threeclass = list(
    X = d012[, preds],
    y = factor(d012$Diabetes_012, levels = c(0, 1, 2),
               labels = c("No diabetes", "Prediabetes", "Diabetes")),
    palette = c("No diabetes" = "#2166AC", "Prediabetes" = "#E69F00", "Diabetes" = "#B2182B"),
    title = "3-class"
  )
)

# stratified subsample preserving class proportions
subsample <- function(y, n = 6000, seed = 42) {
  set.seed(seed)
  idx <- unlist(lapply(split(seq_along(y), y), function(ii) {
    k <- round(n * length(ii) / length(y))
    sample(ii, min(length(ii), max(k, 1)))
  }), use.names = FALSE)
  sort(idx)
}

# 2D faceted figure (shared legend) via ggplot
make_2d <- function(df, palette, title) {
  long <- rbind(
    data.frame(class = df$class, method = "PCA",   D1 = df$PC1,   D2 = df$PC2),
    data.frame(class = df$class, method = "t-SNE", D1 = df$tSNE1, D2 = df$tSNE2),
    data.frame(class = df$class, method = "UMAP",  D1 = df$UMAP1, D2 = df$UMAP2)
  )
  long$method <- factor(long$method, levels = c("PCA", "t-SNE", "UMAP"))
  ggplot(long, aes(D1, D2, colour = class)) +
    geom_point(size = 0.5, alpha = 0.55) +
    facet_wrap(~method, scales = "free", nrow = 1) +
    scale_colour_manual(values = palette, name = NULL) +
    labs(title = paste0(title, ": 2D projections"),
         x = "Dimension 1", y = "Dimension 2") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
}

# 3D figure via scatterplot3d (static, saved to PNG).
# Each axis is scaled to [-1, 1] with fixed limits so the three panels share the
# same box geometry and viewpoint, which makes them easy to compare.
make_3d <- function(df, palette, title, path) {
  cols <- adjustcolor(palette[as.character(df$class)], alpha.f = 0.6)
  nrm <- function(v) { r <- range(v); if (diff(r) == 0) return(rep(0, length(v)))
                       2 * (v - r[1]) / diff(r) - 1 }
  P  <- cbind(nrm(df$PC1),   nrm(df$PC2),   nrm(df$PC3))
  Ts <- cbind(nrm(df$tSNE1), nrm(df$tSNE2), nrm(df$tSNE3))
  U  <- cbind(nrm(df$UMAP1), nrm(df$UMAP2), nrm(df$UMAP3))
  ang <- 50
  lim <- c(-1, 1)
  png(path, width = 2200, height = 820, res = 160)
  par(mfrow = c(1, 3), mar = c(3, 2, 2, 1), oma = c(0, 0, 3, 0), xpd = NA)
  scatterplot3d(P[, 1], P[, 2], P[, 3], xlim = lim, ylim = lim, zlim = lim,
                color = cols, pch = 16, cex.symbols = 0.5, main = "PCA",
                xlab = "PC1", ylab = "PC2", zlab = "PC3",
                angle = ang, box = FALSE, grid = TRUE, col.grid = "grey88")
  legend("topright", legend = levels(df$class), col = palette[levels(df$class)],
         pch = 16, cex = 0.7, bty = "n")
  scatterplot3d(Ts[, 1], Ts[, 2], Ts[, 3], xlim = lim, ylim = lim, zlim = lim,
                color = cols, pch = 16, cex.symbols = 0.5, main = "t-SNE",
                xlab = "t-SNE1", ylab = "t-SNE2", zlab = "t-SNE3",
                angle = ang, box = FALSE, grid = TRUE, col.grid = "grey88")
  legend("topright", legend = levels(df$class), col = palette[levels(df$class)],
         pch = 16, cex = 0.7, bty = "n")
  scatterplot3d(U[, 1], U[, 2], U[, 3], xlim = lim, ylim = lim, zlim = lim,
                color = cols, pch = 16, cex.symbols = 0.5, main = "UMAP",
                xlab = "UMAP1", ylab = "UMAP2", zlab = "UMAP3",
                angle = ang, box = FALSE, grid = TRUE, col.grid = "grey88")
  legend("topright", legend = levels(df$class), col = palette[levels(df$class)],
         pch = 16, cex = 0.7, bty = "n")
  mtext(paste0(title, ": 3D projections (axes scaled to [-1, 1], same viewpoint)"),
        outer = TRUE, cex = 1.05, font = 2, line = 0.5)
  dev.off()
}

embeddings <- list()
overlap <- NULL
counts <- NULL

for (nm in names(datasets)) {
  ds <- datasets[[nm]]
  cat("\n==", nm, "==\n")
  cat("full class counts:\n"); print(table(ds$y))

  idx <- subsample(ds$y, n = 6000, seed = 42)
  Xs <- as.matrix(ds$X[idx, ])
  ys <- ds$y[idx]
  Xz <- scale(Xs)

  t0 <- Sys.time()
  pca <- prcomp(Xz)
  pcs <- pca$x[, 1:3]
  cat("PCA done\n")

  ts <- Rtsne(Xz, dims = 3, perplexity = 30, check_duplicates = FALSE,
              verbose = FALSE)$Y
  cat("t-SNE done:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

  um <- umap(Xz, n_components = 3L, n_neighbors = 15L, min_dist = 0.1,
             metric = "euclidean")$layout
  cat("UMAP done:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

  df <- data.frame(class = ys,
                   PC1 = pcs[, 1], PC2 = pcs[, 2], PC3 = pcs[, 3],
                   tSNE1 = ts[, 1], tSNE2 = ts[, 2], tSNE3 = ts[, 3],
                   UMAP1 = um[, 1], UMAP2 = um[, 2], UMAP3 = um[, 3])

  # 2D and 3D figures
  g2 <- make_2d(df, ds$palette, ds$title)
  ggsave(file.path(outdir, "figures", paste0("dr_2d_", nm, ".png")),
         plot = g2, width = 12, height = 4.8, units = "in", dpi = 150)
  make_3d(df, ds$palette, ds$title,
          file.path(outdir, "figures", paste0("dr_3d_", nm, ".png")))

  # class-overlap metrics in the standardised feature space (kNN purity)
  yc <- as.integer(ys)
  nn <- get.knnx(Xz, Xz, k = 11)$nn.index[, -1]
  same <- rowMeans(matrix(yc[nn], nrow = length(ys)) == yc)
  maj_code <- as.integer(names(which.max(table(yc))))
  majshare <- rowMeans(matrix(yc[nn], nrow = length(ys)) == maj_code)

  ov <- data.frame(dataset = nm, class = levels(ys),
                   n_sub = as.integer(table(ys)),
                   mean_same_class_neighbours = round(as.numeric(tapply(same, ys, mean)), 3),
                   mean_majority_neighbours = round(as.numeric(tapply(majshare, ys, mean)), 3),
                   row.names = NULL)
  overlap <- rbind(overlap, ov)

  cc <- data.frame(dataset = nm, class = levels(ds$y), n = as.integer(table(ds$y)))
  cc$share <- round(cc$n / sum(cc$n), 4)
  counts <- rbind(counts, cc)

  embeddings[[nm]] <- list(df = df, palette = ds$palette, title = ds$title,
                           var_explained = round(summary(pca)$importance[2, 1:3], 4))

  print(ov, row.names = FALSE)
}

saveRDS(embeddings, file.path(outdir, "dr_embeddings.rds"))
write.csv(overlap, file.path(outdir, "tables", "overlap_metrics.csv"), row.names = FALSE)
write.csv(counts, file.path(outdir, "tables", "class_counts.csv"), row.names = FALSE)

cat("\nwrote dr_results/ (embeddings, figures, tables)\n")
cat("PCA variance explained:\n")
for (nm in names(embeddings)) {
  cat(nm, ":", paste(round(100 * embeddings[[nm]]$var_explained, 1), "%", collapse = ", "), "\n")
}
