table_dir <- "tree_results/tables"
figure_dir <- "tree_results/figures"

read_confusion_matrix <- function(filename) {
  matrix_data <- read.csv(file.path(table_dir, filename), row.names = 1, check.names = FALSE)
  as.matrix(matrix_data)
}

plot_confusion_matrix <- function(confusion, title, filename) {
  n <- nrow(confusion)
  colors <- colorRampPalette(c("#EAF2F8", "#0072B2"))(100)
  scaled_values <- confusion / max(confusion)
  display_names <- function(names) {
    gsub("Prediabetes or diabetes", "Prediabetes / diabetes", names)
  }

  png(file.path(figure_dir, filename), width = 1200, height = 1000, res = 180)
  par(mar = c(6, 8, 3, 2), family = "sans")
  image(1:n, 1:n, t(scaled_values[n:1, , drop = FALSE]), col = colors,
        axes = FALSE, xlab = "", ylab = "",
        main = title)
  axis(1, at = 1:n, labels = display_names(colnames(confusion)), las = 1, cex.axis = 0.75)
  axis(2, at = 1:n, labels = rev(display_names(rownames(confusion))), las = 1, cex.axis = 0.75)
  box()

  for (row in 1:n) {
    for (column in 1:n) {
      value <- confusion[row, column]
      text(column, n - row + 1, format(value, big.mark = ","),
           col = ifelse(scaled_values[row, column] > 0.45, "white", "black"), cex = 1.05)
    }
  }
  mtext("Rows = predicted | Columns = actual | Correct predictions are on the diagonal",
        side = 3, line = 0.3, cex = 0.76)
  dev.off()
}

balanced <- read_confusion_matrix("balanced_binary_confusion_matrix.csv")
imbalanced <- read_confusion_matrix("imbalanced_binary_confusion_matrix.csv")
three_class <- read_confusion_matrix("three_class_confusion_matrix.csv")

plot_confusion_matrix(balanced, "Balanced binary tree: test confusion matrix",
                      "balanced_binary_confusion_matrix.png")
plot_confusion_matrix(imbalanced, "Imbalanced binary tree: test confusion matrix",
                      "imbalanced_binary_confusion_matrix.png")
plot_confusion_matrix(three_class, "Three-class tree: test confusion matrix",
                      "three_class_confusion_matrix.png")


writeLines(notes, "tree_results/CONFUSION_MATRIX_NOTES.md")

cat("Compact confusion-matrix plots created.\n")
