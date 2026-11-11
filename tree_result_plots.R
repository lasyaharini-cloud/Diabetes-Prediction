input_file <- "tree_results/tables/imbalanced_binary_variable_importance.csv"
output_file <- "tree_results/figures/imbalanced_binary_variable_importance.png"

importance <- read.csv(input_file)
importance <- importance[order(importance$importance, decreasing = TRUE), ]
top_variables <- head(importance, 8)

top_variables$relative_importance <- 100 * top_variables$importance / max(top_variables$importance)
top_variables <- top_variables[order(top_variables$relative_importance), ]

png(output_file, width = 1800, height = 1200, res = 180)
par(mar = c(5, 12, 3, 2), family = "sans")
bars <- barplot(
  top_variables$relative_importance,
  names.arg = top_variables$variable,
  horiz = TRUE,
  col = "#0072B2",
  border = NA,
  xlim = c(0, 110),
  las = 1,
  cex.names = 0.8,
  xlab = "Relative importance (top variable = 100)",
  main = "Most important variables in the imbalanced binary tree"
)
text(top_variables$relative_importance, bars,
     labels = paste0(round(top_variables$relative_importance), "%"), pos = 4, cex = 0.9)
dev.off()


writeLines(notes, "tree_results/TREE_EDA_AND_MODEL_NOTES.md")

cat("Tree variable-importance plot created.\n")
