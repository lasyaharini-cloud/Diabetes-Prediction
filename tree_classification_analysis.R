# ST563 Diabetes Prediction Project
# Classification Tree Analysis
# Run this file from the project folder using:
# Rscript tree_classification_analysis.R

set.seed(42)

if (!requireNamespace("rpart", quietly = TRUE)) {
  stop("Package 'rpart' is required. Install it with install.packages('rpart').")
}

library(rpart)

# Folder locations for the input data and the files produced by this analysis.
data_dir <- "data"
output_dir <- "tree_results"
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read the two data sets needed for the tree models.
# The full data has 3 outcome groups (0, 1, and 2).
# The 50/50 file is the balanced binary data set.
data_012 <- read.csv(file.path(data_dir, "diabetes_012_health_indicators_BRFSS2021.csv"))
data_bal <- read.csv(file.path(data_dir, "diabetes_binary_5050split_health_indicators_BRFSS2021.csv"))

# These are the 21 health-related variables used to predict diabetes status.
predictors <- setdiff(names(data_012), "Diabetes_012")

# Keep the same numeric coding used in the logistic regression and KNN analysis.
# BMI is continuous. Most other variables are binary or ordered BRFSS variables.
prepare_predictors <- function(data) {
  data[, predictors, drop = FALSE]
}

make_task_data <- function(raw_data, response) {
  data.frame(response = response, prepare_predictors(raw_data), check.names = FALSE)
}

# Make an 80/20 train/test split. Stratifying means each split has about the
# same class proportions as the original data set.
stratified_split <- function(y, train_prop = 0.80, seed = 42) {
  set.seed(seed)
  groups <- split(seq_along(y), y)
  train_ids <- unlist(lapply(groups, function(ids) {
    sample(ids, size = floor(length(ids) * train_prop))
  }), use.names = FALSE)
  train_ids <- sort(train_ids)
  list(train = train_ids, test = setdiff(seq_along(y), train_ids))
}

# This prevents R from giving an error when a metric has a zero denominator.
safe_divide <- function(numerator, denominator) {
  ifelse(denominator == 0, NA_real_, numerator / denominator)
}

# Calculate AUC from predicted probabilities without needing an extra package.
auc_binary <- function(truth, probability, positive_class) {
  is_positive <- truth == positive_class
  n_positive <- sum(is_positive)
  n_negative <- sum(!is_positive)
  if (n_positive == 0 || n_negative == 0) return(NA_real_)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[is_positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

# Calculate the evaluation measures requested for binary classification.
binary_metrics <- function(truth, predicted, probability, positive_class) {
  negative_class <- setdiff(levels(truth), positive_class)
  tp <- sum(truth == positive_class & predicted == positive_class)
  fp <- sum(truth == negative_class & predicted == positive_class)
  fn <- sum(truth == positive_class & predicted == negative_class)
  tn <- sum(truth == negative_class & predicted == negative_class)
  precision <- safe_divide(tp, tp + fp)
  recall <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  f1 <- safe_divide(2 * precision * recall, precision + recall)
  list(
    summary = data.frame(
      accuracy = safe_divide(tp + tn, tp + fp + fn + tn),
      precision = precision,
      recall = recall,
      specificity = specificity,
      f1 = f1,
      auc = auc_binary(truth, probability, positive_class),
      row.names = NULL
    ),
    confusion = table(Predicted = predicted, Actual = truth)
  )
}

# For the 3-class model, calculate a set of metrics for each class and then
# calculate macro averages. Macro averages give every class equal importance.
multiclass_metrics <- function(truth, predicted, probabilities) {
  classes <- levels(truth)
  confusion <- table(Predicted = predicted, Actual = truth)
  per_class <- do.call(rbind, lapply(classes, function(class_name) {
    tp <- sum(truth == class_name & predicted == class_name)
    fp <- sum(truth != class_name & predicted == class_name)
    fn <- sum(truth == class_name & predicted != class_name)
    tn <- sum(truth != class_name & predicted != class_name)
    precision <- safe_divide(tp, tp + fp)
    recall <- safe_divide(tp, tp + fn)
    specificity <- safe_divide(tn, tn + fp)
    # If the model never predicts a class, its precision and F1 are recorded as 0.
    precision <- ifelse(is.na(precision), 0, precision)
    f1 <- safe_divide(2 * precision * recall, precision + recall)
    f1 <- ifelse(is.na(f1), 0, f1)
    data.frame(
      class = class_name,
      support = sum(truth == class_name),
      precision = precision,
      recall = recall,
      specificity = specificity,
      f1 = f1,
      auc_ovr = auc_binary(truth, probabilities[, class_name], class_name),
      row.names = NULL
    )
  }))
  list(
    summary = data.frame(
      accuracy = mean(predicted == truth),
      macro_precision = mean(per_class$precision, na.rm = TRUE),
      macro_recall = mean(per_class$recall, na.rm = TRUE),
      macro_specificity = mean(per_class$specificity, na.rm = TRUE),
      macro_f1 = mean(per_class$f1, na.rm = TRUE),
      macro_auc_ovr = mean(per_class$auc_ovr, na.rm = TRUE),
      row.names = NULL
    ),
    per_class = per_class,
    confusion = confusion
  )
}

# First grow a tree, then use 10-fold cross-validation to choose the pruning
# value. Pruning removes extra splits that do not improve prediction enough.
fit_pruned_tree <- function(train_data) {
  full_tree <- rpart(
    response ~ ., data = train_data, method = "class",
    control = rpart.control(cp = 0.0005, minsplit = 80, maxdepth = 6, xval = 10)
  )
  cp_table <- full_tree$cptable
  best_row <- which.min(cp_table[, "xerror"])
  pruned_tree <- prune(full_tree, cp = cp_table[best_row, "CP"])
  list(full = full_tree, pruned = pruned_tree, cp_table = cp_table, best_row = best_row)
}

# Save a picture of the final tree for the presentation.
save_tree_plot <- function(tree, path, title) {
  png(path, width = 2200, height = 1500, res = 180)
  par(mar = c(1, 1, 3, 1))
  plot(tree, uniform = TRUE, margin = 0.08, main = title)
  text(tree, use.n = TRUE, all = TRUE, cex = 0.62)
  dev.off()
}

# Save an ROC curve for each binary model.
save_roc_plot <- function(truth, probability, positive_class, path, title) {
  thresholds <- seq(0, 1, length.out = 201)
  tpr <- vapply(thresholds, function(t) {
    safe_divide(sum(truth == positive_class & probability >= t), sum(truth == positive_class))
  }, numeric(1))
  fpr <- vapply(thresholds, function(t) {
    safe_divide(sum(truth != positive_class & probability >= t), sum(truth != positive_class))
  }, numeric(1))
  png(path, width = 1400, height = 1200, res = 180)
  plot(fpr, tpr, type = "l", lwd = 3, col = "#0072B2", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False positive rate", ylab = "True positive rate", main = title)
  abline(0, 1, lty = 2, col = "gray45")
  legend("bottomright", legend = paste0("AUC = ", round(auc_binary(truth, probability, positive_class), 3)),
         bty = "n")
  dev.off()
}

# Fit and evaluate one binary classification tree.
run_binary_task <- function(task_name, task_data, positive_class) {
  split <- stratified_split(task_data$response)
  train_data <- task_data[split$train, , drop = FALSE]
  test_data <- task_data[split$test, , drop = FALSE]
  tree_fit <- fit_pruned_tree(train_data)
  # The predicted class is the class with the larger predicted probability.
  probabilities <- predict(tree_fit$pruned, test_data, type = "prob")
  predicted <- factor(colnames(probabilities)[max.col(probabilities)], levels = levels(test_data$response))
  metrics <- binary_metrics(test_data$response, predicted, probabilities[, positive_class], positive_class)

  write.csv(metrics$summary, file.path(table_dir, paste0(task_name, "_metrics.csv")), row.names = FALSE)
  write.csv(as.data.frame.matrix(metrics$confusion), file.path(table_dir, paste0(task_name, "_confusion_matrix.csv")))
  write.csv(as.data.frame(tree_fit$cp_table), file.path(table_dir, paste0(task_name, "_cp_table.csv")), row.names = FALSE)
  write.csv(data.frame(variable = names(tree_fit$pruned$variable.importance),
                       importance = as.numeric(tree_fit$pruned$variable.importance)),
            file.path(table_dir, paste0(task_name, "_variable_importance.csv")), row.names = FALSE)
  save_tree_plot(tree_fit$pruned, file.path(figure_dir, paste0(task_name, "_tree.png")),
                 paste("Pruned classification tree:", task_name))
  save_roc_plot(test_data$response, probabilities[, positive_class], positive_class,
                file.path(figure_dir, paste0(task_name, "_roc.png")), paste("ROC:", task_name))

  list(name = task_name, train = train_data, test = test_data, fit = tree_fit, metrics = metrics)
}

# Fit and evaluate the full 3-class classification tree.
run_multiclass_task <- function(task_name, task_data) {
  split <- stratified_split(task_data$response)
  train_data <- task_data[split$train, , drop = FALSE]
  test_data <- task_data[split$test, , drop = FALSE]
  tree_fit <- fit_pruned_tree(train_data)
  probabilities <- predict(tree_fit$pruned, test_data, type = "prob")
  predicted <- factor(colnames(probabilities)[max.col(probabilities)], levels = levels(test_data$response))
  metrics <- multiclass_metrics(test_data$response, predicted, probabilities)

  write.csv(metrics$summary, file.path(table_dir, paste0(task_name, "_metrics.csv")), row.names = FALSE)
  write.csv(metrics$per_class, file.path(table_dir, paste0(task_name, "_per_class_metrics.csv")), row.names = FALSE)
  write.csv(as.data.frame.matrix(metrics$confusion), file.path(table_dir, paste0(task_name, "_confusion_matrix.csv")))
  write.csv(as.data.frame(tree_fit$cp_table), file.path(table_dir, paste0(task_name, "_cp_table.csv")), row.names = FALSE)
  write.csv(data.frame(variable = names(tree_fit$pruned$variable.importance),
                       importance = as.numeric(tree_fit$pruned$variable.importance)),
            file.path(table_dir, paste0(task_name, "_variable_importance.csv")), row.names = FALSE)
  save_tree_plot(tree_fit$pruned, file.path(figure_dir, paste0(task_name, "_tree.png")),
                 paste("Pruned classification tree:", task_name))

  list(name = task_name, train = train_data, test = test_data, fit = tree_fit, metrics = metrics)
}

# Create the three outcomes used in the project comparison.
# The imbalanced binary outcome combines prediabetes and diabetes as class 1.
balanced_binary <- make_task_data(
  data_bal,
  factor(ifelse(data_bal$Diabetes_binary == 1, "Diabetes", "No diabetes"),
         levels = c("No diabetes", "Diabetes"))
)
imbalanced_binary <- make_task_data(
  data_012,
  factor(ifelse(data_012$Diabetes_012 > 0, "Prediabetes or diabetes", "No diabetes"),
         levels = c("No diabetes", "Prediabetes or diabetes"))
)
three_class <- make_task_data(
  data_012,
  factor(data_012$Diabetes_012, levels = c(0, 1, 2),
         labels = c("No diabetes", "Prediabetes", "Diabetes"))
)

# Run all three tree models.
balanced_result <- run_binary_task("balanced_binary", balanced_binary, "Diabetes")
imbalanced_result <- run_binary_task("imbalanced_binary", imbalanced_binary, "Prediabetes or diabetes")
three_class_result <- run_multiclass_task("three_class", three_class)

# Put the main test-set results into small tables for easy comparison.
binary_comparison <- rbind(
  data.frame(dataset = "Balanced binary", balanced_result$metrics$summary),
  data.frame(dataset = "Imbalanced binary", imbalanced_result$metrics$summary)
)
three_class_comparison <- three_class_result$metrics$summary
three_class_comparison$dataset <- "Three-class imbalanced"
three_class_comparison <- three_class_comparison[, c("dataset", setdiff(names(three_class_comparison), "dataset"))]
write.csv(binary_comparison, file.path(table_dir, "binary_tree_comparison.csv"), row.names = FALSE)
write.csv(three_class_comparison, file.path(table_dir, "three_class_tree_summary.csv"), row.names = FALSE)

# Create short notes using the actual model results. These can be used when
# making slides or explaining the model to the group.
importance <- head(read.csv(file.path(table_dir, "imbalanced_binary_variable_importance.csv")), 5)
round_for_display <- function(data, digits = 3) {
  numeric_columns <- vapply(data, is.numeric, logical(1))
  data[numeric_columns] <- lapply(data[numeric_columns], round, digits = digits)
  data
}
summary_lines <- c(
  "# Classification Tree: Presentation Notes",
  "",
  "## Method",
  "- CART classification trees fit to the same 21 health indicators used by the logistic/KNN analysis.",
  "- BMI remains continuous; documented BRFSS binary/ordinal numeric codes are retained to match the logistic/KNN baseline.",
  "- Each dataset uses a stratified 80/20 split with seed 42.",
  "- The complexity parameter is chosen by minimum 10-fold cross-validation error and the tree is pruned before test evaluation.",
  "",
  "## Test-set results",
  "",
  "### Binary tasks",
  paste(capture.output(print(round_for_display(binary_comparison), row.names = FALSE)), collapse = "\n"),
  "",
  "### Three-class task",
  paste(capture.output(print(round_for_display(three_class_comparison), row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation",
  "- Report recall, precision, F1, and AUC alongside accuracy; accuracy is inflated when the no-diabetes class dominates.",
  "- The three-class per-class table is the key evidence for whether the model recognizes prediabetes rather than only the majority class.",
  "- Variable importance describes predictive splitting value, not causal effect.",
  "",
  "## Top five splitting variables: imbalanced binary tree",
  paste0("- ", importance$variable, ": ", round(importance$importance, 1))
)
writeLines(summary_lines, file.path(output_dir, "presentation_notes.md"))

cat("Classification-tree analysis complete. Results are in ", output_dir, "/\n", sep = "")
cat("Balanced binary test metrics:\n")
print(round(balanced_result$metrics$summary, 3))
cat("Imbalanced binary test metrics:\n")
print(round(imbalanced_result$metrics$summary, 3))
cat("Three-class test metrics:\n")
print(round(three_class_result$metrics$summary, 3))
