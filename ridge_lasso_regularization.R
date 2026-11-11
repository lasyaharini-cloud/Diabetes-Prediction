set.seed(42)
data_dir <- "data"
out_dir <- "regularization_results"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(data_dir, "diabetes_012_health_indicators_BRFSS2021.csv"))
x_names <- setdiff(names(d), "Diabetes_012")
x_raw <- as.matrix(d[, x_names])
y <- as.integer(d$Diabetes_012 > 0)

stratified_split <- function(y, p = .80, seed = 42) {
  set.seed(seed)
  train <- sort(unlist(lapply(split(seq_along(y), y), function(i) sample(i, floor(length(i) * p)))))
  list(train = train, test = setdiff(seq_along(y), train))
}

standardize <- function(x_train, x_new) {
  mu <- colMeans(x_train)
  sd_value <- apply(x_train, 2, sd)
  sd_value[sd_value == 0] <- 1
  list(train = sweep(sweep(x_train, 2, mu, "-"), 2, sd_value, "/"),
       new = sweep(sweep(x_new, 2, mu, "-"), 2, sd_value, "/"))
}

split <- stratified_split(y)
x_train_raw <- x_raw[split$train, , drop = FALSE]
x_test_raw <- x_raw[split$test, , drop = FALSE]
y_train <- y[split$train]
y_test <- y[split$test]
scaled <- standardize(x_train_raw, x_test_raw)
x_train <- scaled$train
x_test <- scaled$new

soft_threshold <- function(z, lambda) sign(z) * pmax(abs(z) - lambda, 0)
predict_prob <- function(model, x) plogis(model$intercept + drop(x %*% model$beta))
log_loss <- function(y, p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

# This fits logistic ridge or lasso. The intercept is never penalized.
fit_penalized <- function(x, y, lambda, type, max_iter = 700, tol = 1e-7, start = NULL) {
  p <- ncol(x)
  if (is.null(start)) {
    b0 <- qlogis(mean(y)); b <- rep(0, p)
  } else {
    b0 <- start$intercept; b <- start$beta
  }
  largest_eigen <- max(eigen(crossprod(x) / nrow(x), symmetric = TRUE, only.values = TRUE)$values)
  step <- 1 / (0.25 * largest_eigen + ifelse(type == "ridge", lambda, 0) + 1e-8)
  for (iteration in seq_len(max_iter)) {
    residual <- predict_prob(list(intercept = b0, beta = b), x) - y
    b0_new <- b0 - step * mean(residual)
    gradient <- drop(crossprod(x, residual)) / nrow(x)
    if (type == "ridge") b_new <- b - step * (gradient + lambda * b)
    if (type == "lasso") b_new <- soft_threshold(b - step * gradient, step * lambda)
    if (max(abs(c(b0_new - b0, b_new - b))) < tol) {
      b0 <- b0_new; b <- b_new; break
    }
    b0 <- b0_new; b <- b_new
  }
  list(intercept = b0, beta = b, iterations = iteration)
}

auc <- function(y, p) {
  pos <- y == 1
  n_pos <- sum(pos); n_neg <- sum(!pos)
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

metrics <- function(y, p, threshold = .5) {
  pred <- as.integer(p >= threshold)
  tp <- sum(y == 1 & pred == 1); fp <- sum(y == 0 & pred == 1)
  fn <- sum(y == 1 & pred == 0); tn <- sum(y == 0 & pred == 0)
  precision <- ifelse(tp + fp == 0, NA, tp / (tp + fp))
  recall <- tp / (tp + fn)
  list(values = c(accuracy = (tp + tn) / length(y), precision = precision,
                  recall = recall, specificity = tn / (tn + fp),
                  f1 = 2 * precision * recall / (precision + recall), auc = auc(y, p)),
       confusion = matrix(c(tn, fn, fp, tp), 2, 2, byrow = TRUE,
                          dimnames = list(Predicted = c("No diabetes", "Prediabetes or diabetes"),
                                          Actual = c("No diabetes", "Prediabetes or diabetes"))))
}

# Select lambda using 5-fold CV on a stratified 30,000-row subset of training data.
set.seed(42)
positive <- sample(which(y_train == 1), round(30000 * mean(y_train == 1)))
negative <- sample(which(y_train == 0), 30000 - length(positive))
tune_id <- sort(c(positive, negative))
x_tune_raw <- x_train_raw[tune_id, , drop = FALSE]
y_tune <- y_train[tune_id]

make_folds <- function(y, k = 5) {
  folds <- vector("list", k)
  for (value in c(0, 1)) {
    ids <- sample(which(y == value))
    assignment <- rep(seq_len(k), length.out = length(ids))
    for (j in seq_len(k)) folds[[j]] <- c(folds[[j]], ids[assignment == j])
  }
  lapply(folds, sort)
}
folds <- make_folds(y_tune)
tune_scaled <- standardize(x_tune_raw, x_tune_raw)$train
lambda_max <- max(abs(drop(crossprod(tune_scaled, y_tune - mean(y_tune))) / nrow(tune_scaled)))
lasso_lambda <- lambda_max * 10^seq(0, -2.25, length.out = 18)
ridge_lambda <- lambda_max * 10^seq(1, -3, length.out = 18)

cv_penalty <- function(type, lambda_grid) {
  fold_results <- matrix(NA_real_, nrow = length(lambda_grid), ncol = length(folds))
  for (i in seq_along(lambda_grid)) {
    loss <- numeric(length(folds))
    for (j in seq_along(folds)) {
      valid <- folds[[j]]
      train <- setdiff(seq_along(y_tune), valid)
      fold_scaled <- standardize(x_tune_raw[train, , drop = FALSE], x_tune_raw[valid, , drop = FALSE])
      model <- fit_penalized(fold_scaled$train, y_tune[train], lambda_grid[i], type)
      loss[j] <- log_loss(y_tune[valid], predict_prob(model, fold_scaled$new))
    }
    fold_results[i, ] <- loss
    cat(type, "lambda", i, "of", length(lambda_grid), "complete\n")
  }
  data.frame(lambda = lambda_grid, log_loss = rowMeans(fold_results),
             standard_error = apply(fold_results, 1, sd) / sqrt(ncol(fold_results)))
}

ridge_cv <- cv_penalty("ridge", ridge_lambda)
lasso_cv <- cv_penalty("lasso", lasso_lambda)
best_ridge <- ridge_cv$lambda[which.min(ridge_cv$log_loss)]
lasso_min_position <- which.min(lasso_cv$log_loss)
lasso_min <- lasso_cv$lambda[lasso_min_position]

lasso_cutoff <- lasso_cv$log_loss[lasso_min_position] + lasso_cv$standard_error[lasso_min_position]
best_lasso <- max(lasso_cv$lambda[lasso_cv$log_loss <= lasso_cutoff])

ridge <- fit_penalized(x_train, y_train, best_ridge, "ridge")
lasso <- fit_penalized(x_train, y_train, best_lasso, "lasso")
baseline <- glm(y ~ ., data = data.frame(y = y_train, x_train_raw, check.names = FALSE), family = binomial())
p_baseline <- predict(baseline, newdata = data.frame(x_test_raw, check.names = FALSE), type = "response")
p_ridge <- predict_prob(ridge, x_test)
p_lasso <- predict_prob(lasso, x_test)
m_baseline <- metrics(y_test, p_baseline)
m_ridge <- metrics(y_test, p_ridge)
m_lasso <- metrics(y_test, p_lasso)
test_results <- rbind(data.frame(model = "Logistic regression", t(m_baseline$values)),
                      data.frame(model = "Ridge logistic regression", t(m_ridge$values)),
                      data.frame(model = "Lasso logistic regression", t(m_lasso$values)))
write.csv(test_results, file.path(tab_dir, "test_metrics_comparison.csv"), row.names = FALSE)
write.csv(ridge_cv, file.path(tab_dir, "ridge_cv_log_loss.csv"), row.names = FALSE)
write.csv(lasso_cv, file.path(tab_dir, "lasso_cv_log_loss.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(m_ridge$confusion), file.path(tab_dir, "ridge_confusion_matrix.csv"))
write.csv(as.data.frame.matrix(m_lasso$confusion), file.path(tab_dir, "lasso_confusion_matrix.csv"))

selected <- data.frame(variable = x_names, standardized_coefficient = lasso$beta)
selected <- selected[abs(selected$standardized_coefficient) > 1e-6, ]
selected <- selected[order(abs(selected$standardized_coefficient), decreasing = TRUE), ]
write.csv(selected, file.path(tab_dir, "lasso_selected_predictors.csv"), row.names = FALSE)

corr_vars <- c("BMI", "GenHlth", "MentHlth", "PhysHlth", "HighBP", "Age", "HighChol")
corr <- cor(d[, corr_vars], method = "spearman")
write.csv(corr, file.path(tab_dir, "selected_predictor_correlations.csv"))
png(file.path(fig_dir, "01_why_regularize_correlations.png"), width = 1800, height = 1350, res = 180)
par(mar = c(7, 7, 3, 2), family = "sans")
n <- ncol(corr)
image(1:n, 1:n, corr[n:1, ], col = colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(100),
      zlim = c(-.5, 1), axes = FALSE, main = "Selected predictor correlations: why test regularization")
axis(1, at = 1:n, labels = colnames(corr), las = 2)
axis(2, at = 1:n, labels = rev(rownames(corr)), las = 2)
for (r in 1:n) for (c in 1:n) text(c, n - r + 1, sprintf("%.2f", corr[r, c]),
                                     col = ifelse(abs(corr[r, c]) > .55, "white", "black"), cex = .85)
dev.off()

# Visual 2: lambda selection.
png(file.path(fig_dir, "02_cross_validated_penalty_selection.png"), width = 1800, height = 900, res = 180)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), family = "sans")
plot(log10(ridge_cv$lambda), ridge_cv$log_loss, type = "b", pch = 19, col = "#0072B2",
     xlab = "log10(lambda)", ylab = "5-fold CV log loss", main = "Ridge penalty selection")
abline(v = log10(best_ridge), col = "#D55E00", lwd = 2, lty = 2)
plot(log10(lasso_cv$lambda), lasso_cv$log_loss, type = "b", pch = 19, col = "#009E73",
     xlab = "log10(lambda)", ylab = "5-fold CV log loss", main = "Lasso penalty selection")
abline(v = log10(lasso_min), col = "#0072B2", lwd = 2, lty = 3)
abline(v = log10(best_lasso), col = "#D55E00", lwd = 2, lty = 2)
legend("topright", legend = c("Minimum loss", "1-SE lasso"), col = c("#0072B2", "#D55E00"),
       lty = c(3, 2), lwd = 2, bty = "n")
par(mfrow = c(1, 1)); dev.off()

# Visual 3: held-out test metrics for the baseline and regularized models.
metric_names <- c("auc", "recall", "f1")
metric_values <- as.matrix(test_results[, metric_names])
rownames(metric_values) <- test_results$model
png(file.path(fig_dir, "03_test_set_model_comparison.png"), width = 1800, height = 1100, res = 180)
par(mar = c(6, 5, 3, 1), family = "sans")
positions <- barplot(metric_values, beside = TRUE, col = c("#0072B2", "#E69F00", "#009E73"),
                    ylim = c(0, 1), names.arg = c("AUC", "Recall", "F1 score"), ylab = "Test-set value",
                    main = "Test-set comparison: logistic, ridge, and lasso")
text(positions, metric_values, labels = sprintf("%.3f", as.vector(metric_values)), pos = 3, cex = .75)
legend("topright", legend = rownames(metric_values), fill = c("#0072B2", "#E69F00", "#009E73"), bty = "n")
dev.off()

# Visual 4: coefficients retained by lasso.
display_selected <- head(selected, 12)
display_selected <- display_selected[order(display_selected$standardized_coefficient), ]
png(file.path(fig_dir, "04_lasso_selected_predictors.png"), width = 1800, height = 1200, res = 180)
par(mar = c(5, 10, 3, 1), family = "sans")
barplot(display_selected$standardized_coefficient, horiz = TRUE, names.arg = display_selected$variable,
        col = ifelse(display_selected$standardized_coefficient > 0, "#D55E00", "#0072B2"),
        border = NA, las = 1, xlab = "Standardized lasso coefficient",
        main = "Predictors retained by lasso at the selected penalty")
abline(v = 0, col = "gray40")
dev.off()

# Visual 5: coefficient paths show what happens as the penalty becomes stronger.
fit_coefficient_path <- function(type, lambda_values) {
  ordered_lambda <- sort(lambda_values, decreasing = TRUE)
  coefficients <- matrix(NA_real_, nrow = length(ordered_lambda), ncol = ncol(x_train),
                         dimnames = list(NULL, colnames(x_train)))
  starting_values <- NULL

  for (i in seq_along(ordered_lambda)) {
    path_fit <- fit_penalized(x_train, y_train, ordered_lambda[i], type, start = starting_values)
    coefficients[i, ] <- path_fit$beta
    starting_values <- path_fit
  }

  list(lambda = ordered_lambda, coefficients = coefficients)
}

ridge_path <- fit_coefficient_path("ridge", ridge_lambda)
lasso_path <- fit_coefficient_path("lasso", lasso_lambda)
path_variables <- head(selected$variable, 8)
path_colors <- grDevices::hcl.colors(length(path_variables), palette = "Dark 3")

draw_coefficient_path <- function(path_data, chosen_lambda, main_title) {
  order_index <- order(path_data$lambda)
  x_values <- log10(path_data$lambda[order_index])
  coefficient_values <- path_data$coefficients[order_index, path_variables, drop = FALSE]
  matplot(x_values, coefficient_values, type = "l", lty = 1, lwd = 2,
          col = path_colors, xlab = "log10(lambda): stronger penalty to the right",
          ylab = "Standardized coefficient", main = main_title)
  abline(h = 0, col = "gray55", lty = 2)
  abline(v = log10(chosen_lambda), col = "#D55E00", lwd = 2, lty = 2)
  legend("topright", legend = path_variables, col = path_colors, lty = 1,
         lwd = 2, cex = 0.72, bty = "n")
}

png(file.path(fig_dir, "05_ridge_lasso_coefficient_paths.png"), width = 2200, height = 1100, res = 180)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), family = "sans")
draw_coefficient_path(ridge_path, best_ridge, "Ridge: all coefficients shrink")
draw_coefficient_path(lasso_path, best_lasso, "Lasso: some coefficients become zero")
par(mfrow = c(1, 1))
dev.off()

show_metric <- function(name) paste0(test_results$model, " = ", sprintf("%.3f", test_results[[name]]))

cat("Ridge/lasso analysis complete. Results are in ", out_dir, "/\n", sep = "")
display_results <- test_results
display_results[vapply(display_results, is.numeric, logical(1))] <-
  lapply(display_results[vapply(display_results, is.numeric, logical(1))], round, digits = 4)
print(display_results)
cat("Lasso retained", nrow(selected), "of", ncol(x_train), "predictors.\n")
