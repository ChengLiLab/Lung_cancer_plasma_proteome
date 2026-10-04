#!/usr/bin/env Rscript

# Random-forest classification of documented bone-metastasis status
# Inputs are the preprocessed/prefiltered protein matrices used for modeling:
# rows = proteins; columns = samples; first column = protein identifier.

suppressPackageStartupMessages({
  library(readxl)
  library(randomForest)
  library(caret)
  library(pROC)
  library(PRROC)
  library(dplyr)
  library(ggplot2)
  library(fastshap)
  library(shapviz)
})

set.seed(123)

data_dir <- file.path("data", "plasma")
out_dir <- file.path("results", "random_forest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bm_file <- file.path(data_dir, "BM_filtered.xlsx")
non_bm_file <- file.path(data_dir, "NO_BM_filtered.xlsx")

train_fraction <- 0.70
top_n <- 20
ntree <- 500
cv_folds <- 5
cv_repeats <- 30
shap_nsim <- 100

read_group <- function(path, label) {
  x <- read_excel(path, .name_repair = "minimal")
  feature <- make.unique(as.character(x[[1]]))
  m <- as.data.frame(x[, -1, drop = FALSE], check.names = FALSE)
  m[] <- lapply(m, function(z) as.numeric(as.character(z)))
  rownames(m) <- feature
  out <- as.data.frame(t(m), check.names = FALSE)
  out$Class <- label
  out
}

bm <- read_group(bm_file, "BM")
non_bm <- read_group(non_bm_file, "NO_BM")

features_bm <- setdiff(names(bm), "Class")
features_non_bm <- setdiff(names(non_bm), "Class")
if (!setequal(features_bm, features_non_bm)) {
  stop("BM and non-BM modeling files must contain the same protein features.")
}

features <- features_bm
non_bm <- non_bm[, c(features, "Class"), drop = FALSE]
bm <- bm[, c(features, "Class"), drop = FALSE]

dat <- bind_rows(non_bm, bm)
dat$Class <- factor(dat$Class, levels = c("NO_BM", "BM"))

if (anyNA(dat[, features, drop = FALSE])) {
  stop("Missing values are present. Use the preprocessed matrices used in the manuscript.")
}

class_weights <- function(y) {
  n <- table(y)
  w <- sum(n) / (length(n) * n)
  as.numeric(w[levels(y)]) |> setNames(levels(y))
}

best_youden_threshold <- function(y, prob_bm) {
  r <- roc(y, prob_bm, levels = c("NO_BM", "BM"), direction = "<", quiet = TRUE)
  as.numeric(coords(r, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1])
}

binary_metrics <- function(y, prob_bm, threshold) {
  pred <- factor(ifelse(prob_bm >= threshold, "BM", "NO_BM"), levels = levels(y))
  cm <- confusionMatrix(pred, y, positive = "BM")
  roc_obj <- roc(y, prob_bm, levels = c("NO_BM", "BM"), direction = "<", quiet = TRUE)
  pr_obj <- pr.curve(
    scores.class0 = prob_bm[y == "BM"],
    scores.class1 = prob_bm[y == "NO_BM"],
    curve = TRUE
  )
  list(
    pred = pred,
    confusion = cm,
    roc = roc_obj,
    pr = pr_obj,
    summary = data.frame(
      ROC_AUC = as.numeric(auc(roc_obj)),
      PR_AUC = pr_obj$auc.integral,
      Accuracy = unname(cm$overall["Accuracy"]),
      Sensitivity = unname(cm$byClass["Sensitivity"]),
      Specificity = unname(cm$byClass["Specificity"]),
      Threshold = threshold
    )
  )
}

train_idx <- createDataPartition(dat$Class, p = train_fraction, list = FALSE)
train <- dat[train_idx, , drop = FALSE]
test <- dat[-train_idx, , drop = FALSE]

rf_rank <- randomForest(
  x = train[, features, drop = FALSE],
  y = train$Class,
  ntree = ntree,
  importance = TRUE,
  classwt = class_weights(train$Class)
)

importance_mat <- importance(rf_rank, type = 1)
importance_df <- data.frame(
  Protein = rownames(importance_mat),
  MeanDecreaseAccuracy = importance_mat[, 1],
  stringsAsFactors = FALSE
) %>% arrange(desc(MeanDecreaseAccuracy))

top_features <- head(importance_df$Protein, min(top_n, nrow(importance_df)))
write.csv(importance_df, file.path(out_dir, "feature_importance_all.csv"), row.names = FALSE)
write.csv(data.frame(Protein = top_features), file.path(out_dir, "top20_features.csv"), row.names = FALSE)

rf_final <- randomForest(
  x = train[, top_features, drop = FALSE],
  y = train$Class,
  ntree = ntree,
  importance = TRUE,
  classwt = class_weights(train$Class),
  keep.inbag = TRUE
)

train_oob_prob <- rf_final$votes[, "BM"]
threshold <- best_youden_threshold(train$Class, train_oob_prob)

test_prob <- predict(rf_final, test[, top_features, drop = FALSE], type = "prob")[, "BM"]
test_eval <- binary_metrics(test$Class, test_prob, threshold)

write.csv(test_eval$summary, file.path(out_dir, "heldout_test_metrics.csv"), row.names = FALSE)
write.csv(
  data.frame(
    Sample = rownames(test),
    Observed = test$Class,
    Probability_BM = test_prob,
    Predicted = test_eval$pred
  ),
  file.path(out_dir, "heldout_test_predictions.csv"),
  row.names = FALSE
)
write.csv(as.data.frame.matrix(test_eval$confusion$table), file.path(out_dir, "heldout_confusion_matrix.csv"))

roc_df <- data.frame(
  Specificity = test_eval$roc$specificities,
  Sensitivity = test_eval$roc$sensitivities
)
p_roc <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = sprintf("Held-out ROC AUC = %.3f", as.numeric(auc(test_eval$roc)))
  ) +
  theme_classic()
ggsave(file.path(out_dir, "heldout_ROC.pdf"), p_roc, width = 5, height = 5)

pr_df <- as.data.frame(test_eval$pr$curve)
names(pr_df)[1:3] <- c("Recall", "Precision", "Threshold")
p_pr <- ggplot(pr_df, aes(Recall, Precision)) +
  geom_line(linewidth = 0.8) +
  labs(
    x = "Recall",
    y = "Precision",
    title = sprintf("Held-out PR AUC = %.3f", test_eval$pr$auc.integral)
  ) +
  theme_classic()
ggsave(file.path(out_dir, "heldout_PR.pdf"), p_pr, width = 5, height = 5)

# Repeated five-fold cross-validation.
# Top-20 RF feature ranking is repeated inside each resampling fold.
folds <- createMultiFolds(dat$Class, k = cv_folds, times = cv_repeats)

cv_predictions <- lapply(seq_along(folds), function(i) {
  tr_idx <- folds[[i]]
  va_idx <- setdiff(seq_len(nrow(dat)), tr_idx)
  tr <- dat[tr_idx, , drop = FALSE]
  va <- dat[va_idx, , drop = FALSE]

  rank_model <- randomForest(
    x = tr[, features, drop = FALSE],
    y = tr$Class,
    ntree = ntree,
    importance = TRUE,
    classwt = class_weights(tr$Class)
  )

  imp <- importance(rank_model, type = 1)
  fold_top <- rownames(imp)[order(imp[, 1], decreasing = TRUE)][seq_len(min(top_n, nrow(imp)))]

  fold_model <- randomForest(
    x = tr[, fold_top, drop = FALSE],
    y = tr$Class,
    ntree = ntree,
    classwt = class_weights(tr$Class)
  )

  fold_threshold <- best_youden_threshold(tr$Class, fold_model$votes[, "BM"])
  va_prob <- predict(fold_model, va[, fold_top, drop = FALSE], type = "prob")[, "BM"]

  data.frame(
    Fold = names(folds)[i],
    Sample = rownames(va),
    Observed = va$Class,
    Probability_BM = va_prob,
    Threshold = fold_threshold,
    Predicted = ifelse(va_prob >= fold_threshold, "BM", "NO_BM"),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(cv_predictions, file.path(out_dir, "repeated_5fold_CV_predictions.csv"), row.names = FALSE)

cv_y <- factor(cv_predictions$Observed, levels = c("NO_BM", "BM"))
cv_eval <- binary_metrics(cv_y, cv_predictions$Probability_BM, 0.5)
cv_summary <- cv_eval$summary
cv_summary$Note <- "AUCs are calculated from pooled out-of-fold probabilities across repeated CV."
write.csv(cv_summary, file.path(out_dir, "repeated_5fold_CV_summary.csv"), row.names = FALSE)

# SHAP interpretation for the final top-20 RF model.
x_all <- dat[, top_features, drop = FALSE]
shap_values <- fastshap::explain(
  rf_final,
  X = train[, top_features, drop = FALSE],
  newdata = x_all,
  pred_wrapper = function(object, newdata) {
    predict(object, newdata = newdata, type = "prob")[, "BM"]
  },
  nsim = shap_nsim,
  adjust = TRUE
)

write.csv(as.data.frame(shap_values), file.path(out_dir, "SHAP_values.csv"))

sv <- shapviz(shap_values, X = x_all)
pdf(file.path(out_dir, "SHAP_beeswarm.pdf"), width = 7, height = 6)
print(sv_importance(sv, kind = "beeswarm", max_display = top_n))
dev.off()

saveRDS(rf_final, file.path(out_dir, "random_forest_top20_model.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
