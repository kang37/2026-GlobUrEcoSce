options(stringsAsFactors = FALSE)

# How many archetypes? A multi-criterion diagnostic.
#
# Average silhouette width alone is a poor guide here: the data are 8.9% dense
# binary codes, 60% of scenario pairs share no category at all, so Jaccard
# distances pile up at 1 and silhouette stays flat around 0.13-0.19 for every k
# with no elbow to read. This script therefore reports five criteria that fail
# in different ways, plus a model-based one that does not use distances at all.
#
#   1. Average silhouette width          - compactness (weak here, reported anyway)
#   2. Share of negatively-silhouetted scenarios - how many sit in the wrong cluster
#   3. Subsample stability (80%, ARI)    - does the solution survive perturbation
#   4. Smallest cluster size             - practical usability
#   5. Interpretability                  - do all clusters have signature categories
#   6. Latent class BIC / AIC            - a probability model for binary indicators,
#                                          which is what this data actually is
#
# Criterion 6 is a Bernoulli mixture (latent class analysis) fitted by EM here
# rather than via poLCA, so the script has no extra dependency.

suppressPackageStartupMessages(library(cluster))

analysis_dir <- Sys.getenv("SCENARIO_ANALYSIS_DIR", unset = "data_proc")
k_values <- 2:12
n_stability <- as.integer(Sys.getenv("SCENARIO_N_STAB", unset = "100"))
n_lca_starts <- as.integer(Sys.getenv("SCENARIO_LCA_STARTS", unset = "20"))
set.seed(20260802)

source_data <- read.csv(file.path(analysis_dir, "source_data_T1_T38.csv"),
                        check.names = FALSE)
category_names <- grep("^T[0-9]+$", names(source_data), value = TRUE)
X <- as.matrix(source_data[category_names])
rownames(X) <- source_data$scenario_id
overall_prevalence <- colMeans(X)
jaccard <- dist(X, method = "binary")
distance_matrix <- as.matrix(jaccard)

# Project-standard signature rule (matches collaborator_full_pipeline_k5.R).
is_signature <- function(prevalence) {
  prevalence >= 0.20 &
    (prevalence - overall_prevalence) >= 0.10 &
    ifelse(overall_prevalence > 0, prevalence / overall_prevalence, 0) >= 1.35
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(z) sum(z * (z - 1) / 2)
  n <- sum(tab)
  index <- comb2(as.vector(tab))
  rows <- comb2(rowSums(tab))
  cols <- comb2(colSums(tab))
  expected <- rows * cols / comb2(n)
  maximum <- (rows + cols) / 2
  if (isTRUE(all.equal(maximum, expected))) return(1)
  (index - expected) / (maximum - expected)
}

# Latent class model ----
# Bernoulli mixture: P(x_i) = sum_c pi_c prod_j p_cj^x_ij (1-p_cj)^(1-x_ij).
# Probabilities are bounded away from 0/1 because sparse columns otherwise
# drive the likelihood to a degenerate boundary solution.
fit_latent_class <- function(X, n_classes, n_starts, max_iter = 500, tol = 1e-7) {
  n <- nrow(X); m <- ncol(X)
  best <- list(loglik = -Inf)
  for (s in seq_len(n_starts)) {
    pi_c <- rep(1 / n_classes, n_classes)
    p <- matrix(runif(n_classes * m, 0.05, 0.95), n_classes, m)
    previous <- -Inf
    for (it in seq_len(max_iter)) {
      # E-step in log space.
      log_lik_mat <- X %*% t(log(p)) + (1 - X) %*% t(log(1 - p))
      log_lik_mat <- sweep(log_lik_mat, 2, log(pi_c), "+")
      row_max <- apply(log_lik_mat, 1, max)
      weights <- exp(log_lik_mat - row_max)
      row_sum <- rowSums(weights)
      loglik <- sum(log(row_sum) + row_max)
      resp <- weights / row_sum
      # M-step.
      nk <- colSums(resp)
      pi_c <- nk / n
      p <- (t(resp) %*% X) / nk
      p <- pmin(pmax(p, 1e-3), 1 - 1e-3)
      if (abs(loglik - previous) < tol * abs(loglik)) break
      previous <- loglik
    }
    if (loglik > best$loglik) best <- list(loglik = loglik, n_par = (n_classes - 1) + n_classes * m)
  }
  best
}

# Diagnostics ----
cat(sprintf("Corpus: %d scenarios x %d categories, density %.1f%%.\n",
            nrow(X), ncol(X), 100 * mean(X)))
cat(sprintf("Running k = %s ...\n", paste(range(k_values), collapse = "-")))

fits <- lapply(k_values, function(k) pam(jaccard, k = k, diss = TRUE))
names(fits) <- as.character(k_values)

results <- do.call(rbind, lapply(seq_along(k_values), function(i) {
  k <- k_values[i]
  fit <- fits[[i]]
  clustering <- fit$clustering
  widths <- fit$silinfo$widths[rownames(X), "sil_width"]

  # Stability: refit on 80% subsamples, compare with the full solution.
  ari <- vapply(seq_len(n_stability), function(b) {
    rows <- sort(sample(seq_len(nrow(X)), floor(0.80 * nrow(X))))
    sub_fit <- pam(as.dist(distance_matrix[rows, rows]), k = k, diss = TRUE)
    adjusted_rand_index(clustering[rows], sub_fit$clustering)
  }, numeric(1))

  signature_sets <- lapply(seq_len(k), function(cl) {
    category_names[is_signature(colMeans(X[clustering == cl, , drop = FALSE]))]
  })

  data.frame(
    k = k,
    avg_silhouette = round(fit$silinfo$avg.width, 3),
    pct_negative_sil = round(100 * mean(widths < 0), 1),
    stability_ari = round(mean(ari), 3),
    stability_sd = round(sd(ari), 3),
    smallest_cluster = min(table(clustering)),
    clusters_without_signature = sum(lengths(signature_sets) == 0),
    distinct_signatures = length(unique(unlist(signature_sets)))
  )
}))

lca <- do.call(rbind, lapply(2:10, function(cc) {
  f <- fit_latent_class(X, cc, n_lca_starts)
  data.frame(k = cc,
             lca_BIC = round(-2 * f$loglik + f$n_par * log(nrow(X)), 1),
             lca_AIC = round(-2 * f$loglik + 2 * f$n_par, 1))
}))

results <- merge(results, lca, by = "k", all.x = TRUE)
write.csv(results, file.path(analysis_dir, "choose_k_diagnostics.csv"),
          row.names = FALSE)

cat("\n=== k diagnostics ===\n")
print(results, row.names = FALSE)

cat("\nBest by each criterion:\n")
cat(sprintf("  highest average silhouette : k = %d\n", results$k[which.max(results$avg_silhouette)]))
cat(sprintf("  fewest misplaced scenarios : k = %d\n", results$k[which.min(results$pct_negative_sil)]))
cat(sprintf("  most stable under resampling: k = %d\n", results$k[which.max(results$stability_ari)]))
cat(sprintf("  lowest latent-class BIC    : k = %d\n", lca$k[which.min(lca$lca_BIC)]))
cat(sprintf("  lowest latent-class AIC    : k = %d\n", lca$k[which.min(lca$lca_AIC)]))

# Figure ----
png(file.path(analysis_dir, "choose_k_diagnostics.png"),
    width = 11, height = 7, units = "in", res = 300, type = "cairo-png")
par(mfrow = c(2, 3), mar = c(4.3, 4.3, 3.0, 1.0))
plot_one <- function(y, main, ylab, highlight = NULL) {
  plot(results$k, y, type = "b", pch = 19, col = "#2077BE", lwd = 1.8,
       xlab = "number of clusters (k)", ylab = ylab, main = main, cex.main = 1.0)
  if (!is.null(highlight)) abline(v = highlight, col = "#C0392B", lty = 2)
  abline(v = 5, col = "#999999", lty = 3)
}
plot_one(results$avg_silhouette, "1. Average silhouette", "width")
plot_one(results$pct_negative_sil, "2. Scenarios in the wrong cluster", "% negative silhouette")
plot_one(results$stability_ari, "3. Stability under 80% resampling", "mean ARI")
plot_one(results$smallest_cluster, "4. Smallest cluster", "n scenarios")
plot_one(results$distinct_signatures, "5. Distinct signature categories", "count")
plot(lca$k, lca$lca_BIC, type = "b", pch = 19, col = "#2077BE", lwd = 1.8,
     xlab = "number of classes", ylab = "BIC", main = "6. Latent class BIC", cex.main = 1.0)
abline(v = lca$k[which.min(lca$lca_BIC)], col = "#C0392B", lty = 2)
abline(v = 5, col = "#999999", lty = 3)
invisible(dev.off())

cat(sprintf("\nWrote diagnostics table and figure to %s/\n", analysis_dir))
