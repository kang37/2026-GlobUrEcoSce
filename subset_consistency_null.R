options(stringsAsFactors = FALSE)

# Is the low agreement between the full-sample and China-removed clusterings a
# China effect, or just a small-sample effect?
#
# The observed number -- ARI = 0.188 between the full-sample labels of the 139
# non-Chinese scenarios and the labels they get when only those 139 are
# clustered -- has no reference point on its own. Any subsample of 614 down to
# ~139 makes clustering unstable and drags ARI down. So this script builds the
# null distribution: draw random subsamples of the SAME SIZE from the whole
# corpus, cluster each independently, compare each back to the full-sample
# solution exactly as the non-Chinese subset was compared, and see where 0.188
# falls.
#
# Sampling unit is the STUDY, never the scenario. Scenarios from one paper are
# near-duplicates of each other; splitting a paper across the in- and
# out-samples would make subsample clusterings artificially easy to reconcile
# with the full solution, inflating the null and manufacturing significance.
# Removing China also removed whole papers, so whole papers are what we draw.

suppressPackageStartupMessages(library(cluster))

analysis_dir <- Sys.getenv("SCENARIO_ANALYSIS_DIR", unset = "data_proc")
n_replicates <- as.integer(Sys.getenv("SCENARIO_N_REP", unset = "1000"))
k_values <- as.integer(strsplit(Sys.getenv("SCENARIO_K_RANGE", unset = "3,4,5,6,7"),
                                ",")[[1]])
final_k <- 5L
set.seed(20260802)

# Data ----
source_data <- read.csv(file.path(analysis_dir, "source_data_T1_T38.csv"),
                        check.names = FALSE)
full_assign <- read.csv(file.path(analysis_dir, "scenario_cluster_assignments.csv"))
no_cn_assign <- read.csv(file.path(analysis_dir, "no_CN",
                                   "scenario_cluster_assignments.csv"))

category_names <- grep("^T[0-9]+$", names(source_data), value = TRUE)
X <- as.matrix(source_data[match(full_assign$scenario_id, source_data$scenario_id),
                           category_names])
rownames(X) <- full_assign$scenario_id
full_labels <- setNames(full_assign$cluster, full_assign$scenario_id)

study <- sub("-[^-]*$", "", rownames(X))
country <- gsub("[^A-Za-z]", "", study)
study_scenarios <- split(rownames(X), study)
study_country <- vapply(study_scenarios, function(ids) country[match(ids[1], rownames(X))],
                        character(1))
study_size <- lengths(study_scenarios)

non_cn_ids <- no_cn_assign$scenario_id
target_n <- length(non_cn_ids)
target_studies <- length(unique(sub("-[^-]*$", "", non_cn_ids)))

cat(sprintf("Corpus: %d scenarios in %d studies.\n", nrow(X), length(study_scenarios)))
cat(sprintf("Target subsample size: %d scenarios / %d studies (the non-Chinese subset).\n",
            target_n, target_studies))

# Helpers ----
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(z) sum(z * (z - 1) / 2)
  n <- sum(tab)
  if (comb2(n) == 0) return(NA_real_)
  index <- comb2(as.vector(tab))
  rows <- comb2(rowSums(tab))
  cols <- comb2(colSums(tab))
  expected <- rows * cols / comb2(n)
  maximum <- (rows + cols) / 2
  if (isTRUE(all.equal(maximum, expected))) return(1)
  (index - expected) / (maximum - expected)
}

cluster_ids <- function(ids, k) {
  M <- X[ids, , drop = FALSE]
  as.integer(pam(dist(M, method = "binary"), k = k, diss = TRUE)$clustering)
}

# How much of each full-sample archetype survives in a subsample solution:
# for every full cluster, the best Jaccard overlap against any subsample cluster.
cluster_recovery <- function(reference, candidate) {
  vapply(sort(unique(reference)), function(r) {
    a <- which(reference == r)
    max(vapply(sort(unique(candidate)), function(cl) {
      b <- which(candidate == cl)
      length(intersect(a, b)) / length(union(a, b))
    }, numeric(1)))
  }, numeric(1))
}

# Draw whole studies until the scenario count first reaches target_n.
draw_by_scenarios <- function(pool, target) {
  shuffled <- sample(pool)
  reached <- which(cumsum(study_size[shuffled]) >= target)[1]
  unlist(study_scenarios[shuffled[seq_len(reached)]], use.names = FALSE)
}

# Draw a fixed number of whole studies.
draw_by_studies <- function(pool, n_studies) {
  unlist(study_scenarios[sample(pool, n_studies)], use.names = FALSE)
}

all_studies <- names(study_scenarios)
cn_studies <- all_studies[study_country == "CN"]

# Observed values ----
observed <- vapply(k_values, function(k) {
  own <- if (k == final_k) no_cn_assign$cluster else cluster_ids(non_cn_ids, k)
  reference <- if (k == final_k) {
    full_labels[non_cn_ids]
  } else {
    cluster_ids(rownames(X), k)[match(non_cn_ids, rownames(X))]
  }
  adjusted_rand_index(reference, own)
}, numeric(1))
names(observed) <- paste0("k", k_values)

cat("\nObserved ARI (non-Chinese subset vs full-sample solution):\n")
print(round(observed, 3))

# Null distributions ----
# For each replicate: draw a subsample, cluster it on its own, and compare with
# the labels those same scenarios carry in the full-sample solution.
run_null <- function(pool, draw, label, k) {
  full_k <- if (k == final_k) full_labels else {
    setNames(cluster_ids(rownames(X), k), rownames(X))
  }
  out <- vapply(seq_len(n_replicates), function(b) {
    ids <- draw(pool)
    own <- cluster_ids(ids, k)
    c(ari = adjusted_rand_index(full_k[ids], own),
      n = length(ids))
  }, numeric(2))
  data.frame(scheme = label, k = k, replicate = seq_len(n_replicates),
             ari = out["ari", ], n_scenarios = out["n", ])
}

cat(sprintf("\nRunning %d replicates per scheme for k = %s ...\n",
            n_replicates, paste(k_values, collapse = ", ")))

null_all <- do.call(rbind, lapply(k_values, function(k) {
  rbind(
    run_null(all_studies, function(p) draw_by_scenarios(p, target_n),
             "random_scenario_matched", k),
    run_null(all_studies, function(p) draw_by_studies(p, target_studies),
             "random_study_matched", k),
    run_null(cn_studies, function(p) draw_by_scenarios(p, target_n),
             "china_only_matched", k)
  )
}))
write.csv(null_all, file.path(analysis_dir, "subset_consistency_null_distribution.csv"),
          row.names = FALSE)

# Summary ----
summarise <- function(df, obs) {
  data.frame(
    scheme = df$scheme[1], k = df$k[1],
    mean_n_scenarios = round(mean(df$n_scenarios), 1),
    null_median = round(median(df$ari), 3),
    null_q025 = round(quantile(df$ari, 0.025), 3),
    null_q975 = round(quantile(df$ari, 0.975), 3),
    observed = round(obs, 3),
    # Left tail: is the non-Chinese subset LESS consistent than a random
    # subsample of the same size?
    p_left = round((1 + sum(df$ari <= obs + 1e-12)) / (nrow(df) + 1), 4)
  )
}

summary_tab <- do.call(rbind, lapply(split(null_all, list(null_all$scheme, null_all$k),
                                            drop = TRUE), function(df) {
  summarise(df, observed[paste0("k", df$k[1])])
}))
summary_tab <- summary_tab[order(summary_tab$k, summary_tab$scheme), ]
rownames(summary_tab) <- NULL
write.csv(summary_tab, file.path(analysis_dir, "subset_consistency_summary.csv"),
          row.names = FALSE)

cat("\n=== Null distribution vs observed ===\n")
print(summary_tab, row.names = FALSE)

# Per-archetype recovery at k = 5 ----
observed_recovery <- cluster_recovery(full_labels[non_cn_ids], no_cn_assign$cluster)
null_recovery <- replicate(n_replicates, {
  ids <- draw_by_scenarios(all_studies, target_n)
  cluster_recovery(full_labels[ids], cluster_ids(ids, final_k))
})
recovery_tab <- data.frame(
  archetype = paste0("A", sort(unique(full_labels))),
  n_non_chinese = as.integer(table(full_labels[non_cn_ids])[
    as.character(sort(unique(full_labels)))]),
  observed_recovery = round(observed_recovery, 3),
  null_median = round(apply(null_recovery, 1, median), 3),
  null_q025 = round(apply(null_recovery, 1, quantile, 0.025), 3),
  null_q975 = round(apply(null_recovery, 1, quantile, 0.975), 3)
)
recovery_tab$p_left <- round(
  (1 + rowSums(null_recovery <= observed_recovery + 1e-12)) / (n_replicates + 1), 4)
write.csv(recovery_tab, file.path(analysis_dir, "subset_consistency_by_archetype.csv"),
          row.names = FALSE)

cat("\n=== Per-archetype recovery at k = 5 (best Jaccard overlap) ===\n")
print(recovery_tab, row.names = FALSE)

# Figure ----
png(file.path(analysis_dir, "subset_consistency_null.png"),
    width = 11, height = 4.6, units = "in", res = 300, type = "cairo-png")
par(mfrow = c(1, 3), mar = c(4.4, 4.2, 3.4, 1.2))
schemes <- c("random_scenario_matched", "random_study_matched", "china_only_matched")
titles <- c("Random subsample\n(scenario-matched, n~139)",
            "Random subsample\n(35 studies)",
            "China-only subsample\n(scenario-matched, n~139)")
for (i in seq_along(schemes)) {
  d <- null_all$ari[null_all$scheme == schemes[i] & null_all$k == final_k]
  hist(d, breaks = 30, col = "#CFE2EF", border = "white", main = titles[i],
       xlab = "Adjusted Rand Index", ylab = if (i == 1) "Replicates" else "",
       xlim = range(c(d, observed["k5"])) + c(-0.02, 0.02), cex.main = 0.98)
  abline(v = observed["k5"], col = "#C0392B", lwd = 2.4)
  abline(v = median(d), col = "#2077BE", lwd = 1.6, lty = 2)
  legend("topright", bty = "n", cex = 0.78,
         legend = c(sprintf("observed = %.3f", observed["k5"]),
                    sprintf("null median = %.3f", median(d))),
         col = c("#C0392B", "#2077BE"), lwd = c(2.4, 1.6), lty = c(1, 2))
}
invisible(dev.off())

cat(sprintf("\nWrote 3 CSVs and 1 figure to %s/\n", analysis_dir))
