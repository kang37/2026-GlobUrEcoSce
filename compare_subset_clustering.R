options(stringsAsFactors = FALSE)

# Does dropping China change the archetype structure, and is the non-Chinese
# solution trustworthy on its own?
#
# Run prepare_cluster_inputs.R three times first:
#   (default)                             -> data_proc
#   SCENARIO_COUNTRY_EXCLUDE=CN           -> data_proc/no_CN
#   SCENARIO_COUNTRY_INCLUDE=CN           -> data_proc/CN_only
#
# Three questions are asked here:
#   1. Are the non-Chinese clusters actually archetypes, or just individual
#      papers? With 35 studies split into 5 groups this is the first thing
#      that has to be ruled out.
#   2. Do the same scenarios group together with and without China present?
#   3. Do the Chinese and non-Chinese solutions describe the same archetypes?

set.seed(20260728)
n_permutations <- as.integer(Sys.getenv("SCENARIO_N_PERM", unset = "9999"))

load_solution <- function(dir) {
  assignments <- read.csv(file.path(dir, "scenario_cluster_assignments.csv"))
  profiles <- read.csv(file.path(dir, "cluster_category_profiles.csv"))
  summary_tab <- read.csv(file.path(dir, "scenario_archetype_summary.csv"))
  assignments$study <- sub("-[^-]*$", "", assignments$scenario_id)
  assignments$country <- gsub("[^A-Za-z]", "", assignments$study)
  list(dir = dir, a = assignments, p = profiles, s = summary_tab)
}

full <- load_solution("data_proc")
no_cn <- load_solution(file.path("data_proc", "no_CN"))
cn <- load_solution(file.path("data_proc", "CN_only"))

# 1. Are clusters just studies? ----
# Share of within-cluster scenario pairs that come from the same paper. If
# clusters were archetypes this should sit near the level you would get by
# chance; if clusters are papers it goes far above it.
same_study_pair_rate <- function(cluster, study) {
  within <- 0
  same <- 0
  for (cl in unique(cluster)) {
    s <- study[cluster == cl]
    n <- length(s)
    if (n < 2) next
    within <- within + choose(n, 2)
    same <- same + sum(choose(table(s), 2))
  }
  same / within
}

pair_diagnostic <- function(sol, label) {
  observed <- same_study_pair_rate(sol$a$cluster, sol$a$study)
  null <- replicate(n_permutations,
                    same_study_pair_rate(sample(sol$a$cluster), sol$a$study))
  p <- (1 + sum(null >= observed - 1e-12)) / (n_permutations + 1)
  composition <- do.call(rbind, lapply(sort(unique(sol$a$cluster)), function(cl) {
    s <- sol$a$study[sol$a$cluster == cl]
    data.frame(
      cluster = cl, n = length(s), n_studies = length(unique(s)),
      n_countries = length(unique(sol$a$country[sol$a$cluster == cl])),
      largest_study_share = round(max(table(s)) / length(s), 3),
      largest_study = names(which.max(table(s)))
    )
  }))
  cat(sprintf("\n--- %s (n=%d scenarios, %d studies) ---\n",
              label, nrow(sol$a), length(unique(sol$a$study))))
  print(composition, row.names = FALSE)
  cat(sprintf("Same-study pair rate: observed %.3f vs %.3f expected by chance (p = %.4f)\n",
              observed, mean(null), p))
  cat(sprintf("Clusters where one paper supplies >40%% of the scenarios: %d of %d\n",
              sum(composition$largest_study_share > 0.40), nrow(composition)))
  invisible(list(observed = observed, expected = mean(null), p = p,
                 composition = composition))
}

cat("=== 1. Are the clusters archetypes, or just papers? ===")
full_pairs <- pair_diagnostic(full, "Full data")
no_cn_pairs <- pair_diagnostic(no_cn, "China removed")
cn_pairs <- pair_diagnostic(cn, "China only")

# 2. Does removing China regroup the remaining scenarios? ----
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  comb2 <- function(x) sum(x * (x - 1) / 2)
  index <- comb2(as.vector(tab))
  rows <- comb2(rowSums(tab))
  cols <- comb2(colSums(tab))
  expected <- rows * cols / comb2(n)
  (index - expected) / ((rows + cols) / 2 - expected)
}

cat("\n\n=== 2. Do the same scenarios stay together without China? ===\n")
shared <- merge(
  full$a[c("scenario_id", "cluster")],
  no_cn$a[c("scenario_id", "cluster")],
  by = "scenario_id", suffixes = c("_full", "_noCN")
)
cat(sprintf("Non-Chinese scenarios present in both solutions: %d\n", nrow(shared)))
cat("Cross-tabulation (rows = full-data archetype, cols = China-removed archetype):\n")
print(table(full = shared$cluster_full, no_CN = shared$cluster_noCN))
cat(sprintf("Adjusted Rand index = %.3f\n",
            adjusted_rand(shared$cluster_full, shared$cluster_noCN)))

shared_cn <- merge(
  full$a[c("scenario_id", "cluster")],
  cn$a[c("scenario_id", "cluster")],
  by = "scenario_id", suffixes = c("_full", "_CN")
)
cat(sprintf("For contrast, Chinese scenarios full vs China-only: ARI = %.3f (n = %d)\n",
            adjusted_rand(shared_cn$cluster_full, shared_cn$cluster_CN),
            nrow(shared_cn)))

# 3. Do the two literatures describe the same archetypes? ----
profile_matrix <- function(sol) {
  wide <- reshape(sol$p[c("cluster", "category", "cluster_prevalence")],
                  idvar = "category", timevar = "cluster", direction = "wide")
  rownames(wide) <- wide$category
  wide <- wide[paste0("T", seq_len(38)), -1, drop = FALSE]
  as.matrix(wide)
}

cn_profiles <- profile_matrix(cn)
no_cn_profiles <- profile_matrix(no_cn)

cat("\n\n=== 3. Do the Chinese and non-Chinese archetypes match? ===\n")
correlations <- cor(no_cn_profiles, cn_profiles)
dimnames(correlations) <- list(
  sprintf("noCN A%d: %s", no_cn$s$cluster, substr(no_cn$s$archetype_name, 1, 34)),
  sprintf("CN A%d", cn$s$cluster)
)
cat("Correlation of 38-category prevalence profiles:\n")
print(round(correlations, 2))

best <- data.frame(
  no_CN_archetype = sprintf("A%d %s", no_cn$s$cluster, no_cn$s$archetype_name),
  best_CN_match = sprintf("A%d %s", cn$s$cluster[apply(correlations, 1, which.max)],
                          cn$s$archetype_name[apply(correlations, 1, which.max)]),
  correlation = round(apply(correlations, 1, max), 3)
)
cat("\nBest match for each non-Chinese archetype:\n")
print(best, row.names = FALSE)

write.csv(best, file.path("data_proc", "subset_archetype_matching.csv"),
          row.names = FALSE)
write.csv(
  rbind(cbind(solution = "full", full_pairs$composition),
        cbind(solution = "no_CN", no_cn_pairs$composition),
        cbind(solution = "CN_only", cn_pairs$composition)),
  file.path("data_proc", "subset_cluster_composition.csv"), row.names = FALSE
)
cat("\nWrote 2 CSVs to data_proc/\n")
