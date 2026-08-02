options(stringsAsFactors = FALSE)

# Is archetype membership related to the country a study is about?
#
# Design note that drives everything below: the 614 scenarios are nested in
# 175 studies, and each study covers exactly one country. Scenarios from the
# same paper share authors, model and framing, so they are not independent
# observations. A chi-square test on 614 rows would treat them as if they
# were and would badly overstate significance. Every p-value here therefore
# comes from permuting COUNTRY LABELS ACROSS STUDIES, which keeps each study's
# internal set of scenarios intact and makes the study the unit of inference.
#
# Second design note: the corpus is 77% Chinese (475/614 scenarios,
# 140/175 studies). No other country has more than 4 studies. So the only
# contrast with real statistical power is China vs everywhere else; the
# per-country breakdown is reported as description, not as inference.

analysis_dir <- Sys.getenv("SCENARIO_ANALYSIS_DIR", unset = "data_proc")
n_permutations <- as.integer(Sys.getenv("SCENARIO_N_PERM", unset = "9999"))
set.seed(20260728)

assignments <- read.csv(file.path(analysis_dir, "scenario_cluster_assignments.csv"))
summary_tab <- read.csv(file.path(analysis_dir, "scenario_archetype_summary.csv"))
X <- read.csv(file.path(analysis_dir, "source_data_T1_T38.csv"), check.names = FALSE)

study_key <- sub("-[^-]*$", "", assignments$scenario_id)   # CN137-1 -> CN137
country <- gsub("[^A-Za-z]", "", study_key)                # CN137   -> CN
cluster <- assignments$cluster
k <- nrow(summary_tab)

# One country per study, so a study-level lookup is well defined.
study_country <- tapply(country, study_key, function(x) x[1])
stopifnot(all(tapply(country, study_key, function(x) length(unique(x))) == 1))

cat(sprintf("%d scenarios, %d studies, %d countries.\n",
            length(cluster), length(study_country), length(unique(country))))

# Permutation machinery ----
# Reassign country labels to whole studies, then rebuild the scenario-level
# table. Study sizes and their internal cluster mixes are preserved.
study_of_scenario <- match(study_key, names(study_country))

permuted_statistic <- function(labels, statistic) {
  shuffled <- sample(labels)
  statistic(shuffled[study_of_scenario], cluster)
}

chisq_statistic <- function(group, clus) {
  tab <- table(group, clus)
  if (nrow(tab) < 2) return(NA_real_)
  suppressWarnings(chisq.test(tab)$statistic)
}

permutation_p <- function(labels, statistic) {
  observed <- statistic(labels[study_of_scenario], cluster)
  null <- replicate(n_permutations, permuted_statistic(labels, statistic))
  # +1 correction: the observed arrangement is itself one of the possibilities.
  list(statistic = as.numeric(observed),
       p = (1 + sum(null >= observed - 1e-9)) / (n_permutations + 1))
}

cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab)$statistic)
  sqrt(as.numeric(chi) / (sum(tab) * (min(dim(tab)) - 1)))
}

# Test 1: all countries ----
full_test <- permutation_p(study_country, chisq_statistic)
full_tab <- table(country, cluster)

# Test 2: China vs rest ----
cn_labels <- ifelse(study_country == "CN", "CN", "Other")
cn_test <- permutation_p(cn_labels, chisq_statistic)
cn_tab <- table(ifelse(country == "CN", "CN", "Other"), cluster)
cn_study_tab <- table(cn_labels)

cat("\n--- Archetype composition, China vs rest (scenario counts) ---\n")
print(cn_tab)
cat("\nRow percentages:\n")
print(round(100 * prop.table(cn_tab, 1), 1))
cat(sprintf("\nStudies: CN = %d, Other = %d\n", cn_study_tab["CN"], cn_study_tab["Other"]))
cat(sprintf("Chi-square = %.2f, Cramer's V = %.3f\n",
            cn_test$statistic, cramers_v(cn_tab)))
cat(sprintf("Study-level permutation p = %.4f  (%d permutations)\n",
            cn_test$p, n_permutations))
naive <- suppressWarnings(chisq.test(cn_tab))
cat(sprintf("For contrast, the naive scenario-level p would be %.2g -- ",
            naive$p.value))
cat("do not report that one.\n")

cat("\nPearson residuals (which archetypes each group leans to):\n")
print(round(suppressWarnings(chisq.test(cn_tab))$residuals, 2))

cat(sprintf("\n--- All %d countries x %d archetypes ---\n", nrow(full_tab), k))
cat(sprintf("Chi-square = %.2f, study-level permutation p = %.4f\n",
            full_test$statistic, full_test$p))

# Per-country description ----
country_table <- as.data.frame.matrix(full_tab)
names(country_table) <- paste0("A", seq_len(k))
country_table$country <- rownames(country_table)
country_table$n_scenarios <- rowSums(full_tab)
country_table$n_studies <- as.integer(table(study_country)[rownames(full_tab)])
country_table <- country_table[order(-country_table$n_scenarios),
                               c("country", "n_studies", "n_scenarios",
                                 paste0("A", seq_len(k)))]
write.csv(country_table, file.path(analysis_dir, "country_cluster_crosstab.csv"),
          row.names = FALSE)

cat("\n--- Per-country archetype counts (description only) ---\n")
print(country_table, row.names = FALSE)

# Category-level mechanism ----
# If China does differ, which of the 38 categories carries the difference?
category_names <- grep("^T[0-9]+$", names(X), value = TRUE)
Xm <- as.matrix(X[match(assignments$scenario_id, X$scenario_id), category_names])
is_cn <- country == "CN"

category_results <- do.call(rbind, lapply(category_names, function(cat_name) {
  column <- Xm[, cat_name]
  observed <- mean(column[is_cn]) - mean(column[!is_cn])
  null <- replicate(n_permutations, {
    shuffled <- sample(cn_labels)[study_of_scenario]
    mean(column[shuffled == "CN"]) - mean(column[shuffled != "CN"])
  })
  data.frame(
    category = cat_name,
    cn_prevalence = round(mean(column[is_cn]), 4),
    other_prevalence = round(mean(column[!is_cn]), 4),
    difference = round(observed, 4),
    p_permutation = (1 + sum(abs(null) >= abs(observed) - 1e-12)) / (n_permutations + 1)
  )
}))
category_results$p_adjusted <- p.adjust(category_results$p_permutation, method = "BH")
category_results <- category_results[order(category_results$p_adjusted,
                                           -abs(category_results$difference)), ]
write.csv(category_results,
          file.path(analysis_dir, "country_category_differences.csv"),
          row.names = FALSE)

cat("\n--- Categories differing most between China and the rest ---\n")
cat("(BH-adjusted permutation p-values; positive difference = more Chinese)\n")
print(head(category_results, 12), row.names = FALSE)
cat(sprintf("\nCategories significant at BH < 0.05: %d of %d\n",
            sum(category_results$p_adjusted < 0.05), length(category_names)))

# Figure ----
cluster_colours <- c("#2077BE", "#3093B1", "#27A0A6", "#6DB73F", "#3C9E39")
plot_countries <- country_table[country_table$n_scenarios >= 5, ]
composition <- t(as.matrix(plot_countries[, paste0("A", seq_len(k))]))
composition <- prop.table(composition, 2)
colnames(composition) <- plot_countries$country

png(file.path(analysis_dir, "country_cluster_composition.png"),
    width = 11, height = 6, units = "in", res = 300, type = "cairo-png")
par(mar = c(5.6, 4.2, 3.2, 12.5), xpd = NA)
bar_x <- barplot(composition, col = cluster_colours[seq_len(k)],
                 border = "white", las = 1, ylab = "Share of scenarios",
                 cex.names = 0.85, ylim = c(0, 1),
                 main = "Archetype composition by country (>= 5 scenarios)")
# Sample size under each bar: the whole point is that only CN is well powered.
mtext(sprintf("%d st.", plot_countries$n_studies), side = 1, at = bar_x,
      line = 1.9, cex = 0.62, col = "#555555")
mtext(sprintf("%d sc.", plot_countries$n_scenarios), side = 1, at = bar_x,
      line = 2.7, cex = 0.62, col = "#555555")
legend("topright", inset = c(-0.30, 0), bty = "n", cex = 0.78,
       fill = cluster_colours[seq_len(k)],
       legend = sprintf("A%d  %s", summary_tab$cluster,
                        summary_tab$archetype_name))
mtext(sprintf(paste("Study-level permutation test, China vs rest:",
                    "p = %.3f (n.s.); Cramer's V = %.3f"),
              cn_test$p, cramers_v(cn_tab)),
      side = 1, line = 4.3, cex = 0.8)
invisible(dev.off())

cat(sprintf("\nWrote 2 CSVs and 1 figure to %s/\n", analysis_dir))
