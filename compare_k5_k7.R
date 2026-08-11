options(stringsAsFactors = FALSE)

# k = 7 as a supplementary solution: how it relates to the published k = 5,
# and whether each of the seven groups is a real archetype rather than the
# fingerprint of one or two papers.

analysis_dir <- "data_proc"
k7_dir <- file.path(analysis_dir, "k7")

k5 <- read.csv(file.path(analysis_dir, "scenario_cluster_assignments.csv"))
k7 <- read.csv(file.path(k7_dir, "scenario_cluster_assignments.csv"))
k5_names <- read.csv(file.path(analysis_dir, "scenario_archetype_summary.csv"))
k7_summary <- read.csv(file.path(k7_dir, "scenario_archetype_summary.csv"))
profiles <- read.csv(file.path(k7_dir, "cluster_category_profiles.csv"))

merged <- merge(k5[c("scenario_id", "cluster")], k7[c("scenario_id", "cluster")],
                by = "scenario_id", suffixes = c("_k5", "_k7"))
merged$study <- sub("-[^-]*$", "", merged$scenario_id)
merged$country <- gsub("[^A-Za-z]", "", merged$study)

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(z) sum(z * (z - 1) / 2)
  index <- comb2(as.vector(tab)); rows <- comb2(rowSums(tab))
  cols <- comb2(colSums(tab)); expected <- rows * cols / comb2(sum(tab))
  (index - expected) / ((rows + cols) / 2 - expected)
}

# Correspondence ----
cross <- table(k5 = merged$cluster_k5, k7 = merged$cluster_k7)
cat("=== k=5 (rows) x k=7 (columns) ===\n")
print(cross)
cat(sprintf("\nAdjusted Rand Index between the two solutions: %.3f\n",
            adjusted_rand_index(merged$cluster_k5, merged$cluster_k7)))

write.csv(cbind(k5_archetype = paste0("A", rownames(cross)),
                as.data.frame.matrix(cross)),
          file.path(analysis_dir, "k5_k7_correspondence.csv"), row.names = FALSE)

# Is each k=7 group carried by more than one paper? ----
composition <- do.call(rbind, lapply(sort(unique(merged$cluster_k7)), function(cl) {
  rows <- merged[merged$cluster_k7 == cl, ]
  studies <- table(rows$study)
  signature <- profiles[profiles$cluster == cl & profiles$signature, ]
  signature <- signature[order(-signature$cluster_prevalence), ]
  data.frame(
    cluster = cl,
    n_scenarios = nrow(rows),
    n_studies = length(studies),
    n_countries = length(unique(rows$country)),
    largest_study = names(which.max(studies)),
    largest_study_share = round(max(studies) / nrow(rows), 3),
    pct_chinese = round(mean(rows$country == "CN"), 3),
    medoid = k7_summary$medoid_scenario[k7_summary$cluster == cl],
    mean_silhouette = k7_summary$mean_silhouette[k7_summary$cluster == cl],
    signature_categories = paste(sprintf("%s(%.0f%%)", signature$category,
                                         100 * signature$cluster_prevalence),
                                 collapse = " "),
    max_lift = round(max(signature$lift), 1)
  )
}))

cat("\n=== Composition of each k=7 group ===\n")
print(composition[c("cluster", "n_scenarios", "n_studies", "n_countries",
                    "largest_study", "largest_study_share", "pct_chinese",
                    "medoid", "mean_silhouette", "max_lift")], row.names = FALSE)
cat("\nSignature categories:\n")
for (i in seq_len(nrow(composition))) {
  cat(sprintf("  B%d: %s\n", composition$cluster[i], composition$signature_categories[i]))
}

write.csv(composition, file.path(analysis_dir, "k7_archetype_composition.csv"),
          row.names = FALSE)

# Where each k=7 group comes from ----
cat("\n=== Origin of each k=7 group in the k=5 solution ===\n")
for (cl in sort(unique(merged$cluster_k7))) {
  origin <- sort(table(merged$cluster_k5[merged$cluster_k7 == cl]), decreasing = TRUE)
  cat(sprintf("  B%d <- %s\n", cl,
              paste(sprintf("A%s:%d", names(origin), origin), collapse = "  ")))
}

cat(sprintf("\nWrote 2 CSVs to %s/\n", analysis_dir))
