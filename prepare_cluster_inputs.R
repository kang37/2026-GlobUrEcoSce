options(stringsAsFactors = FALSE)

# Stage 1 of the scenario-archetype analysis.
#
# Reads the raw literature coding workbook (614 scenarios x T1-T38 binary
# categories) and produces the four files consumed by
# plot_scenario_binary_matrix_dendrogram.R:
#
#   source_data_T1_T38.csv           scenario_id + T1..T38
#   scenario_cluster_assignments.csv scenario_id -> cluster
#   scenario_archetype_summary.csv   one row per archetype, incl. PAM medoid
#   cluster_category_profiles.csv    cluster x category prevalence + signature
#
# Clustering is PAM (partitioning around medoids) on Jaccard dissimilarity.
# Jaccard ignores joint absences, which is the correct choice here: two
# scenarios are not similar merely because neither mentions the same 30
# categories. PAM is used rather than hierarchical cutting because it yields
# an actual medoid scenario per archetype, i.e. a real, citable scenario that
# represents the group.

suppressPackageStartupMessages({
  library(readxl)
  library(cluster)
})

input_xlsx <- Sys.getenv("SCENARIO_INPUT_XLSX",
                         unset = file.path("data_raw", "lit_coding.xlsx"))
output_dir <- Sys.getenv("SCENARIO_ANALYSIS_DIR", unset = "data_proc")
number_of_clusters <- as.integer(Sys.getenv("SCENARIO_K", unset = "5"))

# A category counts as a "signature" of an archetype when it is both common
# within the archetype and disproportionately concentrated there.
signature_min_prevalence <- as.numeric(
  Sys.getenv("SCENARIO_SIGNATURE_MIN_PREVALENCE", unset = "0.25")
)
signature_min_lift <- as.numeric(
  Sys.getenv("SCENARIO_SIGNATURE_MIN_LIFT", unset = "1.5")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260728)

# Short labels for T1-T38, condensed from the long definitions in the coding
# protocol. Used only for human-readable output, never for the clustering.
category_labels <- c(
  T1  = "BAU / trend",              T2  = "Economic growth",
  T3  = "Urban expansion",          T4  = "Growth control",
  T5  = "Compact / infill",         T6  = "Monocentric compact",
  T7  = "Polycentric",              T8  = "Land sparing",
  T9  = "Land sharing",             T10 = "Planning / zoning control",
  T11 = "Cropland protection",      T12 = "Food security",
  T13 = "Production space",         T14 = "Living space",
  T15 = "Ecological space",         T16 = "Ecosystem services",
  T17 = "Urban services",           T18 = "Cultural ES",
  T19 = "Reforestation",            T20 = "Coastal / marine",
  T21 = "Water / blue space",       T22 = "Green infrastructure / NbS",
  T23 = "Low-emission pathway",     T24 = "Middle-path pathway",
  T25 = "High-emission pathway",    T26 = "Low-carbon transition",
  T27 = "Balanced development",     T28 = "Habitat quality",
  T29 = "Ecological network",       T30 = "Patch prioritization",
  T31 = "Ecological degradation",   T32 = "PLE-space integration",
  T33 = "Resilience / adaptation",  T34 = "Equity / well-being",
  T35 = "Demographic pressure",     T36 = "Technology / smart",
  T37 = "Transport / infrastructure", T38 = "Cost / feasibility"
)

# Data preparation ----
# Row 1 holds the T1-T38 headers; column 1 is a sparse study id (one value per
# study, blank for that study's subsequent scenarios); column 2 is the
# scenario id; columns 3-40 are the binary codings.
raw <- as.data.frame(
  read_excel(input_xlsx, sheet = 1, col_names = FALSE, .name_repair = "minimal")
)

category_names <- as.character(unlist(raw[1, 3:40]))
stopifnot(identical(category_names, paste0("T", seq_len(38))))

body <- raw[-1, , drop = FALSE]
scenario_id <- trimws(as.character(body[[2]]))

# Carry the study id down so every scenario knows which paper it came from.
study_id <- as.character(body[[1]])
filled <- study_id
for (i in seq_along(filled)) {
  if (is.na(filled[i]) || !nzchar(trimws(filled[i]))) {
    filled[i] <- if (i == 1) NA_character_ else filled[i - 1]
  }
}
study_id <- filled

X <- vapply(
  body[3:40],
  function(column) as.integer(as.character(column)),
  integer(nrow(body))
)
X[is.na(X)] <- 0L
colnames(X) <- category_names
rownames(X) <- scenario_id

stopifnot(
  !anyDuplicated(scenario_id),
  all(X %in% c(0L, 1L)),
  all(rowSums(X) > 0)
)

# Optional country subset ----
# The country is the alphabetic prefix of the scenario id (CN137-1 -> CN).
# Set SCENARIO_COUNTRY_EXCLUDE=CN to re-run the whole pipeline on the
# non-Chinese literature only, or SCENARIO_COUNTRY_INCLUDE=CN for the mirror.
# Categories are never dropped: all 38 columns are kept even when a subset
# leaves some of them empty, so figures and profiles stay comparable.
country <- gsub("[^A-Za-z]", "", sub("-[^-]*$", "", scenario_id))
include <- Sys.getenv("SCENARIO_COUNTRY_INCLUDE", unset = "")
exclude <- Sys.getenv("SCENARIO_COUNTRY_EXCLUDE", unset = "")
keep <- rep(TRUE, length(scenario_id))
if (nzchar(include)) {
  keep <- keep & country %in% trimws(strsplit(include, ",")[[1]])
}
if (nzchar(exclude)) {
  keep <- keep & !country %in% trimws(strsplit(exclude, ",")[[1]])
}
if (!all(keep)) {
  cat(sprintf("Country filter: keeping %d of %d scenarios (include='%s', exclude='%s').\n",
              sum(keep), length(keep), include, exclude))
  X <- X[keep, , drop = FALSE]
  scenario_id <- scenario_id[keep]
  study_id <- study_id[keep]
  country <- country[keep]
  empty <- colnames(X)[colSums(X) == 0]
  if (length(empty)) {
    cat(sprintf("  %d categories are empty in this subset and carry no information: %s\n",
                length(empty), paste(empty, collapse = ", ")))
  }
  stopifnot(all(rowSums(X) > 0))
}

cat(sprintf(
  "Loaded %d scenarios from %d studies, %d categories, %d present cells.\n",
  nrow(X), length(unique(study_id)), ncol(X), sum(X)
))

write.csv(
  data.frame(scenario_id = scenario_id, study_id = study_id,
             country = country, X, check.names = FALSE),
  file.path(output_dir, "source_data_T1_T38.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# Clustering ----
# dist(method = "binary") is the Jaccard distance for 0/1 data. The same
# metric is reused inside each block by the figure script.
jaccard <- dist(X, method = "binary")

# Report silhouette across a range of k so the choice of k is auditable
# rather than asserted.
diagnostics <- do.call(rbind, lapply(2:10, function(kk) {
  fit <- pam(jaccard, k = kk, diss = TRUE)
  data.frame(k = kk, average_silhouette_width = fit$silinfo$avg.width)
}))
write.csv(diagnostics, file.path(output_dir, "cluster_diagnostics_k2_k10.csv"),
          row.names = FALSE)
cat("Average silhouette width by k:\n")
print(diagnostics, row.names = FALSE)

pam_fit <- pam(jaccard, k = number_of_clusters, diss = TRUE)
clusters <- pam_fit$clustering
medoid_ids <- rownames(X)[pam_fit$id.med]
silhouette_width <- pam_fit$silinfo$widths[rownames(X), "sil_width"]

assignments <- data.frame(
  scenario_id = rownames(X),
  study_id = study_id,
  cluster = as.integer(clusters),
  silhouette_width = round(as.numeric(silhouette_width), 4),
  is_medoid = rownames(X) %in% medoid_ids,
  n_categories = as.integer(rowSums(X))
)
write.csv(assignments, file.path(output_dir, "scenario_cluster_assignments.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Category profiles ----
# Prevalence drives the fill intensity in the figure; the signature flag
# drives the white dots.
overall_prevalence <- colMeans(X)

profiles <- do.call(rbind, lapply(seq_len(number_of_clusters), function(cl) {
  member <- X[clusters == cl, , drop = FALSE]
  prevalence <- colMeans(member)
  lift <- ifelse(overall_prevalence > 0, prevalence / overall_prevalence, 0)
  signature <- prevalence >= signature_min_prevalence &
    lift >= signature_min_lift
  # Never let an archetype end up with no signature at all: fall back to the
  # three most over-represented categories that actually occur in it.
  if (!any(signature)) {
    candidates <- which(prevalence > 0)
    top <- candidates[order(lift[candidates], decreasing = TRUE)][1:3]
    signature[top[!is.na(top)]] <- TRUE
  }
  data.frame(
    cluster = cl,
    category = category_names,
    category_label = category_labels[category_names],
    n_present = as.integer(colSums(member)),
    cluster_size = nrow(member),
    cluster_prevalence = round(prevalence, 4),
    overall_prevalence = round(overall_prevalence, 4),
    lift = round(lift, 3),
    signature = signature
  )
}))
write.csv(profiles, file.path(output_dir, "cluster_category_profiles.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Archetype summary ----
# archetype_name is generated from the two most over-represented signature
# categories. It is a starting point for interpretation, not a finding --
# edit this column by hand once the archetypes have been read properly.
archetype_summary <- do.call(rbind, lapply(seq_len(number_of_clusters), function(cl) {
  profile_cl <- profiles[profiles$cluster == cl, ]
  signature_cl <- profile_cl[profile_cl$signature, ]
  signature_cl <- signature_cl[order(signature_cl$lift, decreasing = TRUE), ]
  data.frame(
    cluster = cl,
    archetype_name = paste(
      utils::head(signature_cl$category_label, 2), collapse = " + "
    ),
    n = sum(clusters == cl),
    share = round(mean(clusters == cl), 4),
    medoid_scenario = medoid_ids[cl],
    mean_silhouette = round(mean(assignments$silhouette_width[clusters == cl]), 4),
    n_signature = nrow(signature_cl),
    signature_categories = paste(signature_cl$category, collapse = ", ")
  )
}))
write.csv(archetype_summary, file.path(output_dir, "scenario_archetype_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nArchetypes:\n")
print(archetype_summary[c("cluster", "n", "medoid_scenario", "mean_silhouette",
                          "archetype_name")], row.names = FALSE)
cat(sprintf("\nWrote 5 files to %s/\n", output_dir))
