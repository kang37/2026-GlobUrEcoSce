# R code for scenario archetype clustering

# Load necessary libraries
# install.packages(c("readr", "dplyr", "ClustOfVar", "factoextra", "cluster", "dendextend"))
library(readr)
library(dplyr)
library(tidyr) # Added for the fill function
library(ClustOfVar) # For hierarchical clustering of variables (or observations)
library(factoextra) # For visualization of clusters
library(cluster)    # For silhouette analysis
library(dendextend) # For dendrogram manipulation

# --- 1. Load Data ---
# Read the CSV file, skipping the first row as it contains long descriptions
# and setting the second row as header.
# We need to manually identify the correct header row if read_csv doesn't handle it well.
# Based on inspection, the true headers (2, 3, ..., 41) are in the second row.
# The first row contains descriptive labels for the futures, we will try to use them as feature names.

# Read the first row separately to get feature descriptions
feature_descriptions <- read_csv("urban_eco_lit_coding.csv", col_names = FALSE, n_max = 1)

# Read the actual data, skipping the first row and using the second row as column names
# We'll use 'col_names = FALSE' for the main data read and then manually assign names.
data_raw <- read_csv("urban_eco_lit_coding.csv", skip = 1, col_names = FALSE)

# --- 2. Data Cleaning/Preparation ---

# Inspect data_raw to determine actual column indices for features and identifiers
# The first few columns seem to be identifiers:
# X1: Study/Group ID (from original Excel's merged cells, maybe from column 'a')
# X2: Scenario Name (from original Excel's 'Future1')
# X3: Scenario Code (from original Excel's 'code')
# X4: BAU status (from original Excel's 'BAU') - This needs to be converted if used as feature.
# X5 onwards are the 0/1 features. Column names are X5, X6, ... X44 based on the original structure.

# Let's try to map the feature descriptions to the correct columns.
# feature_descriptions[1, 5:44] corresponds to the columns X5:X44 in data_raw (if we skip 1 row)
# The actual column numbers for features (2 to 41) from the second row of the CSV are from X5 to X44
# Let's check the number of columns in data_raw, which is 106.
# The features (0/1 values) are from column X5 to X44.
# The original CSV had headers 'True', '2', '3', ... '41' for these columns, meaning 40 features.
# So, the columns of interest for clustering are X5 through X44 (inclusive).

feature_cols_start <- 5
feature_cols_end <- 44 # X44 is the column corresponding to feature '41' from the CSV's second row

# Extract feature matrix (0/1 values)
features_matrix <- data_raw %>%
  select(starts_with("X") & all_of(feature_cols_start:feature_cols_end)) %>%
  as.matrix()

# Convert to numeric (if not already)
features_matrix <- apply(features_matrix, 2, as.numeric)

# Replace NA values with 0 (assuming NA means absence for binary data)
features_matrix[is.na(features_matrix)] <- 0

# Assign descriptive names to features using the first row of the original Excel
# Ensure the number of descriptions matches the number of feature columns
feature_names <- as.character(feature_descriptions[1, feature_cols_start:feature_cols_end])
colnames(features_matrix) <- feature_names

# Prepare scenario metadata (Study ID, Scenario Name, Scenario Code, BAU status)
scenario_metadata <- data_raw %>%
  select(X1, X2, X3, X4) %>%
  rename(Study_ID = X1, Scenario_Name = X2, Scenario_Code = X3, BAU_Status = X4)

# Some Study_ID values are NA due to merged cells in original Excel.
# We need to fill these down for correct study grouping.
scenario_metadata <- scenario_metadata %>%
  fill(Study_ID, .direction = "down")

# --- Address Discrepancies and Clarifications ---
# 1. Number of "Common Futures": The PDF mentioned 59, but we have 40.
#    This script uses the 40 features available in the CSV.
# 2. Number of "Studies/Scenarios": PDF mentioned 37 studies. CSV has 117 scenarios.
#    The 'Study_ID' in scenario_metadata (X1) seems to group these.
#    For now, we will cluster *scenarios*. If clustering studies is required,
#    an aggregation step (e.g., sum of features per study or mode) would be needed.
#    I will add a placeholder for future study-level aggregation if desired.

# --- Optional: Aggregate to Study Level (if each row is a scenario within a study) ---
# If each row is a scenario and multiple scenarios belong to one study,
# we might want to aggregate features per study before clustering.
# For example, sum of '1's for each feature across scenarios within a study.
# study_features_matrix <- scenario_metadata %>%
#   bind_cols(as.data.frame(features_matrix)) %>%
#   group_by(Study_ID) %>%
#   summarise(across(all_of(feature_names), ~sum(cur_column())), .groups = 'drop') %>%
#   tibble::column_to_rownames(var = "Study_ID") %>%
#   as.matrix()
# # Then, normalize or convert to binary presence/absence at study level if sum > 0
# study_features_matrix[study_features_matrix > 0] <- 1

# For this initial implementation, we will cluster individual scenarios.
# If the user intends to cluster 'studies' as in the paper,
# they need to confirm how to aggregate scenarios within a study.

# --- 3. Clustering using Hierarchical Clustering (for scenarios) ---

# Use gower distance for mixed data types (though our features are binary, it's robust)
# For binary data, we can use binary distance from 'dist' function.
# The paper used hclustvar, which is for clustering variables.
# For clustering observations (scenarios), we use hclust on a distance matrix.

# Calculate Gower distance for binary data (can also use 'dist(method = "binary")')
distance_matrix <- dist(features_matrix, method = "binary")

# Perform hierarchical clustering
hc_result <- hclust(distance_matrix, method = "ward.D2") # Ward's method is common for archetypes

# --- 4. Determine Number of Clusters ---

# Elbow method and Silhouette method
# fviz_nbclust(features_matrix, hcut, method = "wss") # Elbow method
# fviz_nbclust(features_matrix, hcut, method = "silhouette") # Silhouette method

# For scree plot visualization of variance explained, typically used with PCA or ClustOfVar for variables.
# For hclust, we primarily rely on dendrogram and elbow/silhouette methods.

# Let's try to find an optimal number of clusters (e.g., 4 as in the paper)
# This part usually requires visual inspection and domain knowledge.
# For demonstration, let's assume 4 clusters, matching the paper's outcome.
k <- 4
clusters <- cutree(hc_result, k = k)

# Add cluster assignments to metadata
scenario_metadata$Cluster <- clusters

# --- 5. Visualization ---

# Dendrogram visualization
plot(hc_result, main = "Dendrogram of Scenario Archetypes",
     xlab = "Scenarios", ylab = "Distance")
rect.hclust(hc_result, k = k, border = 2:(k + 1))

# Heatmap (using heatmap.2 from gplots or a custom approach)
# A simple heatmap using base R
heatmap(features_matrix[hc_result$order, ],
        Colv = NA, # Do not reorder columns
        Rowv = as.dendrogram(hc_result), # Order rows by clustering
        col = c("white", "darkblue"), # Colors for 0 and 1
        labRow = scenario_metadata$Scenario_Name[hc_result$order],
        main = "Heatmap of Scenario Features by Cluster")

# Alternatively, using factoextra for heatmap (if features_matrix is converted to a dataframe)
# fviz_cluster(list(data = features_matrix, cluster = clusters),
#              geom = "point", ggtheme = theme_minimal())

# Calculate frequency of each feature per cluster (for characterization)
feature_frequency_per_cluster <- aggregate(features_matrix, by = list(cluster = clusters), FUN = sum)
# Convert sums to proportions if desired
# feature_frequency_per_cluster_prop <- feature_frequency_per_cluster
# feature_frequency_per_cluster_prop[,-1] <- sweep(feature_frequency_per_cluster_prop[,-1], 1, rowSums(feature_frequency_per_cluster_prop[,-1]), ")

print("Feature Frequencies per Cluster:")
print(feature_frequency_per_cluster)

# Save results
write_csv(scenario_metadata, "scenario_clusters.csv")

message("R script for clustering scenarios has been created as cluster_analysis.R")
message("Output: scenario_clusters.csv with cluster assignments, and feature frequencies per cluster.")
message("Please review the dendrogram and heatmap for visual interpretation of clusters.")

# --- Missing Characterization Data ---
message("\n--- Important Note ---")
message("The original paper characterized clusters using 11 social, economic, and biophysical variables (e.g., GDP per capita, GINI coefficient).")
message("This information was NOT found in 'urban_eco_lit_coding.csv'.")
message("To perform a similar characterization, please provide a separate CSV/Excel file containing these variables, with a column linking them to the 'Study_ID' in this dataset.")
message("Without this data, the characterization step (Fisher's exact tests, Kruskal-Wallis, Dunn tests) cannot be performed as described in the paper.")
