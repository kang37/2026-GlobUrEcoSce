# Scenario archetype clustering at scenario-level.

# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(cluster)
library(dendextend)
library(pheatmap)

# Data preparation ----
# Load data. 
data_raw <- read_csv("urban_eco_lit_coding.csv", skip = 2, col_names = FALSE)

# Extract scenario metadata
scenario_metadata <- data_raw %>%
  select(X1, X2, X3) %>%
  rename(Study_ID = X1, Scenario_Name = X2, Scenario_Code = X3) %>%
  fill(Study_ID, .direction = "down") %>%
  mutate(
    Scenario_Code = coalesce(as.character(Scenario_Code), paste0("Scen_", 1:n())),
    Scenario_Name = ifelse(is.na(Scenario_Name), Scenario_Code, Scenario_Name)
  )

# Extract feature matrix (0/1 values)
feature_cols_start <- 5
feature_cols_end <- 44
features_matrix <- data_raw %>%
  select(all_of(feature_cols_start:feature_cols_end)) %>%
  as.matrix()
features_matrix[is.na(features_matrix)] <- 0
rownames(features_matrix) <- scenario_metadata$Scenario_Code

## Shorten Feature Names ----
short_names <- c(
  "No_Constraints", "Econ_Dev_No_Policy", "Econ_Dev_Max_Benefits", "Econ_Dev_High_Speed",
  "Urban_Expansion_Driven", "Socio_Econ_Pressure", "Urban_Containment", "Reindustrialization",
  "Eco_Protection_Oriented", "Max_Eco_Benefits", "Strict_Land_Control", "Eco_Restore_Pathway_Strict",
  "Eco_Restore_Pathway_Active", "Nature_Based_Solutions", "Ecosystem_Service_Priority", "Sustainable_Dev_Scenario",
  "Demographic_Stabilisation", "Urban_Service_Driven", "Food_Security_Priority", "Production_Space_Priority",
  "Living_Space_Priority", "Cultural_ES_Emphasis", "Open_Population_Policy", "Blue_Infrastructure",
  "Green_Energy_Tech", "Fossil_Fueled_Dev", "Mid_Road_Dev", "Low_Carbon_Transition",
  "LU_Driven_Low_Carbon", "Integrated_Low_Carbon", "Habitat_Security_Priority", "Eco_Connectivity",
  "Eco_Degradation_Mechanism", "Patch_Prioritization", "Compact_Residential_Dev", "Land_Sparing_Mode",
  "Land_Sharing_Mode", "Polycentric_Urban_Dev", "TOD", "Climate_Adaptation"
)
colnames(features_matrix) <- short_names

# Clustering (Scenario-level) ----
distance_matrix <- dist(features_matrix, method = "binary")
hc_result <- hclust(distance_matrix, method = "ward.D2")
k <- 5
clusters <- cutree(hc_result, k = k)

# Visualizations ----
# Save Dendrogram
png("scenario_dendrogram.png", width = 1600, height = 800)
plot(hc_result, main = "Dendrogram of Scenarios",
     xlab = "Scenarios", ylab = "Distance",
     labels = scenario_metadata$Scenario_Code[hc_result$order])
rect.hclust(hc_result, k = k, border = 2:(k + 1))
dev.off()

# Heatmap 1: Scenario x Feature
# Create an annotation data frame for rows
annotation_row <- data.frame(Cluster = factor(clusters))
rownames(annotation_row) <- rownames(features_matrix)

# Save the Scenario x Feature heatmap
png("scenario_feature_heatmap.png", width = 1500, height = 2000, res = 200)
pheatmap(
  features_matrix,
  color = c("#FFFFFF", "#336699"), # White for 0, Blue for 1
  main = "Heatmap of Features per Scenario",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = distance_matrix,
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  fontsize_row = 8,
  annotation_row = annotation_row,
  legend = FALSE
)
dev.off()

# Heatmap 2: Cluster x Feature Frequency
# Calculate feature frequency per cluster
cluster_freq_matrix <- features_matrix %>%
  as.data.frame() %>%
  mutate(Cluster = clusters) %>%
  group_by(Cluster) %>%
  summarise(
    across(all_of(short_names), ~mean(.x, na.rm = TRUE)), .groups = 'drop'
  ) %>%
  select(-Cluster) %>%
  as.matrix()

rownames(cluster_freq_matrix) <- paste("Cluster", 1:k)

# Save the Cluster x Feature Frequency heatmap
png("cluster_frequency_heatmap.png", width = 2500, height = 1300, res = 200)
pheatmap(
  cluster_freq_matrix,
  color = colorRampPalette(c("white", "steelblue"))(50),
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Feature Frequency by Cluster",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  fontsize = 10
)
dev.off()

