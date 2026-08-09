# ============================================================================
# 全部614个scenario：从读取数据、聚类分析到完整大图的全过程R代码
#
# 分析流程：
#   1. 读取614个scenario及T1-T38的0/1数据；
#   2. 填充article/case ID并检查数据；
#   3. 计算scenario之间的Jaccard距离；
#   4. 比较k=2-10的分组结果；
#   5. 使用PAM得到最终五类A1-A5；
#   6. 计算每类规模、代表情景和T1-T38组内出现比例；
#   7. 按统一规则确定典型categories；
#   8. 输出所有结果CSV；
#   9. 绰制k诊断图、五类正文小图和614行完整附件大图。
#
# 最终五类应复现：A1=113、A2=165、A3=76、A4=57、A5=203。
# ============================================================================

options(stringsAsFactors = FALSE)


# ==========================================================================
# 0. 用户设置
# ==========================================================================

# 在“情景原型分析”项目文件夹中运行本代码。
project_dir <- getwd()

# 推荐直接使用已经整理好的CSV，因此无需额外安装Excel读取包。
# 如果需要直接读取原始Excel，把input_mode改成"xlsx"。
input_mode <- "csv"   # 可选："csv"或"xlsx"

input_csv <- file.path(project_dir, "source_data_T1_T38.csv")

# Excel方式：请把文件复制到项目文件夹并修改文件名，
# 或把这里改成原始Excel的完整路径。
input_xlsx <- file.path(project_dir, "614  final code.xlsx")
input_sheet <- "code final code"

output_dir <- file.path(
  project_dir,
  "analysis_outputs",
  "complete_reproducible_pipeline_k5"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 最终论文采用的类别数。
# k诊断仍会比较2-10，但最终结果固定为经过统计和解释共同判断的k=5。
final_k <- 5L

# k稳定性检查：每次随机保留80%的scenario，再与完整结果比较。
# 100次可复现我们此前的诊断；正式运行可提高到500。
k_stability_repetitions <- 100L
set.seed(20260730)

# 完整附件大图尺寸。60英寸高度用于完整呈现614行。
complete_figure_width_in <- 12
complete_figure_height_in <- 60
complete_figure_png_dpi <- 300


# ==========================================================================
# 1. 检查R包
# ==========================================================================

# cluster通常随R安装，用于PAM和轮廓系数。
if (!requireNamespace("cluster", quietly = TRUE)) {
  stop(
    "缺少cluster包。请先运行：install.packages('cluster')"
  )
}

# 只有直接读取Excel时才需要readxl。
if (input_mode == "xlsx" && !requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "读取Excel需要readxl包。请先运行：install.packages('readxl')"
  )
}


# ==========================================================================
# 2. 读取原始数据
# ==========================================================================

if (input_mode == "csv") {
  if (!file.exists(input_csv)) {
    stop("没有找到CSV文件：", input_csv)
  }

  raw <- read.csv(
    input_csv,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )

} else if (input_mode == "xlsx") {
  if (!file.exists(input_xlsx)) {
    stop("没有找到Excel文件：", input_xlsx)
  }

  raw <- as.data.frame(
    readxl::read_excel(
      input_xlsx,
      sheet = input_sheet,
      col_names = TRUE,
      .name_repair = "minimal"
    )
  )

  # 原始Excel前两列分别为article/case ID和scenario ID。
  names(raw)[1:2] <- c("case_id", "scenario_id")

} else {
  stop("input_mode只能是'csv'或'xlsx'。")
}

# 兼容旧CSV中的列名。
if (!"case_id" %in% names(raw)) {
  names(raw)[1] <- "case_id"
}
if (!"scenario_id" %in% names(raw)) {
  names(raw)[2] <- "scenario_id"
}

category_names <- paste0("T", 1:38)

if (!all(c("case_id", "scenario_id", category_names) %in% names(raw))) {
  stop(
    "输入数据必须包含case_id、scenario_id和T1-T38。\n当前列名：",
    paste(names(raw), collapse = ", ")
  )
}

# 只保留本分析需要的40列。
raw <- raw[, c("case_id", "scenario_id", category_names)]


# ==========================================================================
# 3. 整理article ID和0/1矩阵
# ==========================================================================

# 原始表中article ID只写在每篇文章的第一条scenario上。
# 向下填充，使同一篇文章的每个scenario都有完整article ID。
case_id <- raw$case_id

for (i in seq_along(case_id)) {
  if ((is.na(case_id[i]) || trimws(as.character(case_id[i])) == "") && i > 1) {
    case_id[i] <- case_id[i - 1]
  }
}

case_id <- as.integer(case_id)
scenario_id <- trimws(as.character(raw$scenario_id))

# T1-T38全部转换成整数0/1。
X <- as.matrix(
  data.frame(
    lapply(raw[category_names], as.integer),
    check.names = FALSE
  )
)

rownames(X) <- scenario_id
colnames(X) <- category_names


# ==========================================================================
# 4. 数据质量检查
# ==========================================================================

if (nrow(X) != 614) {
  stop("数据应包含614个scenario，目前为：", nrow(X))
}

if (length(unique(case_id)) != 175) {
  stop("数据应包含175篇文章，目前为：", length(unique(case_id)))
}

if (ncol(X) != 38) {
  stop("数据应包含T1-T38共38项，目前为：", ncol(X))
}

if (anyNA(case_id) || anyNA(scenario_id) || anyNA(X)) {
  stop("case ID、scenario ID或T1-T38中存在缺失值。")
}

if (anyDuplicated(scenario_id)) {
  stop("scenario ID存在重复。")
}

if (!all(X %in% c(0L, 1L))) {
  stop("T1-T38必须全部为0或1。")
}

if (any(rowSums(X) == 0)) {
  stop(
    "以下scenario的T1-T38全部为0，请回查编码：",
    paste(scenario_id[rowSums(X) == 0], collapse = ", ")
  )
}

# 国家代码来自scenario ID开头的字母。
country_code <- toupper(sub("^([A-Za-z]+).*", "\\1", scenario_id))
geo_group <- ifelse(country_code == "CN", "China", "non-China")

clean_data <- data.frame(
  case_id = case_id,
  scenario_id = scenario_id,
  country_code = country_code,
  geo_group = geo_group,
  categories_per_scenario = rowSums(X),
  X,
  check.names = FALSE
)

write.csv(
  clean_data,
  file.path(output_dir, "01_clean_614_scenario_data.csv"),
  row.names = FALSE
)

data_audit <- data.frame(
  item = c(
    "number_of_articles",
    "number_of_scenarios",
    "number_of_countries",
    "number_of_categories",
    "China_articles",
    "China_scenarios",
    "non_China_articles",
    "non_China_scenarios"
  ),
  value = c(
    length(unique(case_id)),
    nrow(X),
    length(unique(country_code)),
    ncol(X),
    length(unique(case_id[geo_group == "China"])),
    sum(geo_group == "China"),
    length(unique(case_id[geo_group == "non-China"])),
    sum(geo_group == "non-China")
  )
)

write.csv(
  data_audit,
  file.path(output_dir, "02_data_audit.csv"),
  row.names = FALSE
)


# ==========================================================================
# 5. 计算Jaccard距离
# ==========================================================================

# 对0/1矩阵，stats::dist(method="binary")计算非对称二元距离，
# 即Jaccard distance：只重视共同出现的1，不把共同缺失的0算作相似。
jaccard_distance <- stats::dist(
  X,
  method = "binary"
)

# 保存距离矩阵会产生约614×614个值，文件较大但便于复核。
write.csv(
  as.matrix(jaccard_distance),
  file.path(output_dir, "03_Jaccard_distance_matrix.csv"),
  row.names = TRUE
)


# ==========================================================================
# 6. 辅助函数：Adjusted Rand Index
# ==========================================================================

# Adjusted Rand Index（ARI）比较两次分组的一致程度：
# 1表示完全一致，接近0表示接近随机一致。
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  total_pairs <- choose2(n)

  if (total_pairs == 0) return(NA_real_)

  sum_comb <- sum(choose2(tab))
  row_comb <- sum(choose2(rowSums(tab)))
  col_comb <- sum(choose2(colSums(tab)))
  expected <- row_comb * col_comb / total_pairs
  maximum <- (row_comb + col_comb) / 2

  if (maximum == expected) return(1)
  (sum_comb - expected) / (maximum - expected)
}


# ==========================================================================
# 7. 比较k=2-10
# ==========================================================================

k_values <- 2:10
full_fits <- vector("list", length(k_values))
names(full_fits) <- as.character(k_values)

k_diagnostics <- data.frame(
  k = k_values,
  average_silhouette = NA_real_,
  smallest_cluster_n = NA_integer_,
  smallest_cluster_percent = NA_real_,
  stability_ari_mean = NA_real_,
  stability_ari_sd = NA_real_
)

# 先在全部614个scenario上拟合k=2-10。
for (ii in seq_along(k_values)) {
  k <- k_values[ii]
  message("完整数据PAM：k = ", k)

  fit_k <- cluster::pam(
    jaccard_distance,
    k = k,
    diss = TRUE
  )

  full_fits[[as.character(k)]] <- fit_k
  cluster_sizes_k <- table(fit_k$clustering)

  k_diagnostics$average_silhouette[ii] <- fit_k$silinfo$avg.width
  k_diagnostics$smallest_cluster_n[ii] <- min(cluster_sizes_k)
  k_diagnostics$smallest_cluster_percent[ii] <-
    min(cluster_sizes_k) / nrow(X)
}

# 随机保留80%的scenario，重复检查不同k的分组是否一致。
# 这是对k选择的诊断，不是文章层bootstrap推断。
stability_replicates <- matrix(
  NA_real_,
  nrow = k_stability_repetitions,
  ncol = length(k_values),
  dimnames = list(NULL, paste0("k", k_values))
)

distance_matrix <- as.matrix(jaccard_distance)

for (b in seq_len(k_stability_repetitions)) {
  sampled_rows <- sort(
    sample(
      seq_len(nrow(X)),
      size = floor(0.80 * nrow(X)),
      replace = FALSE
    )
  )

  sampled_distance <- stats::as.dist(
    distance_matrix[sampled_rows, sampled_rows]
  )

  for (ii in seq_along(k_values)) {
    k <- k_values[ii]

    sampled_fit <- cluster::pam(
      sampled_distance,
      k = k,
      diss = TRUE
    )

    stability_replicates[b, ii] <- adjusted_rand_index(
      full_fits[[as.character(k)]]$clustering[sampled_rows],
      sampled_fit$clustering
    )
  }

  if (b %% 10 == 0) {
    message(
      "k稳定性重复：",
      b,
      "/",
      k_stability_repetitions
    )
  }
}

k_diagnostics$stability_ari_mean <- colMeans(
  stability_replicates,
  na.rm = TRUE
)

k_diagnostics$stability_ari_sd <- apply(
  stability_replicates,
  2,
  stats::sd,
  na.rm = TRUE
)

# 指标推荐值：排除最小类低于5%的解，再综合分离度和稳定性排名。
eligible <- k_diagnostics$smallest_cluster_percent >= 0.05
silhouette_rank <- rank(
  -k_diagnostics$average_silhouette,
  ties.method = "average"
)
stability_rank <- rank(
  -k_diagnostics$stability_ari_mean,
  ties.method = "average"
)

k_diagnostics$eligible_minimum_5_percent <- eligible
k_diagnostics$combined_rank <- silhouette_rank + stability_rank

eligible_rows <- which(eligible)
recommended_k_by_metrics <- k_diagnostics$k[
  eligible_rows[
    which.min(k_diagnostics$combined_rank[eligible_rows])
  ]
]

write.csv(
  k_diagnostics,
  file.path(output_dir, "04_k_2_to_10_diagnostics.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    replicate = seq_len(k_stability_repetitions),
    stability_replicates,
    check.names = FALSE
  ),
  file.path(output_dir, "05_k_stability_replicates.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste0("Metric-based recommended k: ", recommended_k_by_metrics),
    paste0("Final manuscript k: ", final_k),
    paste0("Stability repetitions: ", k_stability_repetitions),
    "Final k was judged using diagnostic performance, minimum cluster size,",
    "profile redundancy and substantive interpretability together."
  ),
  file.path(output_dir, "06_k_selection_note.txt")
)


# ==========================================================================
# 8. 绘制k诊断图
# ==========================================================================

draw_k_diagnostics <- function(filename, device = c("pdf", "png")) {
  device <- match.arg(device)

  if (device == "pdf") {
    grDevices::cairo_pdf(
      filename,
      width = 7.2,
      height = 3.8,
      family = "Arial"
    )
  } else {
    grDevices::png(
      filename,
      width = 7.2,
      height = 3.8,
      units = "in",
      res = 600,
      type = "cairo-png"
    )
  }

  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(
    mfrow = c(1, 2),
    mar = c(4.2, 4.2, 1.4, 0.7),
    family = "Arial",
    las = 1
  )

  graphics::plot(
    k_diagnostics$k,
    k_diagnostics$average_silhouette,
    type = "b",
    pch = 21,
    bg = "#2077BE",
    col = "#2077BE",
    lwd = 1.2,
    xlab = "Number of clusters (k)",
    ylab = "Average silhouette",
    xaxt = "n"
  )
  graphics::axis(1, at = k_values)
  graphics::abline(v = final_k, lty = 2, col = "#555555")
  graphics::mtext("a", side = 3, adj = 0, font = 2, line = 0.1)

  graphics::plot(
    k_diagnostics$k,
    k_diagnostics$stability_ari_mean,
    type = "b",
    pch = 21,
    bg = "#3093B1",
    col = "#3093B1",
    lwd = 1.2,
    ylim = c(0, 1),
    xlab = "Number of clusters (k)",
    ylab = "Mean stability (ARI)",
    xaxt = "n"
  )
  graphics::axis(1, at = k_values)
  graphics::arrows(
    k_diagnostics$k,
    pmax(0, k_diagnostics$stability_ari_mean - k_diagnostics$stability_ari_sd),
    k_diagnostics$k,
    pmin(1, k_diagnostics$stability_ari_mean + k_diagnostics$stability_ari_sd),
    angle = 90,
    code = 3,
    length = 0.035,
    col = "#68767F"
  )
  graphics::abline(v = final_k, lty = 2, col = "#555555")
  graphics::mtext("b", side = 3, adj = 0, font = 2, line = 0.1)
}

draw_k_diagnostics(
  file.path(output_dir, "Figure_1_k_diagnostics.pdf"),
  "pdf"
)
draw_k_diagnostics(
  file.path(output_dir, "Figure_1_k_diagnostics.png"),
  "png"
)


# ==========================================================================
# 9. 最终k=5的PAM聚类
# ==========================================================================

final_fit <- cluster::pam(
  jaccard_distance,
  k = final_k,
  diss = TRUE
)

cluster_id <- as.integer(final_fit$clustering)
medoid_row_numbers <- as.integer(final_fit$id.med)

# 每个scenario到其所属cluster代表scenario（medoid）的Jaccard距离。
# 数值越小，说明该scenario越接近本组最有代表性的0/1组合。
jaccard_matrix <- as.matrix(jaccard_distance)
distance_to_assigned_medoid <- vapply(
  seq_len(nrow(X)),
  function(i) {
    assigned_cluster <- cluster_id[i]
    jaccard_matrix[i, medoid_row_numbers[assigned_cluster]]
  },
  numeric(1)
)

archetype_names <- c(
  "1" = "Business-as-Usual and Trend-Based Development",
  "2" = "Economic Growth and Urban Expansion",
  "3" = "Ecosystem-Service-Oriented Conservation and Restoration",
  "4" = "Agricultural Sustainability and Food Security",
  "5" = "Urban Growth Containment for Integrated Landscape Protection and Restoration"
)

cluster_sizes <- data.frame(
  cluster = 1:5,
  n_scenarios = as.integer(table(factor(cluster_id, levels = 1:5)))
)
cluster_sizes$percent_scenarios <-
  cluster_sizes$n_scenarios / nrow(X)
cluster_sizes$archetype <- paste0("A", cluster_sizes$cluster)
cluster_sizes$archetype_name <- unname(
  archetype_names[as.character(cluster_sizes$cluster)]
)

expected_sizes <- c(113, 165, 76, 57, 203)

if (!all(cluster_sizes$n_scenarios == expected_sizes)) {
  warning(
    "当前五类规模为：",
    paste(cluster_sizes$n_scenarios, collapse = ", "),
    "；与此前最终结果113, 165, 76, 57, 203不一致。",
    "请检查输入文件和R/cluster版本。"
  )
}


# ==========================================================================
# 10. 计算每类T1-T38比例和典型categories
# ==========================================================================

category_labels <- c(
  T1 = "Historical trend / BAU / natural development with no additional policy intervention",
  T2 = "Economic and industrial growth-oriented development",
  T3 = "Urban expansion / construction-land growth / build-up growth / urban sprawl",
  T4 = "Urban growth control, containment and policy-regulated construction-land expansion",
  T5 = "Compact / infill / high-density urban development",
  T6 = "Monocentric compact residential development with controlled density",
  T7 = "Polycentric / multi-center urban development",
  T8 = "Land-sparing urban development / human-use concentration to protect land",
  T9 = "Land-sharing urban development / nature embedded within human-use space",
  T10 = "Policy, governance, spatial planning control and regulation",
  T11 = "Agricultural-land / cropland / cultivated-land protection and restoration",
  T12 = "Food security / sustainable food production",
  T13 = "Production-space priority for non-agricultural productive land allocation",
  T14 = "Living-space priority / residential and human-settlement development",
  T15 = "Ecological-space priority / ecological land protection and restoration",
  T16 = "Urban ecosystem-service-priority futures / ecosystem-service benefit maximization",
  T17 = "Urban service-oriented development / public-service and amenity-based urban sustainability",
  T18 = "Cultural ecosystem services: recreation, aesthetics, sense of place and social cohesion",
  T19 = "Reforestation / forest-tree-canopy expansion",
  T20 = "Coastal-marine ecosystem protection, restoration and sea-level-rise adaptation",
  T21 = "Water-resource constraint and inland blue-space protection",
  T22 = "Green infrastructure / nature-based solutions / sponge city / urban green space",
  T23 = "Low-emission / sustainable socio-climate pathway future",
  T24 = "Middle-path / moderate socio-climate pathway future",
  T25 = "High-emission / fossil-fuel-intensive socio-climate pathway future",
  T26 = "Low-carbon / renewable-energy transition / emissions mitigation",
  T27 = "Sustainable and balanced economic, urban/land and agricultural development",
  T28 = "Habitat protection / habitat quality and suitability enhancement",
  T29 = "Ecological network / corridor / connectivity enhancement",
  T30 = "Patch prioritization / structural habitat selection and conservation sequencing",
  T31 = "Habitat loss / ecological degradation / disturbance pressure and network disruption",
  T32 = "Integrated production-living-ecological space coordination and optimization",
  T33 = "Disaster-risk reduction / resilience / climate adaptation",
  T34 = "Social equity, livelihood, public health and human well-being",
  T35 = "Population change / demographic pressure / migration and socioeconomic demand pressure",
  T36 = "Technology, innovation, smart systems and efficiency improvement",
  T37 = "Transport, accessibility and infrastructure-led development",
  T38 = "Cost, land price, budget constraint and economic feasibility"
)

overall_prevalence <- colMeans(X)
profile_list <- vector("list", final_k)

for (cl in 1:final_k) {
  in_cluster <- cluster_id == cl
  within_cluster <- colMeans(X[in_cluster, , drop = FALSE])

  profile_list[[cl]] <- data.frame(
    cluster = cl,
    archetype = paste0("A", cl),
    category = category_names,
    category_label = unname(category_labels[category_names]),
    cluster_prevalence = within_cluster,
    overall_prevalence = overall_prevalence,
    prevalence_difference = within_cluster - overall_prevalence,
    prevalence_ratio = ifelse(
      overall_prevalence > 0,
      within_cluster / overall_prevalence,
      NA_real_
    ),
    medoid_present = as.integer(
      X[medoid_row_numbers[cl], category_names]
    ),
    stringsAsFactors = FALSE
  )
}

profiles <- do.call(rbind, profile_list)

# 典型category判定规则：三个条件必须同时满足。
# 1. 本组至少20%的scenario出现；
# 2. 比全部614个scenario的总体比例高至少10个百分点；
# 3. 至少是总体比例的1.35倍。
profiles$signature <- with(
  profiles,
  cluster_prevalence >= 0.20 &
    prevalence_difference >= 0.10 &
    prevalence_ratio >= 1.35
)


# ==========================================================================
# 11. 生成每个scenario的归属表和五类汇总表
# ==========================================================================

scenario_assignments <- data.frame(
  case_id = case_id,
  scenario_id = scenario_id,
  country_code = country_code,
  geo_group = geo_group,
  cluster = cluster_id,
  archetype = paste0("A", cluster_id),
  archetype_name = unname(
    archetype_names[as.character(cluster_id)]
  ),
  categories_per_scenario = rowSums(X),
  distance_to_medoid = distance_to_assigned_medoid,
  stringsAsFactors = FALSE
)

archetype_summary_list <- vector("list", final_k)

for (cl in 1:final_k) {
  one_profile <- profiles[profiles$cluster == cl, ]
  one_profile <- one_profile[
    order(
      -one_profile$signature,
      -one_profile$prevalence_difference,
      -one_profile$cluster_prevalence
    ),
  ]

  typical_categories <- one_profile$category[
    one_profile$signature
  ]

  if (length(typical_categories) == 0) {
    typical_categories <- head(one_profile$category, 5)
  }

  archetype_summary_list[[cl]] <- data.frame(
    archetype = paste0("A", cl),
    archetype_name = archetype_names[as.character(cl)],
    cluster = cl,
    n_scenarios = sum(cluster_id == cl),
    percent_scenarios = mean(cluster_id == cl),
    n_articles = length(unique(case_id[cluster_id == cl])),
    n_countries = length(unique(country_code[cluster_id == cl])),
    medoid_scenario = scenario_id[medoid_row_numbers[cl]],
    typical_categories = paste(
      typical_categories,
      collapse = ", "
    ),
    mean_categories_per_scenario = mean(
      rowSums(X[cluster_id == cl, , drop = FALSE])
    ),
    stringsAsFactors = FALSE
  )
}

archetype_summary <- do.call(
  rbind,
  archetype_summary_list
)


# ==========================================================================
# 12. 稀有category敏感性检查
# ==========================================================================

keep_categories <- overall_prevalence >= 0.01
removed_categories <- category_names[!keep_categories]

filtered_distance <- stats::dist(
  X[, keep_categories, drop = FALSE],
  method = "binary"
)

filtered_fit <- cluster::pam(
  filtered_distance,
  k = final_k,
  diss = TRUE
)

rare_filter_ari <- adjusted_rand_index(
  cluster_id,
  filtered_fit$clustering
)

analysis_quality <- data.frame(
  metric = c(
    "number_of_scenarios",
    "number_of_articles",
    "number_of_categories",
    "final_k",
    "average_silhouette",
    "k_stability_mean_ari",
    "k_stability_sd_ari",
    "rare_category_rule",
    "removed_rare_categories",
    "rare_filter_solution_ari"
  ),
  value = c(
    nrow(X),
    length(unique(case_id)),
    ncol(X),
    final_k,
    final_fit$silinfo$avg.width,
    k_diagnostics$stability_ari_mean[
      k_diagnostics$k == final_k
    ],
    k_diagnostics$stability_ari_sd[
      k_diagnostics$k == final_k
    ],
    "overall prevalence below 1%",
    paste(removed_categories, collapse = ", "),
    rare_filter_ari
  ),
  stringsAsFactors = FALSE
)


# ==========================================================================
# 13. 输出聚类结果CSV
# ==========================================================================

write.csv(
  scenario_assignments,
  file.path(output_dir, "07_scenario_cluster_assignments.csv"),
  row.names = FALSE
)

write.csv(
  profiles,
  file.path(output_dir, "08_cluster_category_profiles.csv"),
  row.names = FALSE
)

write.csv(
  archetype_summary,
  file.path(output_dir, "09_scenario_archetype_summary.csv"),
  row.names = FALSE
)

write.csv(
  cluster_sizes,
  file.path(output_dir, "10_cluster_sizes.csv"),
  row.names = FALSE
)

write.csv(
  analysis_quality,
  file.path(output_dir, "11_analysis_quality_summary.csv"),
  row.names = FALSE
)


# ==========================================================================
# 14. 五类正文小图
# ==========================================================================

cluster_colours <- c(
  "1" = "#2077BE",
  "2" = "#3093B1",
  "3" = "#27A0A6",
  "4" = "#6DB73F",
  "5" = "#3C9E39"
)

ivory_background <- "#F5F1E8"
grid_colour <- "#AEB4B2"
tree_colour <- "#303434"
outer_frame_colour <- "#8E9392"

mix_colour <- function(background, foreground, amount) {
  amount <- max(0, min(1, amount))
  bg <- grDevices::col2rgb(background) / 255
  fg <- grDevices::col2rgb(foreground) / 255
  value <- bg * (1 - amount) + fg * amount
  grDevices::rgb(value[1], value[2], value[3])
}

prevalence_matrix <- matrix(
  0,
  nrow = 38,
  ncol = 5,
  dimnames = list(category_names, paste0("A", 1:5))
)

signature_matrix <- matrix(
  FALSE,
  nrow = 38,
  ncol = 5,
  dimnames = dimnames(prevalence_matrix)
)

for (cl in 1:5) {
  one_profile <- profiles[profiles$cluster == cl, ]
  one_profile <- one_profile[
    match(category_names, one_profile$category),
  ]
  prevalence_matrix[, cl] <- one_profile$cluster_prevalence
  signature_matrix[, cl] <- one_profile$signature
}

draw_main_profile_figure <- function(
  filename,
  device = c("pdf", "png", "svg")
) {
  device <- match.arg(device)

  if (device == "pdf") {
    grDevices::cairo_pdf(
      filename,
      width = 6.6,
      height = 8.7,
      family = "Arial"
    )
  } else if (device == "png") {
    grDevices::png(
      filename,
      width = 6.6,
      height = 8.7,
      units = "in",
      res = 600,
      type = "cairo-png"
    )
  } else {
    grDevices::svg(
      filename,
      width = 6.6,
      height = 8.7,
      family = "Arial"
    )
  }

  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(
    mar = c(4.0, 4.7, 1.3, 5.9),
    xaxs = "i",
    yaxs = "i",
    family = "Arial"
  )

  graphics::plot(
    NA,
    xlim = c(0.5, 5.5),
    ylim = c(38.5, 0.5),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  for (tt in 1:38) {
    for (cl in 1:5) {
      p <- prevalence_matrix[tt, cl]
      fill_colour <- if (p == 0) {
        ivory_background
      } else {
        mix_colour(
          ivory_background,
          cluster_colours[as.character(cl)],
          0.16 + 0.84 * p
        )
      }

      graphics::rect(
        cl - 0.5,
        tt - 0.5,
        cl + 0.5,
        tt + 0.5,
        col = fill_colour,
        border = grid_colour,
        lwd = 0.95
      )

      if (signature_matrix[tt, cl]) {
        graphics::points(
          cl,
          tt,
          pch = 8,
          col = "#111111",
          cex = 0.98,
          lwd = 1.10
        )
      }
    }
  }

  graphics::axis(
    2,
    at = 1:38,
    labels = category_names,
    las = 1,
    tick = FALSE,
    line = -0.55,
    cex.axis = 0.72,
    font.axis = 2,
    gap.axis = -1
  )

  graphics::axis(
    1,
    at = 1:5,
    labels = paste0("A", 1:5),
    tick = FALSE,
    line = 0.20,
    cex.axis = 0.92,
    font.axis = 2
  )

  graphics::mtext(
    "Future-state category",
    side = 2,
    line = 2.85,
    cex = 0.82,
    font = 2
  )

  graphics::mtext(
    "Scenario archetype",
    side = 1,
    line = 2.20,
    cex = 0.86,
    font = 2
  )

  # 灰色比例尺，与38行主体等高。
  legend_top <- 0.5
  legend_bottom <- 38.5
  legend_left <- 5.78
  legend_right <- 6.06
  legend_grey <- "#4A4F4E"
  legend_steps <- seq(0, 1, length.out = 201)

  for (ii in seq_len(length(legend_steps) - 1)) {
    p0 <- legend_steps[ii]
    p1 <- legend_steps[ii + 1]
    y0 <- legend_bottom - p0 * (legend_bottom - legend_top)
    y1 <- legend_bottom - p1 * (legend_bottom - legend_top)

    graphics::rect(
      legend_left,
      y1,
      legend_right,
      y0,
      col = mix_colour(ivory_background, legend_grey, p1),
      border = NA,
      xpd = NA
    )
  }

  graphics::rect(
    legend_left,
    legend_top,
    legend_right,
    legend_bottom,
    border = grid_colour,
    lwd = 0.90,
    xpd = NA
  )

  legend_values <- c(1.00, 0.75, 0.50, 0.25, 0.00)
  legend_y <- legend_bottom -
    legend_values * (legend_bottom - legend_top)

  graphics::segments(
    legend_right,
    legend_y,
    legend_right + 0.08,
    legend_y,
    col = "#606564",
    lwd = 0.80,
    xpd = NA
  )

  graphics::text(
    legend_right + 0.12,
    legend_y,
    labels = sprintf("%.2f", legend_values),
    adj = c(0, 0.5),
    cex = 0.61,
    xpd = NA
  )

  graphics::text(
    legend_right + 0.63,
    (legend_top + legend_bottom) / 2,
    labels = "Within-cluster prevalence",
    srt = 90,
    cex = 0.61,
    font = 2,
    xpd = NA
  )
}

draw_main_profile_figure(
  file.path(output_dir, "Figure_2_five_archetype_profiles.pdf"),
  "pdf"
)
draw_main_profile_figure(
  file.path(output_dir, "Figure_2_five_archetype_profiles.png"),
  "png"
)
draw_main_profile_figure(
  file.path(output_dir, "Figure_2_five_archetype_profiles.svg"),
  "svg"
)


# ==========================================================================
# 15. 准备614行完整大图的数据
# ==========================================================================

complete_cluster_data <- lapply(1:5, function(cl) {
  scenario_ids_cl <- scenario_id[cluster_id == cl]
  matrix_cl <- X[scenario_ids_cl, , drop = FALSE]

  within_distance <- stats::dist(
    matrix_cl,
    method = "binary"
  )

  within_hclust <- stats::hclust(
    within_distance,
    method = "average"
  )

  ordered_ids <- rownames(matrix_cl)[within_hclust$order]
  one_profile <- profiles[profiles$cluster == cl, ]
  one_profile <- one_profile[
    match(category_names, one_profile$category),
  ]

  list(
    cluster = cl,
    n = nrow(matrix_cl),
    matrix = matrix_cl[ordered_ids, , drop = FALSE],
    hclust = within_hclust,
    ordered_ids = ordered_ids,
    signature_categories = one_profile$category[
      one_profile$signature
    ],
    prevalence = one_profile$cluster_prevalence
  )
})


# ==========================================================================
# 16. 绰制614行完整附件大图
# ==========================================================================

draw_complete_614_figure <- function(
  filename,
  device = c("pdf", "png", "svg"),
  show_scenario_ids = FALSE,
  width_in = complete_figure_width_in,
  height_in = complete_figure_height_in,
  dpi = complete_figure_png_dpi
) {
  device <- match.arg(device)

  if (device == "pdf") {
    grDevices::cairo_pdf(
      filename,
      width = width_in,
      height = height_in,
      family = "Arial"
    )
  } else if (device == "png") {
    grDevices::png(
      filename,
      width = width_in,
      height = height_in,
      units = "in",
      res = dpi,
      type = "cairo-png"
    )
  } else {
    grDevices::svg(
      filename,
      width = width_in,
      height = height_in,
      family = "Arial"
    )
  }

  on.exit(grDevices::dev.off(), add = TRUE)

  block_sizes <- vapply(
    complete_cluster_data,
    `[[`,
    numeric(1),
    "n"
  )

  # 五个数据块之间放置0.55行的空白。
  layout_matrix <- matrix(0, nrow = 9, ncol = 2)
  block_rows <- c(1, 3, 5, 7, 9)
  layout_matrix[block_rows, ] <- matrix(
    1:10,
    nrow = 5,
    ncol = 2,
    byrow = TRUE
  )

  layout_heights <- rep(0.55, 9)
  layout_heights[block_rows] <- block_sizes

  graphics::layout(
    layout_matrix,
    widths = c(3.0, 5.4),
    heights = layout_heights
  )

  for (i in seq_along(complete_cluster_data)) {
    z <- complete_cluster_data[[i]]
    cl <- as.character(z$cluster)
    nr <- nrow(z$matrix)

    circle_cex <- min(
      1.15,
      max(0.62, 0.75 * (height_in / nrow(X)) / 0.12)
    )

    # 左侧：每个archetype内部的完整Jaccard树枝关系。
    graphics::par(
      mar = c(
        if (i == 5) 3.9 else 0.05,
        0.25,
        if (i == 1) 1.30 else 0.05,
        0.03
      ),
      xaxs = "i",
      yaxs = "i"
    )

    graphics::plot(
      stats::as.dendrogram(z$hclust),
      horiz = TRUE,
      leaflab = "none",
      axes = FALSE,
      xlab = "",
      ylab = "",
      edgePar = list(col = tree_colour, lwd = 0.68)
    )

    # 完全相同的0/1结构没有可识别的内部先后关系，
    # 用平行末端和括号表示，不制造假的层级。
    profile_key <- apply(z$matrix, 1, paste0, collapse = "")
    duplicate_runs <- rle(profile_key)
    run_end <- cumsum(duplicate_runs$lengths)
    run_start <- run_end - duplicate_runs$lengths + 1
    tree_height <- max(z$hclust$height)
    bracket_length <- if (tree_height > 0) {
      tree_height * 0.075
    } else {
      0.075
    }

    for (rr in which(duplicate_runs$lengths > 1)) {
      yy <- run_start[rr]:run_end[rr]
      graphics::segments(
        0,
        yy,
        bracket_length,
        yy,
        col = tree_colour,
        lwd = 0.68
      )
      graphics::segments(
        bracket_length,
        min(yy),
        bracket_length,
        max(yy),
        col = tree_colour,
        lwd = 0.68
      )
    }

    # 右侧：T1-T38完整格子。
    graphics::par(
      mar = c(
        if (i == 5) 3.9 else 0.05,
        0.03,
        if (i == 1) 1.30 else 0.05,
        if (show_scenario_ids) 6.0 else 1.55
      ),
      xaxs = "i",
      yaxs = "i",
      family = "Arial"
    )

    graphics::plot(
      NA,
      xlim = c(0.5, 38.5),
      ylim = c(0.5, nr + 0.5),
      axes = FALSE,
      xlab = "",
      ylab = ""
    )

    # 米白色正方形底格。
    for (yy in seq_len(nr)) {
      graphics::rect(
        1:38 - 0.5,
        yy - 0.5,
        1:38 + 0.5,
        yy + 0.5,
        col = ivory_background,
        border = grid_colour,
        lwd = 0.56
      )
    }

    # value=1的格子按本组颜色和组内比例着色。
    for (xx in 1:38) {
      present_y <- which(z$matrix[, xx] == 1)
      if (length(present_y) == 0) next

      intensity <- 0.28 + 0.72 * z$prevalence[xx]

      graphics::rect(
        xx - 0.5,
        present_y - 0.5,
        xx + 0.5,
        present_y + 0.5,
        col = mix_colour(
          ivory_background,
          cluster_colours[cl],
          intensity
        ),
        border = grid_colour,
        lwd = 0.56
      )
    }

    # 典型category在具体scenario中出现时，加白色圆点。
    signature_x <- match(
      z$signature_categories,
      category_names
    )

    for (xx in signature_x[!is.na(signature_x)]) {
      present_y <- which(z$matrix[, xx] == 1)

      if (length(present_y) > 0) {
        graphics::points(
          rep(xx, length(present_y)),
          present_y,
          pch = 21,
          bg = "white",
          col = "white",
          cex = circle_cex
        )
      }
    }

    # 只给矩阵加浅灰色外框；树枝图不加外框。
    usr <- graphics::par("usr")
    graphics::segments(usr[1], usr[3], usr[1], usr[4],
                       col = outer_frame_colour, lwd = 0.95)
    graphics::segments(usr[2], usr[3], usr[2], usr[4],
                       col = outer_frame_colour, lwd = 0.95)
    graphics::segments(usr[1], usr[3], usr[2], usr[3],
                       col = outer_frame_colour, lwd = 0.95)
    graphics::segments(usr[1], usr[4], usr[2], usr[4],
                       col = outer_frame_colour, lwd = 0.95)

    graphics::axis(
      4,
      at = nr / 2,
      labels = sprintf("A%d (n=%d)", z$cluster, z$n),
      las = 1,
      tick = FALSE,
      line = 0.20,
      cex.axis = 0.68,
      font.axis = 2,
      col.axis = cluster_colours[cl]
    )

    if (show_scenario_ids) {
      graphics::axis(
        4,
        at = 1:nr,
        labels = z$ordered_ids,
        las = 1,
        tick = FALSE,
        line = 2.9,
        cex.axis = 0.18,
        col.axis = "#555555"
      )
    }

    if (i == 1) {
      graphics::mtext(
        "All 614 scenarios grouped into five scenario archetypes",
        side = 3,
        line = 0.35,
        adj = 0,
        font = 2,
        cex = 0.98
      )
    }

    if (i == 5) {
      graphics::axis(
        1,
        at = 1:38,
        labels = category_names,
        las = 2,
        tick = FALSE,
        line = -0.15,
        cex.axis = 0.66,
        font.axis = 2
      )
      graphics::mtext(
        "Future-state category",
        side = 1,
        line = 2.90,
        cex = 0.80,
        font = 2
      )
    }
  }
}


# ==========================================================================
# 17. 导出614行完整图和行顺序
# ==========================================================================

complete_prefix <- file.path(
  output_dir,
  "Figure_S1_all_614_scenarios_five_archetypes_complete"
)

draw_complete_614_figure(
  paste0(complete_prefix, ".pdf"),
  "pdf",
  show_scenario_ids = FALSE
)

draw_complete_614_figure(
  paste0(complete_prefix, ".png"),
  "png",
  show_scenario_ids = FALSE
)

draw_complete_614_figure(
  paste0(complete_prefix, ".svg"),
  "svg",
  show_scenario_ids = FALSE
)

draw_complete_614_figure(
  file.path(
    output_dir,
    "Figure_S1_all_614_scenarios_five_archetypes_with_IDs.pdf"
  ),
  "pdf",
  show_scenario_ids = TRUE,
  width_in = 14,
  height_in = 72
)

display_order <- do.call(
  rbind,
  lapply(complete_cluster_data, function(z) {
    data.frame(
      archetype = paste0("A", z$cluster),
      scenario_id = z$ordered_ids,
      stringsAsFactors = FALSE
    )
  })
)

display_order$display_order <- seq_len(nrow(display_order))

write.csv(
  display_order,
  file.path(output_dir, "12_complete_figure_display_order.csv"),
  row.names = FALSE
)


# ==========================================================================
# 18. 最终打印结果
# ==========================================================================

message("============================================================")
message("完整分析已经完成。")
message("输出文件夹：", normalizePath(output_dir))
message("指标推荐k：", recommended_k_by_metrics)
message("论文最终k：", final_k)
message(
  "五类规模：",
  paste(cluster_sizes$n_scenarios, collapse = ", ")
)
message("============================================================")

print(archetype_summary)
print(analysis_quality)
