options(stringsAsFactors = FALSE)

# Supplementary figure: all 614 scenarios (rows) x T1-T38 (columns).
# PAM defines the five archetype blocks. Within each block only, scenarios are
# ordered by average-linkage hierarchical clustering of Jaccard dissimilarity.
# White dots mark every present cell in a signature-category column -- they
# have nothing to do with the medoid scenario. Cell shading encodes the
# category's prevalence within the archetype, so it is constant down a column
# and carries no information about the individual scenario.
#
# Inputs are produced by prepare_cluster_inputs.R -- run that first.

analysis_dir <- Sys.getenv("SCENARIO_ANALYSIS_DIR", unset = "data_proc")
output_dir <- Sys.getenv("SCENARIO_MATRIX_FOLDER", unset = analysis_dir)
figure_prefix <- Sys.getenv(
  "SCENARIO_FIGURE_PREFIX",
  unset = "Figure_S1_614_scenario_binary_matrix"
)
input_file <- file.path(analysis_dir, "source_data_T1_T38.csv")
assignment_file <- file.path(analysis_dir, "scenario_cluster_assignments.csv")
summary_file <- file.path(analysis_dir, "scenario_archetype_summary.csv")
profile_file <- file.path(analysis_dir, "cluster_category_profiles.csv")

missing <- Filter(
  Negate(file.exists),
  c(input_file, assignment_file, summary_file, profile_file)
)
if (length(missing)) {
  stop("Missing input(s):\n  ", paste(missing, collapse = "\n  "),
       "\nRun prepare_cluster_inputs.R first.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(input_file, check.names = FALSE, fileEncoding = "UTF-8-BOM")
category_names <- grep("^T[0-9]+$", names(raw), value = TRUE)
X <- as.matrix(data.frame(lapply(raw[category_names], as.integer)))
rownames(X) <- raw$scenario_id

assignments <- read.csv(assignment_file, check.names = FALSE)
archetype_summary <- read.csv(summary_file, check.names = FALSE)
profiles <- read.csv(profile_file, check.names = FALSE)
number_of_clusters <- nrow(archetype_summary)
X <- X[match(assignments$scenario_id, rownames(X)), , drop = FALSE]
stopifnot(
  nrow(X) == nrow(assignments),
  length(category_names) == 38,
  setequal(rownames(X), assignments$scenario_id)
)
assignments <- assignments[match(rownames(X), assignments$scenario_id), ]

# The published k=5 figure uses the first five of these exact values, so they
# are kept verbatim for any k up to six. Beyond that a ramp is interpolated
# across the same blue-to-green range, which keeps seven or more blocks
# distinguishable instead of crowding three near-identical greens at the end.
default_cluster_colours <- c(
  "#2077BE", "#3093B1", "#27A0A6", "#6DB73F", "#3C9E39", "#2E782C"
)
cluster_colours <- if (number_of_clusters <= length(default_cluster_colours)) {
  setNames(default_cluster_colours[seq_len(number_of_clusters)],
           seq_len(number_of_clusters))
} else {
  setNames(
    colorRampPalette(c("#2077BE", "#27A0A6", "#6DB73F", "#2E782C"))(number_of_clusters),
    seq_len(number_of_clusters)
  )
}
zero_colour <- "#F5F1E8"
grid_colour <- "#AEB4B2"

mix_colour <- function(background, foreground, amount) {
  bg <- col2rgb(background) / 255
  fg <- col2rgb(foreground) / 255
  value <- bg * (1 - amount) + fg * amount
  rgb(value[1], value[2], value[3])
}

parse_categories <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

cluster_data <- lapply(seq_len(number_of_clusters), function(cl) {
  ids <- assignments$scenario_id[assignments$cluster == cl]
  mat <- X[ids, , drop = FALSE]
  jaccard_display <- dist(mat, method = "binary")
  hc <- hclust(jaccard_display, method = "average")
  ordered_ids <- rownames(mat)[hc$order]
  medoid_id <- archetype_summary$medoid_scenario[
    archetype_summary$cluster == cl
  ]
  profile_cl <- profiles[profiles$cluster == cl, ]
  profile_cl <- profile_cl[match(category_names, profile_cl$category), ]
  signature <- profile_cl$category[profile_cl$signature]
  list(
    cluster = cl,
    name = archetype_summary$archetype_name[
      archetype_summary$cluster == cl
    ],
    n = nrow(mat),
    mat = mat[ordered_ids, , drop = FALSE],
    hc = hc,
    ordered_ids = ordered_ids,
    medoid_id = medoid_id,
    signature = signature,
    prevalence = profile_cl$cluster_prevalence
  )
})

draw_binary_matrix <- function(
  file, device = c("pdf", "png"), show_row_labels = FALSE,
  width_in = 12, height_in = 60, dpi = 600
) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::cairo_pdf(file, width = width_in, height = height_in,
                        family = "Arial", onefile = TRUE)
  } else {
    grDevices::png(file, width = width_in, height = height_in, units = "in",
                   res = dpi, type = "cairo-png")
  }
  on.exit(dev.off(), add = TRUE)

  # With row labels the scenario IDs get a column of their own, between the
  # dendrogram and the matrix, so each ID visually links a dendrogram leaf to
  # its matrix row. Without them the figure keeps the plain two-column layout.
  panel_widths <- if (show_row_labels) c(2.4, 0.7, 5.4) else c(2.6, 5.4)
  n_panel <- length(panel_widths)

  block_sizes <- vapply(cluster_data, `[[`, numeric(1), "n")
  layout_matrix <- matrix(
    0, nrow = number_of_clusters * 2 - 1, ncol = n_panel
  )
  block_rows <- seq(1, nrow(layout_matrix), by = 2)
  layout_matrix[block_rows, ] <- matrix(
    seq_len(number_of_clusters * n_panel),
    nrow = number_of_clusters, ncol = n_panel, byrow = TRUE
  )
  layout_heights <- rep(0.5, nrow(layout_matrix))
  layout_heights[block_rows] <- block_sizes
  layout(
    layout_matrix,
    widths = panel_widths,
    heights = layout_heights
  )

  # Size the IDs to the row pitch: about 72% of a row's height, so they are as
  # large as they can be without neighbouring rows colliding.
  row_height_in <- height_in / nrow(X)
  id_cex <- min(1.10, max(0.30, 6.5 * row_height_in))

  for (i in seq_along(cluster_data)) {
    z <- cluster_data[[i]]
    cl <- as.character(z$cluster)
    circle_cex <- min(
      1.10,
      max(0.60, 0.75 * (height_in / nrow(X)) / 0.12)
    )

    # Left: within-archetype Jaccard dendrogram.
    par(
      mar = c(if (i == number_of_clusters) 3.7 else 0.05, 0.2,
              if (i == 1) 1.2 else 0.05, 0),
      xaxs = "i", yaxs = "i"
    )
    plot(
      as.dendrogram(z$hc), horiz = TRUE, leaflab = "none",
      axes = FALSE, xlab = "", ylab = "",
      edgePar = list(col = "#303434", lwd = 0.62)
    )
    # Zero-distance duplicate profiles have no estimable internal hierarchy.
    # Show them honestly as polytomies (parallel terminal ticks plus a bracket)
    # instead of inventing arbitrary binary relationships.
    profile_key <- apply(z$mat, 1, paste0, collapse = "")
    duplicate_runs <- rle(profile_key)
    run_end <- cumsum(duplicate_runs$lengths)
    run_start <- run_end - duplicate_runs$lengths + 1
    tree_height <- max(z$hc$height)
    bracket_length <- if (tree_height > 0) tree_height * 0.075 else 0.075
    for (rr in which(duplicate_runs$lengths > 1)) {
      yy <- run_start[rr]:run_end[rr]
      segments(0, yy, bracket_length, yy, col = "#303434", lwd = 0.62)
      segments(bracket_length, min(yy), bracket_length, max(yy),
               col = "#303434", lwd = 0.62)
    }

    # Middle: scenario IDs, sharing the matrix row grid exactly so each label
    # lines up with both its dendrogram leaf and its matrix row.
    if (show_row_labels) {
      par(
        mar = c(if (i == number_of_clusters) 3.7 else 0.05, 0,
                if (i == 1) 1.2 else 0.05, 0),
        xaxs = "i", yaxs = "i"
      )
      plot(
        NA, xlim = c(0, 1), ylim = c(0.5, nrow(z$mat) + 0.5),
        axes = FALSE, xlab = "", ylab = ""
      )
      text(
        0.5, seq_len(nrow(z$mat)), z$ordered_ids,
        cex = id_cex, col = "#555555", adj = c(0.5, 0.5)
      )
    }

    # Centre/right: binary T1-T38 matrix.
    par(
      mar = c(if (i == number_of_clusters) 3.7 else 0.05, 0,
              if (i == 1) 1.2 else 0.05,
              # The original 1.35 lines clipped the "A1 (n=113)" label off the
              # right edge: the label needs ~0.33in plus the ~0.16in axis
              # offset, i.e. at least ~3.7 lines. 5.0 leaves a clean margin.
              5.0),
      xaxs = "i", yaxs = "i"
    )
    nr <- nrow(z$mat)
    plot(
      NA, xlim = c(0.5, ncol(z$mat) + 0.5), ylim = c(0.5, nr + 0.5),
      axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"
    )
    # Every cell is an ivory square with a visible grey border.
    for (yy in seq_len(nr)) {
      rect(
        seq_len(ncol(z$mat)) - 0.5, yy - 0.5,
        seq_len(ncol(z$mat)) + 0.5, yy + 0.5,
        col = zero_colour, border = grid_colour, lwd = 0.52
      )
    }
    # A present category is shaded by its within-cluster prevalence.
    # The 0.28 floor keeps rare present categories visible in print.
    for (xx in seq_len(ncol(z$mat))) {
      present_y <- which(z$mat[, xx] == 1)
      if (!length(present_y)) next
      intensity <- 0.28 + 0.72 * z$prevalence[xx]
      rect(
        xx - 0.5, present_y - 0.5, xx + 0.5, present_y + 0.5,
        col = mix_colour(zero_colour, cluster_colours[cl], intensity),
        border = grid_colour, lwd = 0.52
      )
    }

    # White circles mark every present cell belonging to a signature category.
    dot_x <- match(z$signature, category_names)
    for (xx in dot_x[!is.na(dot_x)]) {
      present_y <- which(z$mat[, xx] == 1)
      if (length(present_y)) {
        points(
          rep(xx, length(present_y)), present_y,
          pch = 21, bg = "white", col = "white", cex = circle_cex
        )
      }
    }
    # Each cluster's binary matrix has a subtle grey outer frame. The
    # dendrogram intentionally remains open and unframed.
    matrix_usr <- par("usr")
    segments(matrix_usr[1], matrix_usr[3], matrix_usr[1], matrix_usr[4],
             col = "#8E9392", lwd = 0.9)
    segments(matrix_usr[1], matrix_usr[3], matrix_usr[2], matrix_usr[3],
             col = "#8E9392", lwd = 0.9)
    segments(matrix_usr[1], matrix_usr[4], matrix_usr[2], matrix_usr[4],
             col = "#8E9392", lwd = 0.9)
    segments(matrix_usr[2], matrix_usr[3], matrix_usr[2], matrix_usr[4],
             col = "#8E9392", lwd = 0.9)

    axis(
      4, at = nr / 2,
      labels = sprintf("A%d (n=%d)", z$cluster, z$n),
      las = 1, tick = FALSE, line = 0.15, cex.axis = 0.62,
      col.axis = cluster_colours[cl]
    )
    if (i == 1) {
      mtext(
        sprintf(
          "Scenario-category configurations grouped into %d archetypes",
          number_of_clusters
        ),
        side = 3, line = 0.25, adj = 0, font = 2, cex = 0.92
      )
    }
    if (i == number_of_clusters) {
      axis(1, at = seq_along(category_names), labels = category_names,
           las = 2, tick = FALSE, line = -0.2, cex.axis = 0.62)
      mtext("Future-state category", side = 1, line = 2.75, cex = 0.75)
    }
  }
}

figure_width <- as.numeric(Sys.getenv("SCENARIO_FIG_WIDTH", unset = "12"))
figure_height <- as.numeric(Sys.getenv("SCENARIO_FIG_HEIGHT", unset = "60"))

draw_binary_matrix(
  file.path(output_dir, paste0(figure_prefix, "_complete.pdf")),
  device = "pdf", show_row_labels = FALSE,
  width_in = figure_width, height_in = figure_height
)
draw_binary_matrix(
  file.path(output_dir, paste0(figure_prefix, "_complete.png")),
  device = "png", show_row_labels = FALSE,
  width_in = figure_width, height_in = figure_height,
  dpi = 300
)
draw_binary_matrix(
  file.path(output_dir, paste0(figure_prefix, "_with_IDs_complete.pdf")),
  device = "pdf", show_row_labels = TRUE,
  # Wider for the extra ID column, and taller so the row pitch leaves the IDs
  # at a readable ~7pt rather than the hairline size they were before.
  width_in = figure_width + 3, height_in = figure_height * 1.4
)

writeLines(
  c(
    "Figure logic:",
    sprintf("Rows are all %d scenarios; columns are T1-T38.", nrow(X)),
    sprintf("PAM/Jaccard membership defines the %d coloured archetype blocks.",
            number_of_clusters),
    "Within each block, average-linkage clustering on Jaccard distance orders rows.",
    "White circles mark all present cells in signature-category columns.",
    "Zero-distance duplicate profiles are shown as polytomies with terminal",
    "ticks and a bracket; no artificial internal hierarchy is imposed.",
    "The PDF without IDs is the recommended appendix figure.",
    "The tall PDF with IDs is an audit version, not a normal manuscript panel."
  ),
  file.path(output_dir, "README_figure_logic.txt")
)
