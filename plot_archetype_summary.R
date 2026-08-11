options(stringsAsFactors = FALSE)

# A one-page replacement for the 60-inch supplementary matrix.
#
# The tall figure exists so a reader can look up an individual scenario. That is
# an audit job and it belongs in the supplement. A main-text figure has to answer
# something else: what defines each archetype, and how does that compare with the
# corpus as a whole. That is "magnitude in a grid", so the cells use a single-hue
# sequential ramp -- one hue means a value in the top row and a value in the
# bottom row are directly comparable, which per-archetype hues would not allow.
#
# Signature categories are outlined in a warm accent rather than filled with a
# white dot: at 20% prevalence a cell sits near the light end of the ramp and a
# white marker on it would vanish.
#
# Categories carry their short names, not just T-numbers, so the figure reads
# without a lookup table.

analysis_dir <- Sys.getenv("SCENARIO_ANALYSIS_DIR", unset = "data_proc")
output_file <- Sys.getenv("SCENARIO_SUMMARY_FILE",
                          unset = file.path(analysis_dir, "Figure_archetype_summary"))
figure_width <- as.numeric(Sys.getenv("SCENARIO_SUMMARY_WIDTH", unset = "12"))

profiles <- read.csv(file.path(analysis_dir, "cluster_category_profiles.csv"))
summary_tab <- read.csv(file.path(analysis_dir, "scenario_archetype_summary.csv"))
category_names <- paste0("T", seq_len(38))
number_of_clusters <- nrow(summary_tab)

wide <- function(column) {
  m <- matrix(NA_real_, number_of_clusters, length(category_names))
  for (i in seq_len(number_of_clusters)) {
    rows <- profiles[profiles$cluster == summary_tab$cluster[i], ]
    m[i, ] <- as.numeric(rows[[column]][match(category_names, rows$category)])
  }
  m
}
prevalence <- wide("cluster_prevalence")
signature <- wide("signature") == 1
overall <- profiles$overall_prevalence[match(category_names, profiles$category)]
labels <- profiles$category_label[match(category_names, profiles$category)]

ramp <- colorRampPalette(c("#eef4fd", "#cde2fb", "#9ec5f4", "#6da7ec",
                           "#3987e5", "#256abf", "#184f95", "#0d366b"))(101)
accent <- "#EB6834"
ink <- "#2B3137"
muted <- "#7A8288"
cell_of <- function(v) ramp[1 + round(100 * pmin(pmax(v, 0), 1))]

figure_height <- 2.6 + 0.30 * number_of_clusters

row_labels <- sprintf("A%d  %s  (n=%d)", summary_tab$cluster,
                      summary_tab$archetype_name, summary_tab$n)
# Size the left margin to the widest label instead of guessing at it.
pdf(NULL, width = figure_width, height = figure_height)
label_inches <- max(strwidth(row_labels, units = "inches", cex = 0.58))
invisible(dev.off())
left_margin <- ceiling((label_inches + 0.15) / 0.2)

draw_summary <- function() {
  layout(matrix(c(1, 2, 3, 4), nrow = 2),
         widths = c(figure_width - 2.4, 2.4),
         heights = c(0.95, figure_height - 0.95))

  # Top left: corpus-wide prevalence, the baseline archetypes are read against.
  par(mar = c(0.15, left_margin, 1.6, 0.6), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(0.5, 38.5), ylim = c(0, max(overall) * 1.14),
       axes = FALSE, xlab = "", ylab = "")
  rect(seq_along(overall) - 0.36, 0, seq_along(overall) + 0.36, overall,
       col = "#C9D3DB", border = NA)
  axis(2, at = c(0, 0.2), labels = c("0", "20"), las = 1, tick = FALSE,
       line = -0.8, cex.axis = 0.52, col.axis = muted)
  mtext("corpus-wide\nprevalence (%)", side = 2, line = 0.9, cex = 0.5,
        col = muted, las = 1, adj = 1)
  mtext(sprintf("Archetype profiles: %d future-state categories x %d archetypes",
                length(category_names), number_of_clusters),
        side = 3, line = 0.25, adj = 0, font = 2, cex = 0.9, col = ink)

  # Bottom left: archetype x category prevalence.
  par(mar = c(6.8, left_margin, 0.2, 0.6), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(0.5, 38.5), ylim = c(number_of_clusters + 0.5, 0.5),
       axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(number_of_clusters)) {
    rect(seq_len(38) - 0.5, i - 0.5, seq_len(38) + 0.5, i + 0.5,
         col = cell_of(prevalence[i, ]), border = "#FFFFFF", lwd = 0.9)
    marked <- which(signature[i, ])
    if (length(marked)) {
      rect(marked - 0.45, i - 0.45, marked + 0.45, i + 0.45,
           col = NA, border = accent, lwd = 1.9)
    }
  }
  axis(1, at = seq_len(38), labels = sprintf("%s  %s", category_names, labels),
       las = 2, tick = FALSE, line = -0.6, cex.axis = 0.44, col.axis = ink)
  # One combined right-aligned label per row: separate axis + mtext calls
  # collided, because the name is far wider than the "A1 n=83" it sat beside.
  axis(2, at = seq_len(number_of_clusters), labels = row_labels,
       las = 1, tick = FALSE, line = -0.4, cex.axis = 0.58, col.axis = ink)

  # Right column: spacer, then the legend.
  par(mar = c(0.15, 0.4, 1.6, 0.4))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  par(mar = c(6.8, 0.4, 0.2, 0.4))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  bar_y <- seq(0.58, 0.86, length.out = 101)
  rect(0.05, bar_y[-101], 0.17, bar_y[-1], col = ramp[-101], border = NA)
  rect(0.05, 0.58, 0.17, 0.86, col = NA, border = "#C9D3DB", lwd = 0.7)
  text(0.20, c(0.58, 0.72, 0.86), c("0", "50", "100"), cex = 0.52,
       adj = 0, col = muted)
  text(0.05, 0.94, "share of the archetype's scenarios\ncoding the category (%)",
       cex = 0.53, adj = 0, col = ink)
  rect(0.05, 0.40, 0.115, 0.455, col = NA, border = accent, lwd = 1.9)
  text(0.145, 0.428, "signature category", cex = 0.56, adj = 0, col = ink)
  text(0.05, 0.315,
       paste("at least 20% prevalence, 10 points above",
             "the corpus rate, and 1.35x enriched", sep = "\n"),
       cex = 0.48, adj = 0, col = muted)
}

grDevices::cairo_pdf(paste0(output_file, ".pdf"), width = figure_width,
                     height = figure_height, family = "Arial")
draw_summary()
invisible(dev.off())

png(paste0(output_file, ".png"), width = figure_width, height = figure_height,
    units = "in", res = 400, type = "cairo-png")
draw_summary()
invisible(dev.off())

cat(sprintf("Wrote %s.{pdf,png}  (%.1f x %.1f in, %d archetypes)\n",
            output_file, figure_width, figure_height, number_of_clusters))
