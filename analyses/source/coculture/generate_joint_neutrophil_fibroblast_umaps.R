#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrastr)
  library(scales)
})

package_dir <- "input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures"
input_file <- file.path(package_dir, "source_data", "Coculture_UMAP_coordinates_and_metadata.csv")
figure_dir <- file.path(package_dir, "figures", "UMAPs")
source_dir <- file.path(package_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

coords <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
coculture_codes <- c("CF", "UF", "UN", "UA5", "UD")
d <- coords[coords$ConditionCode %in% coculture_codes, , drop = FALSE]

condition_order <- c(
  "Control fibroblast co-culture",
  "UC fibroblast co-culture",
  "UC co-culture + NAMPTi",
  "UC co-culture + a5b1i",
  "UC co-culture + dual blockade"
)
d$condition_label <- factor(d$condition_label, levels = condition_order)
d$cell_compartment <- ifelse(
  d$clean_cell_class == "Neutrophil", "Neutrophil", "Fibroblast candidate"
)
d$cell_compartment <- factor(d$cell_compartment,
                             levels = c("Neutrophil", "Fibroblast candidate"))

cell_colors <- c(
  "Neutrophil" = "#4C78A8",
  "Fibroblast candidate" = "#D81B60"
)

theme_umap <- theme_void(base_family = "Arial", base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0),
    plot.subtitle = element_text(size = 11.5, color = "#333333"),
    strip.text = element_text(face = "bold", size = 10.5, margin = margin(5, 2, 5, 2)),
    strip.background = element_rect(fill = "#F2F2F2", color = NA),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.caption = element_text(size = 9, color = "#666666", hjust = 0),
    plot.margin = margin(10, 12, 10, 12)
  )

axis_arrows <- function(dat) {
  xr <- range(dat$UMAP_1); yr <- range(dat$UMAP_2)
  x0 <- xr[1] + 0.04 * diff(xr); y0 <- yr[1] + 0.04 * diff(yr)
  list(
    annotate("segment", x = x0, xend = x0 + 0.13 * diff(xr), y = y0, yend = y0,
             linewidth = 0.55, arrow = arrow(length = unit(0.10, "in"))),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + 0.13 * diff(yr),
             linewidth = 0.55, arrow = arrow(length = unit(0.10, "in"))),
    annotate("text", x = x0 + 0.065 * diff(xr), y = y0 - 0.035 * diff(yr),
             label = "UMAP 1", size = 3.1),
    annotate("text", x = x0 - 0.025 * diff(xr), y = y0 + 0.065 * diff(yr),
             label = "UMAP 2", angle = 90, size = 3.1)
  )
}

neut <- d[d$cell_compartment == "Neutrophil", , drop = FALSE]
fib <- d[d$cell_compartment == "Fibroblast candidate", , drop = FALSE]

p_joint <- ggplot() +
  geom_point_rast(
    data = neut, aes(UMAP_1, UMAP_2, color = cell_compartment),
    size = 0.20, alpha = 0.42, raster.dpi = 300
  ) +
  geom_point(
    data = fib, aes(UMAP_1, UMAP_2, color = cell_compartment),
    size = 1.05, alpha = 0.95, shape = 21, fill = "#D81B60", stroke = 0.20
  ) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, shape = 16))) +
  labs(
    title = "Neutrophils and fibroblasts in the co-culture atlas",
    subtitle = "Fibroblast candidates are enlarged and layered above neutrophils for visibility",
    color = "Cell compartment",
    caption = paste0(
      "Co-culture conditions only: ", comma(nrow(neut)), " neutrophils and ",
      comma(nrow(fib)), " fibroblast candidates"
    )
  ) + theme_umap + theme(legend.position = "right")
for (layer in axis_arrows(d)) p_joint <- p_joint + layer

p_split <- ggplot() +
  geom_point_rast(
    data = neut, aes(UMAP_1, UMAP_2, color = cell_compartment),
    size = 0.18, alpha = 0.38, raster.dpi = 300
  ) +
  geom_point(
    data = fib, aes(UMAP_1, UMAP_2, color = cell_compartment),
    size = 0.95, alpha = 0.95, shape = 21, fill = "#D81B60", stroke = 0.18
  ) +
  facet_wrap(~condition_label, ncol = 3) +
  scale_color_manual(values = cell_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, shape = 16))) +
  labs(
    title = "Neutrophil-fibroblast occupancy by co-culture condition",
    subtitle = "Every panel uses the same integrated UMAP coordinates and axis limits",
    color = "Cell compartment",
    caption = "Fibroblast candidates are enlarged to remain visible despite their lower abundance."
  ) + theme_umap +
  theme(legend.position = "bottom", panel.spacing = unit(8, "pt"))

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = ragg::agg_png, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}

save_plot(p_joint, "UMAP_4_Joint_neutrophils_and_fibroblasts", 10.5, 7.4)
save_plot(p_split, "UMAP_5_Joint_neutrophils_and_fibroblasts_by_condition", 13.2, 8.6)

write.csv(
  d,
  file.path(source_dir, "Coculture_joint_neutrophil_fibroblast_UMAP_source_data.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(table(Condition = d$condition_label, Cell_compartment = d$cell_compartment)),
  file.path(source_dir, "Coculture_joint_neutrophil_fibroblast_counts.csv"),
  row.names = FALSE
)

cat("Generated joint neutrophil-fibroblast UMAPs in", figure_dir, "\n")
cat("Co-culture cells:", nrow(d), " Neutrophils:", nrow(neut),
    " Fibroblast candidates:", nrow(fib), "\n")
