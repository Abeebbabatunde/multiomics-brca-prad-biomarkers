# Paper 1 figures (black, white and grey) - patient-grouped stability selection
# Run: source("paper1_figures.R")   Requires: ggplot2
# Output: figures/Figure1.pdf ... Figure4.pdf and .tiff (600 dpi)
#
# The values below are the numbers reported in the thesis (Step 12.0).
# To plot directly from the saved analysis outputs instead, replace the data
# frames with read.csv() on the files written by the rebuild scripts
# (<layer>_selection_frequency_rebuild.csv, DIABLO_BRCA_patient_grouped_cv.csv,
#  DIABLO_BRCA_leakage_demonstration.csv,
#  DIABLO_BRCA_selection_leakage_demonstration.csv) and keep the plotting code.

library(ggplot2)
dir.create("figures", showWarnings = FALSE)

theme_pub <- theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.title        = element_text(face = "bold", size = 13),
    plot.title.position = "plot",
    plot.subtitle     = element_text(colour = "grey30", size = 10),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold"),
    axis.text         = element_text(colour = "black"),
    legend.position   = "top",
    plot.caption      = element_text(colour = "grey30", size = 9, hjust = 0)
  )

save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")),  p, width = w, height = h)
  ggsave(file.path("figures", paste0(name, ".tiff")), p, width = w, height = h,
         dpi = 600, compression = "lzw")
}

# ---------------------------------------------------------------- Figure 1
# Workflow diagram drawn with ggplot2 shapes (coordinates are arbitrary units).
boxes <- data.frame(
  x = c(16, 50, 84,  13, 37, 37, 37, 62,  62,  20),
  y = c(62, 62, 62,  40, 46, 40, 34, 40,  14,  14),
  w = c(28, 28, 28,  20, 20, 20, 20, 18,  44,  36),
  h = c(11, 11, 11,  14,  5,  5,  5, 14,  14,  10),
  main = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  label = c(
    "Paired cohorts\n56 breast and 17 prostate\npatient pairs",
    "Label-free inputs\nMethylation M-values,\nvariance-stabilised RNA",
    "Resampling of patients\n50 subsets of 75% of patients;\npairs kept together",
    "Paired limma\nt-test, patient as\nblocking factor\ntop 500 features",
    "LASSO (1-SE rule)",
    "SVM-RFE (20 kept)",
    "Boruta (100 runs)",
    "Consensus\nchosen by at least\n2 of 3 methods",
    "Selection frequency\nShare of the 50 subsets in which a feature\nis chosen; the panel is the features\nat 0.60 or above",
    "Permuted-label control\nLabels permuted within patient;\nsame procedure, 20 subsets"
  )
)
arrows <- data.frame(   # arrowheads
  x = c(30, 64, 84, 25, 25, 25, 50, 62, 38),
  y = c(62, 62, 57, 46, 40, 34, 40, 33, 14),
  xend = c(36, 70, 84, 27, 27, 27, 53, 62, 40),
  yend = c(62, 62, 52, 46, 40, 34, 40, 21, 14)
)
lines_only <- data.frame(  # plain connector segments
  x    = c(23, 25, 47, 47, 47, 50),
  y    = c(40, 34, 46, 40, 34, 34),
  xend = c(25, 25, 50, 50, 50, 50),
  yend = c(40, 46, 46, 40, 34, 46)
)

fig1 <- ggplot() +
  annotate("rect", xmin = 2, xmax = 98, ymin = 27, ymax = 52,
           fill = "grey92", colour = "grey40") +
  annotate("text", x = 4, y = 50, hjust = 0, size = 3, colour = "grey30",
           label = "Inside each subset, the whole selection procedure runs on that subset alone") +
  geom_rect(data = boxes,
            aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
                fill = main, linewidth = main), colour = "black") +
  geom_text(data = boxes, aes(x = x, y = y, label = label), size = 2.9, lineheight = 0.95) +
  geom_segment(data = lines_only, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "black", linewidth = 0.4) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               colour = "black", linewidth = 0.4,
               arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  annotate("text", x = 84, y = 40, size = 2.9, lineheight = 0.95,
           label = "Repeated in each\nof the 50 subsets") +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "grey80"), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 0.9), guide = "none") +
  coord_cartesian(xlim = c(0, 100), ylim = c(6, 70)) +
  labs(title = "Class labels are used only inside each resample of patients") +
  theme_void(base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.margin = margin(8, 8, 8, 8))
save_fig(fig1, "Figure1", 9, 6)

# ---------------------------------------------------------------- Figure 2
freq <- data.frame(
  feature = c("cg25710753", "cg05804846",
              "COL10A1", "SLC24A2", "INHBA", "CAPN11", "CCL11", "ROBO3",
              "LOC107984270", "C8orf88", "CCN4", "MME",
              "NETO2", "PDLIM5", "GMDS"),
  freq = c(0.74, 0.70,
           1.00, 1.00, 0.98, 0.96, 0.84, 0.82, 0.76, 0.74, 0.68, 0.68,
           0.68, 0.66, 0.62)
)
freq$layer <- c(rep("Breast methylation", 2), rep("Breast RNA", 10), rep("Prostate RNA", 3))
freq$layer   <- factor(freq$layer, levels = c("Breast methylation", "Breast RNA", "Prostate RNA"))
freq$feature <- factor(freq$feature, levels = rev(freq$feature))

fig2 <- ggplot(freq, aes(x = freq, y = feature)) +
  annotate("rect", xmin = 0, xmax = 0.10, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_col(fill = "grey35", width = 0.65) +
  geom_vline(xintercept = 0.60, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.2f", freq)), hjust = -0.25, size = 3.2) +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(0, 1.08), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = c(0, 0)) +
  labs(title = "Every rebuilt feature was selected in at least 62% of subsets;\nno permuted-label feature exceeded 10%",
       x = "Selection frequency (50 subsets of 75% of the patients)", y = NULL,
       caption = "Shaded: highest frequency with permuted labels (0.05 to 0.10).\nDashed: panel threshold, 0.60.") +
  theme_pub +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1))
save_fig(fig2, "Figure2", 7, 7)

# ---------------------------------------------------------------- Figure 3
diablo <- data.frame(
  kx = factor(rep(c(1, 2, 3, 5, 10, 20, 30, 50), 2), levels = c(1, 2, 3, 5, 10, 20, 30, 50)),
  series = rep(c("Original labels", "Permuted labels"), each = 8),
  err = c(0.035, 0.020, 0.011, 0.012, 0.018, 0.006, 0.013, 0.014,
          0.498, 0.498, 0.474, 0.504, 0.479, 0.521, 0.507, 0.508)
)

fig3 <- ggplot(diablo, aes(x = kx, y = err, group = series)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40") +
  geom_line(aes(linetype = series), linewidth = 0.7) +
  geom_point(aes(shape = series), size = 2.6, stroke = 0.9) +
  scale_linetype_manual(values = c("Original labels" = "solid", "Permuted labels" = "dashed"), name = NULL) +
  scale_shape_manual(values = c("Original labels" = 16, "Permuted labels" = 1), name = NULL) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2)) +
  annotate("text", x = 1, y = 0.545, hjust = 0, size = 3.2, colour = "grey30",
           label = "0.5 = no information") +
  labs(title = "With patient-grouped folds, DIABLO error was 0.035 or lower;\nwith permuted labels it stayed near 0.5",
       subtitle = "Breast cancer, 56 pairs: classification error by number of features per block",
       x = "Features per block", y = "Classification error") +
  theme_pub
save_fig(fig3, "Figure3", 7, 4.5)

# ---------------------------------------------------------------- Figure 4
leak <- data.frame(
  exp  = factor(rep(c("Fold assignment", "Filter placement"), each = 4),
                levels = c("Fold assignment", "Filter placement")),
  cond = factor(c("By patient", "By patient", "By sample", "By sample",
                  "Inside folds", "Inside folds", "Before CV", "Before CV"),
                levels = c("By patient", "By sample", "Inside folds", "Before CV")),
  kx   = factor(rep(c("5", "20"), 4), levels = c("5", "20")),
  err  = c(0.515, 0.488, 0.484, 0.504, 0.497, 0.511, 0.415, 0.398)
)

fig4 <- ggplot(leak, aes(x = cond, y = err, fill = kx)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey30") +
  geom_text(aes(y = 0.565, label = sprintf("%.3f", err)),
            position = position_dodge(width = 0.8), size = 3.1) +
  facet_grid(. ~ exp, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("5" = "grey80", "20" = "grey20"),
                    name = "Features per block") +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2), expand = c(0, 0)) +
  labs(title = "Only computing the filter before cross-validation\npulled the permuted-label error below 0.5",
       subtitle = "Breast DIABLO, 56 pairs: classification error with labels permuted within patient",
       x = NULL, y = "Classification error",
       caption = "Dashed line: 0.5 = no information.") +
  theme_pub
save_fig(fig4, "Figure4", 7, 4.8)

# Save the R and package versions used for the figures
writeLines(capture.output(sessionInfo()), "figures/sessionInfo_figures.txt")
