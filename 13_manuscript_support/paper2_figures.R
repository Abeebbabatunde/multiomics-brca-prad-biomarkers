# Paper 2 figures (full colour, colour-blind-safe palette) - external validation and clinical associations
# Run: source("paper2_figures.R")   Requires: ggplot2
# Output: figures_paper2/Figure1.pdf, Figure2.pdf and .tiff (600 dpi), and sessionInfo_figures_paper2.txt
# The values are those reported in the thesis (Steps 12.5 and 12.6, Tables 4.19 and 4.20).

library(ggplot2)
dir.create("figures_paper2", showWarnings = FALSE)

navy <- "#1F3A5F"
col_cpg    <- "#9B2D9B"   # breast methylation (purple)
col_brna   <- "#0072B2"   # breast RNA (blue)
col_prna   <- "#D55E00"   # prostate RNA (vermillion)
layer_cols <- c("Breast methylation" = col_cpg, "Breast RNA" = col_brna, "Prostate RNA" = col_prna)

theme_pub <- theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.title          = element_text(face = "bold", size = 13, colour = navy),
    plot.title.position = "plot",
    plot.subtitle       = element_text(colour = "#4A5F7A", size = 10),
    panel.background    = element_rect(fill = "#F2F7FC", colour = NA),
    panel.grid.major.x  = element_line(colour = "#D6E4F2", linewidth = 0.4),
    axis.line           = element_line(colour = navy),
    axis.ticks          = element_line(colour = navy),
    axis.text           = element_text(colour = navy),
    axis.title          = element_text(colour = navy),
    strip.background    = element_blank(),
    strip.text          = element_text(face = "bold", colour = navy),
    legend.position     = "top",
    legend.text         = element_text(colour = navy),
    plot.caption        = element_text(colour = "#4A5F7A", size = 9, hjust = 0),
    plot.caption.position = "plot",
    plot.margin         = margin(8, 14, 8, 8)
  )

save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures_paper2", paste0(name, ".pdf")),  p, width = w, height = h)
  ggsave(file.path("figures_paper2", paste0(name, ".tiff")), p, width = w, height = h,
         dpi = 600, compression = "lzw")
}

# ---------------------------------------------------------------- Figure 1
rows   <- c("Breast methylation (239 / 46)",
            "Breast RNA, 8 of 10 genes (417 / 16)",
            "Breast RNA, CCN4 as WISP1 (417 / 16)",
            "Prostate RNA (47 / 48)")
layers <- c("Breast methylation", "Breast RNA", "Breast RNA", "Prostate RNA")
lab_auc  <- "AUC"
lab_sens <- "Sensitivity at training cut-off, %"

d <- rbind(
  data.frame(panel = lab_auc, row = rows, layer = layers, series = "Rebuilt panel",
             est = c(0.978, 1.000, 1.000, 0.949),
             lo  = c(0.965, 1.000, 1.000, 0.910),
             hi  = c(0.991, 1.000, 1.000, 0.988)),
  data.frame(panel = lab_auc, row = rows, layer = layers, series = "Initial panel",
             est = c(0.950, 0.932, 0.932, 0.971), lo = NA_real_, hi = NA_real_),
  data.frame(panel = lab_sens, row = rows, layer = layers, series = "Rebuilt panel",
             est = c(88.3, 69.1, 66.4, 85.1),
             lo  = c(83.5, 64.4, 61.7, 71.7),
             hi  = c(92.1, 73.5, 70.9, 93.8))
)
d$panel  <- factor(d$panel, levels = c(lab_auc, lab_sens))
d$row    <- factor(d$row, levels = rev(rows))
d$layer  <- factor(d$layer, levels = names(layer_cols))
d$series <- factor(d$series, levels = c("Rebuilt panel", "Initial panel"))
d_ci <- d[!is.na(d$lo) & d$hi > d$lo, ]

fig1 <- ggplot(d, aes(x = est, y = row)) +
  geom_segment(data = d_ci, aes(x = lo, xend = hi, y = row, yend = row, colour = layer),
               linewidth = 1) +
  geom_point(aes(shape = series, colour = layer, fill = series), size = 3.6, stroke = 1.1) +
  facet_wrap(~ panel, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = layer_cols, guide = "none") +
  scale_shape_manual(values = c("Rebuilt panel" = 21, "Initial panel" = 24), name = NULL) +
  scale_fill_manual(values = c("Rebuilt panel" = "#6E7F95", "Initial panel" = "white"), name = NULL,
                    guide = guide_legend(override.aes = list(colour = navy, fill = c("#6E7F95", "white")))) +
  labs(title = "A perfect ranking did not give a usable cut-off:\nbreast RNA reached AUC 1.000 but sensitivity 69.1%",
       x = NULL, y = NULL,
       caption = paste0("Cohort in brackets: tumour / normal samples. Colour marks the layer.\n",
                        "Bars: 95% CI (DeLong for the AUC, Clopper-Pearson for sensitivity). The AUC interval of breast RNA\n",
                        "is degenerate (1.000 to 1.000). Initial-panel AUC (open triangle) is shown without an interval.")) +
  theme_pub +
  theme(panel.spacing.x = unit(1.2, "lines"))
save_fig(fig1, "Figure1", 9.5, 4.4)

# ---------------------------------------------------------------- Figure 2
assoc <- data.frame(
  cancer  = c(rep("Breast: PAM50 subtype", 10), "Prostate: Gleason group"),
  feature = c("CAPN11", "SLC24A2", "ROBO3", "C8orf88", "cg05804846", "cg25710753",
              "COL10A1", "CCN4", "MME", "LOC107984270", "NETO2"),
  layer   = c("Breast RNA", "Breast RNA", "Breast RNA", "Breast RNA", "Breast methylation",
              "Breast methylation", "Breast RNA", "Breast RNA", "Breast RNA", "Breast RNA",
              "Prostate RNA"),
  est = c(0.246, 0.198, 0.135, 0.128, 0.124, 0.107, 0.090, 0.084, 0.078, 0.064, 0.068),
  lo  = c(0.185, 0.150, 0.095, 0.089, 0.083, 0.071, 0.055, 0.047, 0.040, 0.032, 0.019),
  hi  = c(0.316, 0.264, 0.193, 0.179, 0.186, 0.167, 0.144, 0.141, 0.135, 0.118, 0.164)
)
assoc$cancer  <- factor(assoc$cancer, levels = c("Breast: PAM50 subtype", "Prostate: Gleason group"))
assoc$layer   <- factor(assoc$layer, levels = names(layer_cols))
assoc$feature <- factor(assoc$feature, levels = rev(assoc$feature))

fig2 <- ggplot(assoc, aes(x = est, y = feature, colour = layer)) +
  geom_segment(aes(x = lo, xend = hi, y = feature, yend = feature), linewidth = 1) +
  geom_point(size = 3.6) +
  facet_grid(cancer ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = layer_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 0.36), breaks = seq(0, 0.35, 0.05), expand = c(0, 0)) +
  labs(title = "Ten of twelve breast features were associated with PAM50 subtype",
       x = "Epsilon-squared (95% bootstrap interval)", y = NULL,
       caption = "Associations with an FDR below 0.05 (Table 2).\nINHBA and CCL11 had no significant association.") +
  theme_pub +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1))
save_fig(fig2, "Figure2", 7.5, 5)

writeLines(capture.output(sessionInfo()), "figures_paper2/sessionInfo_figures_paper2.txt")
