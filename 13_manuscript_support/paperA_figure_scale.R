# Paper A, Figure 2C: variance explained per view in the prostate MOFA2 model,
# before (count-scale RNA) and after (log2(count + 1) RNA).
# The four values are those reported in the thesis (Step 5.0 and Step 11.3):
#   original model: methylation 72.3 %, RNA 0.27 %
#   re-run model  : methylation 74.0 %, RNA 50.7 %
# If you prefer, replace the numbers below with the values printed by
# 03b_prostate_mofa2_log2.R (r2$r2_total) and the original Step 5.0 run.

library(ggplot2)

d <- data.frame(
  Model = factor(rep(c("Original (count-scale RNA)", "Re-run (log2 RNA)"), each = 2),
                 levels = c("Original (count-scale RNA)", "Re-run (log2 RNA)")),
  View  = factor(rep(c("Methylation", "RNA"), times = 2),
                 levels = c("Methylation", "RNA")),
  Var   = c(72.3, 0.27, 74.0, 50.7)
)
d$Label <- ifelse(d$Var < 1, sprintf("%.2f%%", d$Var), sprintf("%.1f%%", d$Var))

cols <- c("Methylation" = "#9B2D9B", "RNA" = "#0072B2")
navy <- "#1F2A5A"

p <- ggplot(d, aes(x = Model, y = Var, fill = View)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = Label), position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.6, colour = navy) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(limits = c(0, 85), expand = c(0, 0)) +
  labs(x = NULL, y = "Variance explained (%)") +
  theme_classic(base_size = 11) +
  theme(text = element_text(colour = navy),
        axis.text = element_text(colour = navy),
        axis.line = element_line(colour = navy),
        axis.ticks = element_line(colour = navy),
        legend.position = "top",
        plot.background = element_rect(fill = "#EEF3FB", colour = NA),
        panel.background = element_rect(fill = "#EEF3FB", colour = NA),
        legend.background = element_rect(fill = "#EEF3FB", colour = NA))

dir.create("figures_paperA", showWarnings = FALSE)
ggsave("figures_paperA/Figure2C_scale.pdf", p, width = 3.5, height = 3.2, device = cairo_pdf)
ggsave("figures_paperA/Figure2C_scale.tiff", p, width = 3.5, height = 3.2,
       dpi = 600, compression = "lzw")
writeLines(capture.output(sessionInfo()), "figures_paperA/sessionInfo_figure2C.txt")
