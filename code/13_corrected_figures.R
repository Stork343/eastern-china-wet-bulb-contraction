###############################################################################
# Supplementary null-calibration figure used by the current manuscript.
###############################################################################

library(data.table)
library(ggplot2)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))
output_dir <- file.path(project_dir, "output_corrected")

calibration <- fread(file.path(output_dir, "corrected_shift_calibration.csv"))
theme_paper <- theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 11))

figure <- ggplot(calibration, aes(p_value)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 fill = "#1F5A85", colour = "white") +
  geom_vline(xintercept = 0.05, linetype = 2, colour = "#D55E00") +
  labs(x = "Monte Carlo p-value", y = "Null simulations",
       title = "Circular-shift test calibration under the null") +
  theme_paper

ggsave(file.path(output_dir, "figS1_corrected_calibration.pdf"), figure,
       width = 5.8, height = 3.6, device = cairo_pdf)
