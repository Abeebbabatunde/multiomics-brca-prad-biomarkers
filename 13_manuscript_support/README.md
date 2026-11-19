# 13_manuscript_support: figure and check scripts for the manuscripts

Scripts that draw the manuscript figures and run the checks made while the manuscripts were prepared. Each script reads output files of the earlier steps; edit the project folder path at the top before running.

| Script | What it does | Output |
| --- | --- | --- |
| `paper1_figures.R` | Figures 1 to 4 of the combined manuscript (stability selection, leakage controls) in black, white and grey | `figures/` (PDF and 600 dpi TIFF), `sessionInfo` |
| `paper2_figures.R` | Figures 5 and 6 of the combined manuscript (external validation and clinical associations) in colour | `figures_paper2/` (PDF and TIFF), `sessionInfo` |
| `paperA_figure_scale.R` | Figure 2C of the integrative manuscript: variance explained per view before and after log2 transformation, prostate | `figures_paperA/` (PDF and TIFF), `sessionInfo` |
| `check_sex_gender.R` | Counts gender and sex at birth in the cohort patients, for the sex and gender statement | console output |
| `check_breast_dmp_threshold.R` | Recounts breast differentially methylated positions at delta-beta 0.1 and 0.2 from the saved limma table | `breast_DMP_threshold_sensitivity.csv` |

The figure scripts take their values from the result files of Steps 11 and 12; `paperA_figure_scale.R` uses the four variance values reported in the thesis.
