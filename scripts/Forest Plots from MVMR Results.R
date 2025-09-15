# ================================
# Forest Plots from MVMR Results — pQTL/eQTL (log10 + linear) + Runner
# ================================
# Author: Coral Romualdo-Perez
# Notes:
# - Same commenting style as prior scripts (banner + subheaders).
# - Absolutely NO functional logic changes; only clearer comments + pseudo file paths.
# - Keep the ability to make combined plots (pQTL + eQTL) or single-panel plots by
#   commenting/uncommenting the eQTL block exactly as noted below.
# - Color mapping in the log10 version is by `instrument_set` (dataset-specific labels).
#   In-line comment shows how to switch to `instrument_type` (pQTL/eQTL) if desired.
# ================================


# ================================
# CHUNK 0 — Silent library loader
# ================================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(writexl)
})


# ================================
# CHUNK 1 — Forest plot (log10 scale)
# ================================
# make_combined_mvmr_plot()
# - Inputs an Excel workbook with one sheet per instrument source (pQTL / eQTL).
# - By default, only pQTL is activated below (eQTL block is shown but commented).
# - Color mapping is by `instrument_set` (e.g., “PSI IVs”, “eQTL IVs”, “pQTL IVs”);
#   if you prefer a simple pQTL vs eQTL color, switch the aes to `instrument_type`
#   and uncomment the scale noted in-line.
make_combined_mvmr_plot <- function(file_path, pQTL_sheet, eQTL_sheet, output_dir, file_suffix) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- load inputs (comment/uncomment to combine vs single-source) ----
  pqtl_data <- read_excel(file_path, sheet = pQTL_sheet)
  # eqtl_data <- read_excel(file_path, sheet = eQTL_sheet)
  
  # ---- standardize headers (lowercase + safe names) ----
  names(pqtl_data) <- make.names(names(pqtl_data), unique = TRUE)
  names(pqtl_data) <- tolower(names(pqtl_data))
  # names(eqtl_data) <- make.names(names(eqtl_data), unique = TRUE)
  # names(eqtl_data) <- tolower(names(eqtl_data))
  
  # ---- normalize columns expected by plotter ----
  pqtl_data <- pqtl_data %>%
    rename(
      instrument_set = instrument.set,
      adjusted_for   = adjusted.for,
      No_SNPs        = no.snps,
      p_value        = p.value,
      f_stat         = f.stat,
      tier           = tier
    ) %>%
    mutate(instrument_type = "pQTL", No_SNPs = as.character(No_SNPs))
  
  # eqtl_data <- eqtl_data %>%
  #   rename(
  #     instrument_set = instrument.set,
  #     adjusted_for   = adjusted.for,
  #     No_SNPs        = no.snps,
  #     p_value        = p.value,
  #     f_stat         = f.stat
  #   ) %>%
  #   mutate(instrument_type = "eQTL", No_SNPs = as.character(No_SNPs))
  
  # ---- combine (activate eQTL above to include) ----
  mr_data <- bind_rows(pqtl_data) # , eqtl_data)
  
  # ---- compute ORs / CIs + table labels ----
  mr_data <- mr_data %>%
    dplyr::select(outcome, instrument_set, exposure, adjusted_for, No_SNPs,
                  beta, se, p_value, f_stat, instrument_type, tier) %>%
    filter(!is.na(beta) & !is.na(se)) %>%
    mutate(
      p_numeric   = as.numeric(p_value),
      is_significant = p_numeric < 0.05,
      beta = as.numeric(beta),
      se   = as.numeric(se),
      OR        = exp(beta),
      Lower_CI  = exp(beta - 1.96 * se),
      Upper_CI  = exp(beta + 1.96 * se),
      OR_CI_Label = paste0(
        formatC(exp(beta), format = "g", digits = 2), " (",
        formatC(exp(beta - 1.96 * se), format = "g", digits = 2), ", ",
        formatC(exp(beta + 1.96 * se), format = "g", digits = 2), ")"
      ),
      Pval_Label = formatC(p_numeric, format = "g", digits = 2),
      Header = FALSE
    ) %>%
    mutate(RowOrder = row_number())
  
  # ---- add a header row for the table panel ----
  header_row <- mr_data[0, ] %>%
    mutate(
      outcome     = "OUTCOME",
      exposure    = "EXPOSURE",
      adjusted_for= "ADJUSTED FOR",
      No_SNPs     = "NO SNPS",
      OR_CI_Label = "OR (95% CI)",
      Pval_Label  = "IVW P-VALUE",
      tier        = "TIER",
      OR = NA, Lower_CI = NA, Upper_CI = NA,
      RowOrder = 0, Header = TRUE, FontFace = "bold"
    )
  
  plot_data  <- bind_rows(header_row, mr_data) %>%
    mutate(
      Header   = RowOrder == 0,
      Row      = factor(as.character(RowOrder), levels = as.character(unique(RowOrder))),
      FontFace = ifelse(Header, "bold", ifelse(is_significant, "bold", "plain")),
      TextSize = ifelse(Header, 4, 3.2)
    )
  
  strip_data  <- plot_data %>% filter(!Header, RowOrder %% 2 == 0)
  all_levels  <- levels(plot_data$Row)
  header_rows <- plot_data %>% filter(Header)
  data_rows   <- plot_data %>% filter(!Header)
  
  # ---- (A) table panel (left) ----
  table_plot <- ggplot() +
    geom_rect(data = strip_data,
              aes(ymin = as.numeric(Row) - 0.5, ymax = as.numeric(Row) + 0.5),
              xmin = -Inf, xmax = Inf, fill = "gray95") +
    geom_text(data = header_rows, aes(x = 0.5, y = Row, label = outcome,       fontface = "bold"), hjust = 0, size = 13) +
    geom_text(data = header_rows, aes(x = 1.8, y = Row, label = exposure,      fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 2.5, y = Row, label = adjusted_for,  fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 4.9, y = Row, label = No_SNPs,       fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 5.9, y = Row, label = OR_CI_Label,   fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 9.0, y = Row, label = Pval_Label,    fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 8.5, y = Row, label = "TIER",        fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = data_rows,   aes(x = 0.1, y = Row, label = outcome,       fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 2.0, y = Row, label = exposure,      fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 3.5, y = Row, label = adjusted_for,  fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 5.2, y = Row, label = No_SNPs,       fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 6.2, y = Row, label = OR_CI_Label,   fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 8.3, y = Row, label = Pval_Label,    fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 9.5, y = Row, label = tier,          fontface = FontFace), hjust = 0, size = 13) +
    scale_x_continuous(breaks = c(0.1, 2.0, 3.5, 5.2, 6.2, 8.3, 9.5),
                       labels = NULL, limits = c(0.0, 10)) +
    scale_y_discrete(limits = rev(all_levels), drop = FALSE) +
    theme_minimal(base_size = 5) +
    theme(
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
  
  # ---- (B) forest panel (right), log10(OR) ----
  forest_plot <- ggplot(plot_data %>% filter(!Header, !is.na(OR)),
                        aes(x = log10(OR), y = Row, color = instrument_set)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1.2) +
    geom_errorbarh(aes(xmin = log10(pmax(Lower_CI, 1e-08)),
                       xmax = log10(pmin(Upper_CI, 1e+08))),
                   height = 0.3, linewidth = 5) +
    geom_point(shape = 21, fill = "white", size = 2.5) +
    scale_x_continuous(
      name   = expression("log"[10]*" Odds Ratio (CI)"),
      breaks = seq(-7, 7, 1),
      labels = function(x) format(10^x, scientific = TRUE, digits = 2),
      limits = c(-7, 7)
    ) +
    # If you prefer pQTL/eQTL coloring instead of `instrument_set`, switch aes(color = instrument_type)
    # and use the simple palette below:
    # scale_color_manual(values = c("pQTL" = "black", "eQTL" = "gray30"))
    scale_color_manual(values = c(
      "PSI IVs"  = "#0072B2",  # blue
      "eQTL IVs" = "#D55E00",  # orange
      "pQTL IVs" = "#009E73"   # green
    )) +
    scale_y_discrete(limits = rev(all_levels), drop = FALSE) +
    labs(
      x = expression("Log"[10]*" Odds Ratio (CI)\nProtective factor < 1 > Risk factor"),
      y = NULL, color = NULL
    ) +
    theme_minimal(base_size = 26) +
    theme(
      axis.text.x   = element_text(size = 34, angle = 45, hjust = 1),
      axis.title.x  = element_text(size = 38, face = "bold"),
      axis.ticks.y  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.text   = element_text(size = 38),
      legend.title  = element_blank(),
      legend.key    = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
  
  legend <- get_legend(forest_plot + theme(legend.position = "bottom"))
  forest_plot_nolegend <- forest_plot + theme(legend.position = "none")
  
  main_plot  <- plot_grid(table_plot, forest_plot_nolegend, ncol = 1, rel_widths = c(9, 9), align = "h")
  final_plot <- plot_grid(main_plot, legend, ncol = 1, rel_heights = c(1, 0.01))
  
  ggsave(
    filename = file.path(output_dir, paste0("MVMR_ForestPlot_horizontal_", file_suffix, ".png")),
    plot = final_plot, width = 40, height = 34, dpi = 600
  )
  
  write_xlsx(mr_data, file.path(output_dir, paste0("ISo_MVMR_Estimates_", file_suffix, ".xlsx")))
}


# ================================
# CHUNK 2 — Forest plot (linear OR scale)
# ================================
# make_combined_mvmr_plot2()
# - Same behavior as above, but the x-axis is raw OR (not log10).
# - Default color mapping here is by `instrument_type` (pQTL vs eQTL).
make_combined_mvmr_plot2 <- function(file_path, pQTL_sheet, eQTL_sheet, output_dir, file_suffix) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- load inputs (comment/uncomment to combine vs single-source) ----
  pqtl_data <- read_excel(file_path, sheet = pQTL_sheet)
  # eqtl_data <- read_excel(file_path, sheet = eQTL_sheet)
  
  # ---- standardize headers (lowercase + safe names) ----
  names(pqtl_data) <- make.names(names(pqtl_data), unique = TRUE)
  names(pqtl_data) <- tolower(names(pqtl_data))
  # names(eqtl_data) <- make.names(names(eqtl_data), unique = TRUE)
  # names(eqtl_data) <- tolower(names(eqtl_data))
  
  # ---- normalize columns expected by plotter ----
  pqtl_data <- pqtl_data %>%
    rename(
      instrument_set = instrument.set,
      adjusted_for   = adjusted.for,
      No_SNPs        = no.snps,
      p_value        = p.value,
      f_stat         = f.stat,
      tier           = tier
    ) %>%
    mutate(instrument_type = "pQTL", No_SNPs = as.character(No_SNPs))
  
  # eqtl_data <- eqtl_data %>%
  #   rename(
  #     instrument_set = instrument.set,
  #     adjusted_for   = adjusted.for,
  #     No_SNPs        = no.snps,
  #     p_value        = p.value,
  #     f_stat         = f.stat
  #   ) %>%
  #   mutate(instrument_type = "eQTL", No_SNPs = as.character(No_SNPs))
  
  # ---- combine (activate eQTL above to include) ----
  mr_data <- bind_rows(pqtl_data) # , eqtl_data)
  
  # ---- compute ORs / CIs + table labels ----
  mr_data <- mr_data %>%
    dplyr::select(outcome, instrument_set, exposure, adjusted_for, No_SNPs,
                  beta, se, p_value, f_stat, instrument_type, tier) %>%
    filter(!is.na(beta) & !is.na(se)) %>%
    mutate(
      p_numeric   = as.numeric(p_value),
      is_significant = p_numeric < 0.05,
      beta = as.numeric(beta),
      se   = as.numeric(se),
      OR        = exp(beta),
      Lower_CI  = exp(beta - 1.96 * se),
      Upper_CI  = exp(beta + 1.96 * se),
      OR_CI_Label = paste0(
        formatC(exp(beta),               format = "g", digits = 4), " (",
        formatC(exp(beta - 1.96 * se),   format = "g", digits = 2), ", ",
        formatC(exp(beta + 1.96 * se),   format = "g", digits = 2), ")"
      ),
      Pval_Label = formatC(p_numeric, format = "g", digits = 2),
      Header = FALSE
    ) %>%
    mutate(RowOrder = row_number())
  
  # ---- add a header row for the table panel ----
  header_row <- mr_data[0, ] %>%
    mutate(
      outcome     = "OUTCOME",
      exposure    = "EXPOSURE",
      adjusted_for= "ADJUSTED FOR",
      No_SNPs     = "NO SNPS",
      OR_CI_Label = "OR (95% CI)",
      Pval_Label  = "IVW P-VALUE",
      tier        = "TIER",
      OR = NA, Lower_CI = NA, Upper_CI = NA,
      RowOrder = 0, Header = TRUE, FontFace = "bold"
    )
  
  plot_data  <- bind_rows(header_row, mr_data) %>%
    mutate(
      Header   = RowOrder == 0,
      Row      = factor(as.character(RowOrder), levels = as.character(unique(RowOrder))),
      FontFace = ifelse(Header, "bold", ifelse(is_significant, "bold", "plain")),
      TextSize = ifelse(Header, 4, 3.2)
    )
  
  strip_data  <- plot_data %>% filter(!Header, RowOrder %% 2 == 0)
  all_levels  <- levels(plot_data$Row)
  header_rows <- plot_data %>% filter(Header)
  data_rows   <- plot_data %>% filter(!Header)
  
  # ---- (A) table panel (left) ----
  table_plot <- ggplot() +
    geom_rect(data = strip_data,
              aes(ymin = as.numeric(Row) - 0.5, ymax = as.numeric(Row) + 0.5),
              xmin = -Inf, xmax = Inf, fill = "gray95") +
    geom_text(data = header_rows, aes(x = 0.5, y = Row, label = outcome,       fontface = "bold"), hjust = 0, size = 13) +
    geom_text(data = header_rows, aes(x = 1.8, y = Row, label = exposure,      fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 2.5, y = Row, label = adjusted_for,  fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 4.9, y = Row, label = No_SNPs,       fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 5.9, y = Row, label = OR_CI_Label,   fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 9.0, y = Row, label = Pval_Label,    fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = header_rows, aes(x = 8.5, y = Row, label = "TIER",        fontface = "bold"), hjust = 0, size = 12) +
    geom_text(data = data_rows,   aes(x = 0.1, y = Row, label = outcome,       fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 2.0, y = Row, label = exposure,      fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 3.5, y = Row, label = adjusted_for,  fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 5.2, y = Row, label = No_SNPs,       fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 6.2, y = Row, label = OR_CI_Label,   fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 8.3, y = Row, label = Pval_Label,    fontface = FontFace), hjust = 0, size = 13) +
    geom_text(data = data_rows,   aes(x = 9.5, y = Row, label = tier,          fontface = FontFace), hjust = 0, size = 13) +
    scale_x_continuous(breaks = c(0.1, 2.0, 3.5, 5.2, 6.2, 8.3, 9.5),
                       labels = NULL, limits = c(0.0, 10)) +
    scale_y_discrete(limits = rev(all_levels), drop = FALSE) +
    theme_minimal(base_size = 5) +
    theme(
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
  
  # ---- (B) forest panel (right), linear OR ----
  forest_plot <- ggplot(plot_data %>% filter(!Header, !is.na(OR)),
                        aes(x = OR, y = Row, color = instrument_type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.3, linewidth = 7) +
    geom_point(shape = 21, fill = "white", size = 2.5) +
    scale_x_continuous(
      name   = "Odds Ratio (CI)",
      breaks = c(0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15),
      limits = c(0.85, 1.15),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_color_manual(values = c("pQTL" = "black", "eQTL" = "gray30")) +
    scale_y_discrete(limits = rev(all_levels), drop = FALSE) +
    labs(
      x = expression("Odds Ratio (CI)\nProtective factor < 1 > Risk factor"),
      y = NULL, color = NULL
    ) +
    theme_minimal(base_size = 26) +
    theme(
      axis.text.x   = element_text(size = 34, angle = 45, hjust = 1),
      axis.title.x  = element_text(size = 38, face = "bold"),
      axis.ticks.y  = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.text   = element_text(size = 38),
      legend.title  = element_blank(),
      legend.key    = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
  
  legend <- get_legend(forest_plot + theme(legend.position = "bottom"))
  forest_plot_nolegend <- forest_plot + theme(legend.position = "none")
  
  main_plot  <- plot_grid(table_plot, forest_plot_nolegend, ncol = 1, rel_widths = c(9, 9), align = "h")
  final_plot <- plot_grid(main_plot, legend, ncol = 1, rel_heights = c(1, 0.01))
  
  ggsave(
    filename = file.path(output_dir, paste0("MVMR_ForestPlot_horizontal_", file_suffix, ".png")),
    plot = final_plot, width = 40, height = 34, dpi = 600
  )
  
  write_xlsx(mr_data, file.path(output_dir, paste0("ISo_MVMR_Estimates_", file_suffix, ".xlsx")))
}


# ================================
# CHUNK 3 — Runner (examples; edit paths/sheets as needed)
# ================================
# • To generate a combined (pQTL+eQTL) plot:
#   - uncomment the `eqtl_data` lines inside the function(s) above;
#   - pass `eQTL_sheet = "YourSheetName"` here.
# • To generate a single-source plot (e.g., pQTL only), leave `eqtl_data` commented.
# • Pseudo paths are used here for reproducibility.
make_combined_mvmr_plot(
  file_path   = "project/analysis/MR_MVMR_Results.xlsx",
  pQTL_sheet  = "Sheet6",
  # eQTL_sheet  = "MRMV_cyto_isolation_eQTL_sig",
  output_dir  = "project/outputs/figures/mvmr/log10_combined/main_results/Ferkingstad_Folkersen",
  file_suffix = "color"
)

# Example runner for the linear-scale version (disabled by default):
# make_combined_mvmr_plot2(
#   file_path   = "project/analysis/MR_MVMR_Results.xlsx",
#   pQTL_sheet  = "Sheet6",
#   # eQTL_sheet  = "MRMV_cyto_isolation_eQTL_sig",
#   output_dir  = "project/outputs/figures/mvmr/linear_combined/main_results/Ferkingstad_Folkersen",
#   file_suffix = "color_linear"
# )