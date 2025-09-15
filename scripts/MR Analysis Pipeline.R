# ================================
# MR Analysis Pipeline
# ================================
# Author: Coral Romualdo-Perez
# Date: 2025-06-04
# Notes:
# - Pseudo paths are used throughout (replace with your project layout).

# ================================


# ================================
# CHUNK 1 — Load packages (run first)
# ================================
knitr::opts_chunk$set(echo = TRUE)
# Load libraries
library(devtools)
library(MRInstruments)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(genetics.binaRies)
library(MVMR)
library(MRPRESSO)

# PLINK is set if needed
genetics.binaRies::get_plink_binary()
plink_path <- get_plink_binary()


# ================================
# CHUNK 2 — Pipeline function (run second)
# ALPHA = EXPOSURE ; OMEGA = OUTCOME
# ================================
### MR Pipeline Run
run_mr_pipeline <- function(exposure_file, outcome_file, 
                            EXPOSURE_NAME, OUTCOME_NAME, 
                            samplesize_exposure, samplesize_outcome, 
                            ncase_outcome, ncontrol_outcome,
                            output_folder) {
  
  # Create output folder if it doesn't exist, if it exists it will automatically save to folder
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    cat("Created output folder:", output_folder, "\n")
  } else {
    cat("Using existing output folder:", output_folder, "\n")
  }
  
  ######### Debugging #########
  cat("Running MR pipeline for:", EXPOSURE_NAME, "to", OUTCOME_NAME, "\n")
  cat("Exposure sample size:", samplesize_exposure, "\n")
  cat("Outcome sample size:", samplesize_outcome, "\n")
  if (!is.na(ncase_outcome) & !is.na(ncontrol_outcome)) {
    cat("Binary outcome: ncase =", ncase_outcome, ", ncontrol =", ncontrol_outcome, "\n")
  } else {
    cat("Continuous outcome detected.\n")
  }
  ######### END Debugging #########
  
  # Load exposure
  ALPHA_data_exposure <- read.table(exposure_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) # keep as sep = "\t"
  ALPHA_data_exposure$id.exposure <- EXPOSURE_NAME
  
  # Load outcome
  OMEGA_data_outcome <- read.table(outcome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) # keep as sep = "\t"
  OMEGA_data_outcome$id.outcome <- OUTCOME_NAME
  
  # Harmonize
  OMEGA_harmonized_data <- harmonise_data(exposure_dat = ALPHA_data_exposure, outcome_dat = OMEGA_data_outcome, action = 1)
  OMEGA_harmonized_data <- as_tibble(OMEGA_harmonized_data)
  # Save Harmonization Results
  write.csv(
    OMEGA_harmonized_data,
    file.path(output_folder, paste0("Harmonized_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".csv")),
    row.names = FALSE
  )
  
  #  OMEGA_harmonized_data <- OMEGA_harmonized_data %>%
  #    dplyr::filter(
  #      isTRUE(mr_keep.exposure),                              # TwoSampleMR’s keep flag
  #      !is.na(beta.exposure), !is.na(se.exposure),
  #      !is.na(beta.outcome),  !is.na(se.outcome)
  #    )
  
  if (nrow(OMEGA_harmonized_data) == 0) {
    cat("No usable SNPs after harmonisation filtering (mr_keep/NA checks). Skipping MR.\n")
    return(list(harmonized = OMEGA_harmonized_data, mr_results = NULL))
  }
  
  #### WALD RATIO, Single SNP MR ####
  if (nrow(OMEGA_harmonized_data) == 1) {
    cat("\n️ Only one SNP after harmonisation — using Wald ratio only.\n\n")
    
    # Check for NA values
    if (any(is.na(OMEGA_harmonized_data$beta.exposure), 
            is.na(OMEGA_harmonized_data$beta.outcome),
            is.na(OMEGA_harmonized_data$se.exposure),
            is.na(OMEGA_harmonized_data$se.outcome))) {
      cat("Single SNP available but has missing values — skipping Wald ratio.\n")
      return(list(harmonized = OMEGA_harmonized_data, mr_results = NULL))
    }
    
    # Run Wald ratio
    OMEGA_mr_results <- mr(OMEGA_harmonized_data)
    print(OMEGA_mr_results)
    
    # Save MR results
    write.csv(OMEGA_mr_results, file.path(output_folder, paste0(EXPOSURE_NAME, "_to_", OUTCOME_NAME, "_MR_results.csv")), row.names = FALSE)
    
    # Also save F-statistics
    Fstat_exposure <- (OMEGA_harmonized_data$beta.exposure^2) / (OMEGA_harmonized_data$se.exposure^2)
    Fstat_outcome <- (OMEGA_harmonized_data$beta.outcome^2) / (OMEGA_harmonized_data$se.outcome^2)
    mean_F_exposure <- mean(Fstat_exposure, na.rm = TRUE)
    mean_F_outcome <- mean(Fstat_outcome, na.rm = TRUE)
    
    print(paste("Mean F-statistic for", EXPOSURE_NAME, ":", round(mean_F_exposure, 2)))
    print(paste("Mean F-statistic for", OUTCOME_NAME, ":", round(mean_F_outcome, 2)))
    
    fstats_df <- data.frame(
      SNP = OMEGA_harmonized_data$SNP,
      F_stat_exposure = Fstat_exposure,
      F_stat_outcome = Fstat_outcome
    )
    write.csv(fstats_df,
              file.path(output_folder, paste0("Fstats_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".csv")),
              row.names = FALSE)
    
    # Log mean F stats
    writeLines(
      c(
        paste("Mean F-statistic for", EXPOSURE_NAME, ":", round(mean_F_exposure, 2)),
        paste("Mean F-statistic for", OUTCOME_NAME, ":", round(mean_F_outcome, 2))
      ),
      con = file.path(output_folder, paste0("Mean_Fstats_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".txt"))
    )
    
    cat("Single-SNP analysis complete. Skipping multi-SNP diagnostics.\n")
    
    return(list(harmonized = OMEGA_harmonized_data, mr_results = OMEGA_mr_results))
  }
  #### WALD RATIO END ####
  
  cat("\n Multiple SNPs detected — running full MR pipeline.\n\n")
  
  # Basic MR
  OMEGA_mr_results <- mr(OMEGA_harmonized_data)
  print(OMEGA_mr_results)
  
  # Save MR results
  write.csv(OMEGA_mr_results, file.path(output_folder, paste0(EXPOSURE_NAME, "_to_", OUTCOME_NAME, "_MR_results.csv")), row.names = FALSE)
  
  # Scatterplot
  scatterplot_OMEGA <- mr_scatter_plot(OMEGA_mr_results, OMEGA_harmonized_data)
  if (length(scatterplot_OMEGA) > 0) {
    scatterplot_OMEGA_labels <- scatterplot_OMEGA[[1]] +
      geom_text_repel(aes(label = SNP), size = 3, max.overlaps = 200) +
      labs(
        x = paste("SNP effect on", EXPOSURE_NAME),
        y = paste("SNP effect on", OUTCOME_NAME)
      )
    print(scatterplot_OMEGA_labels)
    ggsave(file.path(output_folder, paste0("scatterplot_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".png")), scatterplot_OMEGA_labels, width = 8, height = 6, dpi = 300)
  } else {
    cat("No scatterplot generated.\n")
  }
  
  # Heterogeneity
  OMEGA_heterogeneity <- mr_heterogeneity(OMEGA_harmonized_data)
  print(OMEGA_heterogeneity)
  write.csv(OMEGA_heterogeneity, file.path(output_folder, paste0(OUTCOME_NAME, "_Heterogeneity.csv")), row.names = FALSE)
  
  # Pleiotropy
  OMEGA_pleiotropy_test <- mr_pleiotropy_test(OMEGA_harmonized_data)
  print(OMEGA_pleiotropy_test)
  write.csv(OMEGA_pleiotropy_test, file.path(output_folder, paste0(OUTCOME_NAME, "_Pleiotropy.csv")), row.names = FALSE)
  
  # F Statistics
  Fstat_exposure <- (OMEGA_harmonized_data$beta.exposure^2) / (OMEGA_harmonized_data$se.exposure^2)
  mean_F_exposure <- mean(Fstat_exposure, na.rm = TRUE)
  print(paste("Mean F-statistic for", EXPOSURE_NAME, ":", round(mean_F_exposure, 2)))
  
  Fstat_outcome <- (OMEGA_harmonized_data$beta.outcome^2) / (OMEGA_harmonized_data$se.outcome^2)
  mean_F_outcome <- mean(Fstat_outcome, na.rm = TRUE)
  print(paste("Mean F-statistic for", OUTCOME_NAME, ":", round(mean_F_outcome, 2)))
  
  # Save F-statistics per SNP
  fstats_df <- data.frame(
    SNP = OMEGA_harmonized_data$SNP,
    F_stat_exposure = Fstat_exposure,
    F_stat_outcome = Fstat_outcome
  )
  write.csv(fstats_df,
            file.path(output_folder, paste0("Fstats_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".csv")),
            row.names = FALSE)
  
  # Also log mean F stats
  writeLines(
    c(
      paste("Mean F-statistic for", EXPOSURE_NAME, ":", round(mean_F_exposure, 2)),
      paste("Mean F-statistic for", OUTCOME_NAME, ":", round(mean_F_outcome, 2))
    ),
    con = file.path(output_folder, paste0("Mean_Fstats_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".txt"))
  )
  
  # Ensure samplesize column for leave-one-out
  if (!"samplesize" %in% names(OMEGA_harmonized_data)) {
    if (!is.null(samplesize_outcome) && !is.na(samplesize_outcome)) {
      OMEGA_harmonized_data$samplesize <- samplesize_outcome
      cat("Assigned harmonised data samplesize column:", samplesize_outcome, "\n")
    } else {
      stop("ERROR: samplesize_outcome is missing or NA. Cannot assign required 'samplesize' column for leave-one-out analysis.")
    }
  }
  
  # LOO
  OMEGA_loo <- mr_leaveoneout(OMEGA_harmonized_data)
  plot_loo <- mr_leaveoneout_plot(OMEGA_loo)
  
  write.csv(OMEGA_loo, file.path(output_folder, paste0(OUTCOME_NAME, "_LOO.csv")), row.names = FALSE)
  
  ggsave(filename = file.path(output_folder, paste0("LOO_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".png")),
         plot = plot_loo[[1]],
         width = 6, height = 6, dpi = 300)
  
  ggsave(filename = file.path(output_folder, paste0("LOO_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".pdf")),
         plot = plot_loo[[1]],
         width = 6, height = 6)
  
  # SINGLE SNP ANALYSIS
  OMEGA_singlesnp_results <- mr_singlesnp(OMEGA_harmonized_data)
  write.csv(OMEGA_singlesnp_results, file.path(output_folder, paste0("Single_SNP_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".csv")), row.names = FALSE)
  
  # MR-RAPS
  OMEGA_data_RAPS <- OMEGA_harmonized_data %>%
    dplyr::select(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome) %>%
    na.omit()
  
  if (nrow(OMEGA_data_RAPS) >= 10) {
    res_raps <- tryCatch({
      mr(
        OMEGA_harmonized_data,
        method_list = c("mr_raps"),
        parameters = list(over.dispersion = TRUE, loss.function = "huber")
      )
    }, error = function(e) {
      warning(paste("MR-RAPS failed:", e$message))
      return(NULL)
    })
    
    if (!is.null(res_raps)) {
      write.csv(
        res_raps,
        file.path(output_folder, paste0(EXPOSURE_NAME, "_to_", OUTCOME_NAME, "_MRRAPS.csv")),
        row.names = FALSE
      )
    }
  } else {
    warning("MR-RAPS skipped: fewer than 10 SNPs available.")
  }
  
  # Funnel plot
  OMEGA_singlesnp_results <- mr_singlesnp(OMEGA_harmonized_data)
  OMEGA_funnel_plot <- mr_funnel_plot(OMEGA_singlesnp_results)[[1]]
  OMEGA_funnel_plot_annotated <- OMEGA_funnel_plot +
    geom_text_repel(aes(label = SNP), size = 3)
  print(OMEGA_funnel_plot_annotated)
  ggsave(file.path(output_folder, paste0("funnel_plot_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".png")), plot = OMEGA_funnel_plot_annotated, width = 8, height = 6, dpi = 300)
  
  # MR-PRESSO
  OMEGA_data_MRPRESSO <- data.frame(
    SNP = OMEGA_harmonized_data$SNP,
    BetaOutcome = OMEGA_harmonized_data$beta.outcome,
    BetaExposure = OMEGA_harmonized_data$beta.exposure,
    SdOutcome = OMEGA_harmonized_data$se.outcome,
    SdExposure = OMEGA_harmonized_data$se.exposure
  )
  
  OMEGA_data_MRPRESSO <- na.omit(OMEGA_data_MRPRESSO)
  
  OMEGA_mr_presso_results <- mr_presso(
    BetaOutcome = "BetaOutcome",
    BetaExposure = "BetaExposure",
    SdOutcome = "SdOutcome",
    SdExposure = "SdExposure",
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    data = OMEGA_data_MRPRESSO,
    NbDistribution = 1000,
    SignifThreshold = 0.05
  )
  
  OMEGA_global_test <- OMEGA_mr_presso_results$`Main MR results`
  write.csv(OMEGA_global_test, file.path(output_folder, paste0("MRPRESSO_", OUTCOME_NAME, "_GlobalTest.csv")), row.names = FALSE)
  
  if (!is.null(OMEGA_mr_presso_results$`MR-PRESSO results`$`Outlier Test`)) {
    write.csv(OMEGA_mr_presso_results$`MR-PRESSO results`$`Outlier Test`, file.path(output_folder, paste0("MRPRESSO_", OUTCOME_NAME, "_OutlierTest.csv")), row.names = FALSE)
  }
  
  if (!is.null(OMEGA_mr_presso_results$`MR-PRESSO results`$`Distortion Test`)) {
    write.csv(OMEGA_mr_presso_results$`MR-PRESSO results`$`Distortion Test`, file.path(output_folder, paste0("MRPRESSO_", OUTCOME_NAME, "_DistortionTest.csv")), row.names = FALSE)
  }
  
  ### PATCH: Continuous Outcomes or Binary ###
  is_qtl_outcome <- grepl("pQTL|eQTL", OUTCOME_NAME, ignore.case = TRUE)
  
  if (!is_qtl_outcome) {
    cat("Checking for Steiger filtering eligibility...\n")
    
    # Assign sample sizes safely
    if (!"samplesize.exposure" %in% names(OMEGA_harmonized_data) || all(is.na(OMEGA_harmonized_data$samplesize.exposure))) {
      if (!is.null(samplesize_exposure) && !is.na(samplesize_exposure)) {
        OMEGA_harmonized_data$samplesize.exposure <- samplesize_exposure
        cat("Assigned samplesize.exposure:", samplesize_exposure, "\n")
      } else {
        cat("Warning: samplesize_exposure is missing or NA.\n")
      }
    }
    
    if (!"samplesize.outcome" %in% names(OMEGA_harmonized_data) || all(is.na(OMEGA_harmonized_data$samplesize.outcome))) {
      if (!is.null(samplesize_outcome) && !is.na(samplesize_outcome)) {
        OMEGA_harmonized_data$samplesize.outcome <- samplesize_outcome
        cat("Assigned samplesize.outcome:", samplesize_outcome, "\n")
      } else {
        cat("Warning: samplesize_outcome is missing or NA.\n")
      }
    }
    
    if (!is.na(ncase_outcome) && !is.na(ncontrol_outcome)) {
      OMEGA_harmonized_data$ncase.outcome <- ncase_outcome
      OMEGA_harmonized_data$ncontrol.outcome <- ncontrol_outcome
      cat("Assigned binary outcome case/control:", ncase_outcome, "/", ncontrol_outcome, "\n")
    } else {
      cat("Detected continuous outcome or missing case/control info.\n")
    }
    
    # Final sanity check
    if (nrow(OMEGA_harmonized_data) < 3) {
      cat("Skipping Steiger filtering and plot: too few SNPs (", nrow(OMEGA_harmonized_data), ")\n")
    } else if (any(is.na(OMEGA_harmonized_data$samplesize.exposure)) || any(is.na(OMEGA_harmonized_data$samplesize.outcome))) {
      cat("Skipping Steiger filtering: missing sample size info.\n")
    } else {
      OMEGA_steiger_results <- tryCatch({
        steiger_filtering(OMEGA_harmonized_data)
      }, error = function(e) {
        cat("Steiger filtering failed:", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(OMEGA_steiger_results) &&
          nrow(OMEGA_steiger_results) > 2 &&
          !all(is.na(OMEGA_steiger_results$rsq.exposure)) &&
          !all(is.na(OMEGA_steiger_results$rsq.outcome))) {
        
        # Identify which column holds the directionality
        dir_col <- NULL
        if ("steiger_dir" %in% colnames(OMEGA_steiger_results)) {
          dir_col <- "steiger_dir"
        } else if ("correct_causal_direction" %in% colnames(OMEGA_steiger_results)) {
          dir_col <- "correct_causal_direction"
        } else {
          cat("Warning: No directionality column found in Steiger results.\n")
        }
        
        # Generate plot if column exists
        if (!is.null(dir_col)) {
          steiger_plot <- ggplot(OMEGA_steiger_results, aes_string(x = "rsq.exposure", y = "rsq.outcome", color = dir_col)) +
            geom_point(size = 2) +
            geom_text_repel(aes(label = SNP), size = 3, max.overlaps = 100) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
            theme_minimal() +
            scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "red"),
                               labels = c("TRUE" = "Causal", "FALSE" = "Reverse")) +
            labs(title = paste0("Steiger Filtering: Variance Explained in ", EXPOSURE_NAME, " vs. ", OUTCOME_NAME),
                 x = expression(R^2~"Exposure"),
                 y = expression(R^2~"Outcome"),
                 color = "Direction")
          
          print(steiger_plot)
          ggsave(file.path(output_folder, paste0("steiger_plot_", EXPOSURE_NAME, "_to_", OUTCOME_NAME, ".png")),
                 plot = steiger_plot, width = 8, height = 6, dpi = 300)
        } else {
          cat("Steiger plot not generated — directionality column not found.\n")
        }
        
      } else {
        cat("Steiger plot not generated — insufficient data or error in filtering.\n")
      }
    }
  }
  
  # Return
  return(list(harmonized = OMEGA_harmonized_data, 
              mr_results = OMEGA_mr_results, 
              presso = OMEGA_mr_presso_results))
}


# ================================
# CHUNK 3 — Single analysis convenience wrapper
# ALPHA = EXPOSURE, OMEGA = OUTCOME
# (run this block only if you want ONE pair)
# ================================
# run_mr_pipeline(
#   exposure_file = "project/derived/instruments/DeCODE/Hg37/MR_instruments_FINAL/IL6R/IL6R_MR_instruments_nonremissionMDD.txt",
#   outcome_file  = "project/derived/outcomes/Treatment_Resistant_MDD_Outcome.txt",
#   EXPOSURE_NAME = "IL-6Ra pQTL, cis",
#   OUTCOME_NAME  = "Non-Remission MDD",
#   samplesize_exposure = 35559,     # Exposure total cohort size
#   samplesize_outcome  = 5151,      # Outcome total cohort size
#   ncase_outcome       = 3299,      # Case; put NA for continuous
#   ncontrol_outcome    = 1852,      # Control; put NA for continuous
#   output_folder       = "project/outputs/mr/DeCODE_pQTLs/IL6R_to_NonRemissionMDD"
# )


# ================================
# CHUNK 4 — Batch prep (control center)
# Prepare vectors of exposures/outcomes and metadata
# ================================
# Define exposures (add/remove files as needed)
exposure_files <- c(
  # pQTL (UKB example, active)
  "project/derived/instruments/UKB_PPP/MR_instruments_FINAL/IL6R_MR_instruments_FINAL_percentMDD.txt"
)

# one line for each file, same order as exposures
exposure_names <- c(
  "IL-6R pQTL, cis, UKB"
)

# Define outcomes (add/remove files as needed)
outcome_files <- c(
  "project/derived/outcomes/Remission_MDD_Outcome.txt"
)

# one line for each file, same order as outcomes
outcome_names <- c(
  "Percentage Improvement MDD"
)

# Exposure sample size (length must match exposures)
exposure_sample_size <- c(
  33413   # IL-6Ra pQTL, UKB
)

# Outcome sample sizes (length must match outcomes)
outcome_sample_sizes <- c(
  5218    # Remission MDD
)

# Outcome case/control — NA for continuous
ncase_outcome_list    <- c(NA)  # Remission MDD
ncontrol_outcome_list <- c(NA)  # Remission MDD

# ==== Consistency Check ====
stopifnot(
  length(exposure_files) == length(exposure_names),
  length(exposure_files) == length(exposure_sample_size),
  length(outcome_files)  == length(outcome_names),
  length(outcome_files)  == length(outcome_sample_sizes),
  length(outcome_files)  == length(ncase_outcome_list),
  length(outcome_files)  == length(ncontrol_outcome_list)
)

# Master output folder (rename as needed)
master_output_folder <- "project/outputs/mr/UKB_pQTLs/IL6R_to_PercentageMDD_rerun"
if (!dir.exists(master_output_folder)) dir.create(master_output_folder, recursive = TRUE)


# ================================
# CHUNK 5 — Batch runner (multiple analyses)
# Run after CHUNK 4 is configured
# ================================
######## space saver

# Initialize empty master summary, follows TwoSampleMR Requirements
batch_summary <- data.frame(
  Exposure = character(),
  Outcome  = character(),
  Method   = character(),
  Beta     = numeric(),
  SE       = numeric(),
  Pval     = numeric(),
  stringsAsFactors = FALSE
)

# Start batch run with timer
start_time <- Sys.time()
cat("\n==========================\n")
cat("Starting full MR batch run!\n")
cat("Start time:", format(start_time), "\n")
cat("==========================\n\n")

# Outer loop: exposures
for (i in seq_along(exposure_files)) {
  
  # Inner loop: outcomes
  for (j in seq_along(outcome_files)) {
    
    cat("\n========== [", i, "/", length(exposure_files), "]", exposure_names[i],
        "to", "[", j, "/", length(outcome_files), "]", outcome_names[j], "==========\n")
    
    # Define output subfolder
    output_subfolder <- file.path(master_output_folder, paste0(exposure_names[i], "_to_", outcome_names[j]))
    
    # Try MR run
    result <- tryCatch({
      
      run_result <- run_mr_pipeline(
        exposure_file        = exposure_files[i],
        outcome_file         = outcome_files[j],
        EXPOSURE_NAME        = exposure_names[i],
        OUTCOME_NAME         = outcome_names[j],
        samplesize_exposure  = exposure_sample_size[i],
        samplesize_outcome   = outcome_sample_sizes[j],
        ncase_outcome        = ncase_outcome_list[j],
        ncontrol_outcome     = ncontrol_outcome_list[j],
        output_folder        = output_subfolder
      )
      
      # If MR results exist, add to batch summary
      if (!is.null(run_result$mr_results) && nrow(run_result$mr_results) > 0) {
        batch_summary <- rbind(batch_summary, data.frame(
          Exposure = exposure_names[i],
          Outcome  = outcome_names[j],
          Method   = run_result$mr_results$method,
          Beta     = run_result$mr_results$b,
          SE       = run_result$mr_results$se,
          Pval     = run_result$mr_results$pval
        ))
      }
      
      cat("SUCCESS:", exposure_names[i], "to", outcome_names[j], "\n")
      
    }, error = function(e) {
      cat("ERROR in MR run:", exposure_names[i], "to", outcome_names[j], "\n")
      cat("Reason:", e$message, "\n")
      
      # Optionally log error in batch summary
      batch_summary <<- rbind(batch_summary, data.frame(
        Exposure = exposure_names[i],
        Outcome  = outcome_names[j],
        Method   = "ERROR",
        Beta     = NA,
        SE       = NA,
        Pval     = NA
      ))
      
      return(NULL)
    })
    
    # Small pause if needed to avoid overloading server
    # Sys.sleep(1)
    
  }
}

# End batch run
end_time <- Sys.time()
cat("\n==========================\n")
cat("Full MR batch run complete!\n")
cat("End time:", format(end_time), "\n")
cat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
cat("==========================\n")

# Save batch summary CSV
batch_summary_file <- file.path(master_output_folder, "MR_BATCH_SUMMARY.csv")
write.csv(batch_summary, batch_summary_file, row.names = FALSE)

cat("\n Saved batch summary to:", batch_summary_file, "\n\n")


# ================================
# CHUNK 6 — Quick overlap runner (optional utility)
# ================================
# Load libraries
library(readr)
library(dplyr)

# === Paths ===
exposure_file <- c(
  "project/COJO/final_instruments/pQTLs/IL-1RN_TwoSampleMR_ready.txt"
)

outcome_files <- c(
  "project/derived/outcomes/Isolation_GWAS_Outcome.txt"
)

output_table_path <- "project/outputs/tables/IL1RA_pQTL_to_isolation.csv"

# === Load exposure SNP list ===
exposure_data <- read_tsv(exposure_file, show_col_types = FALSE)
exposure_snps <- exposure_data$SNP

cat("\n Number of exposure SNPs:", length(exposure_snps), "\n\n")

# === Initialize results ===
result_list <- list()

# === Loop through outcomes ===
for (file in outcome_files) {
  outcome_data <- read_tsv(file, show_col_types = FALSE)
  outcome_snps <- outcome_data$SNP
  
  shared_snps <- intersect(exposure_snps, outcome_snps)
  
  # Print to console
  cat("Outcome:", basename(file), "\n")
  cat("  - Total SNPs in outcome:", length(outcome_snps), "\n")
  cat("  - Overlapping SNPs:", length(shared_snps), "\n")
  if (length(shared_snps) > 0) {
    cat("  - Shared SNP IDs:", paste(shared_snps, collapse = ", "), "\n")
  }
  cat("\n")
  
  # Save in results list
  result_list[[basename(file)]] <- data.frame(
    Outcome = basename(file),
    Total_Outcome_SNPs = length(outcome_snps),
    Overlapping_SNPs = length(shared_snps),
    Shared_SNPs = if (length(shared_snps) > 0) paste(shared_snps, collapse = "; ") else ""
  )
}

# === Combine results ===
result_table <- bind_rows(result_list)

# === Save table to CSV ===
write.csv(result_table, output_table_path, row.names = FALSE)

cat("Summary table saved to:", output_table_path, "\n")
