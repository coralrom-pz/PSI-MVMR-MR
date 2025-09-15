# ================================
# MVMR Pipeline
# ================================
# Author: Coral Romualdo-Perez
# Date: 2025-06-16
# Notes:
# - Pseudo paths are used below (swap to your repo layout as needed).

# ================================


# ================================
# CHUNK 1 — Load libraries (run first)
# ================================
knitr::opts_chunk$set(echo = TRUE)
library(TwoSampleMR)
library(MVMR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(readxl)
library(readr)


# ================================
# CHUNK 2 — MVMR pipeline function (run second)
# ================================
run_mvmr_pipeline <- function(exposure1_file, exposure2_file, outcome_file,
                              exposure1_name, exposure2_name, outcome_name,
                              output_dir, subset_exposure2 = TRUE,
                              save_merged_input = TRUE) {
  
  cat("\n=====================================\n")
  cat("Running MVMR pipeline for:", exposure1_name, "+", exposure2_name, "=>", outcome_name, "\n")
  cat("=====================================\n")
  
  # Read exposures and outcome
  exposure1 <- read.table(exposure1_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  exposure2_raw <- read.table(exposure2_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  outcome   <- read.table(outcome_file,   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  cat("[Diagnostics] Raw SNP counts:\n")
  cat("Exposure1 SNPs:", nrow(exposure1), "\n")
  cat("Exposure2 SNPs (pre-subset):", nrow(exposure2_raw), "\n")
  
  # Subset exposure2 to exposure1 SNPs (Burgess-style)
  if (subset_exposure2) {
    exposure2 <- exposure2_raw[exposure2_raw$SNP %in% exposure1$SNP, ]
    cat("Subset", exposure2_name, "to SNPs present in", exposure1_name, "\n")
    cat("Exposure2 SNPs after subsetting:", nrow(exposure2), "\n")
  } else {
    exposure2 <- exposure2_raw
  }
  
  # Merge exposures
  merged_exposures <- merge(exposure1, exposure2, by = "SNP")
  cat("SNPs after merging exposures:", nrow(merged_exposures), "\n")
  
  # Merge with outcome
  merged_mvmr <- merge(merged_exposures, outcome, by = "SNP")
  n_snps <- nrow(merged_mvmr)
  cat("SNPs after merging with outcome:", n_snps, "\n")
  
  # Early exit if not enough SNPs
  if (n_snps < 3) {
    cat("Skipping: insufficient overlapping SNPs (", n_snps, ")\n")
    return(list(
      mvmr_results = NULL, n_snps = n_snps, status = "Skipped",
      fstat = NULL, pleio_q = NA, pleio_df = NA, pleio_p = NA
    ))
  }
  
  # Security checks for missing values
  if (any(is.na(merged_mvmr$beta.exposure.x)) || any(is.na(merged_mvmr$beta.exposure.y))) {
    stop("ERROR: Missing beta values in exposures after merging!")
  }
  if (any(is.na(merged_mvmr$beta.outcome))) {
    stop("ERROR: Missing beta values in outcome after merging!")
  }
  
  # Rename columns for MVMR package
  colnames(merged_mvmr)[colnames(merged_mvmr) == "beta.exposure.x"] <- "beta.exp1"
  colnames(merged_mvmr)[colnames(merged_mvmr) == "se.exposure.x"]   <- "se.exp1"
  colnames(merged_mvmr)[colnames(merged_mvmr) == "beta.exposure.y"] <- "beta.exp2"
  colnames(merged_mvmr)[colnames(merged_mvmr) == "se.exposure.y"]   <- "se.exp2"
  colnames(merged_mvmr)[colnames(merged_mvmr) == "beta.outcome"]    <- "beta.out"
  colnames(merged_mvmr)[colnames(merged_mvmr) == "se.outcome"]      <- "se.out"
  
  # Make output folder
  folder <- file.path(output_dir, paste0(exposure1_name, "_", exposure2_name, "_to_", outcome_name))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  
  # Save merged input file for verification
  if (save_merged_input) {
    write.table(merged_mvmr, file.path(folder, "Merged_Input_Data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Saved merged input data to:", file.path(folder, "Merged_Input_Data.tsv"), "\n")
  }
  
  # Format for MVMR
  BXGs <- cbind(merged_mvmr$beta.exp1, merged_mvmr$beta.exp2)
  seBXGs <- cbind(merged_mvmr$se.exp1, merged_mvmr$se.exp2)
  BYG   <- merged_mvmr$beta.out
  seBYG <- merged_mvmr$se.out
  RSID  <- merged_mvmr$SNP
  
  formatted_data <- format_mvmr(BXGs = BXGs, BYG = BYG, seBXGs = seBXGs, seBYG = seBYG, RSID = RSID)
  
  # Run IVW
  ivw_res <- ivw_mvmr(formatted_data)
  ivw_res <- as.data.frame(ivw_res)
  
  # Guard against empty IVW results
  if (is.null(ivw_res) || nrow(ivw_res) == 0) {
    cat("IVW result is empty or invalid — skipping this analysis.\n")
    return(list(mvmr_results = NULL, n_snps = n_snps, status = "No Results", fstat = NULL))
  }
  
  # Run strength diagnostic
  fstat <- tryCatch({
    out <- strength_mvmr(formatted_data, gencov = 0)
    if (is.null(out) || !is.data.frame(out) || nrow(out) == 0) {
      cat("Warning: strength_mvmr returned empty or invalid output\n")
      NULL
    } else {
      out
    }
  }, error = function(e) {
    cat("Error in strength_mvmr():", e$message, "\n")
    NULL
  })
  
  # Run pleiotropy diagnostic
  pleio <- tryCatch({
    out <- pleiotropy_mvmr(formatted_data, gencov = 0)
    if (is.null(out) || !is.data.frame(out) || nrow(out) == 0) {
      cat("Warning: pleiotropy_mvmr returned empty or invalid output\n")
      NULL
    } else {
      out
    }
  }, error = function(e) {
    cat("Error in pleiotropy_mvmr():", e$message, "\n")
    NULL
  })
  
  # Run heterogeneity diagnostic (qhet_mvmr)
  if (n_snps >= 5) {
    qhet <- tryCatch({
      out <- qhet_mvmr(formatted_data, pcor = 0)
      if (is.null(out) || !is.data.frame(out) || nrow(out) == 0) {
        cat("Warning: qhet_mvmr returned empty or invalid output\n")
        NULL
      } else {
        out
      }
    }, error = function(e) {
      cat("Error in qhet_mvmr():", e$message, "\n")
      NULL
    })
  } else {
    qhet <- NULL
    cat("Skipping qhet_mvmr: too few SNPs (", n_snps, ")\n")
  }
  
  # Save outputs
  if (!is.null(ivw_res) && nrow(ivw_res) > 0) {
    write.csv(ivw_res, file.path(folder, "MVMR_IVW_Results.csv"), row.names = FALSE)
    cat("Saved IVW results.\n")
  }
  
  if (!is.null(fstat) && nrow(fstat) > 0) {
    write.csv(fstat, file.path(folder, "MVMR_Strength.csv"), row.names = FALSE)
    cat("Saved strength diagnostics.\n")
  }
  
  if (!is.null(pleio) && nrow(pleio) > 0) {
    write.csv(pleio, file.path(folder, "MVMR_Pleiotropy.csv"), row.names = FALSE)
    cat("Saved pleiotropy diagnostics.\n")
  } else {
    cat("pleiotropy_mvmr() returned NULL or empty — skipping Pleiotropy output.\n")
  }
  
  cat("Saved all available results for:", exposure1_name, "+", exposure2_name, "to", outcome_name, "| SNPs:", n_snps, "\n")
  
  # Return core outputs
  return(list(
    mvmr_results = ivw_res,
    n_snps = n_snps,
    status = "Success",
    fstat = fstat,
    pleio_q      = if (!is.null(pleio) && "Q_statistic" %in% colnames(pleio)) pleio$Q_statistic[1] else NA,
    pleio_df     = if (!is.null(pleio) && "Q_df" %in% colnames(pleio)) pleio$Q_df[1] else NA,
    pleio_p      = if (!is.null(pleio) && "p-value" %in% colnames(pleio)) pleio[["p-value"]][1] else NA,
    qhet_q   = if (!is.null(qhet) && "Q" %in% colnames(qhet)) qhet$Q[1] else NA,
    qhet_df  = if (!is.null(qhet) && "df" %in% colnames(qhet)) qhet$df[1] else NA,
    qhet_p   = if (!is.null(qhet) && ("p-value" %in% colnames(qhet) || "pvalue" %in% colnames(qhet))) {
      if ("p-value" %in% colnames(qhet)) qhet$`p-value`[1] else qhet$pvalue[1]
    } else NA
  ))
}


# ================================
# CHUNK 3 — Batch MVMR runner
# ================================
# Required headers in the batch CSV:
# Exposure1_File, Exposure2_File, Outcome_File, Exposure1_Name, Exposure2_Name, Outcome_Name
# Specify file path beneath each header: 

# Define input (pseudo paths)
mvmr_df   <- read.csv("project/config/MVMR_input.csv", stringsAsFactors = FALSE)
output_dir <- "project/outputs/mvmr/FOLDER/NAME"

# Optional: Create persistent log file
log_file <- file.path(output_dir, "MVMR_Batch_Run_Log.txt")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("", file = log_file)  # Clear previous log

logcat <- function(...) {
  cat(..., "\n")
  cat(..., "\n", file = log_file, append = TRUE)
}

# Initialize batch summary
batch_summary_mvmr <- data.frame(
  Exposure1 = character(),
  Exposure2 = character(),
  Outcome   = character(),
  Exposure  = character(),
  IVW_Beta  = numeric(),
  IVW_SE    = numeric(),
  IVW_Pval  = numeric(),
  Fstat     = numeric(),
  Pleiotropy_Q = numeric(),
  Pleiotropy_DF = numeric(),
  Pleiotropy_Pval = numeric(),
  qhet_Q = numeric(),
  qhet_DF = numeric(),
  qhet_Pval = numeric(),
  N_SNPs    = numeric(),
  Status    = character(),
  StartTime = character(),
  EndTime   = character(),
  stringsAsFactors = FALSE
)

# Start timer
start_time <- Sys.time()
logcat("\n==========================")
logcat("Starting full MVMR batch run!")
logcat("Start time:", format(start_time))
logcat("==========================\n")

# Loop through rows
for (i in seq_len(nrow(mvmr_df))) {
  exp1 <- mvmr_df$Exposure1_Name[i]
  exp2 <- mvmr_df$Exposure2_Name[i]
  outc <- mvmr_df$Outcome_Name[i]
  
  logcat("\n========== [", i, "/", nrow(mvmr_df), "]", exp1, "+", exp2, "=>", outc, "==========")
  
  time_analysis_start <- Sys.time()
  
  result <- tryCatch({
    # Count SNPs *before* running pipeline for better logging
    exposure1 <- read.table(mvmr_df$Exposure1_File[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    exposure2 <- read.table(mvmr_df$Exposure2_File[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    outcome   <- read.table(mvmr_df$Outcome_File[i],   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    merged_exposures <- merge(exposure1, exposure2, by = "SNP")
    merged_mvmr <- merge(merged_exposures, outcome, by = "SNP")
    n_snps <- nrow(merged_mvmr)
    
    logcat("Initial Exposure1 SNPs:", nrow(exposure1))
    logcat("Initial Exposure2 SNPs:", nrow(exposure2))
    logcat("Merged exposures SNPs:", nrow(merged_exposures))
    logcat("Merged with outcome SNPs:", n_snps)
    if (n_snps < 3) {
      logcat("Skipping due to insufficient overlapping SNPs:", n_snps)
      list(
        mvmr_results = NULL,
        n_snps = n_snps,
        status = "Skipped",
        fstat = NULL,
        pleio_q = NA,
        pleio_df = NA,
        pleio_p = NA,
        qhet_q = NA,
        qhet_df = NA,
        qhet_p = NA
      )
    } else {
      run_mvmr_pipeline(
        exposure1_file = mvmr_df$Exposure1_File[i],
        exposure2_file = mvmr_df$Exposure2_File[i],
        outcome_file   = mvmr_df$Outcome_File[i],
        exposure1_name = exp1,
        exposure2_name = exp2,
        outcome_name   = outc,
        output_dir     = output_dir
      )
    }
    
  }, error = function(e) {
    logcat("ERROR in MVMR run:", exp1, "+", exp2, "=>", outc)
    logcat("Reason:", e$message)
    list(
      mvmr_results = NULL,
      n_snps = NA,
      status = "Error",
      fstat = NULL,
      pleio_q = NA,
      pleio_df = NA,
      pleio_p = NA
    )
  })
  
  time_analysis_end <- Sys.time()
  elapsed_mins <- round(difftime(time_analysis_end, time_analysis_start, units = "mins"), 2)
  logcat("Time elapsed:", elapsed_mins, "minutes")
  
  # Append to batch summary
  if (!is.null(result$mvmr_results) && is.data.frame(result$mvmr_results) && nrow(result$mvmr_results) >= 2) {
    # Print columns for debugging
    logcat("IVW result columns:", paste(colnames(result$mvmr_results), collapse = ", "))
    
    for (j in 1:nrow(result$mvmr_results)) {
      batch_summary_mvmr <- bind_rows(batch_summary_mvmr, data.frame(
        Exposure1 = exp1,
        Exposure2 = exp2,
        Outcome   = outc,
        Exposure  = ifelse(j == 1, exp1, exp2),
        IVW_Beta  = if ("Estimate" %in% colnames(result$mvmr_results) && length(result$mvmr_results$Estimate) >= j) result$mvmr_results$Estimate[j] else NA,
        IVW_SE    = if ("Std. Error" %in% colnames(result$mvmr_results) && length(result$mvmr_results$`Std. Error`) >= j) result$mvmr_results$`Std. Error`[j] else NA,
        IVW_Pval  = if ("Pr(>|t|)" %in% colnames(result$mvmr_results) && length(result$mvmr_results$`Pr(>|t|)`) >= j) {
          result$mvmr_results$`Pr(>|t|)`[j]
        } else NA,
        Fstat     = if (!is.null(result$fstat) && is.data.frame(result$fstat) && ncol(result$fstat) >= j) {
          as.numeric(result$fstat[1, j])
        } else NA,
        Pleiotropy_Q = result$pleio_q,
        Pleiotropy_DF = result$pleio_df,
        Pleiotropy_Pval = result$pleio_p,
        qhet_Q = result$qhet_q,
        qhet_DF = result$qhet_df,
        qhet_Pval = result$qhet_p,
        N_SNPs    = result$n_snps,
        Status    = result$status,
        StartTime = format(time_analysis_start),
        EndTime   = format(time_analysis_end),
        stringsAsFactors = FALSE
      ))
    }
  } else {
    batch_summary_mvmr <- bind_rows(batch_summary_mvmr, data.frame(
      Exposure1 = exp1,
      Exposure2 = exp2,
      Outcome   = outc,
      Exposure  = NA,
      IVW_Beta  = NA,
      IVW_SE    = NA,
      IVW_Pval  = NA,
      Fstat     = NA,
      Pleiotropy_Q = result$pleio_q,
      Pleiotropy_DF = result$pleio_df,
      Pleiotropy_Pval = result$pleio_p,
      qhet_Q = result$qhet_q,
      qhet_DF = result$qhet_df,
      qhet_Pval = result$qhet_p,
      N_SNPs    = result$n_snps,
      Status    = result$status,
      StartTime = format(time_analysis_start),
      EndTime   = format(time_analysis_end),
      stringsAsFactors = FALSE
    ))
  }
}

# End timer
end_time <- Sys.time()
logcat("\n==========================")
logcat("Full MVMR batch run complete!")
logcat("End time:", format(end_time))
logcat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes")
logcat("==========================\n")

# Save batch summary
summary_path <- file.path(output_dir, "MVMR_BATCH_SUMMARY.csv")
tryCatch({
  write.csv(batch_summary_mvmr, summary_path, row.names = FALSE)
  logcat("Saved batch summary to:", summary_path)
}, error = function(e) {
  logcat("ERROR saving batch summary:", e$message)
})