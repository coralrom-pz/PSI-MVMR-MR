# ================================
# SNP Checker + Proxy Replacement (single-job + batch)
# ================================
# Author: Coral Romualdo-Perez
# Date: 2025-06-19
# Notes:
# - Pseudo paths are used throughout (replace with your project layout).

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
})

# ----------------------------------------------------------------
# GLOBAL (pseudo) PATHS (edit to your repo layout)
# ----------------------------------------------------------------
MAP_FILE           <- "project/maps/harmonization_input_map.csv"
VERIFY_OUT_DIR     <- "project/outputs/snp_verification/SUBFOLDER_NAME"

# ================================
# CHUNK 1 — Verify SNVs present in outcome datasets
# ================================
dir.create(VERIFY_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- read job map (must have: Exposure_File, Outcome_File, Exposure_Name, Outcome_Name)
harmonization_map <- read.csv(MAP_FILE, stringsAsFactors = FALSE)

# ---- master summary collector
master_summary <- data.frame(
  Exposure = character(),
  Outcome = character(),
  `Original SNP` = character(),
  `effect_allele (Exposure)` = character(),
  `other_allele (Exposure)` = character(),
  `Present in Outcome (Y/N)` = character(),
  `Effect Allele Mismatch (Y/N)` = character(),
  `effect_allele (Outcome)` = character(),
  `other_allele (Outcome)` = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(harmonization_map))) {
  
  # ---- paths & labels
  exposure_path <- harmonization_map$Exposure_File[i]
  outcome_path  <- harmonization_map$Outcome_File[i]
  exposure_name <- harmonization_map$Exposure_Name[i]
  outcome_name  <- harmonization_map$Outcome_Name[i]
  cat("Running:", exposure_name, "to", outcome_name, "\n")
  
  if (!file.exists(exposure_path) || !file.exists(outcome_path)) {
    warning("Missing required files for ", exposure_name, " or ", outcome_name); next
  }
  
  # ---- load exposure
  exposure_df <- read.table(exposure_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
  colnames(exposure_df) <- gsub("^effect_allele.*$", "effect_allele", colnames(exposure_df))
  colnames(exposure_df) <- gsub("^other_allele.*$",  "other_allele",  colnames(exposure_df))
  
  if (!all(c("SNP","effect_allele","other_allele") %in% colnames(exposure_df))) {
    warning("Missing required columns in exposure file: ", exposure_name); next
  }
  
  exp_snps <- exposure_df %>%
    select(SNP, effect_allele, other_allele) %>%
    rename(exp_SNP = SNP, exp_effect = effect_allele, exp_other = other_allele)
  
  # ---- load outcome
  outcome_df <- read.table(outcome_path, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
  colnames(outcome_df) <- gsub("^effect_allele.*$", "effect_allele", colnames(outcome_df))
  colnames(outcome_df) <- gsub("^other_allele.*$",  "other_allele",  colnames(outcome_df))
  
  if (!all(c("SNP","effect_allele","other_allele") %in% colnames(outcome_df))) {
    warning("Missing required columns in outcome file: ", outcome_name); next
  }
  
  out_snps <- outcome_df %>%
    select(SNP, effect_allele, other_allele) %>%
    rename(out_SNP = SNP, out_effect = effect_allele, out_other = other_allele)
  
  # ---- merge & annotate
  merged <- exp_snps %>%
    left_join(out_snps, by = c("exp_SNP" = "out_SNP")) %>%
    mutate(
      Exposure = exposure_name,
      Outcome  = outcome_name,
      `Present in Outcome (Y/N)` =
        ifelse(is.na(out_effect), "NO", "YES"),
      `Effect Allele Mismatch (Y/N)` =
        ifelse(is.na(out_effect), NA,
               ifelse(exp_effect != out_effect, "YES", "NO"))
    ) %>%
    mutate_all(as.character) %>%
    select(
      Exposure, Outcome, exp_SNP, exp_effect, exp_other,
      `Present in Outcome (Y/N)`, `Effect Allele Mismatch (Y/N)`,
      out_effect, out_other
    ) %>%
    rename(
      `Original SNP`              = exp_SNP,
      `effect_allele (Exposure)`  = exp_effect,
      `other_allele (Exposure)`   = exp_other,
      `effect_allele (Outcome)`   = out_effect,
      `other_allele (Outcome)`    = out_other
    )
  
  # ---- write per-pair result
  per_pair_file <- file.path(
    VERIFY_OUT_DIR,
    paste0(exposure_name, "_to_", outcome_name, "_SNP_verification.csv")
  )
  write.csv(merged, per_pair_file, row.names = FALSE)
  cat("Saved:", per_pair_file, "\n")
  
  # ---- append to master summary
  master_summary <- bind_rows(master_summary, merged %>%
                                select(
                                  Exposure, Outcome, `Original SNP`,
                                  `effect_allele (Exposure)`, `other_allele (Exposure)`,
                                  `Present in Outcome (Y/N)`, `Effect Allele Mismatch (Y/N)`,
                                  `effect_allele (Outcome)`, `other_allele (Outcome)`
                                ))
}

# ---- write master summary
master_file <- file.path(VERIFY_OUT_DIR, "SNP_VERIFICATION_MASTER_SUMMARY.csv")
write.csv(master_summary, master_file, row.names = FALSE)
cat("\nMaster summary saved to:", master_file, "\n")


# ================================
# CHUNK 2 — Replace missing/mismatched SNVs with proxies (single job)
# ================================
# ---- single-job inputs (pseudo paths)
verification_file <- "project/outputs/snp_verification/Isolation_to_Depression_SNP_verification.csv"
raw_file          <- "project/inputs/exposures/isolation_data_exposure.txt"
outcome_file      <- "project/inputs/outcomes/Depression_GWAS_outcome.txt"
clumped_file      <- "project/derived/exposures/isolation_data_clumped_no_rs171697.txt"
proxy_ref_file    <- "project/derived/proxies/500kb_PROXY_REPORT_isolation_data_clumped_no_rs171697.csv"

output_file       <- "project/derived/exposures/isolation_data_with_proxies.txt"
proxy_log_file    <- "project/derived/exposures/isolation_proxy_substitution_log.csv"

# ---- load single-job data
verification <- read.csv(verification_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  rename(
    present  = `Present in Outcome (Y/N)`,
    mismatch = `Effect Allele Mismatch (Y/N)`,
    snp      = `Original SNP`
  )

raw      <- read.table(raw_file,     header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
outcome  <- read.table(outcome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
clumped  <- read.table(clumped_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
proxy_reference <- read.csv(proxy_ref_file, stringsAsFactors = FALSE)

# ---- targets to replace (missing OR allele-mismatch)
to_replace <- verification %>%
  filter(present == "NO" | mismatch == "YES") %>%
  pull(snp)

raw_exposure_snps <- raw     %>% select(SNP, effect_allele.exposure, other_allele.exposure)
outcome_snps      <- outcome %>% select(SNP, effect_allele.outcome,  other_allele.outcome)

# ---- pick best proxies by allele match in both raw exposure and outcome
usable_proxies <- proxy_reference %>%
  filter(original_snp %in% to_replace) %>%
  arrange(original_snp, desc(r2)) %>%
  group_by(original_snp) %>%
  group_modify(~{
    candidates <- .x
    for (j in seq_len(nrow(candidates))) {
      snp <- candidates$proxy_snp[j]
      exp_row <- raw_exposure_snps %>% filter(SNP == snp)
      out_row <- outcome_snps      %>% filter(SNP == snp)
      if (nrow(exp_row) > 0 && nrow(out_row) > 0 &&
          exp_row$effect_allele.exposure == out_row$effect_allele.outcome) {
        return(tibble(proxy_snp = snp, status = "Replaced with proxy"))
      }
    }
    tibble(proxy_snp = NA_character_, status = "Removed (no valid proxy)")
  }) %>%
  ungroup()

# ---- substitution log (one row per SNP needing action)
proxy_log <- tibble(original_snp = to_replace) %>%
  left_join(usable_proxies, by = "original_snp") %>%
  mutate(status = ifelse(is.na(proxy_snp), "Removed (no valid proxy)", status))

# ---- build proxy rows from RAW exposure (keep all available exposure columns)
proxy_rows <- proxy_log %>%
  filter(!is.na(proxy_snp)) %>%
  select(original_snp, proxy_snp) %>%
  left_join(raw, by = c("proxy_snp" = "SNP")) %>%
  mutate(SNP = proxy_snp) %>%
  select(-proxy_snp)

# ---- remove originals that were missing/mismatched; then add proxies
clumped_updated <- clumped %>%
  filter(!(SNP %in% to_replace)) %>%        # drop originals that needed replacement
  filter(!(SNP %in% proxy_log$original_snp))# safety: drop duplicates of originals if present

final_dataset <- bind_rows(clumped_updated, proxy_rows)

# ---- write outputs
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.table(final_dataset, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(proxy_log,       proxy_log_file, row.names = FALSE)
cat("Saved updated dataset to:", output_file, "\n")
cat("Saved proxy substitution log to:", proxy_log_file, "\n")


# ================================
# CHUNK 3 — Batch proxy replacement (strict 1-for-1)
# ================================
BATCH_MAP <- "project/maps/replacement_batch_map.csv"
batch_map <- read.csv(BATCH_MAP, stringsAsFactors = FALSE)

for (i in seq_len(nrow(batch_map))) {
  cat("\n=============================\n")
  cat("Processing batch", i, "\n")
  cat("=============================\n")
  
  tryCatch({
    # ---- per-batch paths
    verification_file <- batch_map$Verification_File[i]
    raw_file          <- batch_map$Raw_Exposure_File[i]
    outcome_file      <- batch_map$Outcome_File[i]
    clumped_file      <- batch_map$Clumped_File[i]
    proxy_ref_file    <- batch_map$Proxy_Ref_File[i]
    output_file       <- batch_map$Output_File[i]
    proxy_log_file    <- batch_map$Proxy_Log_File[i]
    
    # ---- backup original clumped
    clumped_backup_file <- sub("\\.txt$", "_UNMODIFIED.txt", output_file)
    dir.create(dirname(clumped_backup_file), showWarnings = FALSE, recursive = TRUE)
    file.copy(clumped_file, clumped_backup_file, overwrite = TRUE)
    cat("✓ Saved unmodified clumped file to:", clumped_backup_file, "\n")
    
    # ---- load batch data
    verification <- read.csv(verification_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
      rename(
        present  = `Present in Outcome (Y/N)`,
        mismatch = `Effect Allele Mismatch (Y/N)`,
        snp      = `Original SNP`
      )
    
    raw      <- read.table(raw_file,     header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    outcome  <- read.table(outcome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    clumped  <- read.table(clumped_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
    proxy_reference <- read.csv(proxy_ref_file, stringsAsFactors = FALSE)
    
    # ---- targets for replacement (batch version uses present == "NO")
    to_replace <- verification %>% filter(present == "NO") %>% pull(snp)
    
    raw_exposure_snps <- raw     %>% select(SNP, effect_allele.exposure, other_allele.exposure)
    outcome_snps      <- outcome %>% select(SNP, effect_allele.outcome,  other_allele.outcome)
    
    # ---- choose best proxy per original (alleles aligned)
    usable_proxies <- proxy_reference %>%
      filter(original_snp %in% to_replace) %>%
      arrange(original_snp, desc(r2)) %>%
      group_by(original_snp) %>%
      group_modify(~{
        proxies <- .x
        for (j in seq_len(nrow(proxies))) {
          snp <- proxies$proxy_snp[j]
          exp_row <- raw_exposure_snps %>% filter(SNP == snp)
          out_row <- outcome_snps      %>% filter(SNP == snp)
          if (nrow(exp_row) > 0 && nrow(out_row) > 0 &&
              exp_row$effect_allele.exposure == out_row$effect_allele.outcome) {
            return(tibble(proxy_snp = snp, status = "Replaced with proxy", reason = "Alleles aligned"))
          }
        }
        tibble(proxy_snp = NA_character_, status = "Removed (no valid proxy)", reason = "No valid proxy with allele match")
      }, .groups = "drop") %>%
      ungroup()
    
    # ---- detailed log for every SNP in clumped set
    all_log <- tibble(original_snp = clumped$SNP) %>%
      mutate(
        proxy_snp = NA_character_,
        status    = ifelse(original_snp %in% to_replace, "Replacement needed", "Present in outcome, not replaced"),
        reason    = ifelse(original_snp %in% to_replace, "Pending replacement",  "No replacement needed")
      ) %>%
      left_join(usable_proxies, by = "original_snp", suffix = c("", ".y")) %>%
      mutate(
        proxy_snp = coalesce(proxy_snp.y, proxy_snp),
        status    = coalesce(status.y,    status),
        reason    = coalesce(reason.y,    reason)
      ) %>%
      mutate(proxy_snp = ifelse(is.na(proxy_snp), "NA", proxy_snp)) %>%
      select(original_snp, proxy_snp, status, reason)
    
    # ---- build proxy rows from RAW exposure
    replacement_table <- usable_proxies %>% filter(!is.na(proxy_snp)) %>% select(original_snp, proxy_snp)
    
    proxy_rows <- replacement_table %>%
      left_join(raw, by = c("proxy_snp" = "SNP")) %>%
      mutate(SNP = proxy_snp) %>%
      select(-proxy_snp)
    
    # ---- remove originals that needed replacement; combine with proxies
    clumped_clean <- clumped %>% filter(!SNP %in% to_replace)
    final_dataset <- bind_rows(clumped_clean, proxy_rows)
    
    # ---- write outputs
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    write.table(final_dataset, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    write.csv(all_log,       proxy_log_file, row.names = FALSE)
    
    cat("✓ Saved updated dataset to:", output_file, "\n")
    cat("✓ Saved proxy substitution log to:", proxy_log_file, "\n")
    
  }, error = function(e) {
    cat("Error in batch", i, ":", conditionMessage(e), "\n")
  })
}


# ================================
# REQUIRED FILE HEADERS (for replicability)
# ================================
# 1) harmonization_input_map.csv  (MAP_FILE)
#    Columns (exact names):
#      - Exposure_File, Outcome_File, Exposure_Name, Outcome_Name
#
# 2) Exposure files (paths listed in harmonization map)
#    Must include at least:
#      - SNP, effect_allele, other_allele
#    Also accepted (auto-collapsed via regex to base names):
#      - effect_allele.exposure, other_allele.exposure
#    Optional (recommended for MR downstream):
#      - beta.exposure, se.exposure, pval.exposure, eaf.exposure, samplesize.exposure
#
# 3) Outcome files (paths listed in harmonization map)
#    Must include at least:
#      - SNP, effect_allele, other_allele
#    Also accepted:
#      - effect_allele.outcome, other_allele.outcome
#    Optional:
#      - beta.outcome, se.outcome, pval.outcome, samplesize.outcome
#
# 4) Proxy reference file (used in replacement steps)
#    Must include:
#      - original_snp, proxy_snp, r2
#    Optional extras (ignored here): dist, r, etc.
#
# 5) Clumped exposure file (inputs to replacement)
#    Must include:
#      - SNP
#    Recommended for clean merges/logging:
#      - effect_allele.exposure, other_allele.exposure
#      - beta.exposure, se.exposure, pval.exposure


# ================================
# Proxy SNP Pipeline — LDlinkR wrappers + batch (500 kb window)
# ================================
# Author: Coral Romualdo-Perez
# Date: 2025-06-05
# Notes:
# - Pseudo paths are used (replace with your repo layout).
# - Comment style matches prior scripts (banners + subheaders).
# - IMPORTANT: set your LDlink token below
# ================================

suppressPackageStartupMessages({
  library(dplyr); library(jsonlite); library(httr); library(LDlinkR)
})

# ----------------------------------------------------------------
# GLOBAL SETTINGS (edit)
# ----------------------------------------------------------------
# LDlink token (https://ldlink.nih.gov/?tab=apiaccess)
Sys.setenv(LDLINK_TOKEN = "YOUR_TOKEN_ID")

# Output folder for proxy-augmented exposures + reports
PROXY_OUT_DIR <- "project/derived/proxies/500kb"
dir.create(PROXY_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Optional: workaround for occasional NIH SSL updates
httr::set_config(httr::config(ssl_verifypeer = 0))

# ================================
# CHUNK 1 — Function: add_proxies_to_exposure()
# ================================
# Input:
#   exposure_file : path to a tab-delimited exposure file with at least column 'SNP'
#   r2_threshold  : minimum LD r^2 to accept as a proxy (default 0.8)
# Output:
#   - Augmented exposure written to PROXY_OUT_DIR/500kb_<basename>
#   - Proxy report written to PROXY_OUT_DIR/500kb_PROXY_REPORT_<basename>.csv
# Behavior:
#   - For each SNP, queries LDproxy() (EUR, GRCh37, ±500 kb), drops self row,
#     filters by R2 >= threshold, keeps top 10 rows (as in working code),
#     and appends proxy rows to the exposure table (then de-duplicates by SNP).
add_proxies_to_exposure <- function(exposure_file, r2_threshold = 0.8) {
  
  # ---- load exposure
  exposure_data <- read.table(
    exposure_file, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!"SNP" %in% names(exposure_data)) {
    stop("Exposure file missing required column: 'SNP' → ", exposure_file)
  }
  exposure_snps <- exposure_data$SNP
  cat("Found", length(exposure_snps), "exposure SNPs to query for proxies.\n")
  
  # ---- results container
  proxy_results <- data.frame(
    original_snp = character(),
    proxy_snp    = character(),
    r2           = numeric(),
    stringsAsFactors = FALSE
  )
  
  # ---- loop over SNPs
  for (snp in exposure_snps) {
    cat("\n Querying proxies for:", snp, "\n")
    
    proxies <- tryCatch({
      LDlinkR::LDproxy(
        snp          = snp,
        pop          = "EUR",
        r2d          = "r2",
        token        = Sys.getenv("LDLINK_TOKEN"),
        file         = FALSE,
        genome_build = "grch37",
        win_size     = "500000"
      )
    }, error = function(e) {
      cat("LDproxy() error for SNP:", snp, "-", conditionMessage(e), "\n")
      return(NULL)
    })
    
    # ---- process response
    if (!is.null(proxies) && all(c("RS_Number","R2") %in% colnames(proxies))) {
      
      # STEP 1: remove self row
      proxies_no_self <- dplyr::filter(proxies, RS_Number != snp)
      
      # STEP 2: keep R2 ≥ threshold
      proxies_filtered <- dplyr::filter(proxies_no_self, R2 >= r2_threshold)
      
      if (nrow(proxies_filtered) > 0) {
        
        # STEP 3: take top 10 by R2 (original behavior)
        top_proxy <- proxies_filtered[order(-proxies_filtered$R2), ][1:min(10, nrow(proxies_filtered)), ]
        
        # append (original_snp will recycle to match row count)
        proxy_results <- rbind(
          proxy_results,
          data.frame(
            original_snp = snp,
            proxy_snp    = top_proxy$RS_Number,
            r2           = top_proxy$R2,
            stringsAsFactors = FALSE
          )
        )
        
        cat("Top proxy rows for", snp, "→ n =", nrow(top_proxy), "; best R2 =", max(top_proxy$R2), "\n")
        
      } else {
        cat("No proxies ≥ R2 threshold for", snp, "\n")
      }
      
    } else {
      cat("No proxies returned for", snp, "\n")
    }
  } # end SNP loop
  
  # ---- build adjusted exposure (append proxy rows using same columns)
  adjusted_data <- exposure_data
  
  if (nrow(proxy_results) > 0) {
    for (i in seq_len(nrow(proxy_results))) {
      proxy_row    <- proxy_results[i, ]
      proxy_snp    <- proxy_row$proxy_snp
      original_snp <- proxy_row$original_snp
      
      # clone the original SNP row, swap SNP → proxy_snp
      proxy_entry <- dplyr::filter(adjusted_data, SNP == original_snp)
      if (nrow(proxy_entry) > 0) {
        proxy_entry$SNP <- proxy_snp
        adjusted_data   <- rbind(adjusted_data, proxy_entry)
      }
    }
  }
  
  # ---- de-duplicate by SNP (keep first occurrence)
  adjusted_data <- dplyr::distinct(adjusted_data, SNP, .keep_all = TRUE)
  
  # ---- write outputs
  exposure_filename <- basename(exposure_file)
  adjusted_path <- file.path(PROXY_OUT_DIR, paste0("500kb_", exposure_filename))
  report_path   <- file.path(PROXY_OUT_DIR, paste0("500kb_PROXY_REPORT_", exposure_filename, ".csv"))
  
  write.table(adjusted_data, adjusted_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.csv(proxy_results,   report_path,   row.names = FALSE)
  
  cat("\n Adjusted exposure saved:", adjusted_path, "\n")
  cat(" Proxy report saved:",       report_path,   "\n\n")
  
  invisible(list(adjusted = adjusted_path, report = report_path))
}


# ================================
# CHUNK 2 — Batch runner (control center)
# ================================
# Provide a vector of exposure files to process.
# Replace with your own (TwoSampleMR-ready, clumped or FINAL instrument sets).
exposure_files <- c(
  # Examples (comment/uncomment as needed):
  # "project/derived/instruments/Folkersen/Clumped/IL1RA_Cis-Filtered_clumped_p1e5.txt",
  # "project/derived/instruments/Folkersen/Clumped/IL6RA_Cis-Filtered_clumped_p1e5.txt",
  # UKB-PPP pQTL FINAL instruments (examples):
  "project/derived/instruments/UKB_PPP/MR_instruments_FINAL/IL1B_MR_instruments_FINAL.txt",
  "project/derived/instruments/UKB_PPP/MR_instruments_FINAL/IL1RA_MR_instruments_FINAL.txt",
  "project/derived/instruments/UKB_PPP/MR_instruments_FINAL/IL6R_MR_instruments_FINAL.txt",
  "project/derived/instruments/UKB_PPP/MR_instruments_FINAL/TNFA_MR_instruments_FINAL.txt",
  "project/derived/instruments/UKB_PPP/MR_instruments_FINAL/TNFR1_MR_instruments_FINAL.txt",
  "project/derived/instruments/UKB_PPP/MR_instruments_FINAL/TNFR2_MR_instruments_FINAL.txt"
)

for (i in seq_along(exposure_files)) {
  f <- exposure_files[i]
  cat("\n========== [", i, "/", length(exposure_files), "] Processing:", f, "==========\n")
  tryCatch({
    add_proxies_to_exposure(f, r2_threshold = 0.8)
  }, error = function(e) {
    cat("ERROR while processing", f, "→", conditionMessage(e), "\n")
  })
  cat("========== Done:", f, "==========\n\n")
}

cat("Batch complete. All outputs in:", PROXY_OUT_DIR, "\n")


# ================================
# REFERENCE — Original LDlinkR call (single SNP), for record
# ================================
# This mirrors your prior example; use the token from the environment.
# LDlink docs: https://ldlink.nih.gov/
# (Do not hardcode real tokens in shared repos.)

# LDlinkR::LDproxy(
#   snp          = "rs6430286",
#   pop          = "EUR",
#   r2d          = "r2",
#   token        = Sys.getenv("LDLINK_TOKEN"),
#   file         = FALSE,
#   genome_build = "grch37",
#   win_size     = "500000",
#   api_root     = "https://ldlink.nih.gov/LDlinkRest"
# )


# ================================
# REQUIRED INPUT SCHEMA (minimal)
# ================================
# Exposure file(s):
#   - Required: SNP
#   - Recommended (carried through to proxies): effect_allele.exposure, other_allele.exposure,
#     beta.exposure, se.exposure, pval.exposure, eaf.exposure, samplesize.exposure
#
# Output files:
#   - Augmented exposure: tab-delimited, same columns as input, SNP column includes proxies
#   - Proxy report CSV: columns original_snp, proxy_snp, r2 (up to 10 proxies/SNP retained)
