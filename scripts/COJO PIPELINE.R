
# =================================================================================================
# COJO PIPELINE
# -------------------------------------------------------------------------------------------------
# What this contains
# - Clear chunk headers & comments, matching the earlier style.
#
# Inputs you must provide (examples below use pseudo paths):
# - Clumped lead SNPs per exposure (TwoSampleMR format like *_clumped_p1e5.txt)
# - Raw GWAS (e.g., Loneliness_EurRel_Imputed.txt) for region extraction / enrichment
# - COJO outputs (*.jma.cojo / *.cma.cojo / *.given.cojo) generated via terminal runs
#
# Folder layout used here (edit to your repo if needed):
#   project/
#   ├─ data/
#   │  ├─ gwas/
#   │  │  ├─ raw/                     # full GWAS summary stats
#   │  │  └─ clumped/                 # clumped lead SNP files
#   │  ├─ eqtl/
#   │  │  └─ cis/formatted/           # eQTL exposure/outcome formatted files
#   │  └─ pqtl/
#   │     └─ 500KB_Formatted/         # cis-500kb formatted pQTLs (already windowed)
#   └─ outputs/
#      └─ cojo/
#         ├─ cojo_regions/500kb/      # per-SNP region files for COJO
#         │  └─ pQTLs/
#         ├─ cojo_results/500kb/      # COJO terminal outputs (*.jma.cojo etc.)
#         ├─ cojo_cond_results/       # COJO conditional outputs (*.given/cma)
#         ├─ cojo_slct_results/       # COJO -cojo-slct outputs (*.jma.cojo)
#         ├─ final_instruments/       # merged lead + COJO SNPs (per exposure)
#         ├─ final_instruments_cleaned/  # TwoSampleMR-ready exposures
#         └─ combined/                # combined files (txt/xlsx) for MR
#
# Notes
#   • Sample sizes (N) are left as in your working code. Adjust those constants where noted.
#   • Required columns are unchanged from your originals; headers are preserved where needed.
#   • You’ll still run COJO in the terminal (plink/gcta) between chunks; this R code formats IO.
# =================================================================================================

# ================================
# CHUNK 1 — Load packages (silent)
# ================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tools)
  library(openxlsx); library(purrr); library(TwoSampleMR)
})

# ================================
# CHUNK 1A — Select COJO regions (cis window per lead SNP) when CIS filtering not done
# ================================

# File Paths (EDIT as needed)
lead_snps_file <- "project/data/gwas/clumped/isolation_data_clumped_no_rs171697.txt"  # clumped file
gwas_file      <- "project/data/gwas/raw/Loneliness_EurRel_Imputed.txt"              # raw GWAS
output_dir     <- "project/outputs/cojo/cojo_regions/500kb"

# Create output folder
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read in files
lead_snps <- read_tsv(lead_snps_file, col_types = cols())
gwas <- read_table(gwas_file, col_types = cols())

# Prepare GWAS for COJO format
gwas_clean <- gwas %>%
  transmute(
    SNP = SNP,
    CHR = CHR,
    BP  = BP,
    A1  = ALLELE1,
    A2  = ALLELE0,
    freq = A1FREQ,
    b   = BETA,
    se  = SE,
    p   = P_LINREG,
    N   = 452302  # Replace number as needed
  )

# Merge lead SNPs to get their positions
lead_with_pos <- lead_snps %>%
  inner_join(gwas_clean %>% select(SNP, CHR, BP), by = "SNP")

# Extract kb windows, adjust from 300000 to 500000 to 1,000,000 if needed
window <- 1000000

for (i in seq_len(nrow(lead_with_pos))) {
  snp <- lead_with_pos$SNP[i]
  chr <- lead_with_pos$CHR[i]
  bp  <- lead_with_pos$BP[i]
  
  region_data <- gwas_clean %>%
    filter(CHR == chr, BP >= (bp - window), BP <= (bp + window))
  
  outfile <- file.path(output_dir, paste0(snp, "_region_500kb.txt"))
  write.table(region_data, file = outfile, sep = " ", quote = FALSE, row.names = FALSE)
}


# ================================
# CHUNK 1B — Batch-format COJO regions when files are already cis-filtered (pQTLs)
# ================================

# SETTINGS
input_dir   <- "project/data/pqtl/500KB_Formatted"
output_dir  <- "project/outputs/cojo/cojo_regions/500kb/pQTLs"
sample_size <- 21758  # adjust if different

# Create output folder if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Get list of all formatted exposure files
exposure_files <- list.files(
  path = input_dir,
  pattern = "_cis500kb_formatted.txt$",
  full.names = TRUE
)

# Loop through each file
for (exposure_file in exposure_files) {
  cat("\n========== Now processing:", basename(exposure_file), "==========\n")
  
  # Load
  exposure_raw <- read.table(
    exposure_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  # Check required columns
  required_cols <- c(
    "SNP",
    "beta.exposure",
    "se.exposure",
    "effect_allele.exposure",
    "other_allele.exposure",
    "pval.exposure",
    "eaf.exposure"
  )
  
  if (!all(required_cols %in% names(exposure_raw))) {
    cat("Required columns missing in exposure file:", exposure_file, "\n")
    next
  }
  
  # Map to COJO format
  exposure_clean <- exposure_raw %>%
    transmute(
      SNP = SNP,
      A1  = effect_allele.exposure,
      A2  = other_allele.exposure,
      freq = eaf.exposure,
      b    = beta.exposure,
      se   = se.exposure,
      p    = pval.exposure,
      N    = sample_size
    )
  
  # Define output filename
  out_basename <- paste0(
    file_path_sans_ext(basename(exposure_file)),
    "_region_500kb.txt"
  )
  outfile <- file.path(output_dir, out_basename)
  
  # Save COJO-ready region file
  write.table(
    exposure_clean,
    file = outfile,
    sep = " ",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  cat("Saved COJO region file to:", outfile, "\n")
}

cat("\nAll exposure files have been converted to COJO region format!\n")


# ================================
# CHUNK 2 — Remove headers for terminal COJO usage (space/tab formatting)
# ================================

## Set your region folder and sample size
region_dir   <- "project/outputs/cojo/cojo_regions/500kb/pQTLs"
sample_size  <- 21758  # Adjust if different

# List of your region files, adjust as needed
region_files <- c(
  "IL1RA_cis500kb_formatted_region_500kb.txt", "IL6RA_cis500kb_formatted_region_500kb.txt",
  "TNFRSF1A_cis500kb_formatted_region_500kb.txt", "TNFRSF1B_cis500kb_formatted_region_500kb.txt",
  "UPAR_cis500kb_formatted_region_500kb.txt"
)

# Loop through each file and reformat
for (file in region_files) {
  path <- file.path(region_dir, file)
  
  # Read file
  df <- read.table(path, header = TRUE)
  
  # Removes selected column names in the header
  df_cojo <- df[, c("SNP", "A1", "A2", "freq", "b", "se", "p")]
  
  # Add sample size column
  df_cojo$N <- sample_size
  
  # Overwrite the original file (no headers, no row names)
  write.table(df_cojo, path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  message("Reformatted ", file)
}

# ================================
# CHUNK 3.1 — Combine lead SNPs + COJO .jma.cojo signals and export (txt + xlsx)
# ================================

# 1. Load your lead SNPs
lead_snps <- readr::read_tsv("project/data/gwas/clumped/gwas_data_clumped.txt")

# 2. Prepare the lead SNPs in COJO-compatible format
lead_snps_clean <- lead_snps %>%
  transmute(
    SNP,
    A1 = `effect_allele.exposure`,
    A2 = `other_allele.exposure`,
    freq = `eaf.exposure`,
    b = `beta.exposure`,
    se = `se.exposure`,
    p = `pval.exposure`,
    N = 452302,
    source = "lead"
  )

# 2. Define COJO results directory
cojo_dir <- "project/outputs/cojo/cojo_results/500kb" #adjust as needed

# 3. Read .jma.cojo files with explicit types (these are any identified independent signals)
jma_files <- list.files(cojo_dir, pattern = "*.jma.cojo", full.names = TRUE)

read_jma <- function(file) {
  read_table(
    file,
    col_types = cols(
      Chr = col_integer(),
      SNP = col_character(),
      bp = col_integer(),
      refA = col_character(),
      freq = col_double(),
      b = col_double(),
      se = col_double(),
      p = col_double(),
      n = col_double(),
      freq_geno = col_double(),
      bJ = col_double(),
      bJ_se = col_double(),
      pJ = col_double(),
      LD_r = col_double()
    )
  ) %>%
    transmute(
      SNP,
      A1 = refA,
      A2 = NA_character_,
      freq,
      b,
      se,
      p,
      N = n,
      source = "cojo_jma"
    )
}

jma_combined <- purrr::map_dfr(jma_files, read_jma)

# 4. Combine lead SNPs and COJO-selected SNPs only
combined_snps <- bind_rows(lead_snps_clean, jma_combined) %>%
  distinct(SNP, .keep_all = TRUE)  # remove duplicates

# 6. Export for TwoSampleMR
write.table(
  combined_snps %>% select(SNP, A1, A2, freq, b, se, p, N),
  file = "project/outputs/cojo/combined_snps_500kb.txt",
  quote = FALSE, sep = "\t", row.names = FALSE
)

# 7. Export Excel file with full metadata
write.xlsx(
  combined_snps,
  file = "project/outputs/cojo/combined_snps_500kb.xlsx",
  rowNames = FALSE
)

cat("Output saved:\n- project/outputs/cojo/combined_snps_500kb.txt\n- project/outputs/cojo/combined_snps_500kb.xlsx\n")


# ================================
# CHUNK 3.2 — Batch combine clumped lead SNPs + COJO slct SNPs (per exposure; save txt)
# ================================

## For EACH exposure, save its own combined file
##Skip to 3.3 for MR inputs then run Chunk 6

# === CONFIG ===
clumped_dir <- "project/data/eqtl/cis/formatted"
cojo_dir <- "project/outputs/cojo/cojo_slct_results"
out_dir <- "project/outputs/cojo/final_instruments"

# Create output folder if needed
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# === EXPOSURE NAMES ===
exposures <- c("IL-1RN", "IL-1β", "IL-6", "IL-6R", "TNF-α", "TNF-R1", "TNF-R2", "PLAUR")

# === Fix for reading COJO with "nan" in n ===
read_cojo_jma <- function(exposure) {
  cojo_file <- file.path(cojo_dir, paste0(exposure, "_slct.jma.cojo"))
  if (!file.exists(cojo_file)) {
    warning("Missing COJO jma.cojo file for: ", exposure)
    return(NULL)
  }
  
  read_table(
    cojo_file,
    na = c("", "NA", "nan", "NaN"),
    show_col_types = FALSE,
    col_types = cols(
      Chr = col_integer(),
      SNP = col_character(),
      bp = col_integer(),
      refA = col_character(),
      freq = col_double(),
      b = col_double(),
      se = col_double(),
      p = col_double(),
      n = col_number(),
      freq_geno = col_double(),
      bJ = col_double(),
      bJ_se = col_double(),
      pJ = col_double(),
      LD_r = col_double()
    )
  ) %>%
    transmute(
      SNP,
      A1 = refA,
      A2 = NA_character_,
      freq,
      b,
      se,
      p,
      N = n,
      exposure = exposure,
      source = "cojo_jma"
    )
}

# === Read lead SNPs ===
read_lead_snps <- function(exposure) {
  lead_file <- file.path(clumped_dir, paste0(exposure, "_clumped_p1e5.txt"))
  if (!file.exists(lead_file)) {
    warning("Missing lead SNP file for: ", exposure)
    return(NULL)
  }
  
  df <- read_tsv(lead_file, show_col_types = FALSE)
  
  if (!all(c("SNP", "effect_allele.exposure", "other_allele.exposure",
             "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure") %in% names(df))) {
    warning("Required columns missing in lead SNP file for: ", exposure)
    return(NULL)
  }
  
  df %>%
    transmute(
      SNP,
      A1 = effect_allele.exposure,
      A2 = other_allele.exposure,
      freq = eaf.exposure,
      b = beta.exposure,
      se = se.exposure,
      p = pval.exposure,
      N = 452302,  # Replace if exposure sample size differs
      exposure = exposure,
      source = "lead"
    )
}

# === MAIN LOOP ===
for (exposure in exposures) {
  cat("\n========== Processing:", exposure, "==========\n")
  
  lead_snps <- read_lead_snps(exposure)
  cojo_snps <- read_cojo_jma(exposure)
  
  if (is.null(lead_snps) && is.null(cojo_snps)) {
    cat("Skipped - no data found for:", exposure, "\n")
    next
  }
  
  # Combine
  combined <- bind_rows(lead_snps, cojo_snps) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    arrange(SNP)
  
  if (nrow(combined) == 0) {
    cat("No combined SNPs for:", exposure, "\n")
    next
  }
  
  # === SAVE to individual exposure file ===
  out_file <- file.path(out_dir, paste0(exposure, "_final_instruments.txt"))
  
  write.table(
    combined %>% select(SNP, A1, A2, freq, b, se, p, N),
    file = out_file,
    quote = FALSE, sep = "\t", row.names = FALSE
  )
  
  cat("Saved combined instrument file for", exposure, "to:\n", out_file, "\n")
}

# ================================
# CHUNK 3.3 — Combine lead SNPs + new COJO SNPs per exposure (adds only non-lead COJO signals)
# ================================
## Only adds COJO SNPs not already in lead

# === CONFIG ===
clumped_dir    <- "project/data/eqtl/cis/formatted"         # lead/clumped files
cojo_dir       <- "project/outputs/cojo/cojo_slct_results"  # *_slct.jma.cojo
final_dir      <- "project/outputs/cojo/final_instruments"  # intermediate combined
raw_eqtl_dir   <- "project/data/eqtl/cis/formatted"         # raw eQTL reference
out_clean_dir  <- "project/outputs/cojo/final_instruments_cleaned"  # TwoSampleMR-ready

eqtl_sample_size <- 31684   # <-- set your eQTL study N here

# Create output folders if needed
if (!dir.exists(final_dir))     dir.create(final_dir, recursive = TRUE)
if (!dir.exists(out_clean_dir)) dir.create(out_clean_dir, recursive = TRUE)

# === EXPOSURES ===
exposures <- c("IL-1RN", "IL-1β", "IL-6", "IL-6R", "TNF-α", "TNF-R1", "TNF-R2", "PLAUR")

# --- helper: read COJO slct (handle `nan` in n) ---
read_cojo_jma <- function(exposure) {
  cojo_file <- file.path(cojo_dir, paste0(exposure, "_slct.jma.cojo"))
  if (!file.exists(cojo_file)) {
    warning("Missing COJO jma.cojo file for: ", exposure)
    return(NULL)
  }
  readr::read_table(
    cojo_file,
    na = c("", "NA", "nan", "NaN"),     # robust NA handling (n can be 'nan')
    show_col_types = FALSE,
    col_types = readr::cols(
      Chr       = readr::col_integer(),
      SNP       = readr::col_character(),
      bp        = readr::col_integer(),
      refA      = readr::col_character(),
      freq      = readr::col_double(),
      b         = readr::col_double(),
      se        = readr::col_double(),
      p         = readr::col_double(),
      n         = readr::col_number(),
      freq_geno = readr::col_double(),
      bJ        = readr::col_double(),
      bJ_se     = readr::col_double(),
      pJ        = readr::col_double(),
      LD_r      = readr::col_double()
    )
  ) %>%
    dplyr::transmute(
      SNP,
      effect_allele.exposure = refA,
      other_allele.exposure  = NA_character_,
      eaf.exposure           = freq,
      beta.exposure          = b,
      se.exposure            = se,
      pval.exposure          = p,
      samplesize.exposure    = n,          # may be NA; we overwrite to eqtl_sample_size later
      exposure               = exposure,
      id.exposure            = exposure
    )
}

# --- helper: read clumped lead SNPs (guard required columns) ---
read_lead_snps <- function(exposure) {
  lead_file <- file.path(clumped_dir, paste0(exposure, "_clumped_p1e5.txt"))
  if (!file.exists(lead_file)) {
    warning("Missing lead SNP file for: ", exposure)
    return(NULL)
  }
  df <- readr::read_tsv(lead_file, show_col_types = FALSE)
  required_cols <- c(
    "SNP","effect_allele.exposure","other_allele.exposure",
    "eaf.exposure","beta.exposure","se.exposure","pval.exposure"
  )
  if (!all(required_cols %in% names(df))) {
    warning("Required columns missing in lead SNP file for: ", exposure)
    return(NULL)
  }
  df %>%
    dplyr::transmute(
      SNP,
      effect_allele.exposure = effect_allele.exposure,
      other_allele.exposure  = other_allele.exposure,
      eaf.exposure           = eaf.exposure,
      beta.exposure          = beta.exposure,
      se.exposure            = se.exposure,
      pval.exposure          = pval.exposure,
      samplesize.exposure    = 452302,    # lead/clumped N; overwritten to eQTL N in enrichment
      exposure               = exposure,
      id.exposure            = exposure
    )
}

# --- helper: read raw eQTL reference for enrichment ---
read_eqtl_ref <- function(exposure) {
  raw_file <- file.path(raw_eqtl_dir, paste0(exposure, "_EXPOSURE_final.txt"))
  if (!file.exists(raw_file)) {
    warning("Missing raw eQTL file for: ", exposure)
    return(NULL)
  }
  readr::read_tsv(raw_file, show_col_types = FALSE) %>%
    dplyr::select(
      SNP,
      effect_allele_raw = effect_allele.exposure,
      other_allele_raw  = other_allele.exposure,
      eaf_raw           = eaf.exposure,
      beta_raw          = beta.exposure,
      se_raw            = se.exposure,
      pval_raw          = pval.exposure
    ) %>%
    dplyr::distinct(SNP, .keep_all = TRUE)
}

# === MAIN LOOP ===
for (exposure in exposures) {
  cat("\n========== Processing:", exposure, "==========\n")
  
  # 1) Read inputs
  lead_snps <- read_lead_snps(exposure)
  cojo_raw  <- read_cojo_jma(exposure)
  
  if (is.null(lead_snps) && is.null(cojo_raw)) {
    cat("Skipped - no data found for:", exposure, "\n")
    next
  }
  
  # 2) Keep only COJO SNPs not already in lead (if both present)
  if (!is.null(cojo_raw) && !is.null(lead_snps)) {
    cojo_snps <- dplyr::filter(cojo_raw, !SNP %in% lead_snps$SNP)
  } else if (!is.null(cojo_raw)) {
    cojo_snps <- cojo_raw
  } else {
    cojo_snps <- NULL
  }
  
  # 3) Combine lead + (non-lead) COJO and save intermediate file
  combined <- dplyr::bind_rows(lead_snps, cojo_snps) %>%
    dplyr::distinct(SNP, .keep_all = TRUE)
  
  if (nrow(combined) == 0) {
    cat("No combined SNPs for:", exposure, "\n")
    next
  }
  
  out_intermediate <- file.path(final_dir, paste0(exposure, "_final_instruments.txt"))
  write.table(
    combined,
    file = out_intermediate,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("Saved combined instruments (intermediate):", out_intermediate, "\n")
  
  # 4) Enrich with raw eQTL (fills missing alleles/eaf/beta/se/p; set uniform eQTL N)
  eqtl_ref <- read_eqtl_ref(exposure)
  if (is.null(eqtl_ref)) {
    cat("No eQTL reference to enrich; skipping TwoSampleMR-ready write for:", exposure, "\n")
    next
  }
  
  enriched <- combined %>%
    dplyr::left_join(eqtl_ref, by = "SNP") %>%
    dplyr::mutate(
      other_allele.exposure = ifelse(is.na(other_allele.exposure), other_allele_raw, other_allele.exposure),
      eaf.exposure          = ifelse(is.na(eaf.exposure) | eaf.exposure == 0, eaf_raw, eaf.exposure),
      beta.exposure         = ifelse(is.na(beta.exposure), beta_raw, beta.exposure),
      se.exposure           = ifelse(is.na(se.exposure),   se_raw,   se.exposure),
      pval.exposure         = ifelse(is.na(pval.exposure), pval_raw, pval.exposure),
      
      # Unify N to eQTL study size for final files (this mirrored your Chunk C behavior)
      samplesize.exposure   = eqtl_sample_size,
      
      # Keep TwoSampleMR-friendly metadata columns
      exposure              = exposure,
      id.exposure           = exposure,
      mr_keep.exposure      = TRUE,
      pval_origin.exposure  = "reported"
    ) %>%
    dplyr::select(
      SNP,
      effect_allele.exposure,
      other_allele.exposure,
      eaf.exposure,
      beta.exposure,
      se.exposure,
      pval.exposure,
      samplesize.exposure,
      exposure,
      mr_keep.exposure,
      pval_origin.exposure,
      id.exposure
    )
  
  # 5) Drop incomplete rows (essentials)
  cleaned <- enriched %>%
    dplyr::filter(
      !is.na(SNP),
      !is.na(effect_allele.exposure),
      !is.na(other_allele.exposure),
      !is.na(beta.exposure),
      !is.na(se.exposure),
      !is.na(pval.exposure),
      !is.na(samplesize.exposure)
    )
  
  if (nrow(cleaned) == 0) {
    cat("No valid SNPs remaining after enrichment for:", exposure, "\n")
    next
  }
  
  # 6) Write TwoSampleMR-ready per-exposure file
  out_clean <- file.path(out_clean_dir, paste0(exposure, "_TwoSampleMR_ready.txt"))
  write.table(
    cleaned,
    file = out_clean,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("Saved TwoSampleMR-ready file:", out_clean, "\n")
}

# ================================
# CHUNK 4 — Verify combined file; enrich missing fields from original GWAS
# ================================

# Step 1: Load original GWAS file with proper column names
original_gwas <- read_tsv("project/data/gwas/raw/Loneliness_EurRel_Imputed.txt", col_names = FALSE)

# Assign column names manually (based on your data)
colnames(original_gwas) <- c(
  "SNP", "CHR", "BP", "MAF", "A1", "A2", "EAF", "INFO", "OR", "OR_se",
  "beta", "se", "pval", "n_cases", "prop_cases", "n_controls"
)

# Step 3: Remove duplicate SNPs in the original GWAS file
original_gwas_dedup <- original_gwas %>%
  distinct(SNP, .keep_all = TRUE)

# Step 4: Load combined SNP file
combined_snps <- read_tsv("project/outputs/cojo/combined_snps_with_conditional_raw.txt")

# Step 5: Merge and fill missing data
combined_snps_fixed <- combined_snps %>%
  left_join(
    original_gwas_dedup %>%
      select(SNP,
             A1_gwas = A1,
             A2_gwas = A2,
             EAF,
             beta,
             se_gwas = se,
             pval),
    by = "SNP"
  ) %>%
  mutate(
    A2   = ifelse(is.na(A2), A2_gwas, A2),
    freq = ifelse(is.na(freq), EAF, freq),
    b    = ifelse(is.na(b), beta, b),
    se   = ifelse(is.na(se), se_gwas, se),
    p    = ifelse(is.na(p), pval, p)
  ) %>%
  select(SNP, A1, A2, freq, b, se, p, N)

# Step 6: Remove any rows with missing values
combined_snps_final <- combined_snps_fixed %>%
  filter(!is.na(SNP) & !is.na(A1) & !is.na(A2) & !is.na(freq) & !is.na(b) & !is.na(se) & !is.na(p) & !is.na(N))

# Step 7: Save cleaned file
write.table(
  combined_snps_final,
  file = "project/outputs/cojo/combined_snps_with_conditional_cleaned.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ================================
# CHUNK 4.2 — Enrich + reformat final_instruments to TwoSampleMR exposure format (per trait)
# ================================

## One cleaned + properly formatted file per exposure

# === CONFIG ===
final_dir <- "project/outputs/cojo/final_instruments"
eqtl_raw_dir <- "project/data/eqtl/cis/formatted"
out_clean_dir <- "project/outputs/cojo/final_instruments_cleaned"

if (!dir.exists(out_clean_dir)) dir.create(out_clean_dir, recursive = TRUE)

# === EXPOSURE LIST ===
exposures <- c("IL-1RN", "IL-1β", "IL-6", "IL-6R", "TNF-α", "TNF-R1", "TNF-R2", "PLAUR")

for (exposure in exposures) {
  cat("\n========== Processing:", exposure, "==========\n")
  
  # --- 1. Load original raw eQTL summary for enrichment ---
  eqtl_file <- file.path(eqtl_raw_dir, paste0(exposure, "_EXPOSURE_final.txt"))
  if (!file.exists(eqtl_file)) {
    cat("Missing raw eQTL file for:", exposure, "\n")
    next
  }
  
  eqtl_raw <- read_tsv(eqtl_file, show_col_types = FALSE)
  
  eqtl_ref <- eqtl_raw %>%
    select(
      SNP,
      A1_gwas = effect_allele.exposure,
      A2_gwas = other_allele.exposure,
      EAF = eaf.exposure,
      beta_gwas = beta.exposure,
      se_gwas = se.exposure,
      pval_gwas = pval.exposure
    ) %>%
    distinct(SNP, .keep_all = TRUE)
  
  # --- 2. Load combined final_instruments (lead + COJO) ---
  combined_file <- file.path(final_dir, paste0(exposure, "_final_instruments.txt"))
  if (!file.exists(combined_file)) {
    cat("Missing combined instruments for:", exposure, "\n")
    next
  }
  
  combined <- read_tsv(combined_file, show_col_types = FALSE)
  
  # --- 3. Enrich missing values from eQTL reference ---
  enriched <- combined %>%
    left_join(eqtl_ref, by = "SNP") %>%
    mutate(
      effect_allele.exposure = ifelse(is.na(effect_allele.exposure), A1_gwas, effect_allele.exposure),
      other_allele.exposure = ifelse(is.na(other_allele.exposure), A2_gwas, other_allele.exposure),
      eaf.exposure = ifelse(is.na(eaf.exposure), EAF, eaf.exposure),
      beta.exposure = ifelse(is.na(beta.exposure), beta_gwas, beta.exposure),
      se.exposure = ifelse(is.na(se.exposure), se_gwas, se.exposure),
      pval.exposure = ifelse(is.na(pval.exposure), pval_gwas, pval.exposure),
      samplesize.exposure = 452302,
      id.exposure = exposure,
      exposure = exposure
    ) %>%
    select(SNP,
           effect_allele.exposure,
           other_allele.exposure,
           eaf.exposure,
           beta.exposure,
           se.exposure,
           pval.exposure,
           samplesize.exposure,
           exposure,
           id.exposure)
  
  # --- 4. Remove rows with truly essential missing values ---
  cleaned <- enriched %>%
    filter(
      !is.na(SNP),
      !is.na(effect_allele.exposure),
      !is.na(other_allele.exposure),
      !is.na(beta.exposure),
      !is.na(se.exposure),
      !is.na(pval.exposure),
      !is.na(samplesize.exposure)
    )
  
  if (nrow(cleaned) == 0) {
    cat("No valid SNPs remaining after enrichment for:", exposure, "\n")
    next
  }
  
  # --- 5. Format to TwoSampleMR exposure format ---
  formatted <- format_data(
    cleaned,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure",
    phenotype_col = "id.exposure"
  )
  
  # --- 6. Save as TwoSampleMR-ready txt ---
  out_file <- file.path(out_clean_dir, paste0(exposure, "_TwoSampleMR_ready.txt"))
  
  write.table(
    formatted,
    file = out_file,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  cat("Cleaned & saved in TwoSampleMR-ready format:", out_file, "\n")
}

# ================================
# CHUNK 5 — Summarize COJO COND results: secondary signals (per lead SNP; P threshold)
# ================================

# File paths
cond_results_dir <- "project/outputs/cojo/cojo_cond_results/"
region_dir <- "project/outputs/cojo/cojo_regions/500kb/"
output_txt <- "project/outputs/cojo/cojo_conditional_summary.txt"
output_csv <- "project/outputs/cojo/secondary_signals.csv"

# List of lead SNPs, adjust as needed
lead_snps <- c(
  "rs6430286", "rs74338595", "rs4465966", "rs2069117", "rs10456089", "rs7770860",
  "rs773020", "rs13291079", "rs2149351", "rs62085660", "rs613872", "rs1022688"
)

# Result containers
results_list <- list()
summary_lines <- c("COJO Conditional Analysis Summary", "==============================\n")

# Loop through each lead SNP
for (lead in lead_snps) {
  given_file <- file.path(cond_results_dir, paste0(lead, "_cond.given.cojo"))
  cma_file <- file.path(cond_results_dir, paste0(lead, "_cond.cma.cojo"))
  region_file <- file.path(region_dir, paste0(lead, "_region_500kb.txt"))
  
  if (file.exists(given_file) && file.exists(cma_file) && file.exists(region_file)) {
    # Read files
    given <- read.table(given_file, header = TRUE)
    cma <- read.table(cma_file, header = TRUE)
    region <- read.table(region_file, header = FALSE,
                         col.names = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N"))
    
    # Merge region info (alleles) into cma
    cma_annotated <- cma %>%
      left_join(region %>% select(SNP, A1, A2), by = "SNP")
    
    # Filter for genome-wide significant secondary hits, p value can be adjusted from 5e-8 to 1e-5
    secondary_hits <- cma_annotated %>%
      filter(!is.na(pC) & pC < 1e-5) %>% ##Adjust pC here, pC < 5e-8)
      arrange(pC) %>%
      mutate(Lead_SNP = lead) %>%
      select(Lead_SNP, Secondary_SNP = SNP, CHR = Chr, BP = bp, A1, A2, bC, se = bC_se, pC) #extracts necessary MR columns
    
    if (nrow(secondary_hits) > 0) {
      summary_lines <- c(
        summary_lines,
        paste0("Lead SNP: ", lead),
        paste0("  Secondary signals (P < 1e-5): ", nrow(secondary_hits)) ##Adjust p comment here, pC < 5e-8)
      )
      
      for (i in seq_len(nrow(secondary_hits))) {
        row <- secondary_hits[i, ]
        summary_lines <- c(
          summary_lines,
          paste0("    ", row$Secondary_SNP,
                 " (CHR ", row$CHR, ":", row$BP, "), ",
                 row$A1, "/", row$A2, ", bC = ", round(row$bC, 5),
                 ", SE = ", round(row$se, 5),
                 ", P = ", formatC(row$pC, format = "e", digits = 3))
        )
      }
      
      summary_lines <- c(summary_lines, "")
      results_list[[lead]] <- secondary_hits
    }
  } else {
    warning(paste("Missing file(s) for lead SNP:", lead)) #reports if any files are missing for the lead SNPs listed
  }
}

#Combine and export
secondary_signals <- bind_rows(results_list)

# Save outputs
writeLines(summary_lines, con = output_txt)
write.csv(secondary_signals, file = output_csv, row.names = FALSE)

cat("Summary written to:\n", output_txt, "\n", output_csv, "\n")

# ================================
# CHUNK 0 — (Record-keeping) Format the full original GWAS once into COJO-ready columns
# ================================

# EDIT THIS
input_file  <- "project/data/gwas/raw/Loneliness_EurRel_Imputed.txt" # GWAS file
output_file <- "project/outputs/cojo/isolation_for_cojo.txt"         # Output file for COJO
sample_size <- 452302                                                # Replace with your actual N

# Load your raw GWAS file
df <- read_tsv(input_file)

# Reformat to COJO-required format
cojo_df <- df %>%
  transmute(
    SNP = SNP,
    A1 = ALLELE1,  # Effect allele
    A2 = ALLELE0,  # Other allele
    freq = A1FREQ,
    b = BETA,
    se = SE,
    p = P_BOLT_LMM_INF,  
    N = sample_size
  ) %>%
  filter(!is.na(freq) & !is.na(b) & !is.na(se) & !is.na(p))  # Remove any rows with missing data

# Write in space-separated format (not tab-delimited)
write.table(cojo_df, file = output_file, sep = " ", quote = FALSE, row.names = FALSE)

