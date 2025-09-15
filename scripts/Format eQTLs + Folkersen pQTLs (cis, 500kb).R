---
title: "Format eQTLs + Folkersen pQTLs (cis, 500kb)"
author: "Coral Romualdo-Perez"
date: "2025-06-11"
output: html_document
params:
  panel_name: "Folkersen_2020"
window_kb: 500
out_dir: "data/standardized/Folkersen"
ld_bfile: "/Users/coralrom/Desktop/LD Clumping Data/EUR"
---
  
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(glue)
  library(TwoSampleMR); library(ieugwasr); library(genetics.binaRies)
})

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
plink_path <- genetics.binaRies::get_plink_binary()

# ----- Header-agnostic standardizer ----------------------------------------
.pick <- function(dt, keys) {
  nm <- tolower(names(dt))
  keys <- tolower(keys)
  j <- which(nm %in% keys)
  if (!length(j)) return(NULL)
  dt[[j[1]]]
}

standardize_any <- function(path, default_N = NA_real_) {
  raw <- data.table::fread(path, showProgress = FALSE)
  nm  <- tolower(names(raw)); setnames(raw, nm)
  
  get_num <- function(keys) suppressWarnings(as.numeric(.pick(raw, keys)))
  
  out <- tibble(
    snp  = as.character(.pick(raw, c("snp","rsid","rsids","id","hm_rsid","variant_id","hm_variant_id"))),
    chr  = suppressWarnings(as.integer(gsub("^chr","", as.character(.pick(raw, c("chromosome","chr","chrom","chr_hg19","hm_chrom")))))),
    pos  = suppressWarnings(as.integer(.pick(raw, c("base_pair_location","pos","position","bp","pos_hg19","bp_hg19","hm_pos")))),
    ea   = toupper(as.character(.pick(raw, c("effect_allele","ea","a1","allele1","tested_allele","hm_effect_allele")))),
    oa   = toupper(as.character(.pick(raw, c("other_allele","oa","a2","allele2","ref_allele","non_effect_allele","hm_other_allele","allele0")))),
    beta = get_num(c("beta","b","effect","beta.exposure","beta.outcome","hm_beta")),
    se   = get_num(c("standard_error","se","stderr","se_beta","se.exposure","se.outcome")),
    pval = get_num(c("p","pval","p_value","pvalue","p-value","p_ldl","p_cad")),
    eaf  = get_num(c("eaf","effect_allele_frequency","eaf.exposure","eaf.outcome","hm_effect_allele_frequency","maf")),
    n    = get_num(c("n","samplesize","sample_size","samplesize.exposure","samplesize.outcome"))
  )
  
  # If pval missing but beta/se present:
  miss_p <- is.na(out$pval) & !is.na(out$beta) & !is.na(out$se)
  out$pval[miss_p] <- pmin(pmax(2*pnorm(-abs(out$beta[miss_p]/out$se[miss_p])), 1e-300), 1)
  
  # Fill N if known default provided
  if (!all(is.na(out$n)) || !is.na(default_N)) out$n[is.na(out$n)] <- default_N
  
  out %>%
    filter(!is.na(snp) & !is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se)) %>%
    mutate(eaf = ifelse(!is.na(eaf) & (eaf < 0 | eaf > 1), NA_real_, eaf),
           ea = toupper(ea), oa = toupper(oa)) %>%
    arrange(se) %>%
    distinct(snp, .keep_all = TRUE)
}

cis_filter <- function(df, chr, start_bp, end_bp, flank_kb = 500) {
  w <- flank_kb * 1000
  df %>% filter(chr == !!chr, pos >= start_bp - w, pos <= end_bp + w)
}

# ----- Explicit mapping: gene ↔ file ↔ coordinates (GRCh37) -----------------
manifest <- tribble(
  ~gene,   ~chr, ~start,      ~end,        ~path,                                                                 ~panel_N,
  "IL1RA",    2, 113864791, 113891593, here::here("data-raw/folkersen/IL1RA_harmonised.qc.txt"),                 21758,
  "IL1B",     2, 113587328, 113594480, here::here("data-raw/folkersen/IL1B_harmonised.txt"),                     21758,
  "IL6",      7,  22765503,  22771621, here::here("data-raw/folkersen/IL6_harmonised.txt"),                      21758,
  "IL6R",     1, 154377669, 154441926, here::here("data-raw/folkersen/IL6RA_harmonised.qc.txt"),                 21758,
  "TNF",      6,  31543344,  31546113, here::here("data-raw/folkersen/TNF_harmonised.txt"),                      21758,
  "TNFR1",   12,   6437923,   6451280, here::here("data-raw/folkersen/TNFR1_harmonised.qc.txt"),                 21758,
  "TNFR2",    1,  12227060,  12269285, here::here("data-raw/folkersen/TNFR2_harmonised.qc.txt"),                 21758,
  "PLAUR",   19,  44150247,  44174699, here::here("data-raw/folkersen/UPAR_harmonised.qc.txt"),                  21758
)

manifest

# ----- pQTL: standardize → cis-window → write TSV + TwoSampleMR exposure ----
pqtl_out_dir <- file.path(params$out_dir, "pQTL_cis500kb")
dir.create(pqtl_out_dir, showWarnings = FALSE, recursive = TRUE)

pqtl_rows <- pmap_dfr(manifest, function(gene, chr, start, end, path, panel_N) {
  std <- standardize_any(path, default_N = panel_N)
  cis <- cis_filter(std, chr = chr, start_bp = start, end_bp = end, flank_kb = params$window_kb)
  
  # Save raw cis table
  cis_path <- file.path(pqtl_out_dir, glue("{gene}_cis{params$window_kb}kb.tsv"))
  fwrite(cis, cis_path, sep = "\t")
  
  # TwoSampleMR exposure format
  exp <- TwoSampleMR::format_data(
    cis %>% transmute(SNP = snp, beta, se,
                      effect_allele = ea, other_allele = oa,
                      eaf = eaf, pval = pval),
    type = "exposure",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval"
  )
  exp$exposure <- gene; exp$id.exposure <- gene
  
  exp_path <- file.path(pqtl_out_dir, glue("{gene}_EXPOSURE_format.tsv"))
  fwrite(exp, exp_path, sep = "\t")
  
  tibble(gene, chr, start, end,
         n_cis = nrow(cis),
         cis_tsv = cis_path,
         exposure_tsv = exp_path)
})

pqtl_rows

# ----- Clump exposures at two thresholds -----------------------------------
clump_dir <- file.path(pqtl_out_dir, "clumped")
dir.create(clump_dir, showWarnings = FALSE)

clump_once <- function(exp_path, gene, pthresh, r2, kb) {
  d <- fread(exp_path)
  stopifnot(all(c("SNP","pval.exposure") %in% names(d)))
  base <- d %>% select(rsid = SNP, pval = pval.exposure) %>% mutate(id = gene) %>% filter(!is.na(rsid), !is.na(pval))
  base <- base %>% filter(pval < pthresh)
  if (!nrow(base)) return(NULL)
  out <- ieugwasr::ld_clump(
    d = base, clump_kb = kb, clump_r2 = r2, clump_p = 1,
    bfile = params$ld_bfile, plink_bin = plink_path
  )
  if (!nrow(out)) return(NULL)
  joined <- inner_join(d, out, by = c("SNP" = "rsid"))
  out_path <- file.path(clump_dir, glue("{gene}_clumped_p{format(pthresh,scientific=TRUE)}_r2_{r2}_kb{kb}.tsv"))
  fwrite(joined, out_path, sep = "\t")
  out_path
}

clump_summary <- pqtl_rows %>%
  mutate(
    clumped_p5e8 = map_chr(exposure_tsv, ~clump_once(.x, gene=.env$gene, pthresh=5e-8, r2=0.001, kb=1000) %||% NA_character_),
    clumped_p1e5 = map_chr(exposure_tsv, ~clump_once(.x, gene=.env$gene, pthresh=1e-5, r2=0.01,  kb=1000) %||% NA_character_)
  )

clump_summary

# ----- OPTIONAL: cis-eQTL extraction (no MR coercion if only Z is available) -
# Provide your eQTL big file & dictionary (must include Gene, SNPChr, SNPPos, SNP, alleles, p-value).
# Example:
# eqtl_big <- fread("data-raw/eQTLGen_trans.txt.gz") # or cis if you have it

# cytokine reference (GRCh37 coords)
cytokine_genes <- tribble(
  ~gene_name, ~ensembl_id,       ~chr, ~start,      ~end,
  "CRP",      "ENSG00000132693",   1, 159682079, 159684379,
  "TNF",      "ENSG00000232810",   6,  31543344,  31546113,
  "IL1B",     "ENSG00000125538",   2, 113587328, 113594480,
  "IL6",      "ENSG00000136244",   7,  22765503,  22771621,
  "TNFR1",    "ENSG00000067182",  12,   6437923,   6451280,
  "TNFR2",    "ENSG00000028137",   1,  12227060,  12269285,
  "IL1RN",    "ENSG00000136689",   2, 113864791, 113891593,
  "IL6R",     "ENSG00000160712",   1, 154377669, 154441926,
  "PLAUR",    "ENSG00000011422",  19,  44150247,  44174699
)

# If you want to run this step: set `run_eqtl <- TRUE` and provide `eqtl_big`.
run_eqtl <- FALSE
if (run_eqtl) {
  eqtl_out_dir <- file.path(params$out_dir, "eQTL_cis500kb")
  dir.create(eqtl_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  eqtl_rows <- pmap_dfr(cytokine_genes, function(gene_name, ensembl_id, chr, start, end) {
    cis_start <- start - params$window_kb*1000
    cis_end   <- end   + params$window_kb*1000
    # Expect columns: Gene (Ensembl), SNPChr, SNPPos, SNP, AssessedAllele, OtherAllele, Pvalue, Zscore (if present)
    d <- eqtl_big %>%
      filter(Gene == ensembl_id,
             SNPChr == chr,
             SNPPos >= cis_start, SNPPos <= cis_end) %>%
      mutate(chr = as.integer(SNPChr),
             pos = as.integer(SNPPos),
             snp = as.character(SNP))
    if (!nrow(d)) return(tibble(gene = gene_name, n_cis = 0, cis_tsv = NA_character_))
    cis <- d %>%
      transmute(snp, chr, pos,
                ea = toupper(AssessedAllele), oa = toupper(OtherAllele),
                pval = as.numeric(Pvalue),
                z = suppressWarnings(as.numeric(Zscore))) %>%
      arrange(pval)
    out_path <- file.path(eqtl_out_dir, glue("{gene_name}_cis{params$window_kb}kb.tsv"))
    fwrite(cis, out_path, sep = "\t")
    tibble(gene = gene_name, n_cis = nrow(cis), cis_tsv = out_path)
  })
  eqtl_rows
}



# ================================
# Ferkingstad pQTL: Liftover + Format (hg38 -> hg19) + cis windows + MR exposure
# ================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(tidyr)
  library(GenomicRanges); library(rtracklayer); library(TwoSampleMR)
})

# ---- paths (edit) ----
IN_DIR      <- "/Volumes/Marohu/DeCode PQTL"                 # hg38 input files
OUT_DIR     <- "/Volumes/Marohu/DeCode PQTL/Hg37"            # all hg19 outputs
CHAIN_PATH  <- "/Volumes/Marohu/hg38ToHg19.over.chain.gz"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) Helpers
# =========================
to_num <- function(x) suppressWarnings(as.numeric(x))

is_rsid <- function(x) !is.na(x) & grepl("^rs[0-9]+$", x)

# Standardize a Ferkingstad/deCODE-like header (hg38)
standardize_decode_v38 <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  # lowercase for matching; keep originals in dt_raw if you prefer
  names(dt) <- tolower(gsub("\\.", "", names(dt)))
  
  # Minimal required columns in these files:
  need <- c("chrom","pos","effectallele","otherallele","beta","se","pval","n")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing columns in ", basename(path), ": ", paste(miss, collapse=", "))
  
  # Map common extras
  rsids  <- dt[["rsids"]]
  pval   <- dt[["pval"]]
  if (is.null(pval)) pval <- NA_real_
  
  # If p missing but beta/se present, compute two-sided p
  if (all(is.na(pval)) && !all(is.na(dt$beta)) && !all(is.na(dt$se))) {
    pval <- 2*pnorm(-abs(dt$beta/dt$se))
  }
  
  # eaf from ImpMAF if present (deCODE reports imputed MAF); safe for coloc;
  # for MR, it's optional and not always used downstream
  eaf <- dt[["impmaf"]]
  if (!is.null(eaf)) eaf <- to_num(eaf)
  
  tibble(
    chrom38 = as.character(dt$chrom),     # e.g., "chr1"
    pos38   = as.integer(dt$pos),
    rsids   = if ("rsids" %in% names(dt)) as.character(dt$rsids) else NA_character_,
    ea      = toupper(as.character(dt$effectallele)),
    oa      = toupper(as.character(dt$otherallele)),
    beta    = to_num(dt$beta),
    se      = to_num(dt$se),
    pval    = to_num(pval),
    n       = to_num(dt$n),
    eaf     = eaf
  )
}

# Liftover hg38 -> hg19; keep only 1:1 mappings
liftover_to_hg19 <- function(df38, chain_path = CHAIN_PATH) {
  gr38 <- GRanges(seqnames = df38$chrom38, ranges = IRanges(df38$pos38, df38$pos38))
  ch   <- import.chain(chain_path)
  lo   <- liftOver(gr38, ch)
  mapn <- lengths(lo)
  
  keep <- which(mapn == 1L)
  out  <- df38[keep, , drop = FALSE]
  lo1  <- unlist(lo[keep])
  
  out$chr_hg19 <- as.character(seqnames(lo1))
  out$pos_hg19 <- start(lo1)
  
  list(
    mapped  = out,
    summary = list(total = nrow(df38), onetoone = length(keep),
                   multi = sum(mapn > 1L), unmapped = sum(mapn == 0L))
  )
}

# Build “preFormat” frame (prefer rsID; else synthesize chr:pos:EA:OA)
build_preformat <- function(df_mapped) {
  snp <- ifelse(is_rsid(df_mapped$rsids) & df_mapped$rsids != "NA",
                df_mapped$rsids,
                paste0(gsub("^chr","", df_mapped$chr_hg19), ":", df_mapped$pos_hg19, ":",
                       df_mapped$ea, ":", df_mapped$oa))
  tibble(
    SNP           = snp,
    effect_allele = df_mapped$ea,
    other_allele  = df_mapped$oa,
    beta          = df_mapped$beta,
    se            = df_mapped$se,
    pval          = pmin(pmax(df_mapped$pval, 1e-300), 1),  # clamp
    eaf           = df_mapped$eaf,
    samplesize    = df_mapped$n,
    chr_hg19      = gsub("^chr","", df_mapped$chr_hg19),
    pos_hg19      = df_mapped$pos_hg19,
    rsids         = ifelse(df_mapped$rsids %in% c("", "NA"), NA, df_mapped$rsids)
  ) %>%
    # keep best (smallest SE) if duplicate SNPs
    arrange(se) %>% distinct(SNP, .keep_all = TRUE)
}

# =========================
# 2) Run liftover + preFormat on all input files
# =========================
in_files <- list.files(IN_DIR, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
log_rows <- list()

for (f in in_files) {
  message("Processing: ", basename(f))
  std38 <- standardize_decode_v38(f)
  lo    <- liftover_to_hg19(std38)
  
  # write hg19 raw with chr/pos
  raw_out <- file.path(OUT_DIR, sub("\\.txt(\\.gz)?$", "_hg19_raw.tsv", basename(f)))
  fwrite(lo$mapped, raw_out, sep = "\t")
  
  # write preFormat
  pre  <- build_preformat(lo$mapped)
  pre_out <- file.path(OUT_DIR, sub("\\.txt(\\.gz)?$", "_hg19_preFormat.tsv", basename(f)))
  fwrite(pre, pre_out, sep = "\t")
  
  log_rows[[length(log_rows)+1]] <- data.frame(
    file = basename(f),
    out_raw_hg19 = basename(raw_out),
    out_preFormat_hg19 = basename(pre_out),
    total = lo$summary$total,
    mapped_1to1 = lo$summary$onetoone,
    multi = lo$summary$multi,
    unmapped = lo$summary$unmapped,
    success_pct = round(100*lo$summary$onetoone/lo$summary$total, 2),
    stringsAsFactors = FALSE
  )
}
if (length(log_rows)) {
  liftover_summary <- bind_rows(log_rows)
  fwrite(liftover_summary, file.path(OUT_DIR, "liftover_summary.csv"))
  print(liftover_summary)
}

# =========================
# 3) cis-window extraction per gene (±500 kb)
# =========================
cis_dir <- file.path(OUT_DIR, "cis500kb_filtered")
dir.create(cis_dir, showWarnings = FALSE)

# hg19 (GRCh37) coordinates (checked) — consistent names
genes <- tribble(
  ~gene,   ~chr, ~start,      ~end,
  "IL1RA",    2, 113864791, 113891593,
  "IL1B",     2, 113587328, 113594480,
  "IL6",      7,  22765503,  22771621,
  "IL6R",     1, 154377669, 154441926,
  "TNF",      6,  31543344,  31546113,
  "TNFR1",   12,   6437923,   6451280,
  "TNFR2",    1,  12227060,  12269285
  # "PLAUR",  19,  44150247,  44174699
)

# find a preFormat file for a gene (case-insensitive)
find_preformat_for_gene <- function(g) {
  pat <- paste0("^", g, ".*_hg19_preFormat\\.tsv$")
  hits <- list.files(OUT_DIR, pattern = pat, full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) NA_character_ else hits[1]
}

cis_rows <- pmap_dfr(genes, function(gene, chr, start, end) {
  f <- find_preformat_for_gene(gene)
  if (is.na(f)) return(tibble(gene=gene, n_cis=0, cis_file=NA_character_))
  d <- fread(f)
  # ensure types
  d$chr_hg19 <- as.character(d$chr_hg19)
  d$pos_hg19 <- as.integer(d$pos_hg19)
  
  w <- 500*1000L
  s <- max(1L, start - w)
  e <- end + w
  
  cis <- d %>%
    filter(chr_hg19 == as.character(chr),
           !is.na(pos_hg19),
           pos_hg19 >= s, pos_hg19 <= e)
  
  out_path <- file.path(cis_dir, paste0(gene, "_cis500kb_filtered.tsv"))
  fwrite(cis, out_path, sep = "\t")
  tibble(gene=gene, n_cis=nrow(cis), cis_file=out_path)
})

cis_rows
#> per-gene filtered cis tables written to OUT_DIR/cis500kb_filtered

# =========================
# 4) TwoSampleMR exposure format (from cis tables)
# =========================
tsmr_dir <- file.path(OUT_DIR, "TwoSampleMR_exposure")
dir.create(tsmr_dir, showWarnings = FALSE)

exp_rows <- pmap_dfr(cis_rows, function(gene, n_cis, cis_file) {
  if (is.na(cis_file) || !file.exists(cis_file) || n_cis == 0)
    return(tibble(gene=gene, exposure_file=NA_character_, n_kept=0))
  
  d <- fread(cis_file, data.table = FALSE, check.names = FALSE)
  
  need <- c("SNP","beta","se","effect_allele","other_allele","pval")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    warning("Missing columns in ", basename(cis_file), ": ", paste(miss, collapse=", "))
    return(tibble(gene=gene, exposure_file=NA_character_, n_kept=0))
  }
  
  # uppercase alleles (safety)
  d$effect_allele <- toupper(d$effect_allele)
  d$other_allele  <- toupper(d$other_allele)
  
  # format for TwoSampleMR (exposure)
  exp <- TwoSampleMR::format_data(
    d %>% transmute(SNP, beta, se,
                    effect_allele, other_allele,
                    eaf = if ("eaf" %in% names(d)) eaf else NA_real_,
                    pval),
    type = "exposure",
    snp_col="SNP", beta_col="beta", se_col="se",
    effect_allele_col="effect_allele", other_allele_col="other_allele",
    eaf_col="eaf", pval_col="pval"
  )
  exp$exposure <- gene; exp$id.exposure <- gene
  
  # deduplicate by SNP, keep most significant
  if ("pval.exposure" %in% names(exp)) {
    exp <- exp %>% arrange(pval.exposure) %>% distinct(SNP, .keep_all = TRUE)
  } else {
    exp <- exp %>% distinct(SNP, .keep_all = TRUE)
  }
  
  out <- file.path(tsmr_dir, paste0(gene, "_TwoSampleMR_ready_unclumped.tsv"))
  fwrite(exp, out, sep = "\t")
  tibble(gene=gene, exposure_file=out, n_kept=nrow(exp))
})

exp_rows
#> TwoSampleMR-ready exposure files written to OUT_DIR/TwoSampleMR_exposure
